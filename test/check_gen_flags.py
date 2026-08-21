"""Sanity check for the gen-level provenance tags on real samples.

Run it inside the key4hep stack (it needs podio), e.g.::

    python test/check_gen_flags.py -f ztt_2M --n-files 200 --n-workers 60
    python test/check_gen_flags.py --input-list /path/to/one.root

It reports, over the scanned events:

  * which extra neutrals show up (PDG breakdown) and how they distribute across
    decay-mode IDs — the point being that a ``tau -> pi K0_L nu`` keeps ID 0,
  * how many gen taus are flagged as secondary and through which originPDG,
  * the generatorStatus of every tau that has another tau in its ancestry: this
    is the assumption findAllGenTaus relies on (only status==2 is kept), so if
    radiative taus show up with another status they are NOT in the tree at all,
  * the origin breakdown of all generator-level photons.

Extra neutrals and radiative taus are rare, so the scan is parallelised over
files: use ``--n-workers`` to cover enough events for them to appear. Example
``(file, event)`` pairs are printed for each rare category so they can be
inspected afterwards with ``Notebooks/explore_anomalies.py``.
"""

import argparse
import logging
import os
import sys
from collections import Counter, defaultdict
from multiprocessing import Pool

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from podio import root_io

from modules import myutils, tauReco

ORIGIN_NAMES = {
    tauReco.PHOTON_ORIGIN_NOT_A_PHOTON: "not-a-photon",
    tauReco.PHOTON_ORIGIN_PI0: "pi0",
    tauReco.PHOTON_ORIGIN_TAU_FSR: "tau-FSR",
    tauReco.PHOTON_ORIGIN_CHARGED_RAD: "charged-rad",
    tauReco.PHOTON_ORIGIN_OTHER: "other",
    tauReco.PHOTON_ORIGIN_SIMULATION: "simulation",
}

MAX_EXAMPLES = 5


def _pdg_name(pdg):
    try:
        from particle import Particle

        return Particle.from_pdgid(pdg).name
    except Exception:
        return str(pdg)


def _has_tau_ancestor_beyond_copies(mcp):
    """True when a tau appears above the first non-tau ancestor."""
    passed_non_tau = False
    for anc in tauReco._walk_ancestors(mcp):
        abs_pdg = abs(int(anc.getPDG()))
        if not passed_non_tau:
            if abs_pdg == 15:
                continue
            passed_non_tau = True
            continue
        if abs_pdg == 15:
            return True
    return False


def _new_stats():
    return {
        "n_events": 0,
        "n_taus": 0,
        "extra_pdgs": Counter(),
        "id_vs_extra": defaultdict(Counter),
        "secondary_origin": Counter(),
        "n_secondary": 0,
        "n_mother_unresolved": 0,
        "radiative_tau_status": Counter(),
        "photon_origins": Counter(),
        "tau_photon_origins": Counter(),
        "examples": defaultdict(list),
    }


def _merge(into, other):
    into["n_events"] += other["n_events"]
    into["n_taus"] += other["n_taus"]
    into["n_secondary"] += other["n_secondary"]
    into["n_mother_unresolved"] += other["n_mother_unresolved"]
    for key in ("extra_pdgs", "secondary_origin", "radiative_tau_status",
                "photon_origins", "tau_photon_origins"):
        into[key].update(other[key])
    for tau_id, row in other["id_vs_extra"].items():
        into["id_vs_extra"][tau_id].update(row)
    for tag, rows in other["examples"].items():
        room = MAX_EXAMPLES - len(into["examples"][tag])
        if room > 0:
            into["examples"][tag].extend(rows[:room])


def scan_file(task):
    """Scan one ROOT file and return the partial statistics."""
    filename, max_events = task
    stats = _new_stats()
    short = os.path.basename(filename)

    try:
        reader = root_io.Reader([filename])
        events = reader.get("events")
    except Exception as exc:                                    # pragma: no cover
        print(f"WARNING: cannot read {filename}: {exc}")
        return stats

    for event_id, event in enumerate(events):
        if max_events and event_id >= max_events:
            break
        stats["n_events"] += 1
        mc_particles = event.get("MCParticles")

        gen_taus = tauReco.findAllGenTaus(mc_particles)
        for key in gen_taus:
            tau = gen_taus[key]
            stats["n_taus"] += 1
            has_extra = tau.getHasExtraNeutrals()
            stats["id_vs_extra"][tau.getID()]["with" if has_extra else "without"] += 1
            if has_extra:
                pdgs = [int(n.getPDG()) for n in tau.getExtraNeutrals().values()]
                stats["extra_pdgs"].update(pdgs)
                if len(stats["examples"]["extra_neutral"]) < MAX_EXAMPLES:
                    stats["examples"]["extra_neutral"].append(
                        (short, event_id, f"ID={tau.getID()} PDGs={pdgs}")
                    )
            if tau.getIsSecondary():
                stats["n_secondary"] += 1
                stats["secondary_origin"][tau.getOriginPDG()] += 1
                if tau.getMotherTauKey() < 0:
                    stats["n_mother_unresolved"] += 1
                if len(stats["examples"]["secondary"]) < MAX_EXAMPLES:
                    stats["examples"]["secondary"].append(
                        (short, event_id,
                         f"originPDG={tau.getOriginPDG()} motherKey={tau.getMotherTauKey()}")
                    )
            for origin in tau.getConstOrigin().values():
                if origin != tauReco.PHOTON_ORIGIN_NOT_A_PHOTON:
                    stats["tau_photon_origins"][origin] += 1
                    if len(stats["examples"]["const_photon"]) < MAX_EXAMPLES:
                        stats["examples"]["const_photon"].append(
                            (short, event_id,
                             f"ID={tau.getID()} origin={ORIGIN_NAMES.get(origin, origin)}")
                        )

        for part in mc_particles:
            pdg = abs(int(part.getPDG()))
            if pdg == 15:
                if _has_tau_ancestor_beyond_copies(part):
                    status = int(part.getGeneratorStatus())
                    stats["radiative_tau_status"][status] += 1
                    if len(stats["examples"]["radiative_tau"]) < MAX_EXAMPLES:
                        stats["examples"]["radiative_tau"].append(
                            (short, event_id, f"status={status}")
                        )
            elif pdg == 22 and int(part.getGeneratorStatus()) == 1:
                stats["photon_origins"][tauReco.classify_photon_origin(part)] += 1

    return stats


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-f", "--sample", default=None)
    parser.add_argument("--input-list", nargs="+", default=None)
    parser.add_argument("--samples-config", default="config/samples/samples.yaml")
    parser.add_argument("--n-files", type=int, default=2)
    parser.add_argument("--max-events", type=int, default=0,
                        help="Events per file (0 = all)")
    parser.add_argument("--n-workers", type=int, default=1)
    args = parser.parse_args()

    logging.basicConfig(level=logging.ERROR)
    loggers = {k: logging.getLogger(k) for k in ("config", "io", "processing", "pi0mass")}

    filenames, _ = myutils.get_root_trees_path(args.sample, None, loggers, False, args)
    filenames = filenames[: args.n_files]
    if not filenames:
        print("No input files resolved.")
        sys.exit(1)

    print(f"Scanning {len(filenames)} file(s) with {args.n_workers} worker(s)"
          f"{'' if not args.max_events else f', {args.max_events} events/file'}.\n")

    tasks = [(f, args.max_events) for f in filenames]
    total = _new_stats()
    if args.n_workers > 1:
        with Pool(processes=args.n_workers) as pool:
            for partial in pool.imap_unordered(scan_file, tasks):
                _merge(total, partial)
    else:
        for task in tasks:
            _merge(total, scan_file(task))

    print(f"Events scanned: {total['n_events']} | gen taus (status 2): {total['n_taus']}\n")

    print("── Extra neutrals (not counted in the decay-mode ID) ──")
    if not total["extra_pdgs"]:
        print("  none found")
    for pdg, count in total["extra_pdgs"].most_common():
        print(f"  {_pdg_name(pdg):>12} (PDG {pdg:>6}): {count}")

    print("\n── Decay-mode ID vs extra-neutral flag ──")
    for tau_id in sorted(total["id_vs_extra"]):
        row = total["id_vs_extra"][tau_id]
        n_total = row["with"] + row["without"]
        print(f"  ID {tau_id:>4}: {n_total:>8} taus | with extra neutral: {row['with']}")

    print("\n── Secondary taus (tau -> gamma -> tau tau) ──")
    print(f"  flagged secondary: {total['n_secondary']} / {total['n_taus']}")
    for pdg, count in total["secondary_origin"].most_common():
        print(f"    originPDG {pdg} ({_pdg_name(pdg)}): {count}")
    if total["n_mother_unresolved"]:
        print(f"  WARNING: {total['n_mother_unresolved']} secondary tau(s) "
              f"with unresolved mother key")

    print("\n── generatorStatus of taus having a tau ancestor ──")
    print("  (findAllGenTaus only keeps status==2; anything else is invisible to the tree)")
    if not total["radiative_tau_status"]:
        print("  none found in this subsample")
    for status, count in sorted(total["radiative_tau_status"].items()):
        flag = "  <-- kept" if status == 2 else "  <-- NOT kept"
        print(f"    status {status}: {count}{flag}")

    print("\n── Origin of all generator-level photons (status 1) ──")
    for origin, count in total["photon_origins"].most_common():
        print(f"  {ORIGIN_NAMES.get(origin, origin):>12}: {count}")

    print("\n── Origin of photons inside gen tau constituents ──")
    if not total["tau_photon_origins"]:
        print("  none (in this sample the FSR photons hang from the tau copy, "
              "not from the status-2 tau)")
    for origin, count in total["tau_photon_origins"].most_common():
        print(f"  {ORIGIN_NAMES.get(origin, origin):>12}: {count}")

    if total["examples"]:
        print("\n── Examples (file, event, detail) ──")
        for tag, rows in total["examples"].items():
            print(f"  {tag}:")
            for short, event_id, detail in rows:
                print(f"    {short} #{event_id}  {detail}")


if __name__ == "__main__":
    main()
