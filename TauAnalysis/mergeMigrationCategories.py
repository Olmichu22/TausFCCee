#!/usr/bin/env python
"""Build merged migration-efficiency graphs from a tau reconstruction ROOT file.

The per-decay migration plots split the "the gen tau was not reconstructed in
any usable category" outcome into three separate buckets:

* ``Unmatched`` -- no reco tau was matched to the gen tau at all;
* ``NoTau``     -- a reco object was matched but was not identified as a tau
                   (``recoTauId == -1``);
* ``Other``     -- matched and identified, but with an id outside every bucket
                   (in practice ``-20``: one prong plus neutral hadrons).

For plotting purposes the three are usually more useful merged into a single
"lost" category. Efficiencies are stored as ``TGraphAsymmErrors``, which cannot
be added, so the merge is done on the underlying ``TH1`` numerators and the
ratio is recomputed with the same Bayesian interval used upstream.

The output file only holds the merged graphs; it is meant to be declared as an
extra dataset in a CompareAlgs YAML config alongside the original file.
"""

import argparse
import os

import ROOT

# Etiqueta de cada modo gen tal y como aparece en los nombres de histograma.
MODE_TAGS = ["0", "1", "2", "10", "Elec", "Muon"]
# Categorias reco que se funden en una sola: "no reconstruido de forma util".
LOST_CATEGORIES = ["ToUnmatched", "ToNoTau", "ToOther"]
# Mismo intervalo de confianza que usa plotTausLongResults al construir las eficiencias.
DIVIDE_OPTS = "cl=0.683 b(1,1) mode"


def merged_lost_graph(rootfile, variable, tag):
    """Return the merged "lost" efficiency graph for one variable and gen mode.

    Args:
        rootfile (ROOT.TFile): Open tau reconstruction results file.
        variable (str): Either ``"TauP"`` or ``"TauTheta"``.
        tag (str): Gen decay mode tag, e.g. ``"0"``, ``"10"``, ``"Elec"``.

    Returns:
        ROOT.TGraphAsymmErrors: Efficiency named ``hEffiGen<variable><tag>ToLost``,
        or None if the denominator or every numerator is missing.
    """
    denominator = rootfile.Get(f"{variable}{tag}Gen")
    if not denominator:
        print(f"[WARN] missing denominator {variable}{tag}Gen")
        return None

    numerator = None
    for category in LOST_CATEGORIES:
        hist = rootfile.Get(f"{variable}{tag}{category}Matched")
        if not hist:
            print(f"[WARN] missing {variable}{tag}{category}Matched")
            continue
        if numerator is None:
            numerator = hist.Clone(f"{variable}{tag}ToLostMatched")
            numerator.SetDirectory(0)
        else:
            numerator.Add(hist)

    if numerator is None:
        return None

    graph = ROOT.TGraphAsymmErrors()
    graph.Divide(numerator, denominator, DIVIDE_OPTS)
    graph.SetName(f"hEffiGen{variable}{tag}ToLost")
    graph.SetTitle(f"hEffiGen{variable}{tag}ToLost")
    return graph


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-i", "--input", required=True,
                        help="Tau reconstruction results ROOT file")
    parser.add_argument("-o", "--output", default=None,
                        help="Output ROOT file (default: <input dir>/migration_merged.root)")
    args = parser.parse_args()

    output = args.output or os.path.join(os.path.dirname(args.input),
                                         "migration_merged.root")
    
    fin = ROOT.TFile(args.input)
    if not fin or fin.IsZombie():
        raise RuntimeError(f"Could not open {args.input}")

    graphs = []
    for variable in ("TauP", "TauTheta"):
        for tag in MODE_TAGS:
            graph = merged_lost_graph(fin, variable, tag)
            if graph is not None:
                graphs.append(graph)

    fout = ROOT.TFile(output, "RECREATE")
    for graph in graphs:
        graph.Write()
    fout.Close()
    fin.Close()
    print(f"[OK] {len(graphs)} merged graphs written to {output}")


if __name__ == "__main__":
    main()
