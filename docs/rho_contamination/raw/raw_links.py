"""Que le pasa a los dos fotones del pi0 en gen tau->pi pi0 nu (GenTauType 1).

Lee los EDM4hep originales con podio, reproduce la reconstruccion del arbol
(findAllGenTaus / findAllTaus dR<0.4 sin cortes en P, MatchRecoGenTau dR<1) y,
para cada gen tau tipo 1 emparejado con un reco tau de tipo 0/1/2/3+, sigue
cada foton gen del pi0 en sentido MC -> reco con MCTruthRecoLink (peso =
fraccion de la energia/hits del MC que acaba en ese PFO). Si el foton no tiene
links directos se miran sus descendientes de simulacion (conversiones).

Categorias por foton (campo ``cat``):
  a_lost      sin ningun link (ni directo ni via descendientes de sim) o solo links debiles (w_eff<0.1)
  b_in_cone   PFO foton dentro del cono del reco tau (constituyente)
  c_out_cone  PFO foton fuera del cono (dR>0.4 del lider) o en otro tau
  d_charged   PFO cargado (211/11/13): el pion lider, otro cargado, o electron de conversion
  e_neutralh  PFO hadron neutro (2112)
  f_merged    los dos fotones del pi0 -> el MISMO PFO foton
  g_split     un foton -> dos o mas PFO foton con w_eff>=0.1

Uso: python raw_links.py [nfiles] [nevents_por_fichero]   (por defecto 30 ficheros, todos los eventos)
Salida: records.json en este directorio.
"""
import sys, os, math, json
sys.path.insert(0, "/nfs/cms/arqolmo/TausFCCee")
import ROOT
from podio import root_io
from modules import tauReco, electronReco, muonReco, myutils
from multiprocessing import Pool

HERE = os.path.dirname(os.path.abspath(__file__))
SAMPLE = "/pnfs/ciemat.es/data/cms/store/user/cepeda/FCC/FullSim/ZTauTau_PolSM_March24_2M_4"

# Misma config que el arbol (New2MSample_results0.4_tph0.0_tpi0.0_n0.0_g0.0)
dRMax, minPTauPhoton, minPTauPion, PNeutron, generalPCut, dRMatch = 0.4, 0., 0., 0., 0., 1.
BFIELD = 2.0          # T (CLD)
W_SIG = 0.1           # peso efectivo minimo para considerar un link "significativo"
MAX_DESC = 3000       # tope de descendientes de sim a explorar por foton


def mcidx(p):
    return int(p.getObjectID().index)


def p4_of(obj):
    v = ROOT.TLorentzVector()
    m = obj.getMomentum()
    v.SetXYZM(m.x, m.y, m.z, obj.getMass())
    return v


def decode_w(w):
    """peso LCIO: track*1000 + cluster*1000*10000 -> (track, cluster)."""
    enc = int(w)
    return (enc % 10000) / 1000.0, (enc // 10000) / 1000.0


def ancestry(mc, maxd=300):
    out, seen, cur = [], set(), mc
    for _ in range(maxd):
        par = list(cur.getParents())
        if not par:
            break
        cur = par[0]
        i = mcidx(cur)
        if i in seen:
            break
        seen.add(i)
        out.append(cur)
    return out


def is_sim(mc):
    return mc.getGeneratorStatus() == 0 or bool(mc.isCreatedInSimulation())


def gen_ancestor(mc):
    if not is_sim(mc):
        return mc
    for a in ancestry(mc):
        if not is_sim(a):
            return a
    return None


def tau_of(mc, tau_key_by_idx):
    if abs(mc.getPDG()) == 15 and mcidx(mc) in tau_key_by_idx:
        return tau_key_by_idx[mcidx(mc)]
    for a in ancestry(mc):
        if abs(a.getPDG()) == 15 and mcidx(a) in tau_key_by_idx:
            return tau_key_by_idx[mcidx(a)]
    return -1


def is_isr_photon(mc):
    if abs(mc.getPDG()) != 22 or mc.getGeneratorStatus() != 1:
        return False
    for a in ancestry(mc):
        pdg = abs(a.getPDG())
        if pdg == 22:
            continue
        return pdg == 11 and a.getGeneratorStatus() != 1
    return True


def sim_descendants(mc):
    """Descendientes de simulacion (status 0) de mc, BFS con tope."""
    out, seen, stack = [], {mcidx(mc)}, list(mc.getDaughters())
    while stack and len(out) < MAX_DESC:
        d = stack.pop()
        i = mcidx(d)
        if i in seen:
            continue
        seen.add(i)
        if not is_sim(d):
            continue
        out.append(d)
        stack.extend(list(d.getDaughters()))
    return out


def track_p(trk):
    ts = trk.getTrackStates()[0]
    if ts.omega == 0:
        return -1.
    pt = 0.299792458 * BFIELD * 1e-3 / abs(ts.omega)
    return pt * math.sqrt(1. + ts.tanLambda ** 2)


def pfo_info(pf, lead_p4):
    v = p4_of(pf)
    cls = list(pf.getClusters())
    trks = list(pf.getTracks())
    return dict(idx=mcidx(pf), pdg=int(pf.getPDG()), P=v.P(), E=v.E(), theta=v.Theta(),
                dR_lead=myutils.dRAngle(v, lead_p4),
                Ecl=sum(float(c.getEnergy()) for c in cls), ncl=len(cls), ntrk=len(trks),
                Ptrk=(track_p(trks[0]) if trks else -1.))


def process_file(args):
    fname, nmax = args
    fidx = int(fname.split("_")[-1].split(".")[0])
    reader = root_io.Reader([fname])
    records = []
    C = dict(events=0, gen1=0, gen1_extra_neutral=0, gen1_pi0_not_2gamma=0, gen1_unmatched=0,
             gen1_reco_other=0, gen1_reco0=0, gen1_reco1=0, gen1_reco2=0, gen1_reco3p=0)
    for iev, event in enumerate(reader.get("events")):
        if nmax and iev >= nmax:
            break
        C["events"] += 1
        mc = event.get("MCParticles")
        pfos = event.get("PandoraPFOs")
        genTaus = tauReco.findAllGenTaus(mc)
        recoTau_raw = tauReco.findAllTaus(pfos, dRMax, minPTauPhoton, minPTauPion, PNeutron, generalPCut)
        recoEle = electronReco.findAllElectrons(pfos, generalPCut)
        recoMu = muonReco.findAllMuons(pfos, generalPCut)
        recoTaus, k = {}, 0
        for coll in (recoTau_raw, recoEle, recoMu):
            for t in range(len(coll)):
                recoTaus[k] = coll[t]; k += 1
        # PFO idx -> lista de claves de reco tau en los que es constituyente
        pfo_in_tau = {}
        for key, rt in recoTaus.items():
            for _, c in rt.getDaughters().items():
                pfo_in_tau.setdefault(mcidx(c), []).append(key)

        # tau mcidx (status 2 y copias) -> clave
        tau_key_by_idx = {}
        for key, gt in genTaus.items():
            tau_key_by_idx[gt.getMCIdx()] = key
        for p in mc:
            if abs(p.getPDG()) == 15 and p.getGeneratorStatus() == 2 and mcidx(p) in tau_key_by_idx:
                key = tau_key_by_idx[mcidx(p)]
                for a in ancestry(p):
                    if abs(a.getPDG()) != 15:
                        break
                    tau_key_by_idx[mcidx(a)] = key
        for p in mc:
            if abs(p.getPDG()) == 15 and mcidx(p) not in tau_key_by_idx:
                for a in ancestry(p):
                    if abs(a.getPDG()) != 15:
                        break
                    if mcidx(a) in tau_key_by_idx:
                        tau_key_by_idx[mcidx(p)] = tau_key_by_idx[mcidx(a)]
                        break

        # MC -> reco (MCTruthRecoLink: from=reco, to=MC, peso = fraccion del MC en el PFO)
        links_by_mc = {}
        for l in event.get("MCTruthRecoLink"):
            links_by_mc.setdefault(mcidx(l.getTo()), []).append((l.getFrom(), float(l.getWeight())))
        # reco -> MC (RecoMCTruthLink: peso = fraccion del PFO que viene del MC)
        links_by_reco = {}
        for l in event.get("RecoMCTruthLink"):
            links_by_reco.setdefault(mcidx(l.getFrom()), []).append((l.getTo(), float(l.getWeight())))

        nTausType = 0
        for i in range(len(genTaus)):
            gt = genTaus[i]
            if gt.getID() != 1:
                continue
            C["gen1"] += 1
            if gt.getNExtraNeutrals() > 0:
                C["gen1_extra_neutral"] += 1
                continue
            pion_mc, pi0_mc = None, None
            for _, d in gt.getDaughters().items():
                if abs(d.getPDG()) in (211, 321):
                    pion_mc = d
                elif d.getPDG() == 111:
                    pi0_mc = d
            gammas = [d for d in pi0_mc.getDaughters() if d.getGeneratorStatus() == 1 and d.getPDG() == 22] if pi0_mc is not None else []
            if pion_mc is None or len(gammas) != 2:
                C["gen1_pi0_not_2gamma"] += 1
                continue
            findMatch, nTausType = tauReco.MatchRecoGenTau(gt, recoTaus, nTausType, maxDRMatch=dRMatch, selectDecay=-777)
            if findMatch < 0:
                C["gen1_unmatched"] += 1
                continue
            rt = recoTaus[findMatch]
            rid = int(rt.getID())
            if rid == 0: C["gen1_reco0"] += 1
            elif rid == 1: C["gen1_reco1"] += 1
            elif rid == 2: C["gen1_reco2"] += 1
            elif 3 <= rid <= 9: C["gen1_reco3p"] += 1
            else:
                C["gen1_reco_other"] += 1
                C["reco_other_id_%d" % rid] = C.get("reco_other_id_%d" % rid, 0) + 1
                continue
            consts = rt.getDaughters()
            lead = consts[0]
            if abs(lead.getPDG()) != 211:
                C["gen1_reco_other"] += 1
                continue
            lead_p4, pion_p4, pi0_p4 = p4_of(lead), p4_of(pion_mc), p4_of(pi0_mc)
            gam_p4 = [p4_of(g) for g in gammas]
            pion_idx = mcidx(pion_mc)
            gam_idx = {mcidx(g): j for j, g in enumerate(gammas)}
            const_idx = {mcidx(c): kk for kk, c in consts.items()}
            lead_cls = sum(float(c.getEnergy()) for c in lead.getClusters())
            lead_trks = list(lead.getTracks())
            # RecoMCTruthLink del lider: fraccion del cluster del lider que viene de cada MC
            lead_from = {}
            for (m, w) in links_by_reco.get(mcidx(lead), []):
                tw, cw = decode_w(w)
                lead_from[mcidx(m)] = (tw, cw)

            # ---- links de cada foton (directos + via descendientes de sim) ----
            phot = []
            for j, g in enumerate(gammas):
                Eg = gam_p4[j].E()
                links = []
                for (pf, w) in links_by_mc.get(mcidx(g), []):
                    tw, cw = decode_w(w)
                    links.append(dict(pfo=pf, tw=tw, cw=cw, w_eff=max(tw, cw), via="direct", via_pdg=22, via_E=Eg))
                desc = sim_descendants(g)
                n_desc_linked = 0
                for d in desc:
                    for (pf, w) in links_by_mc.get(mcidx(d), []):
                        tw, cw = decode_w(w)
                        Ed = float(d.getEnergy())
                        n_desc_linked += 1
                        links.append(dict(pfo=pf, tw=tw, cw=cw, w_eff=max(tw, cw) * Ed / Eg if Eg > 0 else 0.,
                                          via="sim", via_pdg=int(d.getPDG()), via_E=Ed))
                # agrupa por PFO
                by_pfo = {}
                for L in links:
                    e = by_pfo.setdefault(mcidx(L["pfo"]), dict(pfo=L["pfo"], w_eff=0., tw=0., cw=0., via=set(), via_pdg=set()))
                    e["w_eff"] += L["w_eff"]; e["tw"] = max(e["tw"], L["tw"]); e["cw"] = max(e["cw"], L["cw"])
                    e["via"].add(L["via"]); e["via_pdg"].add(L["via_pdg"])
                lst = sorted(by_pfo.values(), key=lambda e: -e["w_eff"])
                ep = g.getEndpoint()
                phot.append(dict(j=j, mc=g, E=Eg, P=gam_p4[j].P(), theta=gam_p4[j].Theta(),
                                 endpoint_r=math.hypot(ep.x, ep.y), endpoint_z=float(ep.z),
                                 decayed_in_tracker=bool(g.isDecayedInTracker()),
                                 n_desc=len(desc), n_desc_linked=n_desc_linked,
                                 ang_pion=myutils.dRAngle(gam_p4[j], pion_p4),
                                 ang_other=myutils.dRAngle(gam_p4[j], gam_p4[1 - j]),
                                 ang_lead=myutils.dRAngle(gam_p4[j], lead_p4),
                                 links=lst))

            # ---- categoria ----
            main = [p["links"][0] if p["links"] and p["links"][0]["w_eff"] >= W_SIG else None for p in phot]
            merged = (main[0] is not None and main[1] is not None and
                      mcidx(main[0]["pfo"]) == mcidx(main[1]["pfo"]) and main[0]["pfo"].getPDG() == 22)
            out_phot = []
            for j, p in enumerate(phot):
                m = main[j]
                d = dict(E=p["E"], P=p["P"], theta=p["theta"], endpoint_r=p["endpoint_r"], endpoint_z=p["endpoint_z"],
                         decayed_in_tracker=p["decayed_in_tracker"], n_desc=p["n_desc"], n_desc_linked=p["n_desc_linked"],
                         ang_pion=p["ang_pion"], ang_other=p["ang_other"], ang_lead=p["ang_lead"],
                         nlinks=len(p["links"]), n_sig=sum(1 for L in p["links"] if L["w_eff"] >= W_SIG),
                         sum_w=sum(L["w_eff"] for L in p["links"]),
                         links=[dict(pfo_idx=mcidx(L["pfo"]), pdg=int(L["pfo"].getPDG()), P=p4_of(L["pfo"]).P(),
                                     w_eff=L["w_eff"], tw=L["tw"], cw=L["cw"], via=sorted(L["via"]),
                                     via_pdg=sorted(L["via_pdg"]), in_cone=mcidx(L["pfo"]) in const_idx,
                                     is_lead=mcidx(L["pfo"]) == mcidx(lead))
                                for L in p["links"]])
                other_main = mcidx(main[1 - j]["pfo"]) if main[1 - j] is not None else -1
                # PFO foton con peso significativo, sin contar el PFO principal del otro foton
                n_sig_phot = sum(1 for L in p["links"] if L["w_eff"] >= W_SIG and L["pfo"].getPDG() == 22
                                 and mcidx(L["pfo"]) != other_main)
                if m is None:
                    cat = "a_lost"
                    sub = "no_link" if not p["links"] else "weak_links"
                    if p["decayed_in_tracker"] or (p["n_desc"] > 0 and p["endpoint_r"] < 2100):
                        sub = "converted_" + sub
                    pinfo = None
                else:
                    pf = m["pfo"]
                    pinfo = pfo_info(pf, lead_p4)
                    pinfo.update(w_eff=m["w_eff"], tw=m["tw"], cw=m["cw"], via=sorted(m["via"]), via_pdg=sorted(m["via_pdg"]),
                                 in_cone=mcidx(pf) in const_idx, is_lead=mcidx(pf) == mcidx(lead),
                                 in_other_tau=[kk for kk in pfo_in_tau.get(mcidx(pf), []) if kk != findMatch],
                                 purity_cw=lead_from.get(mcidx(p["mc"]), (0., 0.))[1] if mcidx(pf) == mcidx(lead) else None)
                    # fraccion del PFO que viene de ESTE foton (RecoMCTruthLink)
                    pur = 0.
                    for (mm, w) in links_by_reco.get(mcidx(pf), []):
                        if mcidx(mm) == mcidx(p["mc"]):
                            pur = max(pur, decode_w(w)[1])
                    pinfo["frac_pfo_from_gamma"] = pur
                    pdg = abs(pf.getPDG())
                    conv = "sim" in m["via"] and "direct" not in m["via"]
                    if merged:
                        cat, sub = "f_merged", ("in_cone" if pinfo["in_cone"] else "out_cone")
                    elif n_sig_phot >= 2:
                        cat, sub = "g_split", f"{n_sig_phot}_pfo_photons"
                    elif pdg == 22:
                        if pinfo["in_cone"]:
                            cat, sub = "b_in_cone", ("via_conversion" if conv else "direct")
                        else:
                            cat = "c_out_cone"
                            sub = "other_tau" if pinfo["in_other_tau"] else ("dR>0.4" if pinfo["dR_lead"] > dRMax else "in_cone_not_const")
                            if conv: sub += "_conv"
                    elif pdg in (211, 11, 13, 321, 2212):
                        cat = "d_charged"
                        sub = ("lead_pion" if pinfo["is_lead"] else f"pdg{pdg}") + ("_conv" if conv else "")
                    elif pdg == 2112:
                        cat, sub = "e_neutralh", ("in_cone" if pinfo["in_cone"] else "out_cone") + ("_conv" if conv else "")
                    else:
                        cat, sub = "other", f"pdg{pdg}"
                    if merged and n_sig_phot >= 2:
                        sub += "+split"
                d.update(cat=cat, sub=sub, main=pinfo)
                out_phot.append(d)

            # ---- E/p del lider, con y sin el foton absorbido ----
            lead_gam_cw = sum(lead_from.get(mcidx(g), (0., 0.))[1] for g in gammas)
            lead_gam_cw_sim = 0.
            # tambien descendientes de sim de los fotones que caen en el lider
            for p in phot:
                for L in p["links"]:
                    if mcidx(L["pfo"]) == mcidx(lead) and "sim" in L["via"]:
                        pass
            for (m_, w) in links_by_reco.get(mcidx(lead), []):
                if is_sim(m_):
                    ga = gen_ancestor(m_)
                    if ga is not None and mcidx(ga) in gam_idx:
                        lead_gam_cw_sim += decode_w(w)[1]
            lead_info = dict(P=lead_p4.P(), E=lead_p4.E(), theta=lead_p4.Theta(), Ecl=lead_cls,
                             Ptrk=(track_p(lead_trks[0]) if lead_trks else -1.), ntrk=len(lead_trks),
                             ncl=int(lead.clusters_size()),
                             cw_from_gammas=lead_gam_cw, cw_from_gamma_sim=lead_gam_cw_sim,
                             cw_from_pion=lead_from.get(pion_idx, (0., 0.))[1],
                             tw_from_pion=lead_from.get(pion_idx, (0., 0.))[0],
                             EoverP=(lead_cls / lead_p4.P() if lead_p4.P() > 0 else -1.),
                             EoverP_nogamma=(lead_cls * (1. - lead_gam_cw - lead_gam_cw_sim) / lead_p4.P() if lead_p4.P() > 0 else -1.))

            # ---- fotones PFO del reco tau: clasificacion (para 3+) ----
            main_pfo_idx = {mcidx(m["pfo"]): j for j, m in enumerate(main) if m is not None}
            reco_photons = []
            for kk, c in consts.items():
                if kk == 0 or c.getPDG() != 22:
                    continue
                info = pfo_info(c, lead_p4)
                lk = sorted(links_by_reco.get(mcidx(c), []), key=lambda x: -decode_w(x[1])[1])
                if not lk:
                    cls, sub = "no_link", ""
                else:
                    m_, w = lk[0]
                    tw, cw = decode_w(w)
                    ga = gen_ancestor(m_)
                    sim = is_sim(m_)
                    if ga is None:
                        cls, sub = "sim_no_gen", ""
                    elif mcidx(ga) in gam_idx:
                        if not sim:
                            cls = "pi0_photon_main" if main_pfo_idx.get(mcidx(c)) == gam_idx[mcidx(ga)] else "pi0_photon_split"
                        else:
                            cls = "pi0_photon_conv_main" if main_pfo_idx.get(mcidx(c)) == gam_idx[mcidx(ga)] else "pi0_photon_conv_frag"
                        sub = f"gamma{gam_idx[mcidx(ga)]}"
                    elif mcidx(ga) == pion_idx:
                        cls, sub = "shower_pion", ("sim" if sim else "gen_direct")
                    elif abs(ga.getPDG()) == 22:
                        if is_isr_photon(ga):
                            cls, sub = "ISR", ""
                        else:
                            origin, ppdg, _ = tauReco.photon_origin_info(ga)
                            tk = tau_of(ga, tau_key_by_idx)
                            same = tk == i
                            if origin == tauReco.PHOTON_ORIGIN_TAU_FSR:
                                cls = "FSR_tau" if same else "FSR_other_tau"
                            elif origin == tauReco.PHOTON_ORIGIN_CHARGED_RAD:
                                cls = "rad_pion_gen" if same else "rad_other_tau"
                            elif origin == tauReco.PHOTON_ORIGIN_PI0:
                                cls = "pi0_other_tau" if not same else "pi0_same_tau_other"
                            else:
                                cls = "photon_gen_other"
                            sub = f"parent{ppdg}"
                    else:
                        tk = tau_of(ga, tau_key_by_idx)
                        cls = "other_same_tau" if tk == i else ("other_tau" if tk >= 0 else "other")
                        sub = f"pdg{ga.getPDG()}"
                    info.update(cw=cw, tw=tw, mc_pdg=int(m_.getPDG()), mc_E=float(m_.getEnergy()), mc_sim=sim)
                info.update(cls=cls, sub=sub, is_main=mcidx(c) in main_pfo_idx,
                            n_links=len(lk))
                reco_photons.append(info)

            records.append(dict(evt=fidx * 1000 + iev, gen_key=i, reco_id=rid,
                                gen_visP=gt.getvisMomentum().P(), gen_visTheta=gt.getvisMomentum().Theta(),
                                gen_pion_P=pion_p4.P(), gen_pi0_P=pi0_p4.P(),
                                gen_pion_theta=pion_p4.Theta(),
                                ang_gg=myutils.dRAngle(gam_p4[0], gam_p4[1]),
                                ang_pi_pi0=myutils.dRAngle(pion_p4, pi0_p4),
                                reco_P=rt.getMomentum().P(), reco_nconst=len(consts),
                                reco_nphot=sum(1 for kk, c in consts.items() if kk != 0 and c.getPDG() == 22),
                                reco_nneut=sum(1 for kk, c in consts.items() if kk != 0 and c.getPDG() == 2112),
                                lead=lead_info, photons=out_phot, merged=merged, reco_photons=reco_photons))
    return records, C


if __name__ == "__main__":
    nfiles = int(sys.argv[1]) if len(sys.argv) > 1 else 30
    nev = int(sys.argv[2]) if len(sys.argv) > 2 else 0
    files = [os.path.join(SAMPLE, f"out_reco_edm4hep_edm4hep_{i}.root") for i in range(1000, 1000 + nfiles)]
    files = [f for f in files if os.path.exists(f)]
    with Pool(min(8, len(files))) as pool:
        res = pool.map(process_file, [(f, nev) for f in files])
    records, counters = [], {}
    for r, c in res:
        records += r
        for k, v in c.items():
            counters[k] = counters.get(k, 0) + v
    out = os.path.join(HERE, "records.json")
    json.dump(dict(records=records, counters=counters, nfiles=len(files), W_SIG=W_SIG),
              open(out, "w"))
    print(counters, len(records), "->", out)
