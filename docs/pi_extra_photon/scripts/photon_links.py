"""Que es el foton extra en gen tau->pi nu reconstruidos como pi + 1 gamma.

Lee los EDM4hep originales con podio, reproduce la reconstruccion de
TTreesTausLong.py (findAllGenTaus / findAllTaus dR<0.4 sin cortes en P,
MatchRecoGenTau dR<1) y, para cada gen tau tipo 0 cuyo reco tau es tipo 1,
vuelca TODOS los RecoMCTruthLink del PFO foton, la geometria del cluster y la
ascendencia de las MCParticles enlazadas. Salida: records.json.
"""
import sys, os, math, json
sys.path.insert(0, "/nfs/cms/arqolmo/TausFCCee")
import ROOT
from podio import root_io
from modules import tauReco, electronReco, muonReco, myutils
from multiprocessing import Pool

HERE = os.path.dirname(os.path.abspath(__file__))

# Misma config que New2MSample_results0.4_tph0.0_tpi0.0_n0.0_g0.0/config.yaml
dRMax, minPTauPhoton, minPTauPion, PNeutron, generalPCut, dRMatch = 0.4, 0., 0., 0., 0., 1.

NEUTRINOS = {12, 14, 16}


def mcidx(p):
    return int(p.getObjectID().index)


def p4_of(obj):
    v = ROOT.TLorentzVector()
    m = obj.getMomentum()
    v.SetXYZM(m.x, m.y, m.z, obj.getMass())
    return v


def decode_w(w):
    enc = int(w)
    return (enc % 10000) / 1000.0, (enc // 10000) / 1000.0


def ancestry(mc, maxd=300):
    """Lista de ancestros siguiendo el primer padre (hasta la raiz)."""
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
    """Primer ancestro (o el propio) que NO es de simulacion."""
    if not is_sim(mc):
        return mc
    for a in ancestry(mc):
        if not is_sim(a):
            return a
    return None


def tau_of(mc, tau_key_by_idx):
    """Clave del gen tau del que desciende mc (-1 si ninguno)."""
    if abs(mc.getPDG()) == 15 and mcidx(mc) in tau_key_by_idx:
        return tau_key_by_idx[mcidx(mc)]
    for a in ancestry(mc):
        if abs(a.getPDG()) == 15 and mcidx(a) in tau_key_by_idx:
            return tau_key_by_idx[mcidx(a)]
    return -1


def is_isr_photon(mc):
    """Foton status 1 cuyo primer ancestro no-foton es un e+/e- (haz) o no tiene padre."""
    if abs(mc.getPDG()) != 22 or mc.getGeneratorStatus() != 1:
        return False
    for a in ancestry(mc):
        pdg = abs(a.getPDG())
        if pdg == 22:
            continue
        return pdg == 11 and a.getGeneratorStatus() != 1
    return True  # sin padres


def classify_link(mc, this_tau, pion_idx, tau_key_by_idx):
    """Devuelve (categoria, subcategoria) de la MCParticle enlazada."""
    pdg = abs(mc.getPDG())
    st = mc.getGeneratorStatus()
    sim = is_sim(mc)
    if not sim:
        # particula del generador
        if mcidx(mc) == pion_idx:
            return "pion_gen_directo", "link al propio pion cargado"
        if pdg == 22:
            if is_isr_photon(mc):
                return "ISR", "foton ISR gen"
            origin, ppdg, aidx = tauReco.photon_origin_info(mc)
            tk = tau_of(mc, tau_key_by_idx)
            same = "mismo tau" if tk == this_tau else ("otro tau" if tk >= 0 else "sin tau")
            if origin == tauReco.PHOTON_ORIGIN_TAU_FSR:
                return ("FSR_tau" if tk == this_tau else "FSR_otro_tau"), same
            if origin == tauReco.PHOTON_ORIGIN_PI0:
                return ("pi0_mismo_tau" if tk == this_tau else "pi0_otro_tau"), same
            if origin == tauReco.PHOTON_ORIGIN_CHARGED_RAD:
                return ("rad_cargado_gen" if tk == this_tau else "rad_cargado_otro_tau"), f"padre {ppdg} {same}"
            return "foton_gen_otro", f"padre {ppdg} {same}"
        tk = tau_of(mc, tau_key_by_idx)
        return "otra_particula_gen", f"pdg {mc.getPDG()} st {st} " + ("mismo tau" if tk == this_tau else f"tau {tk}")
    # particula de simulacion: sube hasta el generador
    ga = gen_ancestor(mc)
    if ga is None:
        return "sim_sin_ancestro_gen", ""
    gpdg = ga.getPDG()
    if mcidx(ga) == pion_idx:
        return "sim_shower_pion", f"pdg {mc.getPDG()}"
    tk = tau_of(ga, tau_key_by_idx)
    if abs(gpdg) == 22:
        cat, sub = classify_link(ga, this_tau, pion_idx, tau_key_by_idx)
        return "sim_de_" + cat, f"pdg {mc.getPDG()} <- {sub}"
    if tk == this_tau:
        return "sim_de_mismo_tau_no_pion", f"ancestro gen {gpdg}"
    if tk >= 0:
        return "sim_de_otro_tau", f"ancestro gen {gpdg}"
    return "sim_otro", f"ancestro gen {gpdg} st {ga.getGeneratorStatus()}"


def cluster_info(cl):
    pos = cl.getPosition()
    r = math.hypot(pos.x, pos.y)
    sub = list(cl.getSubdetectorEnergies())
    return dict(E=float(cl.getEnergy()), x=pos.x, y=pos.y, z=pos.z, r=r,
                theta=math.atan2(r, pos.z), phi=math.atan2(pos.y, pos.x),
                subE=[float(s) for s in sub], nhits=int(cl.hits_size()),
                type=int(cl.getType()))


def angle_between(a, b):
    """Angulo (rad) entre dos vectores posicion (dicts con x,y,z)."""
    na = math.sqrt(a["x"]**2 + a["y"]**2 + a["z"]**2)
    nb = math.sqrt(b["x"]**2 + b["y"]**2 + b["z"]**2)
    if na == 0 or nb == 0:
        return -1.
    c = (a["x"]*b["x"] + a["y"]*b["y"] + a["z"]*b["z"]) / na / nb
    return math.acos(max(-1., min(1., c)))


def process_file(fname):
    fidx = int(fname.split("_")[-1].split(".")[0])
    reader = root_io.Reader([fname])
    records, isr_records, gen0_records, counters = [], [], [], {"events": 0, "gen0": 0, "gen0_reco1": 0,
                                             "gen0_reco0": 0, "gen0_unmatched": 0,
                                             "gen0_matched": 0}
    for iev, event in enumerate(reader.get("events")):
        counters["events"] += 1
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
        # descendientes tau de copias (por si el status 2 no es la ultima copia)
        for p in mc:
            if abs(p.getPDG()) == 15 and mcidx(p) not in tau_key_by_idx:
                for a in ancestry(p):
                    if abs(a.getPDG()) != 15:
                        break
                    if mcidx(a) in tau_key_by_idx:
                        tau_key_by_idx[mcidx(p)] = tau_key_by_idx[mcidx(a)]
                        break

        # ISR del evento: fotones status 1 con ancestro e+/e- de haz
        isr = []
        for p in mc:
            if is_isr_photon(p):
                v = p4_of(p)
                drs = [myutils.dRAngle(v, genTaus[key].getvisMomentum()) for key in genTaus]
                isr.append(dict(E=v.E(), theta=v.Theta(), mindR_tau=min(drs) if drs else -1.,
                                pdg_parent=[int(a.getPDG()) for a in ancestry(p)[:2]]))
        isr_records.append(dict(evt=fidx*1000+iev, isr=isr,
                                gen_types=[genTaus[k].getID() for k in genTaus]))

        # links reco -> [(mc, w)]
        links_by_reco = {}
        for l in event.get("RecoMCTruthLink"):
            rec, sim = l.getFrom(), l.getTo()
            links_by_reco.setdefault(mcidx(rec), []).append((sim, float(l.getWeight())))

        nTausType = 0
        for i in range(len(genTaus)):
            gt = genTaus[i]
            if gt.getID() != 0:
                continue
            counters["gen0"] += 1
            visP4 = gt.getvisMomentum()
            fsr = []
            for p in mc:
                if abs(p.getPDG()) == 22 and p.getGeneratorStatus() == 1 and tau_of(p, tau_key_by_idx) == i:
                    v = p4_of(p)
                    origin = tauReco.photon_origin_info(p)[0]
                    fsr.append(dict(E=v.E(), dR=myutils.dRAngle(v, visP4), origin=int(origin),
                                    parent=[(int(a.getPDG()), int(a.getGeneratorStatus())) for a in ancestry(p)[:1]]))
            findMatch, nTausType = tauReco.MatchRecoGenTau(gt, recoTaus, nTausType,
                                                           maxDRMatch=dRMatch, selectDecay=-777)
            reco_id = recoTaus[findMatch].getID() if findMatch >= 0 else -99
            gen0_records.append(dict(evt=fidx*1000+iev, gen_key=i, visP=visP4.P(), visTheta=visP4.Theta(),
                                     reco_id=int(reco_id), fsr=fsr, true_mode=gt.getTrueMode(),
                                     n_extra_neutrals=gt.getNExtraNeutrals()))
            if findMatch < 0:
                counters["gen0_unmatched"] += 1
                continue
            counters["gen0_matched"] += 1
            rt = recoTaus[findMatch]
            if rt.getID() == 0:
                counters["gen0_reco0"] += 1
            if rt.getID() != 1:
                continue
            counters["gen0_reco1"] += 1
            # pion cargado del gen tau
            gdaus = gt.getDaughters()
            pion_mc = None
            for _, d in gdaus.items():
                if abs(d.getPDG()) in (211, 321):
                    pion_mc = d
            pion_idx = mcidx(pion_mc) if pion_mc is not None else -1
            consts = rt.getDaughters()
            lead = consts[0]
            gam = [c for kk, c in consts.items() if kk != 0 and abs(c.getPDG()) == 22]
            if len(gam) != 1:
                continue
            gam = gam[0]
            leadP4, gamP4 = p4_of(lead), p4_of(gam)
            lead_cls = [cluster_info(c) for c in lead.getClusters()]
            gam_cls = [cluster_info(c) for c in gam.getClusters()]
            # links del foton
            links = []
            for (m, w) in links_by_reco.get(mcidx(gam), []):
                tw, cw = decode_w(w)
                cat, sub = classify_link(m, i, pion_idx, tau_key_by_idx)
                chain = [(int(a.getPDG()), int(a.getGeneratorStatus())) for a in ancestry(m)[:6]]
                links.append(dict(mc_idx=mcidx(m), pdg=int(m.getPDG()), status=int(m.getGeneratorStatus()),
                                  sim=bool(m.isCreatedInSimulation()), backscatter=bool(m.isBackscatter()),
                                  E=float(m.getEnergy()), w_raw=w, tw=tw, cw=cw,
                                  cat=cat, sub=sub, chain=chain,
                                  vtx_r=math.hypot(m.getVertex().x, m.getVertex().y), vtx_z=float(m.getVertex().z)))
            # links del pion (para saber si comparten MCParticle)
            pion_links = []
            for (m, w) in links_by_reco.get(mcidx(lead), []):
                tw, cw = decode_w(w)
                pion_links.append(dict(mc_idx=mcidx(m), pdg=int(m.getPDG()), status=int(m.getGeneratorStatus()),
                                       sim=bool(m.isCreatedInSimulation()), tw=tw, cw=cw))
            rec = dict(evt=fidx*1000+iev, gen_key=i, gen_visP=gt.getvisMomentum().P(),
                       gen_visTheta=gt.getvisMomentum().Theta(),
                       gen_pion_idx=pion_idx, gen_pion_P=(p4_of(pion_mc).P() if pion_mc is not None else -1),
                       gen_pion_endpoint_r=(math.hypot(pion_mc.getEndpoint().x, pion_mc.getEndpoint().y) if pion_mc is not None else -1),
                       gen_pion_endpoint_z=(float(pion_mc.getEndpoint().z) if pion_mc is not None else 0),
                       gen_pion_decayed_in_tracker=(bool(pion_mc.isDecayedInTracker()) if pion_mc is not None else False),
                       gen_pion_decayed_in_calo=(bool(pion_mc.isDecayedInCalorimeter()) if pion_mc is not None else False),
                       gen_n_extra_neutrals=gt.getNExtraNeutrals(),
                       gen_true_mode=gt.getTrueMode(),
                       other_tau_types=[genTaus[j].getID() for j in genTaus if j != i],
                       pion_P=leadP4.P(), pion_E=leadP4.E(), pion_theta=leadP4.Theta(), pion_ntracks=int(lead.tracks_size()),
                       pion_nclusters=len(lead_cls), pion_clusters=lead_cls,
                       pion_Ecal=sum(c["E"] for c in lead_cls),
                       gam_E=gamP4.E(), gam_theta=gamP4.Theta(), gam_ntracks=int(gam.tracks_size()),
                       gam_nclusters=len(gam_cls), gam_clusters=gam_cls,
                       dR_pfo=myutils.dRAngle(gamP4, leadP4),
                       ang_cluster=(angle_between(gam_cls[0], lead_cls[0]) if gam_cls and lead_cls else -1.),
                       ang_cluster_to_pionP=(angle_between(gam_cls[0], dict(x=leadP4.Px(), y=leadP4.Py(), z=leadP4.Pz())) if gam_cls else -1.),
                       links=links, pion_links=pion_links, n_isr=len(isr),
                       isr_in_cone=sum(1 for x in isr if x["mindR_tau"] < 0.4))
            records.append(rec)
    return records, isr_records, gen0_records, counters


if __name__ == "__main__":
    files = [os.path.join(HERE, f"out_reco_edm4hep_edm4hep_{i}.root") for i in range(1000, 1025) if os.path.exists(os.path.join(HERE, f"out_reco_edm4hep_edm4hep_{i}.root"))]
    with Pool(25) as pool:
        res = pool.map(process_file, files)
    records, isr_records, gen0_records, counters = [], [], [], {}
    for r, ir, g0, c in res:
        records += r; isr_records += ir; gen0_records += g0
        for k, v in c.items():
            counters[k] = counters.get(k, 0) + v
    json.dump(dict(records=records, isr=isr_records, gen0=gen0_records, counters=counters),
              open(os.path.join(HERE, "records.json"), "w"), indent=1)
    print(counters, len(records))
