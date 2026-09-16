"""Tablas y figuras a partir de records.json (salida de photon_links.py)."""
import json, collections, math, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
d = json.load(open(os.path.join(HERE, "records.json")))
recs, g0, isr_ev, cnt = d["records"], d["gen0"], d["isr"], d["counters"]

# ---- categorias condensadas -------------------------------------------------
COND = {
    "FSR_tau": "FSR del propio tau",
    "sim_de_FSR_tau": "FSR del propio tau",          # conversion e+e- / brems del FSR
    "pion_gen_directo": "fragmento shower del pion",
    "sim_shower_pion": "fragmento shower del pion",
    "ISR": "ISR",
    "otra_particula_gen": "K0L del propio tau (tau->K0 pi nu)",
    "sim_de_mismo_tau_no_pion": "otro (mismo tau)",
    "FSR_otro_tau": "FSR del otro tau",
    "pi0_otro_tau": "pi0 del otro tau",
    "pi0_mismo_tau": "pi0 del propio tau",
    "rad_cargado_gen": "radiacion del pion en el generador",
    "sim_de_ISR": "ISR",
}
ORDER = ["FSR del propio tau", "fragmento shower del pion", "ISR",
         "K0L del propio tau (tau->K0 pi nu)", "FSR del otro tau", "pi0 del otro tau",
         "radiacion del pion en el generador", "otro (mismo tau)", "sin link"]
COL = {"FSR del propio tau": "#2a78d6", "fragmento shower del pion": "#eb6834", "ISR": "#1baf7a",
       "K0L del propio tau (tau->K0 pi nu)": "#eda100", "FSR del otro tau": "#e87ba4",
       "pi0 del otro tau": "#008300", "radiacion del pion en el generador": "#4a3aa7",
       "otro (mismo tau)": "#777777", "sin link": "#333333"}
CAT_EN = {"FSR del propio tau": "own tau's FSR", "fragmento shower del pion": "pion shower fragment",
          "ISR": "ISR", "K0L del propio tau (tau->K0 pi nu)": "own tau's K0L (tau->K0 pi nu)",
          "FSR del otro tau": "other tau's FSR", "pi0 del otro tau": "other tau's pi0",
          "radiacion del pion en el generador": "pion radiation in the generator",
          "otro (mismo tau)": "other (same tau)", "sin link": "no link"}

def best_link(r):
    return max(r["links"], key=lambda l: max(l["tw"], l["cw"])) if r["links"] else None

def condense(l):
    if l is None: return "sin link"
    c = COND.get(l["cat"], "otro (mismo tau)")
    if l["cat"] == "sim_de_mismo_tau_no_pion" and "130" in l["sub"]:
        c = "K0L del propio tau (tau->K0 pi nu)"
    if l["cat"] == "otra_particula_gen" and "130" not in l["sub"]:
        c = "otro (mismo tau)"
    return c

for r in recs:
    b = best_link(r)
    r["best"] = b; r["cat"] = condense(b)
    r["cw"] = b["cw"] if b else 0.; r["tw"] = b["tw"] if b else 0.
    r["sum_w_pion"] = sum(l["cw"] for l in r["links"] if COND.get(l["cat"]) == "fragmento shower del pion")
    r["sum_w_fsr"] = sum(l["cw"] for l in r["links"] if COND.get(l["cat"]) == "FSR del propio tau")
    c = r["gam_clusters"][0] if r["gam_clusters"] else None
    r["cl_ecal_frac"] = (c["subE"][0] / c["E"] if c and c["E"] > 0 else -1)
    r["cl_theta"] = c["theta"] if c else -1
    r["cl_region"] = ("barrel" if c and abs(c["z"]) < 2300 else "endcap") if c else "?"

N = len(recs)
cat_cnt = collections.Counter(r["cat"] for r in recs)

lines = []
P = lines.append
P("# Que es el foton extra en gen tau->pi nu reconstruidos como pi + 1 gamma (EDM4hep original)\n")
P(f"Ficheros: out_reco_edm4hep_edm4hep_1000..1024 (25 x 1000 eventos) de Ztt_SM_2M_SigmaZ (ztt_2M_smearing), leidos con podio.")
P("Reconstruccion identica a TTreesTausLong.py / config New2MSample_smearing_results0.4_tph0.0_tpi0.0_n0.0_g0.0:")
P("`findAllGenTaus(mc)`, `findAllTaus(pfos, dRMax=0.4, minP_photon=0, minP_pion=0, PNeutron=0, genminP=0)` (+electrones/muones),")
P("`MatchRecoGenTau(maxDRMatch=1)`. Script: `photon_links.py`; tablas/figuras: `summarize.py`; datos crudos: `records.json`.\n")
P("## Cifras globales\n")
P(f"| eventos | gen tau tipo 0 | matched | reco 0 | reco 1 (pi+1 gamma) | sin match |")
P(f"|---|---|---|---|---|---|")
P(f"| {cnt['events']} | {cnt['gen0']} | {cnt['gen0_matched']} | {cnt['gen0_reco0']} | {cnt['gen0_reco1']} ({100*cnt['gen0_reco1']/cnt['gen0_matched']:.1f} % de los matched) | {cnt['gen0_unmatched']} |\n")
P("Todos los 218 PFO foton tienen exactamente 1 cluster y 0 tracks, y todos tienen al menos un RecoMCTruthLink (0 casos 'sin link').")
P("Todos los links tienen peso de track = 0 (PFO neutro); el peso relevante es el de cluster (w//10000/1000).\n")

P("## Clasificacion por el link de mayor peso del PFO foton\n")
P("| categoria | N | % | <E_gamma> [GeV] | mediana E_gamma | mediana peso cluster | mediana angulo cluster-cluster pion [rad] | mediana dR(PFO gamma, PFO pion) |")
P("|---|---|---|---|---|---|---|---|")
for c in ORDER:
    rr = [r for r in recs if r["cat"] == c]
    if not rr:
        P(f"| {c} | 0 | 0.0 | - | - | - | - | - |"); continue
    E = np.array([r["gam_E"] for r in rr]); w = np.array([r["cw"] for r in rr])
    a = np.array([r["ang_cluster"] for r in rr]); dr = np.array([r["dR_pfo"] for r in rr])
    P(f"| {c} | {len(rr)} | {100*len(rr)/N:.1f} | {E.mean():.2f} | {np.median(E):.2f} | {np.median(w):.2f} | {np.median(a):.3f} | {np.median(dr):.3f} |")
P("")
P("Categorias (segun la MCParticle enlazada, subiendo por parents hasta el generador):")
P("- **FSR del propio tau**: foton status 1 del generador cuyo padre es una copia del propio tau (status 23/51/52, emision QED del shower de Pythia). Incluye los casos en que el link apunta a un e+/e- o gamma de simulacion (conversion/brems) que desciende de ese foton FSR. Estos fotones NO estan entre las hijas del tau status 2, por eso GenTauType=0 (en el tree: GenPhotonOrigin=1, GenPhotonTauKey=-1).")
P("- **fragmento shower del pion**: link directo al pion cargado del tau (status 1; los hits del calorimetro apuntan a la primaria) o a una secundaria de simulacion (generatorStatus 0, isCreatedInSimulation) cuyo ancestro generador es ese pion (gammas de pi0 secundarios, pi+-, protones, electrones del shower hadronico).")
P("- **ISR**: foton status 1 con ancestro e+/e- de haz (status 41/42/4, chain [22,43]/[11,41]).")
P("- **K0L del propio tau**: link al K0L status 1 hijo del K0 del tau (tau->K0 pi nu, GenTauTrueMode 23; extra neutral no contado en GenTauType=0). Pandora reconstruye el cluster del K0L como foton.")
P("- **FSR del otro tau**, **pi0 del otro tau**, **radiacion del pion en el generador** (fotones con padre pi+-/K+- de status 1 en el generador): ver tabla.\n")

# sub-detalle
sub = collections.Counter()
for r in recs:
    b = r["best"]
    if b: sub[(r["cat"], b["cat"], b["pdg"], b["status"], b["sim"])] += 1
P("### Detalle de la MCParticle del mejor link\n")
P("| categoria | cat. fina | PDG | genStatus | creada en sim | N |")
P("|---|---|---|---|---|---|")
for (c, f, pdg, st, sim), n in sorted(sub.items(), key=lambda x: (-x[1])):
    P(f"| {c} | {f} | {pdg} | {st} | {sim} | {n} |")
P("")

# todos los links (no solo el mejor)
nl = np.array([len(r["links"]) for r in recs])
P(f"Numero de links por PFO foton: media {nl.mean():.2f}, max {nl.max()}; con 1 solo link: {(nl==1).sum()} ({100*(nl==1).mean():.0f} %).")
mixed = [r for r in recs if r["cat"] == "FSR del propio tau" and r["sum_w_pion"] > 0.05]
P(f"PFO FSR con contaminacion del shower del pion (peso de cluster del pion > 0.05): {len(mixed)} de {cat_cnt['FSR del propio tau']}; "
  f"mediana de la fraccion de peso que va al pion en esos casos: {np.median([r['sum_w_pion'] for r in mixed]) if mixed else 0:.2f}.")
allcat = collections.Counter()
for r in recs:
    for l in r["links"]:
        allcat[COND.get(l["cat"], l["cat"])] += 1
P(f"Todos los links (no solo el mejor), por categoria: {dict(allcat)}\n")

# ---- geometria ----------------------------------------------------------------
P("## Geometria del cluster foton\n")
P("CLD: ECAL barrel r~2150-2350 mm, |z|<2300; ECAL endcap |z|~2300-2500 (r<2100). Fraccion ECAL = subdetectorEnergies[0]/E.\n")
P("| categoria | N | barrel | endcap | fraccion ECAL >0.95 | 0.30<theta<0.45 o 2.7<theta<2.85 (transicion) | mediana nhits | mediana E_gamma/E_pion PFO |")
P("|---|---|---|---|---|---|---|---|")
for c in ORDER:
    rr = [r for r in recs if r["cat"] == c]
    if not rr: continue
    nb = sum(r["cl_region"] == "barrel" for r in rr); ne = len(rr) - nb
    ecal = sum(r["cl_ecal_frac"] > 0.95 for r in rr)
    th = np.array([r["cl_theta"] for r in rr])
    trans = ((th > 0.30) & (th < 0.45) | (th > 2.70) & (th < 2.85)).sum()
    nh = np.median([r["gam_clusters"][0]["nhits"] for r in rr])
    ratio = np.median([r["gam_E"] / r["pion_E"] for r in rr])
    P(f"| {c} | {len(rr)} | {nb} | {ne} | {ecal} | {trans} | {nh:.0f} | {ratio:.3f} |")
P("")
frag = [r for r in recs if r["cat"] == "fragmento shower del pion"]
fsr = [r for r in recs if r["cat"] == "FSR del propio tau"]
def q(x, p): return np.quantile(np.array(x), p)
P("Distancia angular entre la posicion del cluster foton y la del cluster del pion (rad), cuantiles 10/50/90 %:")
for name, rr in (("fragmento shower", frag), ("FSR propio tau", fsr), ("ISR", [r for r in recs if r["cat"] == "ISR"])):
    a = [r["ang_cluster"] for r in rr]
    P(f"- {name}: {q(a,.1):.3f} / {q(a,.5):.3f} / {q(a,.9):.3f}; dR PFO-PFO: {q([r['dR_pfo'] for r in rr],.1):.3f} / {q([r['dR_pfo'] for r in rr],.5):.3f} / {q([r['dR_pfo'] for r in rr],.9):.3f}")
close = sum(r["ang_cluster"] < 0.05 for r in frag)
P(f"\nFragmentos con cluster a <0.05 rad del cluster del pion (mismo shower partido): {close}/{len(frag)}; "
  f"con cluster a >0.2 rad: {sum(r['ang_cluster']>0.2 for r in frag)}/{len(frag)} (splash lejano del shower hadronico, "
  f"E mediana {np.median([r['gam_E'] for r in frag if r['ang_cluster']>0.2]) if any(r['ang_cluster']>0.2 for r in frag) else 0:.2f} GeV).")
pi_dec = sum(r["gen_pion_decayed_in_tracker"] for r in frag)
P(f"Fragmentos cuyo pion gen interacciono/decayo en el tracker (isDecayedInTracker): {pi_dec}/{len(frag)}; "
  f"E(cluster pion)+E(gamma) vs p(pion): mediana (E_pi+E_gamma)/p = {np.median([(r['pion_Ecal']+r['gam_E'])/r['pion_P'] for r in frag]):.2f}, "
  f"E_pi/p = {np.median([r['pion_Ecal']/r['pion_P'] for r in frag]):.2f}.\n")
fsr_hi = sum(r["gam_E"] > 1 for r in fsr); frag_hi = sum(r["gam_E"] > 1 for r in frag)
P(f"E_gamma > 1 GeV: FSR {fsr_hi}/{len(fsr)}, fragmento {frag_hi}/{len(frag)}, ISR {sum(r['gam_E']>1 for r in recs if r['cat']=='ISR')}/{cat_cnt['ISR']}, "
  f"K0L {sum(r['gam_E']>1 for r in recs if r['cat'].startswith('K0L'))}/{cat_cnt['K0L del propio tau (tau->K0 pi nu)']}.")
P(f"E_gamma > 5 GeV: FSR {sum(r['gam_E']>5 for r in fsr)}/{len(fsr)}, fragmento {sum(r['gam_E']>5 for r in frag)}/{len(frag)}, "
  f"K0L {sum(r['gam_E']>5 for r in recs if r['cat'].startswith('K0L'))}.\n")

# ---- FSR gen en todos los gen0 --------------------------------------------------
P("## Fotones FSR del generador en TODOS los gen tau tipo 0 (no solo los reco 1)\n")
P("Fotones status 1 cuyo ancestro tau es el propio tau (copias status 23/51/52 incluidas).\n")
P("| reco tipo | N gen taus | con >=1 FSR gen (cualquier angulo) | con FSR gen en dR<0.4 | con FSR gen en dR<0.4 y E>0.1 GeV | con FSR en cono y E>1 GeV |")
P("|---|---|---|---|---|---|")
byid = collections.defaultdict(list)
for g in g0: byid[g["reco_id"]].append(g)
for k in sorted(byid, key=lambda k: -len(byid[k])):
    gg = byid[k]
    a = sum(bool(g["fsr"]) for g in gg)
    b = sum(any(f["dR"] < 0.4 for f in g["fsr"]) for g in gg)
    c = sum(any(f["dR"] < 0.4 and f["E"] > 0.1 for f in g["fsr"]) for g in gg)
    e = sum(any(f["dR"] < 0.4 and f["E"] > 1 for f in g["fsr"]) for g in gg)
    lab = {-99: "sin match", -20: "-20 (pi+neutron)", 0: "0 (pi)", 1: "1 (pi+1 gamma)", 2: "2 (pi+2 gamma)"}.get(k, str(k))
    P(f"| {lab} | {len(gg)} | {a} | {b} | {c} | {e} |")
allf = [f for g in g0 for f in g["fsr"]]
E = np.array([f["E"] for f in allf]); dr = np.array([f["dR"] for f in allf])
P(f"\nTotal de fotones FSR gen en {len(g0)} gen taus tipo 0: {len(allf)} ({100*sum(bool(g['fsr']) for g in g0)/len(g0):.1f} % de los taus tienen alguno); "
  f"E cuantiles 10/50/90 %: {q(E,.1):.2e} / {q(E,.5):.3f} / {q(E,.9):.2f} GeV; fraccion en dR<0.4: {(dr<0.4).mean():.2f}.")
# eficiencia de reco del foton FSR en funcion de E (gen0 con 1 FSR en cono, reco 0 vs 1)
incone = [(max(f["E"] for f in g["fsr"] if f["dR"] < 0.4), g["reco_id"]) for g in g0 if any(f["dR"] < 0.4 for f in g["fsr"]) and g["reco_id"] in (0, 1)]
edges = [0, 0.05, 0.1, 0.2, 0.5, 1, 2, 5, 50]
P("\nGen tau tipo 0 con FSR en el cono (reco 0 o 1): fraccion reconstruida como 1 en funcion de la E del FSR mas energetico del cono:\n")
P("| E_FSR [GeV] | N | reco 1 | frac |"); P("|---|---|---|---|")
for lo, hi in zip(edges[:-1], edges[1:]):
    sel = [(e, i) for e, i in incone if lo <= e < hi]
    if sel: P(f"| {lo}-{hi} | {len(sel)} | {sum(i==1 for _, i in sel)} | {sum(i==1 for _, i in sel)/len(sel):.2f} |")
P("")

# ---- ISR --------------------------------------------------------------------------
P("## ISR en el generador (todos los eventos)\n")
n_isr = np.array([len(e["isr"]) for e in isr_ev])
allisr = [x for e in isr_ev for x in e["isr"]]
Ei = np.array([x["E"] for x in allisr]); dri = np.array([x["mindR_tau"] for x in allisr])
P(f"Eventos: {len(isr_ev)}; fotones ISR (status 1, ancestro e+/e- de haz) por evento: media {n_isr.mean():.2f}, con >=1: {100*(n_isr>0).mean():.1f} %.")
P(f"E ISR cuantiles 10/50/90/99 %: {q(Ei,.1):.2e} / {q(Ei,.5):.3f} / {q(Ei,.9):.2f} / {q(Ei,.99):.2f} GeV; E>1 GeV: {100*(Ei>1).mean():.1f} %, E>5 GeV: {100*(Ei>5).mean():.1f} %.")
for ecut in (0., 0.1, 1., 5.):
    m = Ei > ecut
    P(f"- ISR con E>{ecut} GeV: {m.sum()} fotones, {100*(dri[m]<0.4).mean():.1f} % caen a dR<0.4 de algun gen tau (visible); eventos con al menos uno asi: "
      f"{100*np.mean([any(x['E']>ecut and x['mindR_tau']<0.4 for x in e['isr']) for e in isr_ev]):.1f} %.")
P(f"\nEn los 218 casos reco 1: eventos con ISR en cono de algun tau: {sum(r['isr_in_cone']>0 for r in recs)}; mejor link ISR: {cat_cnt['ISR']} "
  f"(E mediana {np.median([r['gam_E'] for r in recs if r['cat']=='ISR']):.2f} GeV).\n")

# ---- conclusiones -------------------------------------------------------------------
P("## Conclusiones\n")
P(f"1. El foton extra es, en {cat_cnt['FSR del propio tau']}/{N} casos ({100*cat_cnt['FSR del propio tau']/N:.0f} %), un foton REAL del generador: FSR emitido por el propio tau en el shower QED de Pythia (padre = copia del tau status 23/51/52). No cuelga del tau status 2, asi que la clasificacion gen (hijas del status 2) lo ignora y da GenTauType=0. La reconstruccion no se equivoca: ve un foton que esta ahi. En el tree se identifica como GenPhotonOrigin==1 con GenPhotonTauKey==-1.")
P(f"2. {cat_cnt['fragmento shower del pion']}/{N} ({100*cat_cnt['fragmento shower del pion']/N:.0f} %) son fragmentos del shower hadronico del pion (link al propio pion o a secundarias de simulacion descendientes de el): cluster ECAL de baja energia, casi siempre a <0.1-0.2 rad del cluster del pion. Son los unicos que en el tree aparecen como RecoPhotonGenMatchIdx==-1 (enlazan al pion, no a un foton gen).")
P(f"3. ISR: {cat_cnt['ISR']}/{N} ({100*cat_cnt['ISR']/N:.0f} %), fotones ISR blandos (mediana {np.median([r['gam_E'] for r in recs if r['cat']=='ISR']):.2f} GeV) que caen en el cono. Solo {100*(Ei>1).mean():.1f} % de los fotones ISR tienen E>1 GeV y {100*(Ei>5).mean():.1f} % E>5 GeV; los de E>5 GeV caen en dR<0.4 de un tau en {100*(dri[Ei>5]<0.4).mean():.1f} % de los casos (ninguno de los 218 fotones extra es un ISR duro).")
P(f"4. K0L de tau->K0 pi nu: {cat_cnt['K0L del propio tau (tau->K0 pi nu)']}/{N}; son los fotones 'extra' mas energeticos (E mediana {np.median([r['gam_E'] for r in recs if r['cat'].startswith('K0L')]) if cat_cnt['K0L del propio tau (tau->K0 pi nu)'] else 0:.1f} GeV) y en realidad el gen tau no es un pi nu puro (GenTauHasExtraNeutrals=1).")
P(f"5. FSR del otro tau / pi0 del otro tau / radiacion del pion en el generador: {cat_cnt['FSR del otro tau']} / {cat_cnt['pi0 del otro tau']} / {cat_cnt['radiacion del pion en el generador']}. Sin link: {cat_cnt['sin link']}.")
P(f"6. Discriminacion: el FSR reconstruido tiene E mediana {np.median([r['gam_E'] for r in fsr]):.2f} GeV ({100*fsr_hi/len(fsr):.0f} % con E>1 GeV; el espectro gen es 1/E pero el reco solo lo ve con eficiencia >0.6 a partir de ~0.1-0.2 GeV) y esta separado del pion (angulo cluster-cluster mediana {np.median([r['ang_cluster'] for r in fsr]):.2f} rad); el fragmento de shower es blando (mediana {np.median([r['gam_E'] for r in frag]):.2f} GeV, {100*frag_hi/len(frag):.0f} % con E>1 GeV) y pegado al cluster del pion (mediana {np.median([r['ang_cluster'] for r in frag]):.3f} rad, 44/74 a <0.05 rad). Un corte en E_gamma (o E_gamma/E_pi) recuperaria la mayoria de los dos, pero hay que medir cuanto pierde en pi pi0 (fotones de pi0 tambien blandos a bajo P del tau). Alternativa mas limpia: en la etiqueta gen, contar como constituyente del tau el FSR del propio tau (GenPhotonOrigin==1 dentro de dR<0.4), con lo que ~55 % de la 'migracion' 0->1 deja de ser migracion.")
P("\n## Figuras\n- fig_link_weight.png: peso de cluster del mejor link por categoria.\n- fig_cluster_distance.png: angulo cluster foton - cluster pion y dR PFO-PFO por categoria.\n- fig_energy.png: energia del PFO foton por categoria y E_gamma/E_pion.\n- fig_fsr_gen.png: espectro de los FSR gen en cono en gen0 y fraccion reconstruida como reco 1.")
open(os.path.join(HERE, "RESULTS.md"), "w").write("\n".join(lines))
print("\n".join(lines))

# ---- figuras ----------------------------------------------------------------------
plt.rcParams.update({"figure.facecolor": "white", "axes.facecolor": "white", "axes.grid": True,
                     "grid.alpha": 0.25, "axes.spines.top": False, "axes.spines.right": False})
cats = [c for c in ORDER if cat_cnt[c] > 0]

def stacked_hist(ax, key, bins, xlabel, log=False):
    data = [np.array([r[key] for r in recs if r["cat"] == c]) for c in cats]
    ax.hist(data, bins=bins, stacked=True, color=[COL[c] for c in cats],
            label=[f"{CAT_EN[c]} ({cat_cnt[c]})" for c in cats], edgecolor="white", linewidth=0.5)
    ax.set_xlabel(xlabel); ax.set_ylabel("photon PFOs / bin")
    if log: ax.set_xscale("log")

fig, ax = plt.subplots(figsize=(8, 4.8), dpi=130)
stacked_hist(ax, "cw", np.linspace(0, 1, 21), "cluster weight of the best RecoMCTruthLink")
ax.set_title("Link weight of the extra photon (gen tau->pi nu, reco pi+1 gamma)")
ax.legend(fontsize=8, loc="upper left"); fig.tight_layout(); fig.savefig(os.path.join(HERE, "fig_link_weight.png")); plt.close(fig)

fig, axs = plt.subplots(1, 2, figsize=(11, 4.6), dpi=130)
stacked_hist(axs[0], "ang_cluster", np.linspace(0, 0.5, 26), "angle photon cluster - pion cluster [rad]")
axs[0].set_title("Distance between clusters")
stacked_hist(axs[1], "dR_pfo", np.linspace(0, 0.4, 21), "dR(photon PFO, pion PFO) [rad]")
axs[1].set_title("Distance between PFOs (theta-phi)")
axs[0].legend(fontsize=7, loc="upper right"); fig.tight_layout(); fig.savefig(os.path.join(HERE, "fig_cluster_distance.png")); plt.close(fig)

for r in recs: r["ratio"] = r["gam_E"] / r["pion_E"]
fig, axs = plt.subplots(1, 2, figsize=(11, 4.6), dpi=130)
stacked_hist(axs[0], "gam_E", np.logspace(-1.3, 1.7, 25), "photon PFO E [GeV]", log=True)
axs[0].set_title("Energy of the extra photon")
stacked_hist(axs[1], "ratio", np.logspace(-3, 0.7, 25), "E_gamma / E_pion (PFO)", log=True)
axs[1].set_title("Energy relative to the pion")
axs[0].legend(fontsize=7, loc="upper right"); fig.tight_layout(); fig.savefig(os.path.join(HERE, "fig_energy.png")); plt.close(fig)

fig, axs = plt.subplots(1, 2, figsize=(11, 4.6), dpi=130)
e0 = np.array([e for e, i in incone if i == 0]); e1 = np.array([e for e, i in incone if i == 1])
bins = np.logspace(-4, 1.7, 30)
axs[0].hist([e0, e1], bins=bins, stacked=True, color=["#2a78d6", "#eb6834"],
            label=[f"reco 0 (pi) ({len(e0)})", f"reco 1 (pi+1 gamma) ({len(e1)})"], edgecolor="white", linewidth=0.5)
axs[0].set_xscale("log"); axs[0].set_xlabel("E of the most energetic gen FSR at dR<0.4 [GeV]"); axs[0].set_ylabel("type-0 gen taus / bin")
axs[0].set_title("Own tau's gen FSR in the cone (gen tau->pi nu)"); axs[0].legend(fontsize=8)
h0, _ = np.histogram(e0, bins); h1, _ = np.histogram(e1, bins)
tot = h0 + h1; frac = np.where(tot > 0, h1 / np.maximum(tot, 1), np.nan)
err = np.where(tot > 0, np.sqrt(frac * (1 - frac) / np.maximum(tot, 1)), np.nan)
cen = np.sqrt(bins[:-1] * bins[1:])
axs[1].errorbar(cen, frac, yerr=err, fmt="o", color="#2a78d6", ms=4, label="reco 1 / (reco 0 + reco 1)")
axs[1].set_xscale("log"); axs[1].set_ylim(0, 1.05); axs[1].set_xlabel("E of the gen FSR in the cone [GeV]"); axs[1].set_ylabel("fraction reconstructed as pi+1 gamma")
axs[1].set_title("Efficiency of seeing the FSR as a photon"); axs[1].legend(fontsize=8)
fig.tight_layout(); fig.savefig(os.path.join(HERE, "fig_fsr_gen.png")); plt.close(fig)
