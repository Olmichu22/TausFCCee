"""Tablas y figuras a partir de records.json (salida de raw_links.py).

Uso: python summary.py [records.json]
Imprime fracciones por categoria (a-g) separadas por RecoTauType 0/1/2/3+ y
por bin de P visible gen del tau, y escribe fig_*.png en este directorio.
"""
import sys, os, json, math, collections
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
fn = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "records.json")
D = json.load(open(fn))
R, C = D["records"], D["counters"]

CATS = ["a_lost", "b_in_cone", "c_out_cone", "d_charged", "e_neutralh", "f_merged", "g_split", "other"]
CAT_LABEL = {"a_lost": "(a) perdido", "b_in_cone": "(b) PFO γ en cono", "c_out_cone": "(c) PFO γ fuera cono",
             "d_charged": "(d) PFO cargado", "e_neutralh": "(e) PFO n", "f_merged": "(f) fusionados",
             "g_split": "(g) partido", "other": "otro"}
# paleta categorica de referencia (orden fijo, un color por categoria)
COL = {"a_lost": "#2a78d6", "b_in_cone": "#eb6834", "c_out_cone": "#1baf7a", "d_charged": "#eda100",
       "e_neutralh": "#e87ba4", "f_merged": "#008300", "g_split": "#4a3aa7", "other": "#e34948"}
GROUPS = ["0", "1", "2", "3+"]
PBINS = [(0, 10), (10, 20), (20, 30), (30, 46)]


def grp(r):
    return "3+" if r["reco_id"] >= 3 else str(r["reco_id"])


def pbin(r):
    for lo, hi in PBINS:
        if lo <= r["gen_visP"] < hi:
            return f"{lo}-{hi}"
    return None


def table(title, rows, cols, counts, fmt="pct"):
    """counts[row][col] -> imprime tabla con N y % por fila."""
    print(f"\n### {title}\n")
    print("| " + " | ".join(["", "N"] + cols) + " |")
    print("|" + "---|" * (len(cols) + 2))
    for rw in rows:
        n = sum(counts[rw].values())
        cells = []
        for c in cols:
            v = counts[rw].get(c, 0)
            cells.append(f"{100.*v/n:.1f}" if n and fmt == "pct" else str(v))
        print("| " + " | ".join([rw, str(n)] + cells) + " |")


# ------------------------------------------------------------------ recuento
print("## Contadores\n")
for k in sorted(C):
    print(f"- {k}: {C[k]}")
nt = collections.Counter(grp(r) for r in R)
print("\n- taus gen tipo 1 con reco 0/1/2/3+ y lider pion:", dict(nt), "total", len(R))

# fotones
ph = [(r, p) for r in R for p in r["photons"]]
cnt = {g: collections.Counter() for g in GROUPS}
for r, p in ph:
    cnt[grp(r)][p["cat"] if p["cat"] in CATS else "other"] += 1
table("Fraccion (%) de fotones del pi0 por categoria y RecoTauType (2 fotones por tau)", GROUPS, CATS, cnt)

# por bin de P visible
for g in GROUPS:
    c2 = {f"{lo}-{hi}": collections.Counter() for lo, hi in PBINS}
    for r, p in ph:
        if grp(r) == g and pbin(r):
            c2[pbin(r)][p["cat"] if p["cat"] in CATS else "other"] += 1
    table(f"RecoTauType {g}: fraccion (%) por categoria y bin de P visible gen (GeV)", list(c2), CATS, c2)

# subcategorias
sub = {g: collections.Counter() for g in GROUPS}
for r, p in ph:
    sub[grp(r)][p["cat"] + ":" + p["sub"]] += 1
subs = sorted({k for g in GROUPS for k in sub[g]})
print("\n### Subcategorias (N) por RecoTauType\n")
print("| subcategoria | " + " | ".join(GROUPS) + " |")
print("|---|" + "---|" * len(GROUPS))
for s in subs:
    print(f"| {s} | " + " | ".join(str(sub[g].get(s, 0)) for g in GROUPS) + " |")

# patron por tau (par de categorias)
pat = {g: collections.Counter() for g in GROUPS}
for r in R:
    k = tuple(sorted(p["cat"][0] for p in r["photons"]))
    pat[grp(r)]["+".join(k)] += 1
pats = sorted({k for g in GROUPS for k in pat[g]})
print("\n### Patron (categoria gamma1 + gamma2) por tau, % por RecoTauType\n")
print("| patron | " + " | ".join(GROUPS) + " |")
print("|---|" + "---|" * len(GROUPS))
for s in pats:
    print(f"| {s} | " + " | ".join(f"{100.*pat[g].get(s,0)/max(1,sum(pat[g].values())):.1f}" for g in GROUPS) + " |")

# propiedades de los perdidos
print("\n### Fotones perdidos (a): propiedades por RecoTauType\n")
print("| reco | N | P mediana (GeV) | P<0.2 GeV (%) | P<0.5 GeV (%) | |cos θ|>0.98 (%) | convertidos (%) | ang γ-π mediana (rad) |")
print("|---|---|---|---|---|---|---|---|")
for g in GROUPS:
    L = [p for r, p in ph if grp(r) == g and p["cat"] == "a_lost"]
    if not L:
        continue
    P = np.array([p["P"] for p in L]); ct = np.abs(np.cos([p["theta"] for p in L]))
    conv = np.array(["converted" in p["sub"] for p in L])
    ang = np.median([p["ang_pion"] for p in L])
    print(f"| {g} | {len(L)} | {np.median(P):.2f} | {100*np.mean(P<0.2):.0f} | {100*np.mean(P<0.5):.0f} | {100*np.mean(ct>0.98):.0f} | {100*np.mean(conv):.0f} | {ang:.3f} |")

# propiedades por categoria (todas las reco juntas)
print("\n### P gen y angulo γ-π por categoria (todos los RecoTauType)\n")
print("| cat | N | P mediana (GeV) | P p10 | P p90 | ang γ-π mediana (rad) | ang γ-γ mediana (rad) |")
print("|---|---|---|---|---|---|---|")
for c in CATS:
    L = [p for r, p in ph if p["cat"] == c]
    if not L:
        continue
    P = np.array([p["P"] for p in L])
    print(f"| {CAT_LABEL[c]} | {len(L)} | {np.median(P):.2f} | {np.percentile(P,10):.2f} | {np.percentile(P,90):.2f} | {np.median([p['ang_pion'] for p in L]):.3f} | {np.median([p['ang_other'] for p in L]):.3f} |")

# (d) detalles
print("\n### (d) PFO cargado: desglose y E/p del pion lider\n")
dch = [(r, p) for r, p in ph if p["cat"] == "d_charged"]
cc = collections.Counter(p["sub"] for r, p in dch)
for k, v in cc.most_common():
    print(f"- {k}: {v}")
absorbed = [r for r in R if any(p["cat"] == "d_charged" and p["sub"].startswith("lead_pion") for p in r["photons"])]
clean = [r for r in R if all(p["cat"] in ("b_in_cone", "a_lost", "c_out_cone") for p in r["photons"]) and r["lead"]["cw_from_gammas"] == 0]
ep_abs = np.array([r["lead"]["EoverP"] for r in absorbed if r["lead"]["EoverP"] > 0])
ep_abs_no = np.array([r["lead"]["EoverP_nogamma"] for r in absorbed if r["lead"]["EoverP"] > 0])
ep_clean = np.array([r["lead"]["EoverP"] for r in clean if r["lead"]["EoverP"] > 0])
if len(ep_abs):
    print(f"- pion absorbe un foton: N={len(ep_abs)}, E/p mediana {np.median(ep_abs):.2f} (media {ep_abs.mean():.2f}); "
          f"restando la parte del foton (RecoMCTruthLink): {np.median(ep_abs_no):.2f} (media {ep_abs_no.mean():.2f})")
    print(f"- E/p>1: {100*np.mean(ep_abs>1):.0f}% absorbiendo vs {100*np.mean(ep_clean>1):.0f}% sin absorber")
    fr = np.array([r["lead"]["cw_from_gammas"] + r["lead"]["cw_from_gamma_sim"] for r in absorbed])
    print(f"- fraccion del cluster del pion que viene del foton: mediana {np.median(fr):.2f}")
    Pg = np.array([p["P"] for r in absorbed for p in r["photons"] if p["cat"] == "d_charged" and p["sub"].startswith("lead_pion")])
    ang = np.array([p["ang_pion"] for r in absorbed for p in r["photons"] if p["cat"] == "d_charged" and p["sub"].startswith("lead_pion")])
    print(f"- P gen del foton absorbido: mediana {np.median(Pg):.2f} GeV; angulo γ-π mediana {np.median(ang):.3f} rad, {100*np.mean(ang<0.05):.0f}% a <0.05 rad")
print(f"- pion sin foton absorbido (control): N={len(ep_clean)}, E/p mediana {np.median(ep_clean):.2f} (media {ep_clean.mean():.2f})")

# (b)/(f)/(g) pureza
print("\n### Pureza del PFO principal (fraccion del PFO que viene del foton, RecoMCTruthLink)\n")
for c in ("b_in_cone", "f_merged", "g_split", "c_out_cone"):
    L = [p["main"]["frac_pfo_from_gamma"] for r, p in ph if p["cat"] == c and p["main"]]
    if L:
        print(f"- {CAT_LABEL[c]}: N={len(L)}, mediana {np.median(L):.2f}, <0.5: {100*np.mean(np.array(L)<0.5):.0f}%")
L = [p["main"]["E"] / p["E"] for r, p in ph if p["cat"] == "b_in_cone" and p["main"]]
print(f"- (b) E(PFO)/E(gen): mediana {np.median(L):.3f}")
L = [p["main"]["E"] / (r["gen_pi0_P"]) for r, p in ph if p["cat"] == "f_merged" and p["main"]]
if L:
    print(f"- (f) E(PFO)/P(pi0 gen): mediana {np.median(L):.3f}")

# (c) fuera del cono: distancia
L = [p["main"]["dR_lead"] for r, p in ph if p["cat"] == "c_out_cone" and p["main"]]
if L:
    print(f"\n- (c) dR PFO-lider: mediana {np.median(L):.2f} rad, >1 rad: {100*np.mean(np.array(L)>1):.0f}%")
    L2 = [p["ang_pion"] for r, p in ph if p["cat"] == "c_out_cone"]
    print(f"- (c) angulo gen γ-π: mediana {np.median(L2):.2f} rad")

# conversiones en el tracker (flag isDecayedInTracker del foton gen)
print("\n### Conversion en el tracker (isDecayedInTracker): % de fotones por categoria y RecoTauType\n")
print("| cat | " + " | ".join(f"{g} (N)" for g in GROUPS) + " | todos |")
print("|---|" + "---|" * (len(GROUPS) + 1))
for c in CATS:
    cells = []
    for g in GROUPS + [None]:
        L = [p["decayed_in_tracker"] for r, p in ph if p["cat"] == c and (g is None or grp(r) == g)]
        cells.append(f"{100*np.mean(L):.0f} ({len(L)})" if L else "-")
    if any(x != "-" for x in cells):
        print(f"| {CAT_LABEL[c]} | " + " | ".join(cells) + " |")
print("\nTaus con >=1 foton del pi0 convertido en el tracker, % por RecoTauType:")
for g in GROUPS:
    L = [any(p["decayed_in_tracker"] for p in r["photons"]) for r in R if grp(r) == g]
    print(f"- reco {g}: {100*np.mean(L):.1f}% (N={len(L)})")
L = [p["endpoint_r"] for r, p in ph if p["decayed_in_tracker"]]
print(f"- radio de conversion (mm): p10 {np.percentile(L,10):.0f}, mediana {np.median(L):.0f}, p90 {np.percentile(L,90):.0f}")
for c, lab in (("d_charged", "pdg11_conv"), ("e_neutralh", None), ("g_split", None)):
    L = [p["endpoint_r"] for r, p in ph if p["cat"] == c and p["decayed_in_tracker"] and (lab is None or p["sub"] == lab)]
    if L:
        print(f"- radio de conversion {CAT_LABEL[c]}{' '+lab if lab else ''}: mediana {np.median(L):.0f} mm, <30 mm: {100*np.mean(np.array(L)<30):.0f}%")

# ------------------------------------------------------------ reco 3+: extras
print("\n### RecoTauType 3+: clasificacion de los PFO foton del cono\n")
r3 = [r for r in R if r["reco_id"] >= 3]
ex = collections.Counter(); allc = collections.Counter()
for r in r3:
    for q in r["reco_photons"]:
        allc[q["cls"]] += 1
        if not q["is_main"]:
            ex[q["cls"]] += 1
print(f"N taus 3+: {len(r3)}, PFO foton: {sum(allc.values())}, extra (no principal de ningun foton del pi0): {sum(ex.values())}\n")
print("| clase | todos | extra | % extra | E mediana extra (GeV) | dR lider mediana |")
print("|---|---|---|---|---|---|")
for k, v in ex.most_common():
    E = np.median([q["E"] for r in r3 for q in r["reco_photons"] if not q["is_main"] and q["cls"] == k])
    dr = np.median([q["dR_lead"] for r in r3 for q in r["reco_photons"] if not q["is_main"] and q["cls"] == k])
    print(f"| {k} | {allc[k]} | {v} | {100.*v/sum(ex.values()):.1f} | {E:.2f} | {dr:.3f} |")
for k in allc:
    if k not in ex:
        print(f"| {k} | {allc[k]} | 0 | 0 | | |")
# patron por tau en 3+: cual es el extra
pat3 = collections.Counter()
for r in r3:
    ks = sorted(q["cls"] for q in r["reco_photons"] if not q["is_main"])
    pat3["+".join(ks) if ks else "(ninguno extra)"] += 1
Ef = [q["E"] / r["photons"][int(q["sub"][-1])]["E"] for r in r3 for q in r["reco_photons"] if q["cls"] == "pi0_photon_conv_frag" and not q["is_main"]]
if Ef:
    print(f"\n- fragmentos de conversion: E(extra)/E(gamma gen) mediana {np.median(Ef):.2f}; el foton gen padre tiene isDecayedInTracker en "
          f"{100*np.mean([r['photons'][int(q['sub'][-1])]['decayed_in_tracker'] for r in r3 for q in r['reco_photons'] if q['cls']=='pi0_photon_conv_frag' and not q['is_main']]):.0f}%")
E = np.array([q["E"] for r in r3 for q in r["reco_photons"] if q["cls"] == "FSR_tau"])
if len(E):
    print(f"- FSR extra: E mediana {np.median(E):.2f} GeV, E>1 GeV {100*np.mean(E>1):.0f}%")
print("\nPatron de extras por tau 3+:")
for k, v in pat3.most_common(12):
    print(f"- {k}: {v} ({100.*v/len(r3):.1f}%)")

# -------------------------------------------------------------------- figuras
plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False,
                     "axes.grid": True, "grid.alpha": 0.25, "grid.linewidth": 0.6})

# 1. P gen del foton por categoria
fig, axs = plt.subplots(1, 2, figsize=(11, 4.2))
bins = np.logspace(-3, np.log10(45), 40)
for c in CATS:
    P = [p["P"] for r, p in ph if p["cat"] == c]
    if len(P) < 5:
        continue
    axs[0].hist(P, bins=bins, histtype="step", lw=1.8, color=COL[c], label=f"{CAT_LABEL[c]} (N={len(P)})", density=True)
axs[0].set_xscale("log"); axs[0].set_xlabel("P gen del fotón del π0 (GeV)"); axs[0].set_ylabel("densidad")
axs[0].set_title("P gen del fotón por categoría (todos los RecoTauType)"); axs[0].legend(fontsize=7.5, frameon=False)
for g, ls in zip(GROUPS, ["-", "--", "-.", ":"]):
    P = [p["P"] for r, p in ph if p["cat"] == "a_lost" and grp(r) == g]
    if len(P) < 5:
        continue
    axs[1].hist(P, bins=bins, histtype="step", lw=1.8, ls=ls, color=COL["a_lost"], label=f"reco {g} (N={len(P)})", density=True)
axs[1].set_xscale("log"); axs[1].set_xlabel("P gen del fotón perdido (GeV)"); axs[1].set_title("(a) fotón perdido, por RecoTauType")
axs[1].legend(fontsize=8, frameon=False)
fig.tight_layout(); fig.savefig(os.path.join(HERE, "fig_photon_P_by_cat.png"), dpi=140); plt.close(fig)

# 2. E/p del pion
fig, ax = plt.subplots(figsize=(6.5, 4.2))
bins = np.linspace(0, 2.0, 41)
ax.hist(ep_clean, bins=bins, histtype="step", lw=1.8, color=COL["b_in_cone"], density=True, label=f"π sin fotón absorbido (N={len(ep_clean)})")
if len(ep_abs):
    ax.hist(ep_abs, bins=bins, histtype="step", lw=1.8, color=COL["d_charged"], density=True, label=f"π absorbe un fotón del π0 (N={len(ep_abs)})")
    ax.hist(ep_abs_no, bins=bins, histtype="step", lw=1.4, ls="--", color=COL["d_charged"], density=True, label="mismo, restando la parte del fotón")
ax.set_xlabel("E(cluster) / p(PFO) del pión líder"); ax.set_ylabel("densidad"); ax.legend(fontsize=8, frameon=False)
ax.set_title("E/p del pión líder")
fig.tight_layout(); fig.savefig(os.path.join(HERE, "fig_EoverP_pion.png"), dpi=140); plt.close(fig)

# 3. angulo foton-pion por categoria
fig, axs = plt.subplots(1, 2, figsize=(11, 4.2))
bins = np.logspace(-3, np.log10(1.5), 36)
for c in CATS:
    A = [p["ang_pion"] for r, p in ph if p["cat"] == c]
    if len(A) < 5:
        continue
    axs[0].hist(A, bins=bins, histtype="step", lw=1.8, color=COL[c], label=f"{CAT_LABEL[c]} (N={len(A)})", density=True)
axs[0].axvline(0.4, color="gray", lw=0.8, ls=":"); axs[0].set_xscale("log")
axs[0].set_xlabel("ángulo gen fotón – pión (rad)"); axs[0].set_ylabel("densidad"); axs[0].legend(fontsize=7.5, frameon=False)
axs[0].set_title("Ángulo γ–π por categoría")
for c in CATS:
    A = [p["ang_other"] for r, p in ph if p["cat"] == c]
    if len(A) < 5:
        continue
    axs[1].hist(A, bins=bins, histtype="step", lw=1.8, color=COL[c], density=True)
axs[1].set_xscale("log"); axs[1].set_xlabel("ángulo gen fotón – fotón (rad)"); axs[1].set_title("Ángulo γ–γ del π0 por categoría")
fig.tight_layout(); fig.savefig(os.path.join(HERE, "fig_angle_by_cat.png"), dpi=140); plt.close(fig)
print("\nfiguras escritas en", HERE)
