# Tablas de la seccion 1.1 (radiativos frente a P visible). Recuperado del codigo inline original; lee ptau_origin.npz.
import numpy as np
d=np.load("ptau_origin.npz"); gP,vP,rt,cat=[d[k] for k in ["gP","gvisP","rtype","cat"]]
rad=gP<44   # taus que han radiado (P total lejos de E_haz)
print("gen tau->pi nu (TrueMode 10), N=",len(gP))
print("frac radiativos (P_tau<44):",round(rad.mean(),4)," (P_tau 20-40):",round(((gP>=20)&(gP<40)).mean(),4))
print()
print("Distribucion en VisP: donde caen los radiativos?")
bins=np.arange(0,50,5)
print(" VisP    N_total  frac_del_bin_radiativa  N_radiativos  frac_de_todos_los_radiativos_en_ese_bin")
for lo,hi in zip(bins[:-1],bins[1:]):
    m=(vP>=lo)&(vP<hi)
    print(f" {lo:2d}-{hi:2d} {m.sum():8d} {rad[m].mean():10.4f} {(m&rad).sum():12d} {(m&rad).sum()/rad.sum():14.4f}")
print()
print("Migracion 0->1 por bin de VisP, separando radiativos y no radiativos:")
print(" VisP    N_norad  mig_norad   N_rad  mig_rad   mig_total")
for lo,hi in zip(bins[:-1],bins[1:]):
    m=(vP>=lo)&(vP<hi)&((rt==0)|(rt==1))
    a=m&~rad; b=m&rad
    print(f" {lo:2d}-{hi:2d} {a.sum():8d} {np.mean(rt[a]==1):9.4f} {b.sum():7d} {np.mean(rt[b]==1) if b.sum() else 0:8.4f} {np.mean(rt[m]==1):10.4f}")
print()
print("Un corte en VisP: que se lleva por delante?")
for cut in [5,10,15,20,25,30]:
    keep=vP>=cut
    print(f" VisP>{cut:2d} GeV: conserva {keep.mean():.3f} de todos los pi-nu, {(keep&rad).sum()/rad.sum():.3f} de los radiativos, pureza radiativa {rad[keep].mean():.4f} (antes {rad.mean():.4f})")
for cut in [40,35,30]:
    keep=vP<=cut
    print(f" VisP<{cut:2d} GeV: conserva {keep.mean():.3f} de todos, {(keep&rad).sum()/rad.sum():.3f} de los radiativos, pureza radiativa {rad[keep].mean():.4f}")
