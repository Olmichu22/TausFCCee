# Figura de la seccion 1.1 (fig_visP_no_separa.png). Recuperado del codigo inline original; lee ptau_origin.npz.
import numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
d=np.load("ptau_origin.npz"); gP,vP,rt=[d[k] for k in ["gP","gvisP","rtype"]]
rad=gP<44; b=np.arange(0,50,2.5); c=0.5*(b[1:]+b[:-1])
fig,ax=plt.subplots(1,3,figsize=(16,4.6))
ax[0].hist(vP[~rad],bins=b,histtype="step",lw=1.8,color="#0072B2",label="no radiativos (P_τ ≥ 44 GeV)")
ax[0].hist(vP[rad],bins=b,histtype="step",lw=1.8,color="#D55E00",label="radiativos (P_τ < 44 GeV)")
ax[0].set_yscale("log"); ax[0].set_xlabel("P visible gen del tau (GeV)"); ax[0].set_ylabel("gen τ→πν / bin")
ax[0].set_title("Los radiativos se reparten por todo el P visible"); ax[0].legend(fontsize=8); ax[0].grid(alpha=.3)
f=[rad[(vP>=lo)&(vP<hi)].mean() for lo,hi in zip(b[:-1],b[1:])]
ax[1].plot(c,f,"-o",color="#D55E00"); ax[1].set_ylim(0,0.12)
ax[1].set_xlabel("P visible gen del tau (GeV)"); ax[1].set_ylabel("fracción radiativa del bin")
ax[1].set_title("Ningún bin de P visible aísla a los radiativos"); ax[1].grid(alpha=.3)
for m,col,lab in [(~rad,"#0072B2","no radiativos"),(rad,"#D55E00","radiativos"),(np.ones_like(rad,bool),"k","todos")]:
    y=[np.mean(rt[m&(vP>=lo)&(vP<hi)&((rt==0)|(rt==1))]==1) for lo,hi in zip(b[:-1],b[1:])]
    ax[2].plot(c,y,"-o",color=col,label=lab,lw=2 if col=="k" else 1.5)
ax[2].set_xlabel("P visible gen del tau (GeV)"); ax[2].set_ylabel("migración reco1/(reco0+reco1)")
ax[2].set_title("El 4-5 % plano es la mezcla de dos poblaciones"); ax[2].legend(fontsize=8); ax[2].grid(alpha=.3)
plt.tight_layout(); plt.savefig("fig_visP_no_separa.png",dpi=130)
print("ok")
