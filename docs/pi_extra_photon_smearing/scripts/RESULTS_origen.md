# Origen del foton extra en gen tau->pi nu reconstruido como pi+gamma

Muestra ztt_2M_smearing, 2000000 eventos procesados (fichero completo). Tau reco: cono dR<0.4, P(foton)>0.5 GeV (genminP).

## 1-2. Clasificacion del foton asignado (gen tau->pi nu, RecoTauType==1)

Categorias: `FSR tau (shower)` = GenPhotonOrigin 1 con GenPhotonTauKey -1 (foton emitido por una copia Pythia del tau antes de la desintegracion; se asigna al tau por angulo, dR<0.4); `FSR desint. (propio)` = origen 1 y TauKey == tau propio (foton hijo directo del tau en la desintegracion, tau->pi nu gamma); `ISR (haz)` = origen 2 con parent |PDG|==11; `rad. pion` = origen 2 con parent 211; `pi0 otro tau` = origen 0 de otro tau; `sin match gen` = RecoPhotonGenMatchIdx==-1 (el PFO no enlaza a ningun foton del generador: fragmento de shower del pion, cluster partido o link al propio pion).


### Global, GenTauType==0, RecoTauType==1

| N | FSR tau (shower) | ISR (haz) | sin match gen | pi0 otro tau | otro (eta...) |
|---|---|---|---|---|---|
| 19932 | 0.540 | 0.030 | 0.430 | 0.000 | 0.000 |

### GenTauType==0, RecoTauType==1, pi nu puro (TrueMode 10)

| N | FSR tau (shower) | ISR (haz) | sin match gen | pi0 otro tau |
|---|---|---|---|---|
| 18020 | 0.563 | 0.031 | 0.406 | 0.000 |

### GenTauType==0, RecoTauType==1, K nu (TrueMode 20)

| N | FSR tau (shower) | ISR (haz) | sin match gen |
|---|---|---|---|
| 1213 | 0.494 | 0.030 | 0.477 |

### GenTauType==0, RecoTauType==1, K0 pi (TrueMode 23)

| N | FSR tau (shower) | sin match gen |
|---|---|---|
| 578 | 0.016 | 0.984 |

Otros TrueMode en type0/reco1: [(np.int32(27), 112), (np.int32(41), 5), (np.int32(28), 4)]

### GenTauType==0, RecoTauType==1, resto de TrueMode

| N | FSR tau (shower) | sin match gen | otro (eta...) |
|---|---|---|---|
| 121 | 0.025 | 0.934 | 0.041 |

### TrueMode 10, RecoTauType==1, por GenVisTauP

| GenVisTauP (GeV) | N | FSR tau (shower) | ISR (haz) | sin match gen | pi0 otro tau |
|---|---|---|---|---|---|
| 0-5 | 1963 | 0.674 | 0.031 | 0.295 | 0.000 |
| 5-10 | 2263 | 0.654 | 0.033 | 0.312 | 0.001 |
| 10-15 | 2273 | 0.613 | 0.029 | 0.358 | 0.000 |
| 15-20 | 2191 | 0.590 | 0.033 | 0.377 | 0.000 |
| 20-25 | 2176 | 0.545 | 0.028 | 0.427 | 0.001 |
| 25-30 | 2045 | 0.545 | 0.024 | 0.431 | 0.000 |
| 30-35 | 1905 | 0.519 | 0.033 | 0.447 | 0.001 |
| 35-40 | 1723 | 0.476 | 0.033 | 0.491 | 0.000 |
| 40-45 | 1371 | 0.392 | 0.036 | 0.572 | 0.000 |
| 45-50 | 110 | 0.118 | 0.018 | 0.864 | 0.000 |

### TrueMode 10, RecoTauType==1, por |cos theta| del tau

| |cos theta| | N | FSR tau (shower) | ISR (haz) | sin match gen | pi0 otro tau |
|---|---|---|---|---|---|
| 0-0.2 | 3131 | 0.642 | 0.023 | 0.334 | 0.001 |
| 0.2-0.4 | 3358 | 0.613 | 0.023 | 0.364 | 0.001 |
| 0.4-0.6 | 3751 | 0.573 | 0.023 | 0.403 | 0.000 |
| 0.6-0.7 | 1909 | 0.506 | 0.019 | 0.475 | 0.000 |
| 0.7-0.75 | 1030 | 0.529 | 0.032 | 0.439 | 0.000 |
| 0.75-0.8 | 1015 | 0.560 | 0.033 | 0.408 | 0.000 |
| 0.8-0.85 | 1010 | 0.521 | 0.041 | 0.439 | 0.000 |
| 0.85-0.9 | 1091 | 0.527 | 0.056 | 0.417 | 0.000 |
| 0.9-0.95 | 1007 | 0.463 | 0.073 | 0.462 | 0.002 |
| 0.95-1 | 718 | 0.400 | 0.052 | 0.547 | 0.001 |

### Fraccion de migracion 0->1 (gen pi nu puro) por |cos theta| del tau (respecto a los reconstruidos con tipo 0 o 1)

| bin cos | N(reco0+reco1) | reco1/(reco0+reco1) | reco1 sin match gen | reco1 FSR shower | reco1 ISR |
|---|---|---|---|---|---|
| 0-0.2 | 57615 | 0.0543 | 0.0182 | 0.0349 | 0.0012 |
| 0.2-0.4 | 61727 | 0.0544 | 0.0198 | 0.0333 | 0.0012 |
| 0.4-0.6 | 67715 | 0.0554 | 0.0223 | 0.0318 | 0.0013 |
| 0.6-0.7 | 33253 | 0.0575 | 0.0272 | 0.0291 | 0.0011 |
| 0.7-0.75 | 19600 | 0.0526 | 0.0231 | 0.0278 | 0.0017 |
| 0.75-0.8 | 21527 | 0.0472 | 0.0192 | 0.0264 | 0.0015 |
| 0.8-0.85 | 22651 | 0.0446 | 0.0196 | 0.0232 | 0.0018 |
| 0.85-0.9 | 24354 | 0.0448 | 0.0187 | 0.0236 | 0.0025 |
| 0.9-0.95 | 25143 | 0.0402 | 0.0185 | 0.0185 | 0.0029 |
| 0.95-1 | 20899 | 0.0345 | 0.0188 | 0.0137 | 0.0018 |

### Fraccion de migracion 0->1 por GenVisTauP (gen pi nu puro)

| P (GeV) | N(reco0+reco1) | reco1/(reco0+reco1) | sin match gen | FSR shower | ISR | pi0 otro tau |
|---|---|---|---|---|---|---|
| 0-5 | 44368 | 0.0443 | 0.0130 | 0.0298 | 0.0014 | 0.0000 |
| 5-10 | 45367 | 0.0499 | 0.0155 | 0.0326 | 0.0016 | 0.0001 |
| 10-15 | 43464 | 0.0524 | 0.0187 | 0.0320 | 0.0015 | 0.0000 |
| 15-20 | 41336 | 0.0531 | 0.0200 | 0.0313 | 0.0018 | 0.0000 |
| 20-25 | 39567 | 0.0550 | 0.0235 | 0.0299 | 0.0015 | 0.0001 |
| 25-30 | 37357 | 0.0548 | 0.0236 | 0.0298 | 0.0013 | 0.0000 |
| 30-35 | 35582 | 0.0535 | 0.0239 | 0.0278 | 0.0018 | 0.0000 |
| 35-40 | 33636 | 0.0513 | 0.0252 | 0.0244 | 0.0017 | 0.0000 |
| 40-45 | 30766 | 0.0446 | 0.0255 | 0.0175 | 0.0016 | 0.0000 |
| 45-50 | 3041 | 0.0362 | 0.0312 | 0.0043 | 0.0007 | 0.0000 |

## 3. Cinematica del foton extra por categoria, comparada con fotones legitimos de pi0 (gen pi pi0, reco 1 o 2)

| categoria | N | P_gamma med (GeV) | P_gamma p10/p90 | P_g/P_pi med | dR(g,pi) med | dR p10/p90 | m(pi+g) med (GeV) | frac dR<0.05 | frac P_g<1 GeV | frac P_g<2 GeV |
|---|---|---|---|---|---|---|---|---|---|---|
| FSR tau (shower) | 10761 | 1.42 | 0.19/12.92 | 0.091 | 0.188 | 0.069/0.346 | 0.720 | 0.051 | 0.425 | 0.564 |
| ISR (haz) | 588 | 0.39 | 0.13/1.69 | 0.023 | 0.285 | 0.128/0.383 | 0.580 | 0.019 | 0.794 | 0.927 |
| sin match gen | 8570 | 0.58 | 0.15/6.22 | 0.031 | 0.069 | 0.020/0.319 | 0.259 | 0.390 | 0.662 | 0.803 |
| pi0 otro tau | 8 | 1.03 | 0.20/4.15 | 0.185 | 0.308 | 0.197/0.353 | 0.600 | 0.000 | 0.375 | 0.750 |
| pi0 propio (gen pi pi0) | 1234679 | 4.72 | 0.80/18.09 | 0.433 | 0.083 | 0.040/0.194 | 0.776 | 0.193 | 0.130 | 0.258 |

Masa reco pi+gamma (solo RecoTauType==1) para pi0 legitimo: med 0.697 GeV (N=104804); para extra en pi nu: med 0.442 GeV

### Cortes sencillos sobre el foton: fraccion que SOBREVIVE en cada categoria (extra en pi nu vs pi0 legitimo en pi pi0)

| corte | FSR tau (shower) | ISR (haz) | sin match gen | pi0 propio (gen pi pi0) | migracion 0->1 residual (pi nu puro) |
|---|---|---|---|---|---|
| P_g > 1 GeV | 0.575 | 0.206 | 0.338 | 0.870 | 0.0226 |
| P_g > 2 GeV | 0.436 | 0.073 | 0.197 | 0.742 | 0.0153 |
| dR < 0.2 | 0.536 | 0.233 | 0.745 | 0.908 | 0.0307 |
| dR < 0.15 | 0.370 | 0.139 | 0.682 | 0.815 | 0.0244 |
| dR > 0.02 | 0.994 | 0.997 | 0.903 | 0.988 | 0.0485 |
| P_g/P_pi > 0.05 | 0.614 | 0.304 | 0.369 | 0.896 | 0.0247 |
| P_g/P_pi > 0.1 | 0.481 | 0.173 | 0.245 | 0.800 | 0.0179 |
| P_g>1 y dR<0.2 | 0.306 | 0.068 | 0.306 | 0.814 | 0.0140 |
| P_g>1 y dR<0.2 y dR>0.02 | 0.302 | 0.066 | 0.260 | 0.804 | 0.0130 |
| P_g/P_pi>0.05 y dR<0.2 | 0.318 | 0.085 | 0.293 | 0.827 | 0.0141 |

(migracion 0->1 sin cortes, pi nu puro: 0.0508; la columna 'pi0 propio' es la fraccion de fotones legitimos de pi0 conservados, evaluada sobre gen pi pi0 con reco 1 o 2)

## 4. Dependencia con P: el foton extra escala con el pion?

| categoria | bin P_tau | N | P_gamma med | P_g/P_pi med | dR med | corr(P_g,P_pi) |
|---|---|---|---|---|---|---|
| FSR tau (shower) | 5-15 | 3066 | 1.92 | 0.197 | 0.188 | -0.05 |
| FSR tau (shower) | 20-30 | 2442 | 1.47 | 0.059 | 0.183 | -0.03 |
| FSR tau (shower) | 35-50 | 1464 | 0.81 | 0.020 | 0.168 | -0.29 |
| ISR (haz) | 5-15 | 147 | 0.38 | 0.040 | 0.299 | 0.00 |
| ISR (haz) | 20-30 | 115 | 0.36 | 0.015 | 0.252 | -0.03 |
| ISR (haz) | 35-50 | 118 | 0.39 | 0.009 | 0.308 | -0.11 |
| sin match gen | 5-15 | 1708 | 0.45 | 0.047 | 0.117 | -0.17 |
| sin match gen | 20-30 | 2171 | 0.64 | 0.026 | 0.051 | -0.66 |
| sin match gen | 35-50 | 2069 | 0.75 | 0.019 | 0.038 | -0.75 |

### Balance de momento: P_pion/P_vis_gen y (P_pion+P_gamma)/P_vis_gen (mediana y p10/p90)

| categoria | bin P_tau | N | P_pi/P_gen med (p10/p90) | (P_pi+P_g)/P_gen med (p10/p90) | frac P_pi/P_gen<0.8 |
|---|---|---|---|---|---|
| FSR tau (shower) | 5-15 | 3066 | 1.000 (0.99/1.01) | 1.194 (1.02/3.05) | 0.002 |
| FSR tau (shower) | 20-30 | 2442 | 1.000 (0.99/1.01) | 1.060 (1.01/1.46) | 0.001 |
| FSR tau (shower) | 35-50 | 1464 | 1.000 (0.99/1.01) | 1.022 (1.00/1.12) | 0.002 |
| ISR (haz) | 5-15 | 147 | 1.001 (0.99/1.01) | 1.042 (1.01/1.17) | 0.000 |
| ISR (haz) | 20-30 | 115 | 1.000 (0.99/1.01) | 1.017 (1.00/1.07) | 0.000 |
| ISR (haz) | 35-50 | 118 | 1.000 (0.99/1.01) | 1.011 (1.00/1.04) | 0.000 |
| sin match gen | 5-15 | 1708 | 0.999 (0.97/1.01) | 1.040 (1.01/1.20) | 0.086 |
| sin match gen | 20-30 | 2171 | 0.999 (0.62/1.01) | 1.019 (0.98/1.11) | 0.117 |
| sin match gen | 35-50 | 2069 | 0.999 (0.88/1.01) | 1.015 (0.99/1.06) | 0.097 |
| pi0 propio (gen pi pi0) | 5-15 | 212632 | 0.485 (0.14/0.83) | 0.810 (0.45/0.98) | 0.869 |
| pi0 propio (gen pi pi0) | 20-30 | 347390 | 0.502 (0.13/0.86) | 0.831 (0.46/0.99) | 0.831 |
| pi0 propio (gen pi pi0) | 35-50 | 330886 | 0.503 (0.10/0.90) | 0.856 (0.44/0.99) | 0.758 |

## 5. RecoTauType==2 y fotones gen legitimos en gen tau->pi nu


### GenTauType==0, RecoTauType==2 (dos fotones, cada foton una fila)

| N | FSR tau (shower) | ISR (haz) | sin match gen | otro (eta...) |
|---|---|---|---|---|
| 3703 | 0.152 | 0.008 | 0.763 | 0.077 |

Combinaciones de origen en reco2 (N taus = 1854):
| combinacion | N | frac |
|---|---|---|
| sin match gen + sin match gen | 1277 | 0.689 |
| FSR tau (shower) + sin match gen | 241 | 0.130 |
| FSR tau (shower) + FSR tau (shower) | 154 | 0.083 |
| otro (eta...) + otro (eta...) | 137 | 0.074 |
| FSR tau (shower) + ISR (haz) | 15 | 0.008 |
| ISR (haz) + sin match gen | 14 | 0.008 |
| otro (eta...) + sin match gen | 10 | 0.005 |
| sin match gen | 5 | 0.003 |

### reco2, fotones pi0 de otro tau (comprobacion)

| N |  |
|---|
| 0 |  |

### Fotones gen fisicos en el cono del gen tau->pi nu (fichero completo, todos los reco)

N gen taus type 0: 480898; TrueMode 10: 431789

| clase | N taus | frac de gen type0 |
|---|---|---|
| con FSR de desintegracion (origen 1, TauKey propio) - cualquier P | 0 | 0.00000 |
| con FSR de shower del tau en dR<0.4, P>0.5 GeV | 11583 | 0.02409 |
|   ... y ese foton reconstruido y asignado al reco tau | 9410 | 0.01957 |
| con ISR (parent e) en dR<0.4, P>0.5 GeV | 545 | 0.00113 |
|   ... y reconstruido y asignado al reco tau | 285 | 0.00059 |
| con radiacion del pion (origen 2, parent 211) en dR<0.4, P>0.5 | 0 | 0.00000 |
| con FSR shower P>0.5 en cono pero reco tipo 0 (foton perdido/no asignado) | 296 | 0.00062 |
| con FSR shower P>0.5 en cono y reco tipo 1 | 7950 | 0.01653 |
| reco tipo 1 (total) | 19948 | 0.04148 |
| reco tipo 1 sin ningun foton gen P>0.5 en el cono | 11763 | 0.02446 |

Eficiencia de que el FSR de shower (en cono, P>0.5) acabe asignado al reco tau, por P del foton gen:

| P_gamma gen (GeV) | N | frac asignada al reco tau | frac reco tau tipo 1 |
|---|---|---|---|
| 0.5-1 | 2377 | 0.791 | 0.664 |
| 1-2 | 2212 | 0.821 | 0.687 |
| 2-5 | 2964 | 0.820 | 0.692 |
| 5-10 | 1902 | 0.823 | 0.700 |
| 10-60 | 2128 | 0.808 | 0.691 |

Spectro de P del FSR de shower (gen, en cono): percentiles 10/50/90 = 0.70/2.85/15.27 GeV

Gen tau->pi nu (TrueMode 10) etiquetados GenTauType==1 por un foton FSR duro de la desintegracion: 0 taus con reco 1/2 en las filas (compara con 697170 gen type1 reco1/2).
