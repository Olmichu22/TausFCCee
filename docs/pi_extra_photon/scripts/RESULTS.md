# Origen del foton extra en gen tau->pi nu reconstruido como pi+gamma

Muestra ztt_2M, 1999994 eventos procesados (fichero completo). Tau reco: cono dR<0.4, P(foton)>0.5 GeV (genminP).

## 1-2. Clasificacion del foton asignado (gen tau->pi nu, RecoTauType==1)

Categorias: `FSR tau (shower)` = GenPhotonOrigin 1 con GenPhotonTauKey -1 (foton emitido por una copia Pythia del tau antes de la desintegracion; se asigna al tau por angulo, dR<0.4); `FSR desint. (propio)` = origen 1 y TauKey == tau propio (foton hijo directo del tau en la desintegracion, tau->pi nu gamma); `ISR (haz)` = origen 2 con parent |PDG|==11; `rad. pion` = origen 2 con parent 211; `pi0 otro tau` = origen 0 de otro tau; `sin match gen` = RecoPhotonGenMatchIdx==-1 (el PFO no enlaza a ningun foton del generador: fragmento de shower del pion, cluster partido o link al propio pion).


### Global, GenTauType==0, RecoTauType==1

| N | FSR tau (shower) | ISR (haz) | sin match gen | pi0 otro tau | otro (eta...) | otro |
|---|---|---|---|---|---|---|
| 18763 | 0.540 | 0.032 | 0.427 | 0.000 | 0.001 | 0.000 |

### GenTauType==0, RecoTauType==1, pi nu puro (TrueMode 10)

| N | FSR tau (shower) | ISR (haz) | sin match gen | pi0 otro tau | otro |
|---|---|---|---|---|---|
| 16961 | 0.560 | 0.033 | 0.406 | 0.000 | 0.000 |

### GenTauType==0, RecoTauType==1, K nu (TrueMode 20)

| N | FSR tau (shower) | ISR (haz) | sin match gen | pi0 otro tau |
|---|---|---|---|---|
| 1168 | 0.531 | 0.033 | 0.435 | 0.001 |

### GenTauType==0, RecoTauType==1, K0 pi (TrueMode 23)

| N | FSR tau (shower) | ISR (haz) | sin match gen |
|---|---|---|---|
| 505 | 0.012 | 0.002 | 0.986 |

Otros TrueMode en type0/reco1: [(np.int32(27), 111), (np.int32(41), 11), (np.int32(28), 7)]

### GenTauType==0, RecoTauType==1, resto de TrueMode

| N | FSR tau (shower) | sin match gen | otro (eta...) |
|---|---|---|---|
| 129 | 0.008 | 0.907 | 0.085 |

### TrueMode 10, RecoTauType==1, por GenVisTauP

| GenVisTauP (GeV) | N | FSR tau (shower) | ISR (haz) | sin match gen | pi0 otro tau | otro |
|---|---|---|---|---|---|---|
| 0-5 | 1643 | 0.721 | 0.043 | 0.235 | 0.001 | 0.001 |
| 5-10 | 2135 | 0.641 | 0.030 | 0.329 | 0.000 | 0.000 |
| 10-15 | 2215 | 0.619 | 0.033 | 0.347 | 0.000 | 0.000 |
| 15-20 | 2214 | 0.578 | 0.035 | 0.386 | 0.001 | 0.000 |
| 20-25 | 2030 | 0.527 | 0.035 | 0.437 | 0.000 | 0.000 |
| 25-30 | 1907 | 0.539 | 0.033 | 0.428 | 0.000 | 0.000 |
| 30-35 | 1753 | 0.500 | 0.028 | 0.472 | 0.000 | 0.000 |
| 35-40 | 1650 | 0.490 | 0.036 | 0.474 | 0.000 | 0.000 |
| 40-45 | 1320 | 0.373 | 0.029 | 0.598 | 0.000 | 0.000 |
| 45-50 | 94 | 0.191 | 0.021 | 0.787 | 0.000 | 0.000 |

### TrueMode 10, RecoTauType==1, por |cos theta| del tau

| |cos theta| | N | FSR tau (shower) | ISR (haz) | sin match gen | pi0 otro tau | otro |
|---|---|---|---|---|---|---|
| 0-0.2 | 2981 | 0.647 | 0.019 | 0.333 | 0.000 | 0.000 |
| 0.2-0.4 | 3098 | 0.613 | 0.022 | 0.365 | 0.000 | 0.000 |
| 0.4-0.6 | 3528 | 0.569 | 0.027 | 0.403 | 0.001 | 0.000 |
| 0.6-0.7 | 1837 | 0.508 | 0.027 | 0.465 | 0.000 | 0.000 |
| 0.7-0.75 | 943 | 0.526 | 0.034 | 0.440 | 0.000 | 0.000 |
| 0.75-0.8 | 976 | 0.517 | 0.048 | 0.434 | 0.000 | 0.000 |
| 0.8-0.85 | 1010 | 0.524 | 0.052 | 0.424 | 0.000 | 0.000 |
| 0.85-0.9 | 1045 | 0.509 | 0.062 | 0.428 | 0.001 | 0.000 |
| 0.9-0.95 | 901 | 0.489 | 0.071 | 0.438 | 0.001 | 0.000 |
| 0.95-1 | 642 | 0.354 | 0.059 | 0.586 | 0.002 | 0.000 |

### Fraccion de migracion 0->1 (gen pi nu puro) por |cos theta| del tau (respecto a los reconstruidos con tipo 0 o 1)

| bin cos | N(reco0+reco1) | reco1/(reco0+reco1) | reco1 sin match gen | reco1 FSR shower | reco1 ISR |
|---|---|---|---|---|---|
| 0-0.2 | 54754 | 0.0545 | 0.0182 | 0.0352 | 0.0010 |
| 0.2-0.4 | 58410 | 0.0531 | 0.0194 | 0.0325 | 0.0012 |
| 0.4-0.6 | 64155 | 0.0550 | 0.0222 | 0.0313 | 0.0015 |
| 0.6-0.7 | 31309 | 0.0587 | 0.0273 | 0.0298 | 0.0016 |
| 0.7-0.75 | 18560 | 0.0510 | 0.0224 | 0.0267 | 0.0017 |
| 0.75-0.8 | 20251 | 0.0483 | 0.0209 | 0.0249 | 0.0023 |
| 0.8-0.85 | 21730 | 0.0465 | 0.0197 | 0.0243 | 0.0024 |
| 0.85-0.9 | 23081 | 0.0453 | 0.0194 | 0.0230 | 0.0028 |
| 0.9-0.95 | 23851 | 0.0379 | 0.0166 | 0.0185 | 0.0027 |
| 0.95-1 | 19028 | 0.0338 | 0.0198 | 0.0119 | 0.0020 |

### Fraccion de migracion 0->1 por GenVisTauP (gen pi nu puro)

| P (GeV) | N(reco0+reco1) | reco1/(reco0+reco1) | sin match gen | FSR shower | ISR | pi0 otro tau |
|---|---|---|---|---|---|---|
| 0-5 | 39354 | 0.0418 | 0.0098 | 0.0301 | 0.0018 | 0.0001 |
| 5-10 | 42885 | 0.0499 | 0.0164 | 0.0319 | 0.0015 | 0.0000 |
| 10-15 | 41752 | 0.0531 | 0.0184 | 0.0329 | 0.0017 | 0.0000 |
| 15-20 | 39451 | 0.0562 | 0.0216 | 0.0324 | 0.0020 | 0.0001 |
| 20-25 | 37634 | 0.0540 | 0.0236 | 0.0284 | 0.0019 | 0.0000 |
| 25-30 | 35503 | 0.0537 | 0.0230 | 0.0289 | 0.0018 | 0.0000 |
| 30-35 | 33944 | 0.0517 | 0.0244 | 0.0258 | 0.0014 | 0.0000 |
| 35-40 | 32015 | 0.0516 | 0.0244 | 0.0252 | 0.0019 | 0.0000 |
| 40-45 | 29666 | 0.0445 | 0.0266 | 0.0166 | 0.0013 | 0.0000 |
| 45-50 | 2925 | 0.0321 | 0.0253 | 0.0062 | 0.0007 | 0.0000 |

## 3. Cinematica del foton extra por categoria, comparada con fotones legitimos de pi0 (gen pi pi0, reco 1 o 2)

| categoria | N | P_gamma med (GeV) | P_gamma p10/p90 | P_g/P_pi med | dR(g,pi) med | dR p10/p90 | m(pi+g) med (GeV) | frac dR<0.05 | frac P_g<1 GeV | frac P_g<2 GeV |
|---|---|---|---|---|---|---|---|---|---|---|
| FSR tau (shower) | 10124 | 1.46 | 0.20/12.95 | 0.093 | 0.184 | 0.067/0.346 | 0.723 | 0.056 | 0.417 | 0.566 |
| ISR (haz) | 608 | 0.37 | 0.13/1.34 | 0.022 | 0.282 | 0.132/0.380 | 0.558 | 0.015 | 0.844 | 0.938 |
| sin match gen | 8012 | 0.57 | 0.15/6.07 | 0.029 | 0.070 | 0.021/0.314 | 0.263 | 0.379 | 0.667 | 0.801 |
| pi0 otro tau | 7 | 0.40 | 0.16/4.21 | 0.037 | 0.350 | 0.138/0.367 | 0.549 | 0.000 | 0.714 | 0.714 |
| pi0 propio (gen pi pi0) | 1182518 | 4.70 | 0.79/17.99 | 0.423 | 0.084 | 0.040/0.192 | 0.781 | 0.188 | 0.130 | 0.260 |

Masa reco pi+gamma (solo RecoTauType==1) para pi0 legitimo: med 0.690 GeV (N=105976); para extra en pi nu: med 0.447 GeV

### Cortes sencillos sobre el foton: fraccion que SOBREVIVE en cada categoria (extra en pi nu vs pi0 legitimo en pi pi0)

| corte | FSR tau (shower) | ISR (haz) | sin match gen | pi0 propio (gen pi pi0) | migracion 0->1 residual (pi nu puro) |
|---|---|---|---|---|---|
| P_g > 1 GeV | 0.583 | 0.156 | 0.333 | 0.870 | 0.0225 |
| P_g > 2 GeV | 0.434 | 0.062 | 0.199 | 0.740 | 0.0151 |
| dR < 0.2 | 0.548 | 0.235 | 0.749 | 0.911 | 0.0310 |
| dR < 0.15 | 0.391 | 0.118 | 0.685 | 0.817 | 0.0249 |
| dR > 0.02 | 0.990 | 0.997 | 0.913 | 0.984 | 0.0485 |
| P_g/P_pi > 0.05 | 0.617 | 0.260 | 0.360 | 0.895 | 0.0243 |
| P_g/P_pi > 0.1 | 0.486 | 0.127 | 0.232 | 0.797 | 0.0175 |
| P_g>1 y dR<0.2 | 0.315 | 0.044 | 0.300 | 0.816 | 0.0140 |
| P_g>1 y dR<0.2 y dR>0.02 | 0.309 | 0.044 | 0.262 | 0.802 | 0.0131 |
| P_g/P_pi>0.05 y dR<0.2 | 0.327 | 0.051 | 0.292 | 0.828 | 0.0142 |

(migracion 0->1 sin cortes, pi nu puro: 0.0506; la columna 'pi0 propio' es la fraccion de fotones legitimos de pi0 conservados, evaluada sobre gen pi pi0 con reco 1 o 2)

## 4. Dependencia con P: el foton extra escala con el pion?

| categoria | bin P_tau | N | P_gamma med | P_g/P_pi med | dR med | corr(P_g,P_pi) |
|---|---|---|---|---|---|---|
| FSR tau (shower) | 5-15 | 2920 | 1.88 | 0.201 | 0.185 | -0.01 |
| FSR tau (shower) | 20-30 | 2260 | 1.44 | 0.058 | 0.170 | -0.10 |
| FSR tau (shower) | 35-50 | 1405 | 0.85 | 0.021 | 0.170 | -0.28 |
| ISR (haz) | 5-15 | 150 | 0.36 | 0.037 | 0.296 | 0.05 |
| ISR (haz) | 20-30 | 142 | 0.42 | 0.017 | 0.280 | -0.14 |
| ISR (haz) | 35-50 | 107 | 0.37 | 0.009 | 0.271 | -0.04 |
| sin match gen | 5-15 | 1652 | 0.44 | 0.052 | 0.120 | -0.07 |
| sin match gen | 20-30 | 2014 | 0.63 | 0.025 | 0.055 | -0.58 |
| sin match gen | 35-50 | 1942 | 0.71 | 0.018 | 0.041 | -0.56 |

### Balance de momento: P_pion/P_vis_gen y (P_pion+P_gamma)/P_vis_gen (mediana y p10/p90)

| categoria | bin P_tau | N | P_pi/P_gen med (p10/p90) | (P_pi+P_g)/P_gen med (p10/p90) | frac P_pi/P_gen<0.8 |
|---|---|---|---|---|---|
| FSR tau (shower) | 5-15 | 2920 | 1.000 (0.99/1.01) | 1.203 (1.02/2.87) | 0.003 |
| FSR tau (shower) | 20-30 | 2260 | 1.000 (0.99/1.01) | 1.059 (1.01/1.50) | 0.000 |
| FSR tau (shower) | 35-50 | 1405 | 1.000 (0.99/1.01) | 1.024 (1.00/1.12) | 0.003 |
| ISR (haz) | 5-15 | 150 | 1.000 (0.99/1.01) | 1.037 (1.01/1.14) | 0.007 |
| ISR (haz) | 20-30 | 142 | 1.000 (0.99/1.01) | 1.018 (1.00/1.06) | 0.000 |
| ISR (haz) | 35-50 | 107 | 1.001 (0.99/1.01) | 1.012 (1.00/1.03) | 0.000 |
| sin match gen | 5-15 | 1652 | 1.000 (0.97/1.01) | 1.041 (1.01/1.22) | 0.087 |
| sin match gen | 20-30 | 2014 | 0.999 (0.62/1.01) | 1.020 (0.99/1.10) | 0.115 |
| sin match gen | 35-50 | 1942 | 0.999 (0.96/1.01) | 1.014 (0.99/1.07) | 0.091 |
| pi0 propio (gen pi pi0) | 5-15 | 199528 | 0.496 (0.16/0.83) | 0.813 (0.46/0.98) | 0.864 |
| pi0 propio (gen pi pi0) | 20-30 | 335339 | 0.505 (0.13/0.87) | 0.832 (0.46/0.99) | 0.830 |
| pi0 propio (gen pi pi0) | 35-50 | 316075 | 0.508 (0.10/0.91) | 0.858 (0.44/0.99) | 0.755 |

## 5. RecoTauType==2 y fotones gen legitimos en gen tau->pi nu


### GenTauType==0, RecoTauType==2 (dos fotones, cada foton una fila)

| N | FSR tau (shower) | ISR (haz) | sin match gen | pi0 otro tau | otro (eta...) |
|---|---|---|---|---|---|
| 3331 | 0.170 | 0.008 | 0.728 | 0.003 | 0.091 |

Combinaciones de origen en reco2 (N taus = 1668):
| combinacion | N | frac |
|---|---|---|
| sin match gen + sin match gen | 1088 | 0.652 |
| FSR tau (shower) + sin match gen | 221 | 0.132 |
| FSR tau (shower) + FSR tau (shower) | 166 | 0.100 |
| otro (eta...) + otro (eta...) | 147 | 0.088 |
| ISR (haz) + sin match gen | 12 | 0.007 |
| FSR tau (shower) + ISR (haz) | 11 | 0.007 |
| otro (eta...) + sin match gen | 9 | 0.005 |
| sin match gen | 5 | 0.003 |

### reco2, fotones pi0 de otro tau (comprobacion)

| N | pi0 otro tau |
|---|---|
| 11 | 1.000 |

### Fotones gen fisicos en el cono del gen tau->pi nu (fichero completo, todos los reco)

N gen taus type 0: 479425; TrueMode 10: 430112

| clase | N taus | frac de gen type0 |
|---|---|---|
| con FSR de desintegracion (origen 1, TauKey propio) - cualquier P | 0 | 0.00000 |
| con FSR de shower del tau en dR<0.4, P>0.5 GeV | 11595 | 0.02419 |
|   ... y ese foton reconstruido y asignado al reco tau | 8944 | 0.01866 |
| con ISR (parent e) en dR<0.4, P>0.5 GeV | 539 | 0.00112 |
|   ... y reconstruido y asignado al reco tau | 262 | 0.00055 |
| con radiacion del pion (origen 2, parent 211) en dR<0.4, P>0.5 | 0 | 0.00000 |
| con FSR shower P>0.5 en cono pero reco tipo 0 (foton perdido/no asignado) | 313 | 0.00065 |
| con FSR shower P>0.5 en cono y reco tipo 1 | 7600 | 0.01585 |
| reco tipo 1 (total) | 18779 | 0.03917 |
| reco tipo 1 sin ningun foton gen P>0.5 en el cono | 10949 | 0.02284 |

Eficiencia de que el FSR de shower (en cono, P>0.5) acabe asignado al reco tau, por P del foton gen:

| P_gamma gen (GeV) | N | frac asignada al reco tau | frac reco tau tipo 1 |
|---|---|---|---|
| 0.5-1 | 2341 | 0.755 | 0.626 |
| 1-2 | 2376 | 0.794 | 0.673 |
| 2-5 | 2880 | 0.777 | 0.662 |
| 5-10 | 1823 | 0.770 | 0.654 |
| 10-60 | 2175 | 0.758 | 0.661 |

Spectro de P del FSR de shower (gen, en cono): percentiles 10/50/90 = 0.70/2.78/15.45 GeV

Gen tau->pi nu (TrueMode 10) etiquetados GenTauType==1 por un foton FSR duro de la desintegracion: 0 taus con reco 1/2 en las filas (compara con 672351 gen type1 reco1/2).
