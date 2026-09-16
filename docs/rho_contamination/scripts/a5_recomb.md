## Clase reco 3g+: mejor pareja de fotones y resto, por tipo gen
| gen | N | frac |m_par-m_pi0|<0.05 | <0.03 | P_resto/P_par p10/50/90 | frac P_resto/P_par<0.2 | <0.3 | <0.5 | m(all)<0.25 | m(all)<0.3 |
|---|---|---|---|---|---|---|---|---|---|
| 1 rho | 93733 | 0.779 | 0.674 | 0.03/0.37/3.28 | 0.365 | 0.455 | 0.570 | 0.247 | 0.309 |
| 2 pi2pi0 | 249084 | 0.975 | 0.947 | 0.16/0.90/4.21 | 0.129 | 0.200 | 0.323 | 0.036 | 0.077 |
| 3 pi3pi0 | 30846 | 0.982 | 0.965 | 0.60/2.07/7.21 | 0.012 | 0.027 | 0.072 | 0.003 | 0.007 |
| 0 pi | 1099 | 0.781 | 0.663 | 0.14/0.83/5.24 | 0.152 | 0.226 | 0.356 | 0.338 | 0.402 |
| e | 5256 | 0.678 | 0.475 | 0.10/0.53/2.80 | 0.211 | 0.307 | 0.474 | 0.264 | 0.360 |

## Reclasificar 3g+ -> 2g con cada criterio (fichero completo)
Ganancia = gen1 3g+ que pasa a 2g / todos los gen1 ; contaminacion = gen2+gen3+otros 3g+ que pasan a 2g, relativo al tamano actual de la clase 2g (589866)
| criterio | gen1 recuperados | Δeff(gen1->2g) | gen2 -> 2g | gen3 -> 2g | otros -> 2g | contaminacion añadida / clase 2g | pureza gen1 de la clase 2g resultante |
|---|---|---|---|---|---|---|---|
| baseline | 0 | 0 | 0 | 0 | 0 | 0 | 0.953 |
| A: m(all)<0.25 | 23180 | +0.022 | 9003 | 104 | 1955 | 0.019 | 0.938 |
| A: m(all)<0.30 | 28971 | +0.028 | 19084 | 226 | 2577 | 0.037 | 0.922 |
| A: m(all)<0.40 | 37505 | +0.036 | 57141 | 754 | 3648 | 0.104 | 0.870 |
| B: par<0.05 & Prest/Ppar<0.1 | 19462 | +0.019 | 12814 | 84 | 549 | 0.023 | 0.934 |
| B: par<0.05 & Prest/Ppar<0.2 | 28287 | +0.027 | 30898 | 305 | 1071 | 0.055 | 0.907 |
| B: par<0.05 & Prest/Ppar<0.3 | 33997 | +0.032 | 48009 | 732 | 1522 | 0.085 | 0.884 |
| B: par<0.05 & Prest/Ppar<0.5 | 41799 | +0.040 | 77927 | 2055 | 2346 | 0.140 | 0.846 |
| B': par<0.05 & Prest_max/Ppar<0.2 | 29062 | +0.028 | 41408 | 1202 | 1223 | 0.074 | 0.892 |
| C: par<0.05 & m(all)<0.40 | 35267 | +0.034 | 56173 | 702 | 2787 | 0.101 | 0.872 |
| D: m(all)<0.3 | (par<0.05 & Prest/Ppar<0.3) | 45024 | +0.043 | 59496 | 898 | 3216 | 0.108 | 0.869 |

## Calidad del pi0 en gen1 3g+ (criterio D) : P_est/P_gen(pi0)
- suma de todos los fotones: p10/50/90 = 0.945/1.043/1.227 ; frac |r-1|<0.1: 0.668
- solo mejor pareja: p10/50/90 = 0.560/0.958/1.070 ; frac |r-1|<0.1: 0.556
- referencia 2g: p10/50/90 = 0.936/1.018/1.101 ; frac |r-1|<0.1: 0.847

## 1g: P_reco/P_gen del foton de pi0 conservado (¿absorbe al perdido?)
- 1g: p10/50/90 = 0.917/1.065/1.447 ; frac r>1.15: 0.316 ; frac r>1.3: 0.162
- 2g: p10/50/90 = 0.886/1.018/1.167 ; frac r>1.15: 0.117 ; frac r>1.3: 0.039
- 1g con dR(gg) gen en [0,0.02): frac 0.272, P_reco/P_gen del conservado p50 1.261, frac >1.15: 0.755, P gen del perdido p50 4.88
- 1g con dR(gg) gen en [0.02,0.05): frac 0.157, P_reco/P_gen del conservado p50 1.054, frac >1.15: 0.220, P gen del perdido p50 1.86
- 1g con dR(gg) gen en [0.05,0.2): frac 0.298, P_reco/P_gen del conservado p50 1.023, frac >1.15: 0.113, P gen del perdido p50 0.44
- 1g con dR(gg) gen en [0.2,2): frac 0.271, P_reco/P_gen del conservado p50 1.016, frac >1.15: 0.172, P gen del perdido p50 0.08
