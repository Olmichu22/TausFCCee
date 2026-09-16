# Migración de τ→π π0 ν (GenTauType 1): qué le pasa a cada fotón del π0

Muestra: 30 ficheros `ZTauTau_PolSM_March24_2M_4/out_reco_edm4hep_edm4hep_1000-1029.root`
(30 000 eventos). Reconstrucción idéntica al árbol (`findAllGenTaus` / `findAllTaus`
dR<0.4 sin cortes en P, `MatchRecoGenTau` dR<1). Cada fotón gen del π0 (hijos status 1)
se sigue en sentido MC → reco con `MCTruthRecoLink` (existe; from=PFO, to=MC, peso = fracción
de la energía/hits del MC que cae en ese PFO, codificado LCIO: track en las unidades bajas,
cluster en las altas). Si el fotón no tiene link directo se miran sus descendientes de
simulación (conversiones e+e-), escalando el peso por E(descendiente)/E(γ). El PFO "principal"
es el de mayor peso efectivo (mínimo 0.1). `RecoMCTruthLink` (reco → MC) se usa para la pureza
del PFO, para la parte del cluster del pión que viene del fotón, y para clasificar los PFO
fotón de los reco 3+. Scripts: `raw_links.py` (→ `records.json`), `summary.py` (→ `summary.log`, figuras).

## Estadística

- gen tau tipo 1: 15 926. Excluidos: 205 con neutro extra (K0…), 201 con π0 que no va a 2γ (Dalitz).
- Sin reco tau a dR<1: 1 454 (9.4 %). Reco con otro ID: 2 488 (16 %): 1 737 son el misID π→n
  de Pandora (ID −20), 570 electrón (−11), 142 muón (−13), 31 tres prongs/otros.
- Analizados (líder π±): **reco 0: 130, reco 1: 1 610, reco 2: 8 445, reco 3+: 1 393** (11 578).
  Del total emparejado con líder pión: 13.9 % migran a 1γ, 12.0 % a 3+γ, 1.1 % a 0γ.

## Categoría de cada fotón (2 por tau), % por RecoTauType


|  | N | a_lost | b_in_cone | c_out_cone | d_charged | e_neutralh | f_merged | g_split | other |
|---|---|---|---|---|---|---|---|---|---|
| 0 | 260 | 27.7 | 0.0 | 54.6 | 5.4 | 9.2 | 2.3 | 0.8 | 0.0 |
| 1 | 3220 | 16.9 | 34.5 | 7.0 | 7.0 | 4.3 | 29.9 | 0.4 | 0.0 |
| 2 | 16890 | 0.3 | 97.6 | 0.1 | 0.5 | 0.2 | 0.6 | 0.6 | 0.0 |
| 3+ | 2786 | 0.3 | 71.1 | 0.0 | 1.1 | 0.5 | 0.6 | 26.3 | 0.0 |


Por bin de P visible gen del tau (GeV):

**RecoTauType 0**


|  | N | a_lost | b_in_cone | c_out_cone | d_charged | e_neutralh | f_merged | g_split | other |
|---|---|---|---|---|---|---|---|---|---|
| 0-10 | 74 | 27.0 | 0.0 | 58.1 | 4.1 | 9.5 | 0.0 | 1.4 | 0.0 |
| 10-20 | 118 | 28.0 | 0.0 | 58.5 | 2.5 | 10.2 | 0.0 | 0.8 | 0.0 |
| 20-30 | 48 | 29.2 | 0.0 | 41.7 | 12.5 | 4.2 | 12.5 | 0.0 | 0.0 |
| 30-46 | 20 | 25.0 | 0.0 | 50.0 | 10.0 | 15.0 | 0.0 | 0.0 | 0.0 |

**RecoTauType 1**


|  | N | a_lost | b_in_cone | c_out_cone | d_charged | e_neutralh | f_merged | g_split | other |
|---|---|---|---|---|---|---|---|---|---|
| 0-10 | 258 | 26.4 | 48.4 | 14.7 | 3.5 | 5.8 | 0.8 | 0.4 | 0.0 |
| 10-20 | 854 | 23.1 | 45.2 | 11.9 | 6.8 | 4.8 | 7.5 | 0.7 | 0.0 |
| 20-30 | 870 | 15.7 | 33.3 | 6.0 | 7.5 | 4.5 | 32.4 | 0.6 | 0.0 |
| 30-46 | 1238 | 11.4 | 25.0 | 2.7 | 7.6 | 3.4 | 49.8 | 0.2 | 0.0 |

**RecoTauType 2**


|  | N | a_lost | b_in_cone | c_out_cone | d_charged | e_neutralh | f_merged | g_split | other |
|---|---|---|---|---|---|---|---|---|---|
| 0-10 | 648 | 0.8 | 96.1 | 1.4 | 0.3 | 0.2 | 0.0 | 1.2 | 0.0 |
| 10-20 | 4614 | 0.4 | 97.7 | 0.1 | 0.3 | 0.3 | 0.2 | 0.9 | 0.0 |
| 20-30 | 4790 | 0.3 | 97.7 | 0.1 | 0.7 | 0.3 | 0.5 | 0.5 | 0.0 |
| 30-46 | 6834 | 0.2 | 97.7 | 0.0 | 0.5 | 0.1 | 1.0 | 0.5 | 0.0 |

**RecoTauType 3+**


|  | N | a_lost | b_in_cone | c_out_cone | d_charged | e_neutralh | f_merged | g_split | other |
|---|---|---|---|---|---|---|---|---|---|
| 0-10 | 60 | 0.0 | 90.0 | 0.0 | 0.0 | 0.0 | 0.0 | 10.0 | 0.0 |
| 10-20 | 610 | 0.3 | 75.2 | 0.0 | 0.8 | 1.1 | 0.3 | 22.1 | 0.0 |
| 20-30 | 808 | 0.5 | 70.3 | 0.1 | 1.1 | 0.1 | 0.5 | 27.4 | 0.0 |
| 30-46 | 1308 | 0.2 | 68.9 | 0.0 | 1.2 | 0.4 | 0.9 | 28.4 | 0.0 |

Patrón por tau (par de categorías de los dos fotones), % por RecoTauType (solo los >1 %):

| patrón | reco 0 | reco 1 | reco 2 | reco 3+ |
|---|---|---|---|---|
| a+b (uno perdido, otro contado) | 0.0 | **33.3** | 0.4 | 0.4 |
| f+f (los dos fusionados en un PFO) | 2.3 | **29.9** | 0.6 | 0.6 |
| b+d (uno absorbido por un PFO cargado) | 0.0 | **13.6** | 0.9 | 1.5 |
| b+c (uno fuera del cono) | 0.0 | **13.2** | 0.2 | 0.0 |
| b+e (uno como hadrón neutro) | 0.0 | **8.3** | 0.4 | 0.6 |
| b+g (uno partido en ≥2 PFO fotón) | 0.0 | 0.7 | 0.8 | **47.9** |
| b+b | 0.0 | 0.0 | **96.3** | 45.9 |
| c+c | **37.7** | 0.2 | 0.0 | 0.0 |
| a+c | **23.1** | 0.1 | 0.0 | 0.0 |
| a+a | **10.0** | 0.1 | 0.0 | 0.0 |
| a+d / c+d / c+e / e+e | 6.9 / 3.8 / 6.2 / 3.1 | | | |
| g+g | 0.0 | 0.0 | 0.0 | 1.8 |

## Subcategorías (N)


| subcategoria | 0 | 1 | 2 | 3+ |
|---|---|---|---|---|
| a_lost:converted_no_link | 29 | 170 | 13 | 2 |
| a_lost:converted_weak_links | 5 | 6 | 0 | 0 |
| a_lost:no_link | 38 | 366 | 39 | 7 |
| a_lost:weak_links | 0 | 1 | 0 | 0 |
| b_in_cone:direct | 0 | 1068 | 15907 | 1828 |
| b_in_cone:via_conversion | 0 | 43 | 580 | 154 |
| c_out_cone:dR>0.4 | 113 | 154 | 12 | 0 |
| c_out_cone:dR>0.4_conv | 20 | 67 | 10 | 0 |
| c_out_cone:in_cone_not_const | 0 | 3 | 2 | 0 |
| c_out_cone:other_tau | 7 | 1 | 0 | 1 |
| c_out_cone:other_tau_conv | 2 | 0 | 0 | 0 |
| d_charged:lead_pion | 4 | 100 | 19 | 5 |
| d_charged:lead_pion_conv | 3 | 4 | 4 | 1 |
| d_charged:pdg11_conv | 6 | 119 | 59 | 24 |
| d_charged:pdg211_conv | 1 | 3 | 0 | 0 |
| e_neutralh:out_cone | 9 | 2 | 0 | 0 |
| e_neutralh:out_cone_conv | 15 | 135 | 40 | 13 |
| f_merged:in_cone | 0 | 964 | 100 | 18 |
| f_merged:out_cone | 6 | 0 | 0 | 0 |
| g_split:2_pfo_photons | 2 | 13 | 98 | 633 |
| g_split:3_pfo_photons | 0 | 1 | 5 | 94 |
| g_split:4_pfo_photons | 0 | 0 | 2 | 6 |


`_conv` = el link llega a través de descendientes de simulación (el fotón gen no tiene link directo);
`converted_` en (a) = el fotón tiene `isDecayedInTracker` o descendientes de sim con endpoint dentro del tracker.

## Conversión en el tracker


| cat | 0 (N) | 1 (N) | 2 (N) | 3+ (N) | todos |
|---|---|---|---|---|---|
| (a) perdido | 44 (72) | 31 (543) | 25 (52) | 22 (9) | 32 (676) |
| (b) PFO γ en cono | - | 4 (1111) | 4 (16487) | 8 (1982) | 4 (19580) |
| (c) PFO γ fuera cono | 15 (142) | 29 (225) | 38 (24) | 0 (1) | 25 (392) |
| (d) PFO cargado | 71 (14) | 56 (226) | 77 (82) | 83 (30) | 64 (352) |
| (e) PFO n | 58 (24) | 99 (137) | 100 (40) | 100 (13) | 94 (214) |
| (f) fusionados | 17 (6) | 4 (964) | 26 (100) | 33 (18) | 7 (1088) |
| (g) partido | 100 (2) | 100 (14) | 98 (105) | 98 (733) | 98 (854) |

Taus con >=1 foton del pi0 convertido en el tracker, % por RecoTauType:
- reco 0: 50.0% (N=130)
- reco 1: 35.7% (N=1610)
- reco 2: 9.6% (N=8445)
- reco 3+: 61.1% (N=1393)
- radio de conversion (mm): p10 35, mediana 436, p90 1607
- radio de conversion (d) PFO cargado pdg11_conv: mediana 15 mm, <30 mm: 54%
- radio de conversion (e) PFO n: mediana 341 mm, <30 mm: 5%
- radio de conversion (g) partido: mediana 521 mm, <30 mm: 4%


## Propiedades por categoría


| reco | N | P mediana (GeV) | P<0.2 GeV (%) | P<0.5 GeV (%) | |cos θ|>0.98 (%) | convertidos (%) | ang γ-π mediana (rad) |
|---|---|---|---|---|---|---|---|
| 0 | 72 | 0.17 | 57 | 78 | 32 | 47 | 0.570 |
| 1 | 543 | 0.09 | 77 | 93 | 5 | 32 | 0.321 |
| 2 | 52 | 0.05 | 79 | 94 | 4 | 25 | 0.395 |
| 3+ | 9 | 0.10 | 78 | 89 | 0 | 22 | 0.274 |

### P gen y angulo γ-π por categoria (todos los RecoTauType)

| cat | N | P mediana (GeV) | P p10 | P p90 | ang γ-π mediana (rad) | ang γ-γ mediana (rad) |
|---|---|---|---|---|---|---|
| (a) perdido | 676 | 0.09 | 0.02 | 0.46 | 0.341 | 0.350 |
| (b) PFO γ en cono | 19580 | 4.36 | 0.78 | 15.63 | 0.084 | 0.042 |
| (c) PFO γ fuera cono | 392 | 0.75 | 0.16 | 5.41 | 0.467 | 0.287 |
| (d) PFO cargado | 352 | 3.04 | 0.36 | 10.10 | 0.083 | 0.057 |
| (e) PFO n | 214 | 1.27 | 0.57 | 2.58 | 0.115 | 0.084 |
| (f) fusionados | 1088 | 11.78 | 2.10 | 29.44 | 0.063 | 0.017 |
| (g) partido | 854 | 8.01 | 2.23 | 18.50 | 0.076 | 0.029 |


- pdg11_conv: 208
- lead_pion: 128
- lead_pion_conv: 12
- pdg211_conv: 4
- pion absorbe un foton: N=140, E/p mediana 1.15 (media 1.14); restando la parte del foton (RecoMCTruthLink): 0.99 (media 0.99)
- E/p>1: 77% absorbiendo vs 55% sin absorber
- fraccion del cluster del pion que viene del foton: mediana 0.07
- P gen del foton absorbido: mediana 1.28 GeV; angulo γ-π mediana 0.054 rad, 46% a <0.05 rad
- pion sin foton absorbido (control): N=6016, E/p mediana 1.02 (media 1.01)

### Pureza del PFO principal (fraccion del PFO que viene del foton, RecoMCTruthLink)

- (b) PFO γ en cono: N=19580, mediana 0.98, <0.5: 4%
- (f) fusionados: N=1088, mediana 0.42, <0.5: 54%
- (g) partido: N=854, mediana 0.00, <0.5: 98%
- (c) PFO γ fuera cono: N=392, mediana 1.00, <0.5: 25%
- (b) E(PFO)/E(gen): mediana 1.016
- (f) E(PFO)/P(pi0 gen): mediana 1.016

- (c) dR PFO-lider: mediana 0.53 rad, >1 rad: 10%
- (c) angulo gen γ-π: mediana 0.47 rad

### Conversion en el tracker (isDecayedInTracker): % de fotones por categoria y RecoTauType

| cat | 0 (N) | 1 (N) | 2 (N) | 3+ (N) | todos |
|---|---|---|---|---|---|
| (a) perdido | 44 (72) | 31 (543) | 25 (52) | 22 (9) | 32 (676) |
| (b) PFO γ en cono | - | 4 (1111) | 4 (16487) | 8 (1982) | 4 (19580) |
| (c) PFO γ fuera cono | 15 (142) | 29 (225) | 38 (24) | 0 (1) | 25 (392) |
| (d) PFO cargado | 71 (14) | 56 (226) | 77 (82) | 83 (30) | 64 (352) |
| (e) PFO n | 58 (24) | 99 (137) | 100 (40) | 100 (13) | 94 (214) |
| (f) fusionados | 17 (6) | 4 (964) | 26 (100) | 33 (18) | 7 (1088) |
| (g) partido | 100 (2) | 100 (14) | 98 (105) | 98 (733) | 98 (854) |

Taus con >=1 foton del pi0 convertido en el tracker, % por RecoTauType:
- reco 0: 50.0% (N=130)
- reco 1: 35.7% (N=1610)
- reco 2: 9.6% (N=8445)
- reco 3+: 61.1% (N=1393)
- radio de conversion (mm): p10 35, mediana 436, p90 1607
- radio de conversion (d) PFO cargado pdg11_conv: mediana 15 mm, <30 mm: 54%
- radio de conversion (e) PFO n: mediana 341 mm, <30 mm: 5%
- radio de conversion (g) partido: mediana 521 mm, <30 mm: 4%


## Reco 3+: qué son los fotones extra


N taus 3+: 1393, PFO foton: 4571, extra (no principal de ningun foton del pi0): 1875

| clase | todos | extra | % extra | E mediana extra (GeV) | dR lider mediana |
|---|---|---|---|---|---|
| pi0_photon_conv_frag | 1035 | 1028 | 54.8 | 1.89 | 0.122 |
| FSR_tau | 422 | 422 | 22.5 | 2.19 | 0.162 |
| shower_pion | 323 | 320 | 17.1 | 0.62 | 0.146 |
| photon_gen_other | 58 | 55 | 2.9 | 6.02 | 0.084 |
| pi0_photon_split | 29 | 28 | 1.5 | 0.31 | 0.073 |
| ISR | 18 | 18 | 1.0 | 0.46 | 0.288 |
| FSR_other_tau | 2 | 2 | 0.1 | 0.76 | 0.082 |
| other_same_tau | 1 | 1 | 0.1 | 1.24 | 0.244 |
| rad_pion_gen | 1 | 1 | 0.1 | 0.16 | 0.160 |
| pi0_photon_main | 1835 | 0 | 0 | | |
| pi0_photon_conv_main | 847 | 0 | 0 | | |

- fragmentos de conversion: E(extra)/E(gamma gen) mediana 0.25; el foton gen padre tiene isDecayedInTracker en 100%
- FSR extra: E mediana 2.19 GeV, E>1 GeV 69%

Patron de extras por tau 3+:
- pi0_photon_conv_frag: 530 (38.0%)
- FSR_tau: 319 (22.9%)
- pi0_photon_conv_frag+pi0_photon_conv_frag: 159 (11.4%)
- shower_pion: 144 (10.3%)
- FSR_tau+pi0_photon_conv_frag: 33 (2.4%)
- pi0_photon_conv_frag+pi0_photon_conv_frag+pi0_photon_conv_frag: 24 (1.7%)
- pi0_photon_split: 24 (1.7%)
- photon_gen_other+photon_gen_other: 22 (1.6%)
- shower_pion+shower_pion: 20 (1.4%)
- pi0_photon_conv_frag+shower_pion: 16 (1.1%)
- ISR: 13 (0.9%)
- FSR_tau+FSR_tau: 13 (0.9%)


`pi0_photon_conv_frag` = PFO fotón cuyo mejor link es un e±/γ de simulación descendiente de un fotón del π0
que se convirtió en el tracker (no es el PFO principal de ese fotón). `pi0_photon_split` = PFO extra
enlazado directamente al fotón gen (cluster partido sin conversión).

## Figuras

- `fig_photon_P_by_cat.png`: P gen del fotón por categoría (izq.) y del fotón perdido por RecoTauType (der.).
- `fig_EoverP_pion.png`: E(cluster)/p del pión líder cuando absorbe un fotón del π0 frente a cuando no, y restando la parte del fotón.
- `fig_angle_by_cat.png`: ángulo gen γ–π y γ–γ por categoría.

## Conclusiones

1. **Reco 1 (π + 1γ, 13.9 % de los emparejados) tiene tres causas del mismo tamaño y una cuarta menor.**
   - 33 % un fotón **perdido** (a+b): fotón muy blando, P mediana 0.09 GeV, 77 % por debajo de 0.2 GeV y 93 % de 0.5 GeV;
     sólo el 5 % está fuera de aceptancia (|cos θ|>0.98) y un 31 % convirtió en el tracker. Domina a P visible bajo
     (26 % de los fotones a P_vis<10 GeV frente a 11 % a 30-46 GeV).
   - 30 % los dos fotones **fusionados** en un solo PFO fotón (f+f): π0 duro (P mediana del fotón 11.8 GeV, ángulo γγ
     mediano 0.017 rad); el PFO lleva toda la energía del π0 (E(PFO)/P(π0) = 1.016). Es la causa dominante a P visible alto:
     50 % de los fotones a P_vis 30-46 GeV, 7.5 % a 10-20 GeV. Sólo el 4 % convirtió: es puramente resolución angular del ECAL.
   - 13.6 % un fotón **absorbido por un PFO cargado** (b+d): la mitad por el propio pión líder (100+4 casos, fotón a 0.054 rad
     del pión, P mediana 1.3 GeV), la otra mitad por un electrón de conversión (119 casos: el fotón convirtió en el beam pipe /
     vértice, radio mediano 15 mm, y el e± se reconstruye como PFO electrón que la reconstrucción ignora).
   - 13.2 % un fotón **fuera del cono** (b+c): fotón de P mediana 0.75 GeV a 0.47 rad del pión; el PFO está a dR mediano 0.53 del
     líder. Sólo relevante a P visible bajo (15 % a <10 GeV, 3 % a 30-46 GeV).
   - 8.3 % un fotón reconstruido como **hadrón neutro** (b+e): 99 % son conversiones en el tracker (radio mediano 340 mm) cuyo
     e+e- acaba en un PFO 2112 que además está fuera del cono (135/137 casos).
2. **Reco 0 (130 taus, 1.1 %) son taus blandos con el cono demasiado estrecho:** 38 % los dos fotones fuera del cono, 23 % uno
   fuera y otro perdido, 10 % los dos perdidos; el fotón perdido aquí sí es de aceptancia en un 32 %. No hay caso b.
3. **Reco 3+ (12 %) es, en primer lugar, conversión:** el 48 % de los taus tiene un fotón "partido" en ≥2 PFO fotón y el 98 %
   de esos fotones partidos convirtió en el tracker (radio mediano 520 mm). El 61 % de los taus 3+ tiene al menos un fotón
   convertido, frente al 9.6 % de los reco 2. De los 1 875 PFO fotón extra: 55 % fragmentos de conversión (E/E_γ mediana 0.25,
   a 0.12 rad del líder), 22.5 % FSR del tau (E mediana 2.2 GeV, 69 % >1 GeV), 17 % fragmentos del shower del pión
   (E mediana 0.6 GeV), 3 % otros fotones gen, 1.5 % cluster partido sin conversión, 1 % ISR.
4. **El pión que absorbe un fotón se ve en E/p:** E/p mediano 1.15 (77 % > 1) frente a 1.02 (55 % > 1) en el control; restando
   la parte del cluster que el `RecoMCTruthLink` atribuye al fotón vuelve a 0.99. La fracción absorbida es pequeña
   (mediana 7 % del cluster) y el fotón está a <0.05 rad en el 46 % de los casos. Es un efecto del 6 % de la migración a reco 1.
5. **La transición reco 2 → reco 1 con la energía es la fusión, no la pérdida:** la fracción de fotones perdidos en reco 1 baja
   de 26 % a 11 % entre P_vis<10 y 30-46 GeV mientras la fusión sube de 1 % a 50 %. Recuperar el reco 1 a P alto pasa por
   reconocer el PFO fotón único con la masa/energía del π0 (E(PFO)/P(π0)≈1), no por buscar un segundo fotón.
