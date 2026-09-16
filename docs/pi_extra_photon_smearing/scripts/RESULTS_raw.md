# Que es el foton extra en gen tau->pi nu reconstruidos como pi + 1 gamma (EDM4hep original)

Ficheros: out_reco_edm4hep_edm4hep_1000..1024 (25 x 1000 eventos) de Ztt_SM_2M_SigmaZ (ztt_2M_smearing), leidos con podio.
Reconstruccion identica a TTreesTausLong.py / config New2MSample_smearing_results0.4_tph0.0_tpi0.0_n0.0_g0.0:
`findAllGenTaus(mc)`, `findAllTaus(pfos, dRMax=0.4, minP_photon=0, minP_pion=0, PNeutron=0, genminP=0)` (+electrones/muones),
`MatchRecoGenTau(maxDRMatch=1)`. Script: `photon_links.py`; tablas/figuras: `summarize.py`; datos crudos: `records.json`.

## Cifras globales

| eventos | gen tau tipo 0 | matched | reco 0 | reco 1 (pi+1 gamma) | sin match |
|---|---|---|---|---|---|
| 25000 | 5996 | 5754 | 4487 | 222 (3.9 % de los matched) | 242 |

Todos los 218 PFO foton tienen exactamente 1 cluster y 0 tracks, y todos tienen al menos un RecoMCTruthLink (0 casos 'sin link').
Todos los links tienen peso de track = 0 (PFO neutro); el peso relevante es el de cluster (w//10000/1000).

## Clasificacion por el link de mayor peso del PFO foton

| categoria | N | % | <E_gamma> [GeV] | mediana E_gamma | mediana peso cluster | mediana angulo cluster-cluster pion [rad] | mediana dR(PFO gamma, PFO pion) |
|---|---|---|---|---|---|---|---|
| FSR del propio tau | 141 | 63.5 | 3.94 | 1.02 | 1.00 | 0.133 | 0.161 |
| fragmento shower del pion | 63 | 28.4 | 2.22 | 0.42 | 1.00 | 0.031 | 0.075 |
| ISR | 12 | 5.4 | 0.84 | 0.31 | 1.00 | 0.243 | 0.351 |
| K0L del propio tau (tau->K0 pi nu) | 6 | 2.7 | 15.81 | 17.19 | 1.00 | 0.090 | 0.059 |
| FSR del otro tau | 0 | 0.0 | - | - | - | - | - |
| pi0 del otro tau | 0 | 0.0 | - | - | - | - | - |
| radiacion del pion en el generador | 0 | 0.0 | - | - | - | - | - |
| otro (mismo tau) | 0 | 0.0 | - | - | - | - | - |
| sin link | 0 | 0.0 | - | - | - | - | - |

Categorias (segun la MCParticle enlazada, subiendo por parents hasta el generador):
- **FSR del propio tau**: foton status 1 del generador cuyo padre es una copia del propio tau (status 23/51/52, emision QED del shower de Pythia). Incluye los casos en que el link apunta a un e+/e- o gamma de simulacion (conversion/brems) que desciende de ese foton FSR. Estos fotones NO estan entre las hijas del tau status 2, por eso GenTauType=0 (en el tree: GenPhotonOrigin=1, GenPhotonTauKey=-1).
- **fragmento shower del pion**: link directo al pion cargado del tau (status 1; los hits del calorimetro apuntan a la primaria) o a una secundaria de simulacion (generatorStatus 0, isCreatedInSimulation) cuyo ancestro generador es ese pion (gammas de pi0 secundarios, pi+-, protones, electrones del shower hadronico).
- **ISR**: foton status 1 con ancestro e+/e- de haz (status 41/42/4, chain [22,43]/[11,41]).
- **K0L del propio tau**: link al K0L status 1 hijo del K0 del tau (tau->K0 pi nu, GenTauTrueMode 23; extra neutral no contado en GenTauType=0). Pandora reconstruye el cluster del K0L como foton.
- **FSR del otro tau**, **pi0 del otro tau**, **radiacion del pion en el generador** (fotones con padre pi+-/K+- de status 1 en el generador): ver tabla.

### Detalle de la MCParticle del mejor link

| categoria | cat. fina | PDG | genStatus | creada en sim | N |
|---|---|---|---|---|---|
| FSR del propio tau | FSR_tau | 22 | 1 | False | 137 |
| fragmento shower del pion | pion_gen_directo | -211 | 1 | False | 27 |
| fragmento shower del pion | pion_gen_directo | 211 | 1 | False | 21 |
| ISR | ISR | 22 | 1 | False | 12 |
| K0L del propio tau (tau->K0 pi nu) | otra_particula_gen | 130 | 1 | False | 6 |
| fragmento shower del pion | sim_shower_pion | 22 | 0 | True | 4 |
| FSR del propio tau | sim_de_FSR_tau | -11 | 0 | True | 3 |
| fragmento shower del pion | sim_shower_pion | 2212 | 0 | True | 3 |
| fragmento shower del pion | pion_gen_directo | -321 | 1 | False | 3 |
| fragmento shower del pion | sim_shower_pion | 11 | 0 | True | 2 |
| fragmento shower del pion | sim_shower_pion | 211 | 0 | True | 1 |
| fragmento shower del pion | pion_gen_directo | 321 | 1 | False | 1 |
| fragmento shower del pion | sim_shower_pion | 13 | 0 | True | 1 |
| FSR del propio tau | sim_de_FSR_tau | 11 | 0 | True | 1 |

Numero de links por PFO foton: media 1.27, max 7; con 1 solo link: 169 (76 %).
PFO FSR con contaminacion del shower del pion (peso de cluster del pion > 0.05): 17 de 141; mediana de la fraccion de peso que va al pion en esos casos: 0.18.
Todos los links (no solo el mejor), por categoria: {'fragmento shower del pion': 117, 'FSR del propio tau': 148, 'ISR': 12, 'K0L del propio tau (tau->K0 pi nu)': 6}

## Geometria del cluster foton

CLD: ECAL barrel r~2150-2350 mm, |z|<2300; ECAL endcap |z|~2300-2500 (r<2100). Fraccion ECAL = subdetectorEnergies[0]/E.

| categoria | N | barrel | endcap | fraccion ECAL >0.95 | 0.30<theta<0.45 o 2.7<theta<2.85 (transicion) | mediana nhits | mediana E_gamma/E_pion PFO |
|---|---|---|---|---|---|---|---|
| FSR del propio tau | 141 | 100 | 41 | 139 | 7 | 45 | 0.073 |
| fragmento shower del pion | 63 | 34 | 29 | 62 | 4 | 14 | 0.028 |
| ISR | 12 | 6 | 6 | 12 | 2 | 16 | 0.035 |
| K0L del propio tau (tau->K0 pi nu) | 6 | 3 | 3 | 2 | 0 | 506 | 2.436 |

Distancia angular entre la posicion del cluster foton y la del cluster del pion (rad), cuantiles 10/50/90 %:
- fragmento shower: 0.009 / 0.031 / 0.262; dR PFO-PFO: 0.023 / 0.075 / 0.306
- FSR propio tau: 0.060 / 0.133 / 0.320; dR PFO-PFO: 0.068 / 0.161 / 0.340
- ISR: 0.123 / 0.243 / 0.355; dR PFO-PFO: 0.154 / 0.351 / 0.382

Fragmentos con cluster a <0.05 rad del cluster del pion (mismo shower partido): 36/63; con cluster a >0.2 rad: 12/63 (splash lejano del shower hadronico, E mediana 0.18 GeV).
Fragmentos cuyo pion gen interacciono/decayo en el tracker (isDecayedInTracker): 8/63; E(cluster pion)+E(gamma) vs p(pion): mediana (E_pi+E_gamma)/p = 1.03, E_pi/p = 0.98.

E_gamma > 1 GeV: FSR 71/141, fragmento 15/63, ISR 3/12, K0L 6/6.
E_gamma > 5 GeV: FSR 32/141, fragmento 7/63, K0L 6.

## Fotones FSR del generador en TODOS los gen tau tipo 0 (no solo los reco 1)

Fotones status 1 cuyo ancestro tau es el propio tau (copias status 23/51/52 incluidas).

| reco tipo | N gen taus | con >=1 FSR gen (cualquier angulo) | con FSR gen en dR<0.4 | con FSR gen en dR<0.4 y E>0.1 GeV | con FSR en cono y E>1 GeV |
|---|---|---|---|---|---|
| 0 (pi) | 4487 | 707 | 228 | 15 | 2 |
| -20 (pi+neutron) | 824 | 149 | 64 | 27 | 18 |
| sin match | 242 | 47 | 17 | 10 | 6 |
| 1 (pi+1 gamma) | 222 | 157 | 145 | 139 | 74 |
| -11 | 106 | 17 | 6 | 3 | 2 |
| -13 | 52 | 12 | 6 | 1 | 0 |
| 2 (pi+2 gamma) | 25 | 10 | 8 | 8 | 7 |
| -1 | 11 | 1 | 1 | 0 | 0 |
| -21 | 6 | 1 | 1 | 0 | 0 |
| 4 | 5 | 2 | 1 | 0 | 0 |
| 7 | 4 | 0 | 0 | 0 | 0 |
| 3 | 4 | 2 | 2 | 2 | 2 |
| 15 | 3 | 0 | 0 | 0 | 0 |
| 12 | 2 | 0 | 0 | 0 | 0 |
| 10 | 1 | 0 | 0 | 0 | 0 |
| 11 | 1 | 1 | 1 | 1 | 0 |
| 14 | 1 | 0 | 0 | 0 | 0 |

Total de fotones FSR gen en 5996 gen taus tipo 0: 1242 (18.4 % de los taus tienen alguno); E cuantiles 10/50/90 %: 1.05e-04 / 0.025 / 5.54 GeV; fraccion en dR<0.4: 0.41.

Gen tau tipo 0 con FSR en el cono (reco 0 o 1): fraccion reconstruida como 1 en funcion de la E del FSR mas energetico del cono:

| E_FSR [GeV] | N | reco 1 | frac |
|---|---|---|---|
| 0-0.05 | 194 | 3 | 0.02 |
| 0.05-0.1 | 25 | 3 | 0.12 |
| 0.1-0.2 | 22 | 12 | 0.55 |
| 0.2-0.5 | 41 | 39 | 0.95 |
| 0.5-1 | 15 | 14 | 0.93 |
| 1-2 | 24 | 23 | 0.96 |
| 2-5 | 18 | 18 | 1.00 |
| 5-50 | 34 | 33 | 0.97 |

## ISR en el generador (todos los eventos)

Eventos: 25000; fotones ISR (status 1, ancestro e+/e- de haz) por evento: media 2.74, con >=1: 100.0 %.
E ISR cuantiles 10/50/90/99 %: 1.58e-08 / 0.000 / 0.09 / 1.65 GeV; E>1 GeV: 1.9 %, E>5 GeV: 0.2 %.
- ISR con E>0.0 GeV: 68540 fotones, 1.8 % caen a dR<0.4 de algun gen tau (visible); eventos con al menos uno asi: 4.8 %.
- ISR con E>0.1 GeV: 6590 fotones, 2.7 % caen a dR<0.4 de algun gen tau (visible); eventos con al menos uno asi: 0.7 %.
- ISR con E>1.0 GeV: 1331 fotones, 2.0 % caen a dR<0.4 de algun gen tau (visible); eventos con al menos uno asi: 0.1 %.
- ISR con E>5.0 GeV: 159 fotones, 0.6 % caen a dR<0.4 de algun gen tau (visible); eventos con al menos uno asi: 0.0 %.

En los 218 casos reco 1: eventos con ISR en cono de algun tau: 15; mejor link ISR: 12 (E mediana 0.31 GeV).

## Conclusiones

1. El foton extra es, en 141/222 casos (64 %), un foton REAL del generador: FSR emitido por el propio tau en el shower QED de Pythia (padre = copia del tau status 23/51/52). No cuelga del tau status 2, asi que la clasificacion gen (hijas del status 2) lo ignora y da GenTauType=0. La reconstruccion no se equivoca: ve un foton que esta ahi. En el tree se identifica como GenPhotonOrigin==1 con GenPhotonTauKey==-1.
2. 63/222 (28 %) son fragmentos del shower hadronico del pion (link al propio pion o a secundarias de simulacion descendientes de el): cluster ECAL de baja energia, casi siempre a <0.1-0.2 rad del cluster del pion. Son los unicos que en el tree aparecen como RecoPhotonGenMatchIdx==-1 (enlazan al pion, no a un foton gen).
3. ISR: 12/222 (5 %), fotones ISR blandos (mediana 0.31 GeV) que caen en el cono. Solo 1.9 % de los fotones ISR tienen E>1 GeV y 0.2 % E>5 GeV; los de E>5 GeV caen en dR<0.4 de un tau en 0.6 % de los casos (ninguno de los 218 fotones extra es un ISR duro).
4. K0L de tau->K0 pi nu: 6/222; son los fotones 'extra' mas energeticos (E mediana 17.2 GeV) y en realidad el gen tau no es un pi nu puro (GenTauHasExtraNeutrals=1).
5. FSR del otro tau / pi0 del otro tau / radiacion del pion en el generador: 0 / 0 / 0. Sin link: 0.
6. Discriminacion: el FSR reconstruido tiene E mediana 1.02 GeV (50 % con E>1 GeV; el espectro gen es 1/E pero el reco solo lo ve con eficiencia >0.6 a partir de ~0.1-0.2 GeV) y esta separado del pion (angulo cluster-cluster mediana 0.13 rad); el fragmento de shower es blando (mediana 0.42 GeV, 24 % con E>1 GeV) y pegado al cluster del pion (mediana 0.031 rad, 44/74 a <0.05 rad). Un corte en E_gamma (o E_gamma/E_pi) recuperaria la mayoria de los dos, pero hay que medir cuanto pierde en pi pi0 (fotones de pi0 tambien blandos a bajo P del tau). Alternativa mas limpia: en la etiqueta gen, contar como constituyente del tau el FSR del propio tau (GenPhotonOrigin==1 dentro de dR<0.4), con lo que ~55 % de la 'migracion' 0->1 deja de ser migracion.

## Figuras
- fig_link_weight.png: peso de cluster del mejor link por categoria.
- fig_cluster_distance.png: angulo cluster foton - cluster pion y dR PFO-PFO por categoria.
- fig_energy.png: energia del PFO foton por categoria y E_gamma/E_pion.
- fig_fsr_gen.png: espectro de los FSR gen en cono en gen0 y fraccion reconstruida como reco 1.