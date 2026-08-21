# Plan de implementación: convenio $z=\cos\theta_{\tau^-}$

Estado: **fases 1–3 implementadas** (11 ago 2026). Pendiente: fases 0 y 4, que son
ejecuciones sobre datos (validación previa, regeneración de árboles y test `--sin-eff`).

Notas de implementación:

- `tau_pdg` es **keyword-only obligatorio** en todas las funciones de peso: un call
  site olvidado da `TypeError` (§3.3).
- `RhoHistFromTree_MDecs_parallel.py` **aborta** si el árbol de entrada no trae
  `tau{1,2}_tauPDG` / `tau{1,2}_recoCharge` (`_require_charge_branches`): la lectura
  de ramas cae a 0.0 por defecto, lo que marcaría todo hemisferio como τ⁺ en
  silencio. Los árboles anteriores hay que regenerarlos.
- Hemispheres **with no matched gen tau** carry the sentinels `tauPDG = -999` and
  `decayID = -999` instead of 0. This is not an edge case: with fewer than 2 gen
  taus in the event `_select_decay_modes` returns `gen_idx = -1` for both
  hemispheres, which is every event of any sample declared `has_gen_taus: false`
  in `config/samples/samples.yaml` (today: `bhabha`). It also covers ττ events
  where fewer than 2 taus pass the `getGeneratorStatus() == 2` filter of
  `findAllGenTaus`.
  - `tauPDG`: `_gen_tau_pdg` **raises `ValueError`** on anything other than ±15,
    and its five call sites are guarded with `_has_gen_tau` — gen weights are left
    as the producer wrote them, `optimalVar` = −999, and the joint weight falls
    back to the per-tau product. A 0 there read silently as τ⁺, the very failure
    mode this convention removes.
  - `decayID`: a 0 means "gen pion", so these hemispheres were classified as
    SIGNAL in the pion channel by `_classify_hemisphere`. With −999 they fall to
    `BGOther`. Side effect: the pion `vism_cut` (which only applies to
    `decayID == 0`) no longer touches them — a no-op in practice, since the gen
    `visM` of such a hemisphere is 0.0 and `0.0 < vism_cut` passed them anyway.
  - `weight_P1` / `weight_M1`: **1.0**, not the generic 0.0 of the other gen
    branches. They are multiplicative, so the neutral element is 1: with no gen
    tau there is nothing to reweight and the event must stay as it is, not
    disappear. This is not cosmetic — `_compute_joint_weights` falls back to the
    product of the two per-tau weights when either hemisphere lacks a gen tau, so
    a 0.0 here silently empties the `corr_P1`/`corr_M1` variants of `ALL_ALL` for
    any sample made of such events (bhabha). The old code hid this: it called
    `newAtauJoint` with a null 4-vector, which returns exactly 1.0 (`Theta()` of a
    null `TLorentzVector` is 0, `cosThetaStar` clamps to 0, `sumHH = 0` →
    numerator = denominator), so `corr` kept weight 1 while `P1`/`M1` were already
    being zeroed. Nothing in the fit reads these (`makeCosBins_MDecs.py` uses the
    `reco_*` variants over `SIGNAL_SIGNAL`, and the reco weights are computed
    normally for these hemispheres from `recoCharge`), but a template that should
    equal the nominal must not come out empty.
- El stem de salida lleva la marca `zTaum_`, **después** de los ids del par para no
  romper `_parse_pair_ids` de `makeCosBins_MDecs.py` (§4.6).
- La validación 4.2 (invarianza del joint) se ha comprobado analíticamente sobre
  taus back-to-back sintéticos: `newAtauJoint(p4(dec0),H,Hp,·,tau_pdg=15)` ==
  `newAtauJoint(p4(dec1),Hp,H,·,tau_pdg=-15)` con diferencia 0.

Fecha original del plan: 11 ago 2026.

Documento autocontenido — puede implementarse sin el contexto de la conversación
que lo originó.

## Referencias

- `../ExamplesFCCFullSim/docs/signo_costheta_carga_tau.md` — diagnóstico completo,
  derivación y test experimental del 10 ago 2026. **Leer antes de implementar.**
- `../ExamplesFCCFullSim/docs/reweighting_tau_polarization.md` — fórmulas del repesado.
- J. Alcaraz, *Reweight_tautau.pdf*, mayo 2026, eqs. (1)–(3).

---

## 1. El problema en una línea

Los observables de spin ($\omega$, $H_V$, $H_\ell$) y el término cruzado $h_1h_2$
están **bien**. Lo que está mal es el signo de $z=\cos\theta$: debe ser siempre el
del $\tau^-$ y en el código es el del hemisferio analizado.

$$P(z)=-\frac{A_\tau(1+z^2)+2A_e z}{(1+z^2)+2A_eA_\tau z}$$

El término $A_\tau(1+z^2)$ es par (inmune); $2A_ez$ y $2A_eA_\tau z$ son impares.
Como los taus son back-to-back, $\cos\theta_{\tau^+}=-\cos\theta_{\tau^-}$, así que
$z\to-z \iff A_e\to-A_e$. Usar el hemisferio en vez del $\tau^-$ equivale a invertir
$A_e$ en ~50% de los eventos.

**Consecuencia:** $A_\tau$ sobrevive; $A_e$ y $\sin^2\theta_\text{eff}$ **no son
actualmente medidas independientes de los datos** — el fit devuelve la asunción de
entrada ($P_{fit}=\tilde P$, §5 del doc de diagnóstico; verificado 1:1 en §7).

## 2. Qué NO hay que tocar

Crítico para no romper nada por exceso de celo:

- **Los observables no llevan signo de carga.** $\omega$ (`wVariab` usa
  `genTauP4.Vect()` como eje, natural para cada carga), $z_R$ (fracción de energía,
  ciega a la carga), $H_\ell$. Ni `optimalVar` ni `omega` cambian.
- La corrección por helicidad del antineutrino ya está absorbida: para el $\tau^+$,
  $d\Gamma\propto 1-h^+\omega^+ = 1+h^-\omega^+$. Ambos hemisferios se rigen por la
  **misma** $P(z)$ del $\tau^-$.
- El binado del YAML: ya es `[-1,1] × 20`. Cambia el significado de los bins, no el
  rango.
- `makeCosBins_MDecs.py`, `fitPolAssym.py`: sin cambios.
- Las 4 funciones `_depc` de `weightsPol.py` (L298, 344, 428, 517).
- `modules/rhoTreeUtils.py:222-236`: tercer camino de pesos, pero `get_entry_vars`
  solo se usa en L296 del propio módulo y nada externo la importa → **código muerto**.
  Si se revive, necesita `genTauPDG`, que no está en `SCALAR_BRANCHES_WEIGHTS`.
- `RhoAnalysis/genHelicityHistos.py`: **ya es correcto** (L189-194 identifica `tau_m`
  por `tauPDG == 15` y pasa `tauMinusP4`). Es la implementación de referencia. Solo
  necesita el `tau_pdg=15` explícito si se hace el parámetro obligatorio (§4.1).

## 3. Decisiones ya cerradas

| # | Cuestión | Resolución |
|---|---|---|
| 1 | ¿El eje $+Z$ coincide con la dirección del $e^-$ incidente? | **Sí.** Todas las cards de `../SimProdScripts/CLDFCC_sim/` usan `Beams:idA = 11` / `Beams:idB = -11`, sin `frameType`/`pzA`/`pzB`. Pythia8 manda el haz A por $+z$. **No hace falta signo global.** |
| 2 | Fallback si $\lvert q_\text{reco}\rvert \neq 1$ | **Ignorar el caso** (decisión del usuario, 11 ago 2026): en principio siempre hay carga. Ver §4.3 para el tratamiento mínimo. |
| 3 | ¿`tau_pdg` obligatorio o con default? | **Obligatorio**, sin valor por defecto (decisión del usuario, 11 ago 2026). Un call site olvidado debe dar `TypeError`, no calcular mal en silencio — que es exactamente el modo de fallo que estamos arreglando. Implica migrar los ~20 call sites de golpe, incluido el `tau_pdg=15` explícito en `genHelicityHistos.py`. |

**Nota de segundo orden (no bloquea):** `CLDConfig/Sim/cld_steer.py:28` aplica
`SIM.crossingAngleBoost = 0.015` (15 mrad), activo vía `run_sequence_CLD.py:196`.
Inclina el eje de haces ~0.86° respecto a $+z$. No afecta al signo; queda como
sistemático conocido sobre la definición de $z$.

---

## 4. Fases

### Fase 0 — Validación previa (sin tocar código, sin regenerar)

Sobre árboles gen-only existentes (ya tienen `tauPDG`), comparar `Omega_plus` vs
`Omega_minus` **dentro de bins de $\cos\theta$**.

- Deben **diferir claramente ahora**.
- Integrados en $\cos\theta$ coinciden ya (el promedio de carga de $P(z)$ es
  $\approx-A_\tau$, plano; lo que sobrevive es $\mathcal{O}(A^3)\sim3\times10^{-3}$).
  **Por eso el bug pasó desapercibido — esa comparación no vale como test.**

Confirma el diagnóstico en estos datos antes de invertir en la regeneración.

Aprovechar para **contar cuántos τ reco tienen $\lvert q\rvert\neq1$** (informativo,
ver §4.3).

### Fase 1 — Núcleo (`modules/`)

**`modules/weightsPol.py`** — helper nuevo:

```python
def _z_taum(p4, tau_pdg):
    """cos(theta) del tau-. tau_pdg: 15 (tau-) o -15 (tau+)."""
    return (1.0 if int(tau_pdg) == 15 else -1.0) * math.cos(p4.Theta())
```

Parámetro `tau_pdg` (obligatorio, §3.3) en las 6 funciones vivas:

| Función | Línea del `costheta` |
|---|---|
| `newAtau` | 117 |
| `newAtauLep` | 141 |
| `newAtauFromH` | 168 |
| `newAtauJoint` | 196 |
| `newAtauJoint_had_had` | 223 |
| `newAtauJoint_had_lep` | 254 |

Más `newAtauRhoOmega` (alias de `newAtauFromH`, L176-178): propaga el argumento.

**`modules/optimalVariabRho.py`** — `wVariab` tiene una **copia inline** de la
fórmula de $P(z)$ (L66-80), no llama a `_compute_Ptau`. Es un sitio independiente:

- Añadir `tau_pdg` y aplicar el signo en L74 (`costheta_tau`).
- Sustituir la copia inline (L76-80) por llamadas a `weightsPol._compute_Ptau` para
  que no vuelvan a divergir.
- `optimal_var` (L178) solo usa el elemento `[3]` ($\omega$, sin pesos): propagar el
  argumento, el observable no cambia.

### Fase 2 — Productores de árbol (ramas nuevas + regeneración)

La carga **sí se reconstruye**, con redundancia; simplemente se descarta al llenar:

| Objeto | Fuente | Dónde |
|---|---|---|
| τ hadrónico reco | `getCharge()`, y `setPDG(±15)` ya derivado de la carga | `tauReco.py:640-643` |
| e/μ reco | `pf.getCharge()`, `PDGID=pf.getPDG()` | `electronReco.py:29-31`, `muonReco.py:30-32` |
| π cargado líder | PDG firmado del constituyente (el `abs()` es solo del test) | `_extract_pion_p4` |
| τ gen | `setPDG(±15)` | `tauReco.py:544-547` |

**`RhoAnalysis/analysisRHOTree_MDecs_parallel.py`** (reco completo):

1. `_build_reco_candidates` (L119-155): añadir `"charge": tau.getCharge()` en las 3
   ramas (τ hadrónico, μ, e). **Aquí es donde se pierde hoy.**
2. `_TAU_SCALAR_SUFFIXES` (L57-69): añadir `tauPDG` y `recoCharge`.
3. `_fill_tau_branches_mdecs`: llenar `recoCharge` desde el candidato y `tauPDG`
   desde `gen_tau_obj.getPDG()` (ya es ±15).
4. Pesos reco (L292-300): pasar el signo desde `recoCharge`.
   Pesos gen (L364-380): desde `tauPDG`.
5. `recoLepPDG` (L262) se deja como está — hoy vale `float(abs(reco_id))`, que ni es
   un PDG ni tiene signo, pero tiene consumidores. `recoCharge` es la fuente única.

**`RhoAnalysis/genOnlyRHOTree_MDecs_parallel.py`**:

1. `recoCharge` como espejo del gen (`-1 if tauPDG==15 else +1`). `tauPDG` ya existe
   (L142).
2. Pesos (L205-237): pasar `tau_pdg`.

**Regeneración de árboles.** Es el grueso del coste. Los gen-only se pueden validar
sin regenerar; los de reco completo no.

### Fase 3 — Stage de histogramas y fit

**`RhoAnalysis/RhoHistFromTree_MDecs_parallel.py`**:

| Sitio | Líneas | Fuente del signo |
|---|---|---|
| `_recompute_weights` | 350-363 | `tau_vars["tauPDG"]` |
| `_recompute_reco_weights` | 425-435 | `recoCharge` |
| Joint gen | 493, 497, 504 | `tauPDG` |
| Joint reco | 644, 649, 655 | `recoCharge` |
| Eje del fit `OptimalReco_vs_CosThetaVis` | 171 | `recoCharge` |
| Eje del fit `OptimalNN_vs_CosThetaVis` | 178 | `recoCharge` |

Los dos ejes del fit pasan a:

```python
lambda v, sh: _reco_tau_sign(v) * math.cos(v["recoVisTheta"])
```

**Bug adicional que esto arregla de paso.** `_recompute_joint_weights` llama a
`newAtauJoint` con `p4(vars_dec0)` o `p4(vars_dec1)` según la rama (L493 vs L497) y
escribe el **mismo** `w` en ambos hemisferios (L509-510). Sin signo, esos dos
$\cos\theta$ son opuestos: el peso `corr` de un par ρ-leptón depende del orden en que
vinieran los taus en el árbol. Con el signo aplicado por hemisferio ambas ramas
convergen al mismo valor. La validación 4.2 lo verifica.

### Fase 4 — Validación

1. **Back-to-back** (gen): $\lvert\cos\theta_{\tau^-}+\cos\theta_{\tau^+}\rvert<\epsilon$.
2. **Invarianza del joint**: `newAtauJoint(p4(dec0),...)` == `newAtauJoint(p4(dec1),...)`
   evento a evento.
3. **`Omega_plus` vs `Omega_minus` por bin de $\cos\theta$**: deben coincidir tras
   plegar por signo (contraste con Fase 0, donde difieren).
4. **Test `--sin-eff` — criterio de éxito principal.** Repetir §7 del doc de
   diagnóstico con 0.2312 y 0.2300, misma muestra:
   - **Antes:** el fit sigue la asunción 1:1 (0.231223 → 0.230024).
   - **Después:** debe devolver ≈0.2312 (el valor de generación) en **ambos** casos.
   - Si sigue siguiendo la asunción, **el arreglo no ha funcionado**.
5. **$A_\tau$ se mueve poco** (parte par, ya estaba a salvo). Si se mueve mucho, hay
   un bug en la propagación del signo.
6. Versionar el stem de salida para no pisar `Binned_histograms_MDecs` existentes.
   **Los resultados previos dejan de ser comparables.**

---

## 4.3 Tratamiento de la carga anómala

Por decisión (§3.2) se asume $\lvert q\rvert=1$ siempre. Tratamiento mínimo, sin
rama de fallback:

```python
def _reco_tau_sign(v):
    """+1 si el hemisferio es tau- (carga<0), -1 si es tau+."""
    return 1.0 if v.get("recoCharge", -1.0) < 0 else -1.0
```

Añadir un **contador en el log** de casos con $\lvert q\rvert\neq1$ (posible cuando
`charge_condition=False`; `buildTauFromPion` los acepta, `tauReco.py:448`). No
cambia el comportamiento, pero evita que el caso pase inadvertido si algún día
aparece con estadística no despreciable.

---

## 5. Resumen de ficheros

| Fichero | Cambio |
|---|---|
| `modules/weightsPol.py` | Helper `_z_taum` + `tau_pdg` en 6 funciones + alias |
| `modules/optimalVariabRho.py` | `tau_pdg` en `wVariab` (L74) + deduplicar $P(z)$ |
| `RhoAnalysis/analysisRHOTree_MDecs_parallel.py` | Propagar carga + 2 ramas nuevas + 2 call sites de pesos |
| `RhoAnalysis/genOnlyRHOTree_MDecs_parallel.py` | Rama `recoCharge` + call sites de pesos |
| `RhoAnalysis/RhoHistFromTree_MDecs_parallel.py` | 4 grupos de call sites + 2 ejes del fit |
| `RhoAnalysis/genHelicityHistos.py` | Solo `tau_pdg=15` explícito (ya correcto) |

Sin cambios: `config/histograms/rho_analysis_config_mdecs.yml`,
`makeCosBins_MDecs.py`, `fitPolAssym.py`.
