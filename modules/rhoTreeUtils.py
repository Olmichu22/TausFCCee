import ROOT
import math
from modules import optimalVariabRho
from modules import weightsPol
from modules import myutils


FILL_RULES = [
    # ── Gen ───────────────────────────────────────────────────────────────────
    ("Gen",     "Omega_GEN",            lambda v: v["gen_w"],                          None),
    ("Gen",     "OmegaCosThetaTau_GEN", lambda v: v["gen_w"],                          lambda v: v["gen_cos_theta_tau"]),
    # 2D Omega vs Z-mass (renombrado: Omega_GEN_ZGenMass_SIGNAL en lugar de Omega_GEN_SIGNAL_ZGenMass)
    ("Gen",     "Omega_GEN_ZGenMass",   lambda v: v["gen_w"],                          lambda v: v["ZGenMass"]),
    ("Gen",     "Omega_GEN_ZVisMass",   lambda v: v["gen_w"],                          lambda v: v["ZGenVisMass"]),
    ("Gen",     "CosTheta_GEN",         lambda v: v["gen_cos_theta"],                  None),
    ("Gen",     "CosPsi_GEN",           lambda v: v["gen_cos_psi"],                    None),
    ("Gen",     "CosThetaTau_GEN",      lambda v: v["gen_cos_theta_tau"],              None),
    ("Gen",     "CosThetaRho_GEN",      lambda v: math.cos(v["genMesonTheta"]),        None),
    ("Gen",     "CosThetaRho",          lambda v: math.cos(v["recoMesonTheta"]),       None),
    ("Gen", "MesonPvsLeptonP_GEN", lambda v: v["genMesonP"],                     lambda v: v["genLepP"]),
    ("Gen", "Omega_GEN_Gen_LepP", lambda v: v["gen_w"],                          lambda v: v["genLepP"]),
    # Solo categoría SIGNAL en el YAML (renombrado: GenVisZMass_SIGNAL)
    ("Gen",     "GenVisZMass",          lambda v: v["ZGenVisMass"],                    None),
    # Categorías SIGNAL y BG (renombrado: MesonType_SIGNAL y MesonType_BG)
    ("Gen",     "MesonType",            lambda v: v["genTauID"],                       None),
    ("Gen", "GenZMass", lambda v: v["ZGenMass"], None),
    ("Gen", "GenMuOmegaVar", lambda v: v["genOptimalvarLep"], None, lambda v: abs(v["genLepPDG"]) == 13),
    ("Gen", "GenElOmegaVar",  lambda v: v["genOptimalvarLep"], None, lambda v: abs(v["genLepPDG"]) == 11),
    ("Gen", "GenLepP_omega_max_0", lambda v: v["genLepP"], None, lambda v: v["gen_w"] < 0),
    ("Gen", "GenLepTheta_omega_max_0", lambda v: v["genLepTheta"], None, lambda v: v["gen_w"] < 0),
    ("Gen", "GenLepPhi_omega_max_0", lambda v: v["genLepPhi"], None, lambda v: v["gen_w"] < 0),
    ("Gen", "GenTauP_omega_max_0", lambda v: v["genTauP"], None, lambda v: v["gen_w"] < 0),
    ("Gen", "GenTauTheta_omega_max_0", lambda v: v["genTauTheta"], None, lambda v: v["gen_w"] < 0),
    ("Gen", "GenTauPhi_omega_max_0", lambda v: v["genTauPhi"], None, lambda v: v["gen_w"] < 0),
    ("Gen", "GenWeightP1LepvsGenWeightRhoP1", lambda v: v["weight_P1"], lambda v: v["weight_lep_P1"], lambda v: v["gen_w"] * v["weight_P1"]),
    ("Gen", "GenWeightM1LepvsGenWeightRhoM1", lambda v: v["weight_M1"], lambda v: v["weight_lep_M1"], lambda v: v["gen_w"] * v["weight_M1"]),
    # ── Reco ──────────────────────────────────────────────────────────────────
    ("Reco",    "Omega",                lambda v: v["w"],                              None),
    ("Reco",    "OmegaCosTheta",        lambda v: v["w"],                              lambda v: v["cos_theta_rho"]),
    ("Reco",    "CosTheta",             lambda v: v["cos_theta"],                      None),
    ("Reco",    "CosPsi",               lambda v: v["cos_psi"],                        None),
    ("Reco",    "RecoMesonEOverBeamE",  lambda v: v["recoMesonE"] / v["beamE"],        None),
    ("Reco",    "RecoMesonCosTheta",    lambda v: math.cos(v["recoMesonTheta"]),       None),
    # Solo categoría BGMuon en el YAML
    ("Reco",    "RecoMesonE",           lambda v: v["recoMesonE"],                     None),
    # 2D ángulo-phi y pT-theta (BGMuon y BGEle en el YAML; renombrado: RecoMesonPhiTheta_BGMuon, etc.)
    ("Reco",    "RecoMesonPhiTheta",    lambda v: v["recoMesonTheta"],                 lambda v: v["recoMesonPhi"]),
    ("Reco",    "RecoMesonPtTheta",     lambda v: v["recoMesonTheta"],                 lambda v: v["recoMesonPt"]),
    ("Reco",    "RecoMeson_X",          lambda v: v["Optimal_var_x"],                  None),
    ("Reco",    "RecoMeson_P",          lambda v: v["recoMesonP"],                     None),
    # Solo categoría SIGNAL en el YAML (renombrado: RecoZMass_SIGNAL)
    ("Reco",    "RecoZMass",            lambda v: v["zmass"],                          None),
    # Categorías SIGNAL y BG (renombrado: RecoMesonType_SIGNAL y RecoMesonType_BG)
    ("Reco",    "RecoMesonType",        lambda v: v["recoTauID"],                      None),
    ("Reco",    "DeltaR_LepMeson",      lambda v: v["dR_lep_meson"],                   None),
    ("Reco",    "MesonP",               lambda v: v["recoMesonP"],                     None),
    ("Reco",    "LeptonP",              lambda v: v["leptonP"],                        None),
    ("Reco",    "CosTheta_Tau",         lambda v: v["cos_theta_rho"],                  None),
    # ── Matched ───────────────────────────────────────────────────────────────
    ("Matched", "MesonEOverBeamE",      lambda v: v["genMesonE"] / v["beamE"],         None),
    ("Matched", "MesonCosTheta",        lambda v: math.cos(v["genMesonTheta"]),        None),
]


WEIGHT_VALUES = {
    "nominal":     lambda v, w: w,
    "P1":          lambda v, w: w * v["weight_P1"],
    "M1":          lambda v, w: w * v["weight_M1"],
    "corr_P1":     lambda v, w: w * v["weight_P1"] * v["weight_lep_P1"],
    "corr_M1":     lambda v, w: w * v["weight_M1"] * v["weight_lep_M1"],
    "lep_P1":      lambda v, w: w * v["weight_lep_P1"],
    "lep_M1":      lambda v, w: w * v["weight_lep_M1"],
    "corr_LEP_P1": lambda v, w: w * v["weight_lep_P1"] * v["weight_P1"],
    "corr_LEP_M1": lambda v, w: w * v["weight_lep_M1"] * v["weight_M1"],
}

def fill_category(hists, vars_dict, category: str, base_weight: float):
    for rule in FILL_RULES:
        level, var_name, x_fn, y_fn = rule[:4]
        condition_fn = rule[4] if len(rule) > 4 else None
        if condition_fn is not None and not condition_fn(vars_dict):
            continue
        if var_name not in hists[level]:
            continue
        if category not in hists[level][var_name]:
            continue
        for w_name, w_hist in hists[level][var_name][category].items():
            effective_weight = WEIGHT_VALUES[w_name](vars_dict, base_weight)
            x_val = x_fn(vars_dict)
            if y_fn is None:
                w_hist.Fill(x_val, effective_weight)
            else:
                w_hist.Fill(x_val, y_fn(vars_dict), effective_weight)

SCALAR_BRANCHES = [
    ("weight_P1",     "weight_P1",     float),
    ("weight_M1",     "weight_M1",     float),
    ("weight_lep_P1", "weight_lep_P1", float),
    ("weight_lep_M1", "weight_lep_M1", float),
    ("genTauP",         "genTauP",       float),
    ("genTauTheta",     "genTauTheta",   float),
    ("genTauPhi",       "genTauPhi",     float),
    ("beamE",            "beamE",         float),
    ("recoMesonE",       "recoMesonE",    float),
    ("recoMesonTheta",   "recoMesonTheta",float),
    ("recoMesonPhi",     "recoMesonPhi",  float),
    ("recoMesonP",       "recoMesonP",    float),
    ("genMesonE",        "genMesonE",     float),
    ("genMesonP",        "genMesonP",     float),
    ("genMesonTheta",    "genMesonTheta", float),
    ("cos_theta",        "cos_theta",     float),
    ("cos_psi",          "cos_psi",       float),
    ("cos_theta_rho",    "cos_theta_rho", float),
    ("gen_w",            "genOmega",      float),
    ("w",                "omega",         float),
    ("gen_cos_theta",    "gen_cos_theta", float),
    ("gen_cos_psi",      "gen_cos_psi",   float),
    ("gen_cos_theta_tau","gen_cos_theta_tau", float),
    ("genTauID",         "genTauID",      int),
    ("recoTauID",        "recoTauID",     int),
    ("leptonP",          "lepP",          float),
    ("leptonE",          "lepE",          float),
    ("leptonPhi",        "lepPhi",        float),
    ("leptonTheta",      "lepTheta",      float),
    ("ZGenMass",         "GenZMass",      float),
    ("ZGenVisMass",      "GenZVisMass",   float),
    ("ZRecoMass",        "ZMass",         float),
    ("genLepP",          "genLepP",       float),
    ("genLepE",          "genLepE",       float),
    ("genLepTheta",      "genLepTheta",   float),
    ("genLepPhi",        "genLepPhi",     float),
    ("genLepPDG",        "genLepPDG",     float),
    ("genOptimalvarLep",  "genOptimalvarLep", float),
    ("genOptimalvarPi", "genOptimalvarPi", float),
]

SCALAR_BRANCHES_LEPTAU = [
    ("genLepTauP",     "genLepTauP",     float),
    ("genLepTauE",     "genLepTauE",     float),
    ("genLepTauTheta", "genLepTauTheta", float),
    ("genLepTauPhi",   "genLepTauPhi",   float),
    ("genLepTauM",     "genLepTauM",     float),
]

SCALAR_BRANCHES_WEIGHTS = [
    ("gentauP",    "genTauP",    float),
    ("gentauTheta","genTauTheta",float),
    ("gentauPhi",  "genTauPhi",  float),
    ("gentauE",    "genTauE",    float),
    ("genMesonP",  "genMesonP",  float),
    ("genMesonPhi","genMesonPhi",float),
    ("genMesonTheta","genMesonTheta",float),
    ("genPionP",   "genPionP",   float),
    ("genPionPhi", "genPionPhi", float),
    ("genPionTheta","genPionTheta",float),
    ("genPionE",   "genPionE",   float),
    ("genLepE",          "genLepE",       float),
    
    # ("leptonE",          "lepE",          float),
]

SCALAR_BRANCHES_STORED_WEIGHTS = [
    ("weight_P1",     "weight_P1",     float),
    ("weight_M1",     "weight_M1",     float),
    ("weight_lep_P1", "weight_lep_P1", float),
    ("weight_lep_M1", "weight_lep_M1", float),
]

def fill_special_signal_histograms(special_histograms, vars_dict, weight):
    """Histogramas de señal que requieren lógica de evento más allá de un simple lookup.
    Solo deben ir aquí los que NO pueden expresarse como (variable, categoría, peso).
    Actualmente: bins de Z-mass visible y Pi0Mass (con smearing).
    """
    # Bins de ZVisMass (lógica condicional de qué bin rellenar)
    z_vis = vars_dict["ZGenVisMass"]
    bin_ranges = [(0, 40), (40, 70), (70, 100)]  # ajustar según config
    for idx, (lo, hi) in enumerate(bin_ranges, start=1):
        if lo <= z_vis < hi:
            special_histograms["Gen"]["Events"][f"Omega_GEN_SIGNAL_ZVisMass_Bin{idx}"].Fill(
                vars_dict["gen_w"], weight)
            break
  
def extract_scalars(entry, branch_defs):
    """Extrae ramas escalares de un TTree entry en un dict."""
    return {key: typ(getattr(entry, branch)) for key, branch, typ in branch_defs}


def extract_scalars_optional(entry, branch_defs, default=0.0):
    """Like extract_scalars but silently uses `default` for missing branches (backward compat)."""
    return {key: typ(getattr(entry, branch, default)) for key, branch, typ in branch_defs}


def make_p4(P, theta, phi, E):
    """Construye un TLorentzVector desde coordenadas esféricas."""
    v = ROOT.TLorentzVector()
    v.SetPxPyPzE(
        P * math.sin(theta) * math.cos(phi),
        P * math.sin(theta) * math.sin(phi),
        P * math.cos(theta),
        E,
    )
    return v


def get_entry_vars(entry, proccesing_cfg, logger_process):
    sin_eff = proccesing_cfg.get("sin_eff", None)
    compute_weights = proccesing_cfg.get("compute_weights", False)
    decay_mode = proccesing_cfg.get("decay_mode", None)
    try:
        v = extract_scalars(entry, SCALAR_BRANCHES)

        v["mesonp4"] = make_p4(v["recoMesonP"], v["recoMesonTheta"], v["recoMesonPhi"], v["recoMesonE"])
        v["leptonp4"] = make_p4(v["leptonP"],    v["leptonTheta"],    v["leptonPhi"],    v["leptonE"])

        if compute_weights:
            v.update(extract_scalars(entry, SCALAR_BRANCHES_WEIGHTS))
            v["genTauP4"]   = make_p4(v["gentauP"],   v["gentauTheta"],   v["gentauPhi"],   v["gentauE"])
            v["genMesonP4"] = make_p4(v["genMesonP"],  v["genMesonTheta"], v["genMesonPhi"],  v["genMesonE"])
            v["genPionP4"]  = make_p4(v["genPionP"],   v["genPionTheta"],  v["genPionPhi"],   v["genPionE"])
            # Dispatch on decay type: each decay channel has its own weight formula
            if v["genTauID"] == 1:  # rho: full formula (needs pion in rho RF)
                (_, _, _, _, v["weight_P1"], v["weight_M1"]) = optimalVariabRho.wVariab(
                    v["genTauP4"], v["genMesonP4"], v["genPionP4"], v["beamE"], sin_eff=sin_eff
                )
            elif v["genTauID"] in (0, 10):  # pion (alpha=1) or a1 (alpha=0.12)
                v["weight_P1"] = weightsPol.newAtau(v["genTauP4"], v["genMesonP4"], v["genTauID"], +1, sin_eff = sin_eff)
                v["weight_M1"] = weightsPol.newAtau(v["genTauP4"], v["genMesonP4"], v["genTauID"], -1, sin_eff = sin_eff)
            else:  # lepton backgrounds: no polarization sensitivity
                v["weight_P1"] = 1.0
                v["weight_M1"] = 1.0
            v["genLepP4"] = make_p4(v["genLepP"], v["genLepTheta"], v["genLepPhi"], v["genLepE"])
            v.update(extract_scalars_optional(entry, SCALAR_BRANCHES_LEPTAU, default=0.0))
            v["genLepTauP4"] = make_p4(v["genLepTauP"], v["genLepTauTheta"], v["genLepTauPhi"], v["genLepTauE"])
            if int(v["genLepPDG"]) in [11, 13] and v["genLepTauE"] > 0:
                v["weight_lep_P1"] = weightsPol.newAtauLep(v["genLepP4"], v["genLepTauP4"], v["beamE"], +1, sin_eff=sin_eff)
                v["weight_lep_M1"] = weightsPol.newAtauLep(v["genLepP4"], v["genLepTauP4"], v["beamE"], -1, sin_eff=sin_eff)
            else:
                v["weight_lep_P1"] = 1.0
                v["weight_lep_M1"] = 1.0
        else:
            v.update(extract_scalars_optional(entry, SCALAR_BRANCHES_STORED_WEIGHTS, default=1.0))

    except AttributeError as e:
        logger_process.error("Missing branch in tree '%s'", e)
        return None  # o `continue` si estás en un loop

    return v

def process_tree(rho_vars_extremes_trees,
                root_histograms_super,
                special_histograms_super,
                weight,
                selectGEN,
                cuts_cfg,
                proccesing_cfg,
                logger_process,
                other_BG_id,):
  
  # Counters (only meaningful for the "original" tree)
  totalEvents = 0
  selectedEvents = 0
  sumWeights = 0.0
  sumWeightsP1 = 0.0
  sumWeightsM1 = 0.0

  # Cuts extraction
  tauPCut = cuts_cfg.get("tauPCut", 0)
  meson_cut = cuts_cfg.get("meson_cut", 0)
  lepton_cut = cuts_cfg.get("lepton_cut", 0)
  zmass_cut = cuts_cfg.get("zmass_cut", 0)
  angle_sep = cuts_cfg.get("angle_sep", 0)
  
  # Iterate throught original, max and min extremes (if present)
  for tree_key, tree in rho_vars_extremes_trees.items():
          # tree key case histograms and resolution 
          root_histograms = root_histograms_super[tree_key]
          
          logger_process.info(
              "Refilling histograms for tree key '%s' with %d entries",
              tree_key,
              tree.GetEntries(),
          )          

          n_entries = tree.GetEntries()
          
          # Extract entries from tree
          for i in range(n_entries):
              tree.GetEntry(i)
              entry = tree  # para que sea más corto escribir

              # Contadores globales solo para el árbol "original"
              if tree_key == "original":
                  totalEvents += 1

              # Extraer variables del árbol (mismas ramas que definiste)
              vars_dict = get_entry_vars(entry, proccesing_cfg, logger_process)

              # Variable x definida en el análisis original
              if vars_dict["beamE"] != 0:
                  Optimal_var_x = 2.0 * vars_dict["recoMesonE"] / vars_dict["beamE"] - 1.0
              else:
                  Optimal_var_x = 0.0
              
              vars_dict["Optimal_var_x"] = Optimal_var_x
              
              if vars_dict["recoMesonP"] < tauPCut:
                  continue

              z_p4 = vars_dict["mesonp4"] + vars_dict["leptonp4"]
              zmass = z_p4.M()
              vars_dict["zmass"] = zmass

              if vars_dict["recoMesonP"] < meson_cut[0] or vars_dict["recoMesonP"] > meson_cut[1]:
                  continue # Ignoring this event to get the bk
              if vars_dict["leptonP"] < lepton_cut[0] or vars_dict["leptonP"] > lepton_cut[1]:
                  continue # Ignoring this event to get the bk
              if zmass < zmass_cut[0] or zmass > zmass_cut[1]:
                  continue # Ignoring this event to get the bk
              
              dR_between = myutils.dRAngle(vars_dict["mesonp4"], vars_dict["leptonp4"])
              vars_dict["dR_lep_meson"] = dR_between
              if dR_between < angle_sep[0] or dR_between > angle_sep[1]:
                  continue # Ignoring this event to get the bk
                  
              
              
              # pt a partir de P y theta
              vars_dict["recoMesonPt"] = vars_dict["recoMesonP"] * math.sin(vars_dict["recoMesonTheta"])

              # -----------------------------------------------------------------
              # Histogramas "ALL"
              # -----------------------------------------------------------------
              # fill_histograms_all_events(root_histograms, vars_dict, weight)
              fill_category(root_histograms, vars_dict, "ALL", weight)

              # -----------------------------------------------------------------
              # Clasificación SIGNAL vs BG
              # (usando genTauID y selectGEN)
              # -----------------------------------------------------------------
              if vars_dict["genTauID"] == selectGEN:
                  # Consideramos este evento como señal
                  if tree_key == "original":
                      selectedEvents += 1
                      sumWeights += weight
                      sumWeightsP1 += weight * vars_dict["weight_P1"]
                      sumWeightsM1 += weight * vars_dict["weight_M1"]
                  fill_category(root_histograms, vars_dict, "SIGNAL", weight)
                  fill_special_signal_histograms(special_histograms_super[tree_key], vars_dict, weight)
              else:
                  # fill_general_histograms_BG_events(root_histograms, vars_dict, weight)
                  fill_category(root_histograms, vars_dict, "BG", weight)

                  # Sub-categorías de fondo según genTauID
                  if vars_dict["genTauID"] == -13:  # muones
                    fill_category(root_histograms, vars_dict, "BGMuon", weight)

                  elif vars_dict["genTauID"] == -11:  # electrones
                    fill_category(root_histograms, vars_dict, "BGEle", weight)
                      

                  elif vars_dict["genTauID"] == 0:  # piones
                    fill_category(root_histograms, vars_dict, "BGPion", weight)

                  elif vars_dict["genTauID"] == 1:  # rho
                    fill_category(root_histograms, vars_dict, "BGRho", weight)

                  elif vars_dict["genTauID"] == 10:  # a1
                    fill_category(root_histograms, vars_dict, "BGA1", weight)

                  else:  # other BG
                    other_BG_id[vars_dict["genTauID"]] = other_BG_id.get(vars_dict["genTauID"], 0) + 1 # Pythonic!
                    fill_category(root_histograms, vars_dict, "BGOther", weight)

  return {
      "totalEvents": totalEvents,
      "selectedEvents": selectedEvents,
      "sumWeights": sumWeights,
      "sumWeightsP1": sumWeightsP1,
      "sumWeightsM1": sumWeightsM1,
  }
