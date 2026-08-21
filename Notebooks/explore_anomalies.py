import sys
from pathlib import Path

# Sube un nivel desde Notebooks/ hasta la raíz del repo
repo_root = Path.cwd().parent  
sys.path.append(str(repo_root))

from TauAnalysis.TTreesTausLong import get_final_state_constituents
from modules.tauReco import findAllGenTaus
from modules import myutils
from particle import Particle
from modules.ParticleObjects import GenParticle, RecoParticle


from dataclasses import dataclass
import logging
from tqdm import tqdm
import pandas as pd


import sys
from pathlib import Path
from multiprocessing import Pool, cpu_count
import argparse
repo_root = Path.cwd().parent
sys.path.append(str(repo_root))

from modules.tauReco import findAllGenTaus
from podio import root_io
loggers = {
    "config": logging.getLogger("config"),
    "io": logging.getLogger("io"),
    "processing": logging.getLogger("processing"),
    "pi0mass": logging.getLogger("pi0mass")
}
@dataclass
class arguments:
    samples_config: str
    input_list: list

sample = "ztt_2M"
gatr_path = None
samples_config="/nfs/cms/arqolmo/TausFCCee/config/samples/samples.yaml"
input_list = []
args = arguments(samples_config=samples_config, input_list=input_list)
# loggers = None
args.samples_config
filenames, mlpf_results = myutils.get_root_trees_path(
        sample, gatr_path, loggers, False, args
    )

def search_file(filename, print_every=1000):
    """Busca en un único archivo el primer evento con más de 2 genTaus."""
    cases = {}
    file_reader = root_io.Reader(filename)
    for eventid, event in enumerate(file_reader.get("events")):
        # if eventid % print_every == 0:
            # print(f"[{filename}] id {eventid}")
        mc_particles = event.get("MCParticles")
        genTaus: dict[int, GenParticle] = findAllGenTaus(mc_particles)
        if len(genTaus)!=2:
            continue
        if genTaus[0].getID() != 0 or genTaus[1].getID() != 0:
            continue 
        for genTau in genTaus:
            daugs = genTaus[genTau].getDaughters()
            if len(daugs)>1:
                n_charged = 0
                for d in daugs:
                    if daugs[d].getCharge()!=0:
                        n_charged+=1
                # print(f"Más de un hijo para desintegraciones a un pion {filename}, evento {eventid}")
                if n_charged>1:
                    if not cases.get(filename):
                        cases[filename]=set()
                    cases[filename].add(eventid)  
    return cases

def collapse_chain(path):
    """Colapsa PDGs consecutivos repetidos en una cadena."""
    collapsed = [path[0]]
    for pdg in path[1:]:
        if pdg[0] != collapsed[-1][0]:
            collapsed.append(pdg)
    return collapsed

def print_final_state_tree(particle, depth=0, path=None, stop_at_gen_status1=True):
    if path is None:
        path = [(particle.getPDG(), particle.getObjectID().index)]

    daughters = list(particle.getDaughters())
    is_generator_final = particle.getGeneratorStatus() == 1

    # Estado final de generador: para. Si no hay hijas tampoco sigas (evita crash).
    if (stop_at_gen_status1 and is_generator_final) or len(daughters) == 0:
        chain = " -> ".join(Particle.from_pdgid(pdg).name+" id:"+str(ID_) for pdg, ID_ in collapse_chain(path))
        tag = "GEN-FINAL" if is_generator_final else "LEAF(sim/other)"
        print(f"{'  ' * depth}[{tag}] PDG={particle.getPDG()} "
              f"(status={particle.getGeneratorStatus()}) | chain: {chain}")
        return

    for daughter in daughters:
        print_final_state_tree(daughter, depth + 1, path + [(daughter.getPDG(), daughter.getObjectID().index)],
                                stop_at_gen_status1=stop_at_gen_status1)
        
if __name__ == "__main__":
    # filenames viene de tu myutils.get_root_trees_path(...)
    parser = argparse.ArgumentParser()
    parser.add_argument("--explore", action="store_true")
    
    args = parser.parse_args()
    if args.explore:
        n_workers = 120
        cases = {}
        found = None
        with Pool(processes=n_workers) as pool:
            for result in tqdm(pool.imap_unordered(search_file, filenames), total=len(filenames)):
                cases.update(result)
                # if result is not None:
                    # found = result
                    # pool.terminate()  # corta el resto de workers en cuanto hay un match
                    # break
    # Formato largo: una fila por (filename, eventid)
        rows = [
            {"filename": fname, "eventid": eid}
            for fname, eids in cases.items()
            for eid in eids
        ]
        df = pd.DataFrame(rows, columns=["filename", "eventid"])
        df.to_csv("Cases_more_mesons.csv", index=False)
        print(f"Total de casos encontrados: {len(df)}")
    else:
        df = pd.read_csv("/nfs/cms/arqolmo/TausFCCee/Cases_more_mesons.csv")

        # Agrupamos los eventid por archivo, para no reabrir el reader por cada evento
        events_by_file = df.groupby("filename")["eventid"].apply(list).to_dict()

        stop = False
        for filename, eventids in events_by_file.items():
            if stop:
                break

            wanted_events = set(eventids)
            reader = root_io.Reader([filename])

            for eventid, event in enumerate(reader.get("events")):
                if eventid not in wanted_events:
                    continue

                print(f"\n{'='*80}")
                print(f"File: {filename}")
                print(f"Event: {eventid}")
                print(f"{'='*80}")

                mc_p = event.get("MCParticles")
                genTaus = findAllGenTaus(mc_p)

                print("#### GEN TAUS ###")
                for gen in genTaus:
                    tau = genTaus[gen]
                    print(f"\n-- GenTau id={gen} -- {tau}")
                    daugs = tau.getDaughters()
                    print(f"   nDaughters={len(daugs)}")
                    for d in daugs:
                        try:
                            name = Particle.from_pdgid(daugs[d].getPDG()).name
                        except Exception:
                            name = "?"
                        print(f"     -> {name} (PDG={daugs[d].getPDG()}, "
                              f"status={daugs[d].getGeneratorStatus()}, "
                              f"charge={daugs[d].getCharge()})")

                print("\n#### MC PARTICLES (Tree from e-) ###")
                for particle in mc_p:
                    if particle.getPDG() == 11 and len(list(particle.getParents())) == 0:
                        print_final_state_tree(particle)

                try:
                    user_input = input("\n[Enter] siguiente evento | 'q' salir: ")
                except EOFError:
                    user_input = "q"

                if user_input.strip().lower() == "q":
                    stop = True
                    break