#!/usr/bin/env python3
import os
import numpy as np
from Bio.PDB import PDBParser, Superimposer
import csv
import itertools
import sys
# ---------------------------- CONFIG ----------------------------
TARGET_DIR   = "fpocket_results/target"    # contient A20 outward (_POS) et conformations négatives (_NEG)
TEMPLATE_DIR = "fpocket_results/template"  # A19 outward

N_USER = 52 # top N poches à sélectionner
DIST_CUTOFF = 4.0
SO_MIN = 0.8
SO_MAX = 0.95


PDB_TEMPLATE = sys.argv[1]
PDB_TARGET   = sys.argv[2]


parser = PDBParser(QUIET=True)

# ---------------------------- FONCTIONS ----------------------------

# --- 1. ALIGNEMENT GLOBAL ---
def get_matching_ca_atoms(struct1, struct2):
    CA1, CA2 = [], []
    for ch1, ch2 in zip(struct1.get_chains(), struct2.get_chains()):
        for r1, r2 in zip(ch1, ch2):
            if r1.id[0] == " " and r2.id[0] == " ":
                if "CA" in r1 and "CA" in r2:
                    CA1.append(r1["CA"])
                    CA2.append(r2["CA"])
    return CA1, CA2

def compute_global_alignment():
    print("Computing global structure alignment...")
    struct_temp = parser.get_structure("template", PDB_TEMPLATE)
    struct_target = parser.get_structure("target", PDB_TARGET)
    CA_temp, CA_target = get_matching_ca_atoms(struct_temp, struct_target)
    sup = Superimposer()
    sup.set_atoms(CA_temp, CA_target)
    R, t = sup.rotran
    print("✔ Global alignment computed.")
    return R, t

# --- 2. TRANSFORMATION DES COORDONNÉES ---
def transform_coords(atom, R, t):
    return np.dot(atom.coord, R) + t

def get_ca_atoms_aligned(pocket_file, R=None, t=None):
    struct = parser.get_structure("pocket", pocket_file)
    atoms, residues = [], []
    for model in struct:
        for chain in model:
            for res in chain:
                if res.id[0] != " ":
                    continue
                if "CA" in res:
                    atom = res["CA"]
                    if R is not None:
                        atom.coord = transform_coords(atom, R, t)
                    atoms.append(atom)
                    residues.append(res)
    return atoms, residues

# Extraction des numéros de résidus d'une poche
def extract_residue_numbers(pocket_file):
    struct = parser.get_structure("atm", pocket_file)
    residue_numbers = set()
    for model in struct:
        for chain in model:
            for res in chain:
                if res.id[0] == " ":
                    residue_numbers.add(str(res.id[1]))
    return ",".join(sorted(residue_numbers, key=lambda x: int(x)))

# --- 3. SCORE D’OVERLAP SO ---
def compute_SO(temp_pocket_file, target_pocket_file, R, t):
    atoms_temp, residues_temp = get_ca_atoms_aligned(temp_pocket_file, None, None)
    atoms_target, residues_target = get_ca_atoms_aligned(target_pocket_file, R, t)

    matched = 0
    coords_temp = np.array([a.coord for a in atoms_temp])
    coords_target = np.array([a.coord for a in atoms_target])

    for coord_target in coords_target:
        for coord_temp in coords_temp:
            if np.linalg.norm(coord_target - coord_temp) < DIST_CUTOFF:
                matched += 1
                break

    union = len(coords_temp) + len(coords_target)
    SO = 2 * matched / union if union > 0 else 0
    return SO, matched, union

# --- 4. SCORE DE COMPATIBILITÉ SC ---
def compute_SC(target_pocket_file, negative_pockets, R, t):
    atoms_target, residues_target = get_ca_atoms_aligned(target_pocket_file, None, None)
    total_residues = len(residues_target)
    compatible_count = 0

    for neg_pocket_file in negative_pockets:
        atoms_neg, residues_neg = get_ca_atoms_aligned(neg_pocket_file, R, t)
        for resT, atomT in zip(residues_target, atoms_target):
            for resN, atomN in zip(residues_neg, atoms_neg):
                if np.linalg.norm(atomT.coord - atomN.coord) < DIST_CUTOFF:
                    if resT.get_resname()[0] == resN.get_resname()[0]:
                        compatible_count += 1
                        break
    SC = compatible_count / total_residues if total_residues > 0 else 0
    return SC

# --- 5. SCORE GLOBAL ---
def compute_global_score(SO, SC):
    return SO * (1 - SC)

# --- 6. FPocket UTIL ---
def get_pos_pocket_dir(base_dir):
    for d in os.listdir(base_dir):
        if d.endswith("_POS"):
            pockets = os.path.join(base_dir, d, "pockets")
            if os.path.isdir(pockets):
                return pockets
    raise FileNotFoundError(f"No _POS folder in {base_dir}")

def get_neg_pockets(base_dir):
    neg_files = []
    for d in os.listdir(base_dir):
        if d.endswith("_NEG"):
            pockets_dir = os.path.join(base_dir, d, "pockets")
            if os.path.isdir(pockets_dir):
                files = [
                    os.path.join(pockets_dir, f)
                    for f in os.listdir(pockets_dir)
                    if f.endswith("_atm.pdb")
                ]
                neg_files.extend(files)
    return neg_files

def get_top_n_pockets(pocket_dir, n):
    files = [f for f in os.listdir(pocket_dir) if f.endswith("_atm.pdb")]
    files.sort(key=lambda x: int(x.split("pocket")[1].split("_")[0]))
    return [os.path.join(pocket_dir, f) for f in files[:n]]

# ---------------------------- MAIN ----------------------------
def main():
    R, t = compute_global_alignment()

    template_pocket_dir = get_pos_pocket_dir(TEMPLATE_DIR)
    target_pocket_dir   = get_pos_pocket_dir(TARGET_DIR)

    total_template = len([f for f in os.listdir(template_pocket_dir) if f.endswith("_atm.pdb")])
    total_target   = len([f for f in os.listdir(target_pocket_dir) if f.endswith("_atm.pdb")])

    template_pockets = get_top_n_pockets(template_pocket_dir, N_USER)
    target_pockets   = get_top_n_pockets(target_pocket_dir, N_USER)
    negative_pockets = get_neg_pockets(TARGET_DIR)

    results = []

    total_jobs = len(template_pockets) * len(target_pockets)
    job_id = 0

    print(f"\n=== Début du calcul des scores ({total_jobs} comparaisons) ===\n")

    for temp_pocket_file, target_pocket_file in itertools.product(template_pockets, target_pockets):

        job_id += 1
        print(f"[{job_id}/{total_jobs}] Recherche en cours : "
              f"{os.path.basename(temp_pocket_file)} vs {os.path.basename(target_pocket_file)} ...")

        SO, matched, union = compute_SO(temp_pocket_file, target_pocket_file, R, t)
        SC = compute_SC(target_pocket_file, negative_pockets, R, t)
        score_global = compute_global_score(SO, SC)

        residues_list = extract_residue_numbers(target_pocket_file)

        results.append({
            "target": os.path.basename(target_pocket_file),
            "SO": SO,
            "SC": SC,
            "score_global": score_global,
            "residues": residues_list
        })

    results_sorted = sorted(results, key=lambda x: x["score_global"], reverse=True)[:N_USER]

    out_csv = "pocket_scores_selected.csv"
    with open(out_csv, "w", newline="") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=["target","SO","SC","score_global","residues"]
        )
        writer.writeheader()
        for r in results_sorted:
            writer.writerow(r)

    print(f"✔ Done! Top {N_USER} poches écrites dans {out_csv}")

if __name__ == "__main__":
    main()
