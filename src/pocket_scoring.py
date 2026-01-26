#!/usr/bin/env python3

"""
Pocket scoring script to select target-positive-specific pockets for selective docking
"""

import os
import numpy as np
from Bio.PDB import PDBParser, Superimposer
import csv
import itertools
import sys


# =========================
#   CONFIGURATION
# =========================

# Directories containing fpocket results
FPOCKET_BASE = sys.argv[3]

TARGET_DIR   = os.path.join(FPOCKET_BASE, "target")     # Target protein: positive (_POS) and negative (_NEG) fpocket results
TEMPLATE_DIR = os.path.join(FPOCKET_BASE, "template")   # Template protein: positive (_POS) fpocket results


# Pocket filtering parameters
N_USER = 52        # Number of top pockets to keep
DIST_CUTOFF = 4.0  # Distance cutoff (Å) for residue overlap
SO_MIN = 0.8       # Minimum geometric overlap score
SO_MAX = 0.95      # Maximum geometric overlap score

# Protein structures (positive conformations) for alignment
PDB_TEMPLATE = sys.argv[1] # Target positive conformation
PDB_TARGET   = sys.argv[2] # Template positive conformation

parser = PDBParser(QUIET=True)


# =========================
#   GLOBAL ALIGNMENT
# =========================

def get_matching_ca_atoms(struct1, struct2):
    """
    Extracts corresponding Cα atoms from two aligned protein structures.

    Parameters :
    ------------
    struct1 : Bio.PDB.Structure.Structure
        First protein structure object
    struct2 : Bio.PDB.Structure.Structure
        Second protein structure object (must have same chain/residue numbering)

    Returns :
    ---------
    tuple[list, list]
        A tuple containing two lists:
        - First list: Cα atom objects from struct1
        - Second list: Corresponding Cα atom objects from struct2

    Notes :
    -------
    - Structures must be aligned and share identical chain/residue numbering
    - Only standard residues with a Cα atom are considered
    """
    CA1, CA2 = [], []

    for ch1, ch2 in zip(struct1.get_chains(), struct2.get_chains()):
        for r1, r2 in zip(ch1, ch2):
            if r1.id[0] == " " and r2.id[0] == " ":
                if "CA" in r1 and "CA" in r2:
                    CA1.append(r1["CA"])
                    CA2.append(r2["CA"])

    return CA1, CA2


def compute_global_alignment():
    """
    Computes global structural alignment between template and target proteins (positive conformations).

    Parameters :
    ------------
    None (uses module-level constants PDB_TEMPLATE and PDB_TARGET)

    Returns :
    ---------
    tuple[numpy.ndarray, numpy.ndarray]
        A tuple containing:
        - R : numpy.ndarray (3x3)
            Rotation matrix for optimal alignment
        - t : numpy.ndarray (1x3)
            Translation vector for optimal alignment

    Process :
    ---------
    1. Parses template and target structures from PDB files
    2. Extracts corresponding Cα atoms using get_matching_ca_atoms()
    3. Computes optimal superimposition using BioPython's Superimposer
    4. Returns rotation matrix and translation vector

    Notes :
    -------
    - Uses Cα atoms for alignment
    - Structures must share matching chain/residue numbering
    """
    print("Computing global structure alignment...")

    struct_temp = parser.get_structure("template", PDB_TEMPLATE)
    struct_target = parser.get_structure("target", PDB_TARGET)

    CA_temp, CA_target = get_matching_ca_atoms(struct_temp, struct_target)
    sup = Superimposer()
    sup.set_atoms(CA_temp, CA_target)
    R, t = sup.rotran
    print("Global alignment computed.")

    return R, t


# =========================
#   COORDINATE TRANSFORMATION
# =========================

def transform_coords(atom, R, t):
    """
    Applies a rigid-body transformation to atomic coordinates.

    Parameters :
    ------------
    atom : Bio.PDB.Atom.Atom
        Atom object containing coordinates to transform
    R : numpy.ndarray (3x3)
        Rotation matrix
    t : numpy.ndarray (1x3) or (3,)
        Translation vector

    Returns :
    ---------
    numpy.ndarray
        Transformed coordinates as a 3-element numpy array

    Transformation :
    ----------------
    Performs the linear transformation:
        new_coords = dot(atom.coord, R) + t

    Where:
        - atom.coord is the original (x, y, z) coordinates
        - R is the 3x3 rotation matrix (applied first)
        - t is the translation vector (applied after rotation)

    Notes :
    -------
    - Applies the global structural alignment (rotation + translation) 
    so that the pocket coordinates of the target are in the same frame 
    as the template
    - Necessary for meaningful comparison of pockets
    """
    return np.dot(atom.coord, R) + t


def get_ca_atoms_aligned(pocket_file, R=None, t=None):
    """
    Extracts Cα atoms from a pocket PDB file with optional structural alignment.

    Parameters :
    ------------
    pocket_file : str
        Path to the PDB file containing the pocket structure
    R : numpy.ndarray, optional
        3x3 rotation matrix for structural alignment
    t : numpy.ndarray, optional
        Translation vector for structural alignment (1x3 or (3,))

    Returns :
    ---------
    tuple[list, list]
        A tuple containing:
        - atoms : list[Bio.PDB.Atom.Atom]
            List of Cα atom objects (coordinates are transformed if R/t provided)
        - residues : list[Bio.PDB.Residue.Residue]
            List of corresponding residue objects

    Process :
    ---------
    1. Parses the pocket structure from the PDB file
    2. Iterates through all models, chains, and residues
    3. Selects only standard residues (ignores hetero residues)
    4. Extracts Cα atoms and their parent residues
    5. If R and t are provided, transforms atom coordinates in-place

    Notes :
    -------
    - Extracts only standard residues with Cα atoms
    - Applies global alignment (R, t) if provided, to compare pockets in the same frame
    """
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


# =========================
#   POCKET RESIDUE NUMBERS
# =========================

def extract_residue_numbers(pocket_file):
    """
    Extracts residue numbers from a pocket PDB file.

    Parameters :
    ------------
    pocket_file : str
        Path to the PDB file containing the pocket structure

    Returns :
    ---------
    str
        Comma-separated string of residue numbers sorted numerically

    Notes :
    -------
    - Useful for docking and pocket comparison after structural alignment
    """
    struct = parser.get_structure("atm", pocket_file)
    residue_numbers = set()

    for model in struct:
        for chain in model:
            for res in chain:
                if res.id[0] == " ":
                    residue_numbers.add(str(res.id[1]))

    return ",".join(sorted(residue_numbers, key=lambda x: int(x)))


# =========================
#   SO : OVERLAP SCORE
# =========================

def compute_SO(temp_pocket_file, target_pocket_file, R, t):
    """
    Computes Shape Overlap (SO) score between two pockets after structural alignment.

    Parameters :
    ------------
    temp_pocket_file : str
        Path to template pocket PDB file
    target_pocket_file : str
        Path to target pocket PDB file
    R : numpy.ndarray (3x3)
        Rotation matrix from global structure alignment
    t : numpy.ndarray (1x3)
        Translation vector from global structure alignment

    Returns :
    ---------
    tuple[float, int, int]
        A tuple containing:
        - SO : float
            Shape Overlap score (0 to 1), where 1 indicates perfect spatial overlap
        - matched : int
            Number of Cα atoms in target pocket that have a counterpart in template pocket
        - union : int
            Total number of unique Cα atoms in both pockets (sum of counts)

    Calculation :
    -------------
    1. Extracts Cα atoms from both pockets (target atoms are transformed by R/t)
    2. For each target Cα, checks if any template Cα is within DIST_CUTOFF distance
    3. SO = 2 × matched / (len(template) + len(target))

    Notes :
    -------
    - Useful for choosing pockets that resemble each other but are not exactly the 
    same between template and target positive conformations
    """
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


# =========================
#   SC : COMPATIBILITY SCORE
# =========================

def compute_SC(target_pocket_file, negative_pockets, R, t):
    """
    Computes Compatibility Score (SC) for a target pocket.

    Parameters :
    ------------
    target_pocket_file : str
        Path to target pocket PDB file
    negative_pockets : list[str]
        List of file paths to negative/control pocket PDB files
    R : numpy.ndarray (3x3)
        Rotation matrix from global structure alignment
    t : numpy.ndarray (1x3)
        Translation vector from global structure alignment

    Returns :
    ---------
    float
        Compatibility Score (0 to 1), where 1 indicates perfect conservation

    Calculation :
    -------------
    1. Extracts Cα atoms and residues from target pocket
    2. For each negative pocket (aligned using R/t):
    a. For each target residue, checks if any negative pocket residue is:
        - Within DIST_CUTOFF distance (spatial proximity)
        - Has the same first letter of residue name (e.g., 'ALA' vs 'ARG' = different)
    3. SC = (compatible_count / total_target_residues)

    Notes :
    -------
    - Useful for choosing pockets that are specific to the target positive conformation
    and not present in negative/control conformations
    """
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


# =========================
#   GLOBAL SCORE
# =========================

def compute_global_score(SO, SC):
    """
    Computes global similarity score from Shape Overlap and Compatibility Score.

    Parameters :
    ------------
    SO : float
        Shape Overlap score (0 to 1)
    SC : float
        Compatibility Score (0 to 1)

    Returns :
    ---------
    float
        Global similarity score (0 to 1)

    Notes :
    -------
    - Global score = SO * (1 - SC)
    - Encourages pockets that:
        Are similar between target positive and template positive (high SO, but <1)
        Are specific to target positive vs negative conformations (low SC)
    - Higher global score = pockets most promising for selective targeting
    """
    return SO * (1 - SC)


# =========================
#   POCKET SELECTION
# =========================

def get_pos_pocket_dir(base_dir):
    """
    Locates the directory containing positive pocket structures.

    Parameters :
    ------------
    base_dir : str
        Base directory containing subdirectories with pocket data

    Returns :
    ---------
    str
        Path to the 'pockets' subdirectory within the first found '*_POS' folder

    Process :
    ---------
    1. Scans all entries in base_dir
    2. Finds the first directory whose name ends with '_POS'
    3. Checks if 'pockets' subdirectory exists within it
    4. Returns the full path to this 'pockets' directory

    Notes :
    -------
    - Assumes a specific directory structure: base_dir/*_POS/pockets/
    """
    for d in os.listdir(base_dir):
        if d.endswith("_POS"):
            pockets = os.path.join(base_dir, d, "pockets")
            if os.path.isdir(pockets):
                return pockets
            
    raise FileNotFoundError(f"No _POS folder in {base_dir}")


def get_neg_pockets(base_dir):
    """
    Collects all negative pocket structure files from a directory hierarchy.

    Parameters :
    ------------
    base_dir : str
        Base directory containing subdirectories with negative pocket data

    Returns :
    ---------
    list[str]
        List of full paths to negative pocket PDB files

    Process :
    ---------
    1. Scans all entries in base_dir
    2. For each directory ending with '_NEG':
    a. Checks for 'pockets' subdirectory
    b. Collects all files ending with '_atm.pdb' within it
    3. Returns combined list of all found negative pocket files

    Notes :
    -------
    - Assumes a specific directory structure: base_dir/*_NEG/pockets/*_atm.pdb
    - Files are expected to have '_atm.pdb' suffix (fpocket outputs)
    """
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
    """
    Selects the top N pocket files from a directory based on pocket numbering.

    Parameters :
    ------------
    pocket_dir : str
        Directory containing pocket PDB files
    n : int
        Number of top pocket files to return

    Returns :
    ---------
    list[str]
        List of full paths to the first N pocket files, sorted by pocket number

    Notes :
    -------
    - Assumes pocket numbering in filename format: *pocketX_* where X is an integer
    - Files are sorted numerically by pocket number
    - Typical use: get top 5 pockets for docking or analysis
    """
    files = [f for f in os.listdir(pocket_dir) if f.endswith("_atm.pdb")]
    files.sort(key=lambda x: int(x.split("pocket")[1].split("_")[0]))
    return [os.path.join(pocket_dir, f) for f in files[:n]]


# =========================
#          MAIN
# =========================

def main():
    """
    Main execution function for pocket similarity scoring analysis.

    Parameters :
    ------------
    Uses constants: TEMPLATE_DIR, TARGET_DIR, N_USER, DIST_CUTOFF

    Outputs :
    ---------
    - CSV file 'pocket_scores_ranked.csv' with columns:
    * template_pocket, target_pocket: Pocket filenames
    * SO, SC, global_score: Computed similarity scores
    * template_residues, target_residues: Comma-separated residue numbers

    Process Flow :
    --------------
    1. Computes global structural alignment between template and target proteins
    2. Locates positive pocket directories for both proteins
    3. Checks that N_USER does not exceed available positive pockets
    4. Selects top N_USER pockets from each protein
    5. Gathers negative pockets from target directory for compatibility scoring
    6. Performs all-vs-all pocket comparisons (N_USER * N_USER combinations)
    7. For each pocket pair:
        a. Computes Shape Overlap (SO) score
        b. Computes Compatibility Score (SC) 
        c. Calculates global similarity score
        d. Extracts residue numbers for both pockets
    8. Sorts results by global_score (descending) and selects top N_USER
    9. Writes results to CSV file

    Notes :
    -------
    - Assumes pocket files follow naming convention: *_atm.pdb
    """
    R, t = compute_global_alignment()

    template_pocket_dir = get_pos_pocket_dir(TEMPLATE_DIR)
    target_pocket_dir   = get_pos_pocket_dir(TARGET_DIR)

    total_template = len([f for f in os.listdir(template_pocket_dir) if f.endswith("_atm.pdb")])
    total_target   = len([f for f in os.listdir(target_pocket_dir) if f.endswith("_atm.pdb")])

    # N_USER Verification
    if N_USER > total_template:
        print(f"Error : N_USER = {N_USER} but only {total_template} pockets POS found in TEMPLATE.")
        return
    if N_USER > total_target:
        print(f"Error : N_USER = {N_USER} but only {total_target} pockets POS found in TARGET.")
        return

    template_pockets = get_top_n_pockets(template_pocket_dir, N_USER)
    target_pockets   = get_top_n_pockets(target_pocket_dir, N_USER)
    negative_pockets = get_neg_pockets(TARGET_DIR)

    results = []

    total_jobs = len(template_pockets) * len(target_pockets)
    job_id = 0

    print(f"\n=== Start of the score calculation ({total_jobs} comparisons) ===\n")

    for temp_pocket_file, target_pocket_file in itertools.product(template_pockets, target_pockets):

        job_id += 1
        print(f"[{job_id}/{total_jobs}] Research in progress : "
              f"{os.path.basename(temp_pocket_file)} vs {os.path.basename(target_pocket_file)} ...")

        SO, matched, union = compute_SO(temp_pocket_file, target_pocket_file, R, t)
        SC = compute_SC(target_pocket_file, negative_pockets, R, t)
        global_score = compute_global_score(SO, SC)

        results.append({
            "template_pocket": os.path.basename(temp_pocket_file),
            "target_pocket": os.path.basename(target_pocket_file),
            "SO": SO,
            "SC": SC,
            "global_score": global_score,
            "template_residues": extract_residue_numbers(temp_pocket_file),
            "target_residues": extract_residue_numbers(target_pocket_file)
        })

    results_sorted = sorted(results, key=lambda x: x["global_score"], reverse=True)[:N_USER]

    out_csv = "pocket_scores_ranked.csv"
    with open(out_csv, "w", newline="") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=["template_pocket", "target_pocket", "SO", "SC", "global_score", "template_residues", "target_residues"]

        )
        writer.writeheader()
        for r in results_sorted:
            writer.writerow(r)

    print(f"Done! Top {N_USER} pockets written in {out_csv}")


if __name__ == "__main__":
    main()
