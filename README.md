# Paris_Team_2

# MEET-EU PROJECT: SLC6A20 Selective Inhibitor Discovery Pipeline

**Date:** Janvier 2026

## Authors

**MEET-EU Project team from Sorbonne University**

- **Christine ANTON JESUTHASAN**
- **Ndeye Fatou DIENG**
- **Cheikh LO**
- **Adi THOHA**
- **Badiallo DRAMÉ**


## Project Overview

This project was developed as part of **Meet-EU**, an international team-based course organized by six 4EU+ member universities (Copenhagen, Heidelberg, Milan, Paris, Prague, Warsaw) during the academic year 2025-2026. Meet-EU promotes interdisciplinary training through research, bringing together students from computer science, physics, mathematics, chemistry, biology, and biotechnologies to address ambitious and open research problems in groups of 4-5 students.

This project implements an automated computational pipeline for discovering **selective inhibitors of SLC6A20 (SIT1)**, a human transporter involved in glycine and proline uptake and implicated in diabetes and COVID-19 severity.

A key challenge is achieving selectivity against its close homolog SLC6A19 (B0AT1), which shares the same LeuT fold and exhibits high structural similarity.

SLC6A20 belongs to the SLC6 transporter family, whose members universally operate through a gated-pore mechanism. This mechanism involves cycling through three conserved conformational states:

- **Outward-open conformation** — extracellular side open; can expose allosteric sites.
- **Occluded conformation** — substrate (when present) is trapped between both gates.
- **Inward-open conformation** — intracellular side open; often stabilized by inhibitors in related SLC6 transporters.

This pipeline focuses on identifying an **allosteric pocket in the outward-open conformation of SLC6A20**, which provides a promising opportunity to achieve selectivity over SLC6A19.

## Pipeline Architecture

The pipeline is organized into three major components:

1. **Structure-based pocket identification and scoring**
   - Identifies and scores binding pockets across template and target conformations to identify target-specific sites for selective docking.

2. **Ligand library generation and preparation**
   - Generates and prepares a focused chemical library through ChEMBL similarity search from known inhibitors.

3. **Comparative molecular docking and selectivity ranking**
   - Docks all ligands against both transporters and computes selectivity scores to prioritize target-preferential binders.

All workflow steps are fully orchestrated and automated using **Nextflow**.

## Pipeline Structure

```
root_dir/
├── data/
│   ├── target/
│   │   ├── positive/              # Target protein positive conformation
│   │   │   └── SIT1_model_oo.pdb
│   │   └── negative/              # Target protein negative conformations
│   │       └── 8I91_chainB.pdb
│   │       └── 8WM3_chainC.pdb
│   ├── template/
│   │   └── positive/              # Template protein positive conformation
│   │       └── 8WBY_chainA.pdb
│   └── known_inhibitors.tsv
├── src/
│   ├── fpocket.sh                 # Pocket detection wrapper
│   ├── pocket_scoring.py          # Pocket similarity scoring
│   ├── resolve_seeds_chembl.py
│   ├── chembl_homolog_search.py
│   ├── ligands_to_3d.py
│   ├── minimize_ligands.py
│   ├── sdf_to_pdbqt_obabel.py
│   └── docking.py
├── main.nf                     # Nextflow pipeline
├── envs/
│   ├── environment.yml
│   ├── slc6env.yml                   
├── workflow_diagram.png
├── README.md
└── CODE_OF_CONDUCT.md
```

## Requirements

### 1. Required files

#### Protein structures

- **Target protein (SLC6A20)** in PDB format:
  - Positive conformations (`data/target/positive/*.pdb`): conformations of interest where we aim to find selective inhibitors.
  - Negative conformations (`data/target/negative/*.pdb`): other possible conformations used to filter out non-specific pockets.
- **Template protein (SLC6A19)** positive conformation: `data/template/positive/*.pdb`

#### Known Inhibitors

- TSV file: `data/known_inhibitors.tsv`
  - A tab separated file with 3 columns:
    - `seed_name` (Unique ID)
    - `query_type` (`name`, `chembl_id`, `smiles`)
    - `query` (the value to search)

### 2. User configurable parameters

The user has control over the following parameters in `src/pocket_scoring.py`:

- `N_USER` (number of top pockets to analyze)
- `DIST_CUTOFF` (Distance (Å) for residue overlap)
- `SO_MIN` (Minimum shape overlap score)
- `SO_MAX` (Maximum shape overlap score)

The user can choose to include or exclude negative conformations for analysis. Default logic for this implementation:

- **Template (SLC6A19)**: Only positive conformation is used
- **Target (SLC6A20)**: Both positive and negative conformations are analyzed to ensure selected pockets are specific to the conformation of interest (positive) and not shared with other conformations (negative)

## Installation and dependencies

To run this pipeline, you need the following installed:

### 1. Conda Environment

a. If the `slc6env` virtual environment already exists, activate it using:
```bash
conda activate envs/slc6env
```

b. If no environment exists, create it (it contains all required dependencies):
- First load Miniconda: `module load miniconda`
- Then navigate to your preferred path: `cd your_path/envs`
- And create the environment: `conda env create -f environment.yml -p your_env_name`

c. If a library is missing, add it to `envs/environment.yml` and recreate or update the environment.

### 2. Nextflow (≥21.04)

Install it with the official script:
```bash
curl -s https://get.nextflow.io | bash
```

### 3. Java (≥11)

Install Java 11 or newer, which is required to run Nextflow.

### 4. PyMOL (3.1.6.1)

Install PyMOL for the vizualization. 

## Usage

To launch the full pipeline run the following command from the project root:
```bash
nextflow main.nf
```

## Pipeline Workflow

The pipeline executes the following steps:

### Part 1: Pocket-Based Analysis

1. **FPOCKET Detection** (`RUN_FPOCKET`)
   - Detects binding pockets in all protein conformations
   - Output: `results/fpocket_results/{target,template}/`

2. **Pocket Scoring** (`POCKET_SCORING`)
   - Aligns template and target structures
   - Computes Shape Overlap (SO) and Compatibility (SC) scores
   - Selects top N_USER pockets based on global score: SO × (1 - SC)
   - Output: `results/pocket_scoring/pocket_scores_ranked.csv`

### Part 2: Ligand Preparation

3. **Seed Resolution** (`RESOLVE_SEEDS`)
   - Resolves compound names/IDs to chemical structures via ChEMBL API
   - Output: `results/ligands/seeds/seeds.sdf`

4. **ChEMBL Similarity Search** (`CHEMBL_HOMOLOGS`)
   - Expands seed library using Tanimoto similarity (threshold: 40%)
   - Output: `results/ligands/homologs/chembl_homologs.sdf`

5. **3D Conformation Generation** (`LIGANDS_3D`)
   - Generates 3D coordinates using ETKDGv3 algorithm
   - Output: `results/ligands/3d/ligands_3d.sdf`

6. **Energy Minimization** (`MINIMIZE_LIGANDS`)
   - Refines structures using UFF force field
   - Output: `results/ligands/minimized/ligands_minimized.sdf`

7. **PDBQT Conversion** (`SDF_TO_PDBQT`)
   - Prepares ligands for docking (adds hydrogens and charges)
   - Output: `results/ligands/pdbqt/*.pdbqt`

### Part 3: Molecular Docking & Selectivity Analysis

8. **Comparative Docking** (`DOCKING`)
   - Docks all ligands against template (SLC6A19) and target (SLC6A20) pockets
   - Uses AutoDock Vina with parallel processing
   - Calculates binding affinity scores for each ligand-pocket pair

9. **Selectivity Scoring**
   - Computes selectivity score for each ligand
   - Ranks ligands by preference for target over template
   - Identifies top 10 most selective compounds

## Outputs

```
root_dir/
├── results/
│   ├── fpocket_results/
│   │   └── Output directories from fpocket for each protein conformation (see fpocket documentation)
│   ├── ligands/
│   │   └── PDBQT files of all ligands retrieved and processed from the ChEMBL API
│   ├── docking_results/
│   │   ├── docking_selectivity_scores.tsv      # Ranked docking and selectivity scores for each ligand-pocket pair
│   │   └── top10_ligands/
│   │       └── PDBQT files of the 10 most selective ligands (both template and target pockets)
│   └── pocket_scores_ranked.csv   # List of binding site residues with specificity scores.
```

## References

- **Fpocket User Manual**: [manual_fpocket2.pdf](https://fpocket.sourceforge.net/manual_fpocket2.pdf)

- Xu, J., Hu, Z., Dai, L. *et al.* **Molecular basis of inhibition of the amino acid transporter B0AT1 (SLC6A19)**. *Nat Commun* 15, 7224 (2024). https://doi.org/10.1038/s41467-024-51748-1

- Li, Y., Chen, Y., Zhang, Y. *et al.* **Structural insight into the substrate recognition and transport mechanism of amino acid transporter complex ACE2-B0AT1 and ACE2-SIT1**. *Cell Discov* 9, 93 (2023). https://doi.org/10.1038/s41421-023-00596-2

- Claire Colas, Elodie Laine, **Targeting Solute Carrier Transporters through Functional Mapping**, Trends in Pharmacological Sciences, https://doi.org/10.1016/j.tips.2020.11.005

- J. Eberhardt, D. Santos-Martins, A. F. Tillack, and S. Forli, **AutoDock Vina 1.2.0: New Docking Methods, Expanded Force Field, and Python Bindings**, J. Chem. Inf. Model. (2021). DOI 10.1021/acs.jcim.1c00203

- O. Trott, A. J. Olson, **AutoDock Vina: improving the speed and accuracy of docking with a new scoring function, efficient optimization and multithreading**, J. Comp. Chem. (2010). DOI 10.1002/jcc.21334. Please see https://github.com/ccsb-scripps/AutoDock-Vina for more information.

- Zayats, Vasilina, et al. (2023), **Conservation of knotted and slipknotted topology in transmembrane transporters**. *Biophysical Journal* 122.23: 4528-4541

- **ChEMBL**: Mendez, D *et al* (2019) ChEMBL: **towards direct deposition of bioassay data**. Nucl. Acids Res, 47, D930-40, https://www.ebi.ac.uk/chembl/

- **RDKit**: https://www.rdkit.org/docs/api-docs.html

- **Nextflow**: https://www.nextflow.io/

## How to Contribute

We welcome contributions to improve and extend this pipeline! Here are several ways you can help:

- Reporting issues
- Suggesting enhancements
- Contributing code
- Improving documentation
- Adding tests
- Fixing known issues

If you have any suggestions, questions, or would like to discuss potential contributions, feel free to contact us on Github! We appreciate all contributions, big or small!








