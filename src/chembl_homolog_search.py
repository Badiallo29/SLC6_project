#!/usr/bin/env python3
"""
ChEMBL similarity expansion from seed structures.

Input :
  seeds/seeds.sdf

Output :
  homologs/chembl_homologs.sdf
  homologs/chembl_edges.tsv
"""

import argparse
import time
import requests
from pathlib import Path
from rdkit import Chem
from rdkit.Chem import AllChem

# Base URL for the ChEMBL REST API
CHEMBL_BASE = "https://www.ebi.ac.uk/chembl/api/data"


def chembl_similarity(chembl_id, threshold=40, page_size=1000):
    # Perform a similarity search in ChEMBL for a given molecule
    # chembl_id : reference ChEMBL compound ID
    # threshold : minimum similarity score (Tanimoto, in %)
    # page_size : maximum number of results per API page

    # Initial similarity search endpoint
    url = f"{CHEMBL_BASE}/similarity/{chembl_id}/{threshold}"
    params = {"format": "json", "limit": page_size}

    results = []

    # Iterate over paginated API results
    while url:
        r = requests.get(url, params=params, timeout=60)
        r.raise_for_status()  # Fail on HTTP errors
        data = r.json()

        # Accumulate similarity hits
        results.extend(data.get("molecules", []))

        # Follow pagination link if available
        url = data.get("page_meta", {}).get("next")

    return results


def main(in_sdf, out_sdf, out_edges, threshold):
    # Load seed molecules from the input SDF file
    suppl = Chem.SDMolSupplier(in_sdf, removeHs=False)

    # Dictionary to store unique homolog molecules (deduplicated by ChEMBL ID)
    out_mols = {}

    # List of edges linking seeds to homologs with similarity scores
    edges = []

    # Iterate over all molecules in the seed SDF
    for mol in suppl:
        if mol is None:
            # Skip unreadable molecules
            continue

        # Only process molecules explicitly marked as seeds
        if not mol.HasProp("is_seed"):
            continue

        seed_name = mol.GetProp("seed_name")
        chembl_id = mol.GetProp("chembl_id") if mol.HasProp("chembl_id") else None

        # Skip seeds that cannot be queried in ChEMBL
        if not chembl_id:
            print(f"[INFO] Skipping similarity for {seed_name} (no ChEMBL ID)")
            continue

        print(f"[INFO] Similarity search for {seed_name} ({chembl_id})")

        # Query ChEMBL for similar compounds
        hits = chembl_similarity(chembl_id, threshold)

        # Process each similarity hit
        for h in hits:
            h_id = h.get("molecule_chembl_id")
            sim = h.get("similarity")

            # Extract molecular structure information
            structs = h.get("molecule_structures") or {}
            smiles = structs.get("canonical_smiles")

            # Skip incomplete entries
            if not h_id or not smiles:
                continue

            # Parse homolog SMILES into an RDKit molecule
            mol_h = Chem.MolFromSmiles(smiles)
            if mol_h is None:
                continue

            # Annotate homolog molecule with metadata
            mol_h.SetProp("_Name", h_id)
            mol_h.SetProp("chembl_id", h_id)
            mol_h.SetProp("canonical_smiles", smiles)
            mol_h.SetProp("is_seed", "false")
            mol_h.SetProp("best_seed", seed_name)
            mol_h.SetProp("similarity", str(sim))

            # Deduplicate homologs using the ChEMBL ID
            if h_id not in out_mols:
                out_mols[h_id] = mol_h

            # Record the seed–homolog relationship
            edges.append((seed_name, chembl_id, h_id, sim))

        # Sleep briefly to avoid overloading the ChEMBL API
        time.sleep(0.2)

    # Write all unique homolog molecules to an SDF file
    writer = Chem.SDWriter(out_sdf)
    for mol in out_mols.values():
        writer.write(mol)
    writer.close()

    # Write the seed–homolog edge list as a TSV file
    with open(out_edges, "w") as f:
        f.write("seed_name\tseed_chembl_id\thomolog_chembl_id\tsimilarity\n")
        for e in edges:
            f.write("\t".join(map(str, e)) + "\n")


if __name__ == "__main__":
    # Parse command-line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("--in_sdf", required=True)
    parser.add_argument("--out_sdf", required=True)
    parser.add_argument("--out_edges", required=True)
    parser.add_argument("--threshold", type=int, default=40)
    args = parser.parse_args()

    # Run the ChEMBL similarity expansion workflow
    main(args.in_sdf, args.out_sdf, args.out_edges, args.threshold)

