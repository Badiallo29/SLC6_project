#!/usr/bin/env python3
"""
Resolve known inhibitors into chemical structures.

Input  : seeds/known_inhibitors.tsv
Output : seeds/seeds.sdf
         seeds/seeds_resolve.tsv

Supports:
- name       -> resolved via ChEMBL API
- chembl_id  -> resolved via ChEMBL API
- smiles     -> parsed directly with RDKit
"""

import argparse
import csv
import requests
from rdkit import Chem
from rdkit.Chem import AllChem

# Base URL for the ChEMBL REST API
CHEMBL_BASE = "https://www.ebi.ac.uk/chembl/api/data"


# ChEMBL helpers

def chembl_search_by_name(name, timeout=30):
    # Query ChEMBL by compound name
    r = requests.get(
        f"{CHEMBL_BASE}/molecule/search",
        params={"q": name, "format": "json"},
        timeout=timeout,
    )
    r.raise_for_status()  # Raise an exception for HTTP errors
    data = r.json()

    # Extract list of candidate molecules
    mols = data.get("molecules", [])
    if not mols:
        # No hit found
        return None, None

    # Take the top hit for deterministic behavior
    m = mols[0]
    chembl_id = m.get("molecule_chembl_id")
    smiles = (m.get("molecule_structures") or {}).get("canonical_smiles")
    return chembl_id, smiles


def chembl_get_by_id(chembl_id, timeout=30):
    # Query ChEMBL by a specific ChEMBL identifier
    r = requests.get(
        f"{CHEMBL_BASE}/molecule/{chembl_id}.json",
        timeout=timeout,
    )
    r.raise_for_status()  # Raise an exception for HTTP errors
    m = r.json()

    # Extract canonical SMILES representation
    smiles = (m.get("molecule_structures") or {}).get("canonical_smiles")
    return chembl_id, smiles



# RDKit helper

def rdkit_from_smiles(smiles, seed_name):
    # Parse a SMILES string into an RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        # RDKit failed to parse the SMILES
        return None

    # Set the molecule name for identification in the SDF
    mol.SetProp("_Name", seed_name)

    # Compute 2D coordinates for visualization and consistency
    AllChem.Compute2DCoords(mol)
    return mol


# Main logic
def main(in_tsv, out_sdf, out_report):
    # Store resolution results for the report
    resolved_rows = []

    # Store successfully resolved RDKit molecules
    mols = []

    # Read input TSV containing known inhibitors
    with open(in_tsv, newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            # Extract and normalize input fields
            seed_name = row["seed_name"].strip()
            qtype = row["query_type"].strip().lower()
            query = row["query"].strip()

            chembl_id = ""
            smiles = ""
            status = "OK"

            try:
                # Resolve the query depending on its type
                if qtype == "name":
                    # Resolve by compound name via ChEMBL
                    chembl_id, smiles = chembl_search_by_name(query)
                    if chembl_id is None or smiles is None:
                        status = "NOT_FOUND"

                elif qtype == "chembl_id":
                    # Resolve by ChEMBL identifier
                    chembl_id, smiles = chembl_get_by_id(query)
                    if smiles is None:
                        status = "MISSING_SMILES"

                elif qtype == "smiles":
                    # Directly use the provided SMILES string
                    smiles = query
                    chembl_id = ""

                else:
                    # Unsupported query type
                    status = f"BAD_QUERY_TYPE:{qtype}"

                # If resolution was successful, try to build an RDKit molecule
                if status == "OK":
                    mol = rdkit_from_smiles(smiles, seed_name)
                    if mol is None:
                        status = "RDKIT_PARSE_FAIL"
                    else:
                        # Annotate molecule with metadata
                        mol.SetProp("seed_name", seed_name)
                        mol.SetProp("is_seed", "true")
                        mol.SetProp("canonical_smiles", smiles)
                        if chembl_id:
                            mol.SetProp("chembl_id", chembl_id)
                        mols.append(mol)

            except Exception as e:
                # Catch any unexpected error during resolution
                status = f"ERROR:{type(e).__name__}"

            # Record the resolution outcome for this seed
            resolved_rows.append({
                "seed_name": seed_name,
                "query_type": qtype,
                "query": query,
                "chembl_id": chembl_id,
                "canonical_smiles": smiles,
                "status": status,
            })

    # Write successfully resolved molecules to an SDF file
    writer = Chem.SDWriter(out_sdf)
    for mol in mols:
        writer.write(mol)
    writer.close()

    # Write resolution report as a TSV file
    with open(out_report, "w", newline="", encoding="utf-8") as f:
        fieldnames = [
            "seed_name",
            "query_type",
            "query",
            "chembl_id",
            "canonical_smiles",
            "status",
        ]
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(resolved_rows)



if __name__ == "__main__":
    # Parse command-line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("--in_tsv", required=True, help="known_inhibitors.tsv")
    parser.add_argument("--out_sdf", required=True, help="Output SDF (seeds)")
    parser.add_argument("--out_report_tsv", required=True, help="Resolution report TSV")
    args = parser.parse_args()

    # Run the seed resolution pipeline
    main(args.in_tsv, args.out_sdf, args.out_report_tsv)

