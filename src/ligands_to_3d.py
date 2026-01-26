#!/usr/bin/env python3
import argparse
from rdkit import Chem
from rdkit.Chem import AllChem


def main(in_sdf, out_sdf, seed=0):
    # Load molecules from the input SDF file (keep existing hydrogens if any)
    suppl = Chem.SDMolSupplier(in_sdf, removeHs=False)

    # Prepare an SDF writer for 3D-embedded molecules
    writer = Chem.SDWriter(out_sdf)

    # Configure RDKit ETKDG v3 parameters for 3D conformation generation
    params = AllChem.ETKDGv3()
    params.randomSeed = seed      # Ensure reproducibility

    # Counters for successful and failed embeddings
    ok = 0
    fail = 0

    # Iterate over molecules from the input SDF
    for mol in suppl:
        if mol is None:
            # Skip invalid or unreadable molecules
            fail += 1
            continue

        # Add explicit hydrogens and preserve coordinates when possible
        mol = Chem.AddHs(mol, addCoords=True)

        # Generate a 3D conformation using distance geometry
        cid = AllChem.EmbedMolecule(mol, params)
        if cid != 0:
            # Embedding failed
            fail += 1
            continue

        # Write successfully embedded molecule to the output SDF
        writer.write(mol)
        ok += 1

    # Close the SDF writer
    writer.close()

    # Print a summary of the 3D embedding process
    print(f"[3D] OK={ok} FAIL={fail}")


if __name__ == "__main__":
    # Parse command-line arguments
    ap = argparse.ArgumentParser()
    ap.add_argument("--in_sdf", required=True)   # Input SDF file
    ap.add_argument("--out_sdf", required=True) # Output SDF with 3D coordinates
    ap.add_argument("--seed", type=int, default=0)  # Random seed for reproducibility
    args = ap.parse_args()

    # Run the 3D embedding workflow
    main(args.in_sdf, args.out_sdf, args.seed)

