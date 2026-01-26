#!/usr/bin/env python3
import argparse
import subprocess
from pathlib import Path
from rdkit import Chem

def safe_name(mol, i):
    # Generate a safe and consistent filename for a molecule

    # First, try to use the RDKit internal molecule name (_Name)
    if mol.HasProp("_Name"):
        n = mol.GetProp("_Name").strip()  # Remove leading/trailing spaces
        if n:
            # Replace characters that may cause filesystem issues
            return n.replace(" ", "_").replace("/", "_")

    # If _Name is not available or empty, try a custom property (seed_name)
    if mol.HasProp("seed_name"):
        n = mol.GetProp("seed_name").strip()
        if n:
            # Apply the same sanitization rules
            return n.replace(" ", "_").replace("/", "_")

    # Fallback: generate a unique default ligand name using the index
    return f"ligand_{i:05d}"


def main(in_sdf, out_dir):
    # Create output directory if it does not already exist
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    # Load molecules from the input SDF file (keep hydrogens if present)
    suppl = Chem.SDMolSupplier(in_sdf, removeHs=False)

    # Counters for successful and failed conversions
    ok = 0
    fail = 0

    # Iterate over molecules in the SDF supplier
    for i, mol in enumerate(suppl):
        # Skip invalid or unreadable molecules
        if mol is None:
            fail += 1
            continue

        # Generate a safe base name for the current molecule
        name = safe_name(mol, i)

        # Temporary SDF file and final PDBQT output path
        tmp_sdf = out_dir / f"{name}.sdf"
        out_pdbqt = out_dir / f"{name}.pdbqt"

        # Write the molecule to a temporary SDF file
        w = Chem.SDWriter(str(tmp_sdf))
        w.write(mol)
        w.close()

        try:
            # Convert SDF to PDBQT using Open Babel
            # -h: add hydrogens
            # --partialcharge gasteiger: compute Gasteiger partial charges
            subprocess.run(
                ["obabel", str(tmp_sdf), "-O", str(out_pdbqt), "-h", "--partialcharge", "gasteiger"],
                check=True,
                stdout=subprocess.DEVNULL,  # Suppress standard output
                stderr=subprocess.DEVNULL,  # Suppress error output
            )
            ok += 1  # Conversion succeeded
        except Exception:
            fail += 1  # Conversion failed
        finally:
            # Remove the temporary SDF file to keep the directory clean
            if tmp_sdf.exists():
                tmp_sdf.unlink()

    # Print a summary of the conversion process
    print(f"[PDBQT] OK={ok} FAIL={fail} OUT={out_dir}")


if __name__ == "__main__":
    # Parse command-line arguments
    ap = argparse.ArgumentParser()
    ap.add_argument("--in_sdf", required=True)   # Input SDF file
    ap.add_argument("--out_dir", required=True) # Output directory for PDBQT files
    args = ap.parse_args()

    # Run the main conversion workflow
    main(args.in_sdf, args.out_dir)

