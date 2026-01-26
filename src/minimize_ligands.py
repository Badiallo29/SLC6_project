#!/usr/bin/env python3
import argparse
from rdkit import Chem
from rdkit.Chem import AllChem


def minimize_uff(mol, max_iters=200):
    # Perform energy minimization using the UFF force field
    try:
        # Initialize the UFF force field for the molecule
        ff = AllChem.UFFGetMoleculeForceField(mol)
        ff.Initialize()

        # Run the minimization for a maximum number of iterations
        ff.Minimize(maxIts=max_iters)
        return True
    except Exception:
        # Return False if force-field setup or minimization fails
        return False


def main(in_sdf, out_sdf, out_fail_sdf=None, max_iters=200):
    # Load molecules from the input SDF file (keep hydrogens if present)
    suppl = Chem.SDMolSupplier(in_sdf, removeHs=False)

    # Writer for successfully minimized molecules
    w_ok = Chem.SDWriter(out_sdf)

    # Optional writer for molecules that fail minimization
    w_fail = Chem.SDWriter(out_fail_sdf) if out_fail_sdf else None

    # Counters for successful and failed minimizations
    ok = 0
    fail = 0

    # Iterate over molecules in the input SDF
    for mol in suppl:
        if mol is None:
            # Skip invalid or unreadable molecules
            fail += 1
            continue

        # Ensure explicit hydrogens are present before minimization
        if mol.GetNumAtoms() == Chem.RemoveHs(mol).GetNumAtoms():
            mol = Chem.AddHs(mol, addCoords=True)

        # Run UFF energy minimization
        success = minimize_uff(mol, max_iters=max_iters)

        if success:
            # Write successfully minimized molecule
            w_ok.write(mol)
            ok += 1
        else:
            # Handle minimization failure
            fail += 1
            if w_fail:
                w_fail.write(mol)

    # Close output writers
    w_ok.close()
    if w_fail:
        w_fail.close()

    # Print a summary of the minimization step
    print(f"[MIN] OK={ok} FAIL={fail} (max_iters={max_iters})")


if __name__ == "__main__":
    # Parse command-line arguments
    ap = argparse.ArgumentParser()
    ap.add_argument("--in_sdf", required=True)         # Input SDF file
    ap.add_argument("--out_sdf", required=True)        # Output SDF for minimized molecules
    ap.add_argument("--out_fail_sdf", default=None)    # Optional SDF for failed molecules
    ap.add_argument("--max_iters", type=int, default=200)  # Max minimization iterations
    args = ap.parse_args()

    # Run the minimization workflow
    main(args.in_sdf, args.out_sdf, args.out_fail_sdf, args.max_iters)

