#!/usr/bin/env bash
set -euo pipefail

# -------------------------------
# Config
# -------------------------------
INPUT_ROOT=${1:-data}          # Premier argument ou 'data' par défaut
OUTPUT_ROOT=fpocket_results

# Boucle sur target/template
for TYPE in target template; do
  for CONF in positive negative; do

    # Skip negative for template
    if [[ "$TYPE" == "template" && "$CONF" == "negative" ]]; then
      continue
    fi

    # Parcours des fichiers PDB
    for PDB in "$INPUT_ROOT/$TYPE/$CONF"/*.pdb; do

      NAME=$(basename "$PDB" .pdb)

      if [[ "$CONF" == "positive" ]]; then
        TAG="POS"
      else
        TAG="NEG"
      fi

      echo "-----------------------------------------"
      echo "Running fpocket on $PDB ($TYPE, $TAG)"

      # Aller dans le dossier du PDB pour que fpocket crée le _out ici
      PDB_DIR=$(dirname "$PDB")
      pushd "$PDB_DIR" > /dev/null

      fpocket -f $(basename "$PDB")

      OUT_DIR="${NAME}_out"
      FINAL_DIR="$OLDPWD/$OUTPUT_ROOT/$TYPE/${OUT_DIR}_${TAG}"

      # Si fpocket a généré des poches
      if [[ -d "$OUT_DIR" ]]; then
        mkdir -p "$OLDPWD/$OUTPUT_ROOT/$TYPE"
        mv "$OUT_DIR" "$FINAL_DIR"
        echo "Pockets saved to $FINAL_DIR"
      else
        echo "No pockets detected for $PDB, skipping."
      fi

      popd > /dev/null

    done
  done
done

echo "-----------------------------------------"
echo "Fpocket analysis completed."

