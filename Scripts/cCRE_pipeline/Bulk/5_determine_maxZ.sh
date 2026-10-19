#!/bin/bash

mode=$1
signalDir=/outputs/signal-output
scriptDir=/Scripts
outDir=/outputs/maxZ

cd "${signalDir}/${mode}/"

txt_files=(*.txt)

echo "Determining maxZ across ${#txt_files[@]} files for ${mode}..."
python "${scriptDir}/select-max-zscore.py" "${txt_files[@]}" \
    > "${outDir}/colon-${mode}-maxZ.txt"

echo "Done. Output: ${outDir}/colon-${mode}-maxZ.txt"
