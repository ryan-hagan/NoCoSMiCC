#!/bin/bash

scriptDir=/Scripts
signalDir=/outputs/sc-signal-output
#outDir=/outputs

for mode in scATAC H3K4me3 H3K27ac CTCF; do
    echo "Determining maxZ for $mode..."
    cd "$signalDir/$mode"
    python "$scriptDir/select-max-zscore.py" *.txt \
        > "$signalDir/$mode-maxZ.txt"
    echo "  Done: $signalDir/$mode-maxZ.txt"
done

echo "All complete."

# Verify
wc -l "$signalDir"/*-maxZ.txt
