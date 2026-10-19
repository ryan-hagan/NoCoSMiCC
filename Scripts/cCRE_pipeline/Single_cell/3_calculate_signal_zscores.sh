#!/bin/bash

mode=$1
files=/files/${mode}-List_colon.txt
output=/outputs/sc-signal-output
scriptDir=/Scripts
peaks=/outputs/colon-epithelial-rCARs/Filtered-CARs/scATAC-rCAR-labelled.bed
bigwigDir=/files/Bigwigs

mkdir -p "$output"

# Set widths
if [ "$mode" == "CTCF" ] || [ "$mode" == "scATAC" ]; then
    width=0
else
    width=500
fi

echo "Mode: $mode"

while IFS=$'\t' read -r dset dsig _; do

    echo "  Processing $dset / $dsig..."

    # Temp working directory
    tmp_dir=$(mktemp -d -t nocosmicc-XXXXXX)
    cd "$tmp_dir" || exit 1

    awk -F "\t" -v w="$width" '{printf "%s\t%.0f\t%.0f\t%s\n", $1, $2-w, $3+w, $4}' "$peaks" \
        | awk '{if ($2 < 0) print $1 "\t" 0 "\t" $3 "\t" $4; else print $0}' \
        | sort -u > little

    # Locate bigWig
    if [ -f "${bigwigDir}/${mode}/${dsig}.bw" ]; then
    bw="${bigwigDir}/${mode}/${dsig}.bw"
    elif [ -f "${bigwigDir}/${mode}/${dsig}.bigWig" ]; then
    bw="${bigwigDir}/${mode}/${dsig}.bigWig"
    fi
    
    # Extract signal
    bigWigAverageOverBed -bedOut=out2.bed "$bw" little out2

    # Log zscore normalisation
    python "$scriptDir/log.zscore.normalization.py" out2 > l

    # Sort and write output
    sort -k2,2rg l \
        | awk 'BEGIN {rank=0; before=0; running=1} {
            if ($2 != before) rank = running
            print $1 "\t" $2 "\t" $3 "\t" rank
            before=$2
            running += 1}' | sort -k1,1 > "${output}/$mode/${dset}-${dsig}.txt"

    echo "  Done -> ${output}/$mode/${dset}-${dsig}.txt"

    cd /tmp && rm -rf "$tmp_dir"

done < "$files"
