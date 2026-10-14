#!/bin/bash

# Set base dir
base_dir=. 

# Set paths
dir=$base_dir/outputs/colon-epithelial-rCARs
outdir=$dir/Filtered-CARs

mkdir -p "$outdir"
cd "$dir/Processed-CARs"

# Combine all samples
echo "Combining CARs..."
cat output.sc.* > "$dir/scATAC-CAR-All.bed"

# 10th percentile cutoff
cutoff=0.00256986

# Filter
echo "Filtering CARs..."
awk -v cutoff="$cutoff" 'BEGIN{OFS="\t"} $6 >= cutoff' \
    "$dir/scATAC-CAR-All.bed" | \
sort -k1,1 -k2,2n > "$outdir/scATAC-CAR-Filtered-All.bed"

echo "Selecting best peak..."
sort -k4,4 -k6,6rn "$outdir/scATAC-CAR-Filtered-All.bed" | \
    sort -k4,4 -u | \
    sort -k1,1 -k2,2n \
    > "$outdir/scATAC-rCAR-Filtered.bed"

echo "Adding labels..."
awk 'BEGIN{OFS="\t"} {printf "%s\t%s\t%s\tColoncCRE%d\n", $1, $2, $3, NR}' \
    "$outdir/scATAC-rCAR-Filtered.bed" \
    > "$outdir/scATAC-rCAR-labelled.bed"

echo "Done."
wc -l "$outdir/scATAC-rCAR-Filtered.bed"
wc -l "$outdir/scATAC-rCAR-labelled.bed"
