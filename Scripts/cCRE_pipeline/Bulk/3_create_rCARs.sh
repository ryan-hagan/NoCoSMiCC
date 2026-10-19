#!/bin/bash

# Define directories
dir=/outputs/Bulk-colon-rCARs
scriptDir=/Scripts

cd $dir/Processed-CARs

# Calculate the 10th percentile of average signal in DNase regions
paste output.ENCFF291LZE output.ENCFF299OOV output.ENCFF405NTZ output.ENCFF452PQB output.ENCFF613BDA output.ENCFF753MXL | awk '{ print $6}' > merged_DNase_signal.txt
sort -n merged_DNase_signal.txt  | awk '{all[NR] = $0} END{print all[int(NR*0.1 - 0.5)]}' # 0.111376

# Calculate the 10th percentile of average signal in ATAC regions
paste output.ENCFF668GUI output.ENCFF033RPN output.ENCFF811ERB output.ENCFF509NRA output.ENCFF057BIJ output.ENCFF961XDO \
output.ENCFF784HME output.ENCFF049WJI output.ENCFF796DRU | awk '{ print $6}' > merged_ATAC_signal.txt
sort -n merged_ATAC_signal.txt  | awk '{all[NR] = $0} END{print all[int(NR*0.1 - 0.5)]}' # 2.2947

## Process the DNase data
echo -e "Combining DNase peaks..."
cat output.ENCFF291LZE output.ENCFF299OOV output.ENCFF405NTZ output.ENCFF452PQB output.ENCFF613BDA output.ENCFF753MXL > tmp
mv tmp $dir/DNase-All.bed

# Set cutoff
dnase_cutoff=0.111376

# Filter peaks
echo -e "Filtering DNase peaks..."
cd $dir
awk '{if ($1 !~ /_/ && $3-$2 >= 50 && $6 >= '$dnase_cutoff') print $0}' \
    DNase-All.bed | grep -v "chrEBV" | grep -v "chrM" > DNase-Filtered.bed

mkdir scratch
cp DNase-Filtered.bed scratch/tmp.bed
cd scratch

# Sort peaks
echo -e "Sorting DNase peaks..."
sort -k1,1 -k2,2n tmp.bed > sorted
rm -f rPeaks
num=$(wc -l sorted | awk '{print $1}')

# Merge peaks
echo -e "Merging DNase peaks..."
while [ $num -gt 0 ]
do
    echo -e "\t" $num
    bedtools merge -i sorted -c 4,6 -o collapse,collapse > merge
    python $scriptDir/pick.best.peak.py merge > peak-list
    awk 'FNR==NR {x[$1];next} ($4 in x)' peak-list sorted >> rPeaks
    bedtools intersect -v -a sorted -b rPeaks > remaining
    mv remaining sorted
    num=$(wc -l sorted | awk '{print $1}')
done

# Tidy
mv rPeaks ../tmp.bed
cd ../
rm -r scratch

sort -k1,1 -k2,2n tmp.bed > sorted.dnase.bed

awk '{print $1 "\t" $2 "\t" $3 "\t" $4}' sorted.dnase.bed> rDNase-Filtered.bed

## Process the ATAC data
echo -e "Combining ATAC peaks..."
cd $dir/Processed-CARs
cat output.ENCFF668GUI output.ENCFF033RPN output.ENCFF811ERB output.ENCFF509NRA output.ENCFF057BIJ output.ENCFF961XDO \
output.ENCFF784HME output.ENCFF049WJI output.ENCFF796DRU > tmp
mv tmp $dir/BulkATAC-All.bed

# Set cutoff
atac_cutoff=2.2947 #custom cutoff

# Filter peaks
echo -e "Filtering ATAC peaks..."
cd $dir
awk '{if ($1 !~ /_/ && $3-$2 >= 50 && $6 >= '$atac_cutoff') print $0}' \
    BulkATAC-All.bed | grep -v "chrEBV" | grep -v "chrM" > BulkATAC-Filtered.bed

mkdir scratch
cp BulkATAC-Filtered.bed scratch/tmp.bed
cd scratch

# Sort peaks
echo -e "Sorting ATAC peaks..."
sort -k1,1 -k2,2n tmp.bed > sorted
rm -f rPeaks
num=$(wc -l sorted | awk '{print $1}')

# Merge peaks
echo -e "Merging ATAC peaks..."
while [ $num -gt 0 ]
do
    echo -e "\t" $num
    bedtools merge -i sorted -c 4,6 -o collapse,collapse > merge
    python $scriptDir/pick.best.peak.py merge > peak-list
    awk 'FNR==NR {x[$1];next} ($4 in x)' peak-list sorted >> rPeaks
    bedtools intersect -v -a sorted -b rPeaks > remaining
    mv remaining sorted
    num=$(wc -l sorted | awk '{print $1}')
done

mv rPeaks ../tmp.bed
cd ../
rm -r scratch

sort -k1,1 -k2,2n tmp.bed > sorted.atac.bed

awk '{print $1 "\t" $2 "\t" $3 "\t" $4}' sorted.atac.bed > Bulk-rATAC-Filtered.bed

## Combine DNase and ATAC peaks 

cat rDNase-Filtered.bed Bulk-rATAC-Filtered.bed | sort -k1,1 -k2,2n | mergeBed -i stdin > Bulk-rCARs.bed

## Add identifier label to the cCREs
awk '{printf "%s\t%s\t%s\tBulkcCRE%d\n", $1, $2, $3, NR}' Bulk-rCARs.bed > Bulk-rCARs-labelled.bed
