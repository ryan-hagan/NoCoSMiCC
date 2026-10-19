#!/bin/bash

dataDir=/outputs
rdhs=/outputs/Bulk-colon-rCARs/Bulk-rCARs-labelled.bed
scriptDir=/Scripts
filesDir=/files

tss=$filesDir/hg38_tss.bed
prox=$filesDir/hg38_tss_4k.bed
chromSizes=$filesDir/hg38.chrom.sizes

# Re-sort rCAR BED
rdhs_sorted=$(dirname $rdhs)/rdhs_sorted.bed
bedtools sort -g $chromSizes -i $rdhs > $rdhs_sorted
rdhs=$rdhs_sorted

cd $dataDir/maxZ

# ---- Build TSS proximity files ----
echo "Building TSS files..."
bedtools sort -i $tss -g $chromSizes > tss_sorted.bed
bedtools slop -i tss_sorted.bed -g $chromSizes -b 2000 > $prox
echo "Done."

# ---- Split cCREs into groups ----
echo "Splitting cCREs into groups..."
awk '$2 > 1.28 {print $1}' colon-BulkATAC-maxZ.txt colon-DNase-maxZ.txt | sort -u > list
awk 'FNR==NR {x[$1];next} ($4 in x)' list $rdhs > bed

bedtools intersect -u -a bed -b tss_sorted.bed > tss
bedtools intersect -v -a bed -b tss_sorted.bed > a1
bedtools intersect -u -a a1 -b $prox | bedtools sort -g $chromSizes -i - > prox
bedtools intersect -v -a bed -b $prox > distal

bedtools closest -d -g $chromSizes -a prox -b tss_sorted.bed > tmp
python $scriptDir/calculate-center-distance.py tmp agnostic > new
awk '{if ($2 >= -200 && $2 <= 200) print $0}' new > center-distance
awk '{if ($2 < -2000 || $2 > 2000) print $0}' new > far
awk 'FNR==NR {x[$1];next} ($4 in x)' center-distance prox >> tss
awk 'FNR==NR {x[$1];next} ($4 in x)' far prox >> distal
cat center-distance far > new
awk 'FNR==NR {x[$1];next} !($4 in x)' new prox > tmp
mv tmp prox
rm new

# ---- TSS elements ----
awk 'FNR==NR {x[$4];next} ($1 in x)' tss colon-H3K4me3-maxZ.txt | \
    awk '{if ($2 > 1.28) print $0}' > PLS

awk 'FNR==NR {x[$4];next} ($1 in x)' tss colon-H3K4me3-maxZ.txt | \
    awk '{if ($2 <= 1.28) print $0}' > no1

awk 'FNR==NR {x[$1];next} ($1 in x)' no1 colon-H3K27ac-maxZ.txt | \
    awk '{if ($2 > 1.28) print $0}' > pELS

awk 'FNR==NR {x[$1];next} ($1 in x)' no1 colon-H3K27ac-maxZ.txt | \
    awk '{if ($2 <= 1.28) print $0}' > no2

awk 'FNR==NR {x[$1];next} ($1 in x)' no2 colon-CTCF-maxZ.txt | \
    awk '{if ($2 > 1.28) print $0}' > CTCFonly

awk 'FNR==NR {x[$1];next} ($1 in x)' no2 colon-CTCF-maxZ.txt | \
    awk '{if ($2 <= 1.28) print $0}' > CAonly

# ---- Proximal elements ----
awk 'FNR==NR {x[$4];next} ($1 in x)' prox colon-H3K27ac-maxZ.txt | \
    awk '{if ($2 > 1.28) print $0}' >> pELS

awk 'FNR==NR {x[$4];next} ($1 in x)' prox colon-H3K27ac-maxZ.txt | \
    awk '{if ($2 <= 1.28) print $0}' > no1

awk 'FNR==NR {x[$1];next} ($1 in x)' no1 colon-H3K4me3-maxZ.txt | \
    awk '{if ($2 > 1.28) print $0}' > DNaseK4

awk 'FNR==NR {x[$1];next} ($1 in x)' no1 colon-H3K4me3-maxZ.txt | \
    awk '{if ($2 <= 1.28) print $0}' > no2

awk 'FNR==NR {x[$1];next} ($1 in x)' no2 colon-CTCF-maxZ.txt | \
    awk '{if ($2 > 1.28) print $0}' >> CTCFonly

awk 'FNR==NR {x[$1];next} ($1 in x)' no2 colon-CTCF-maxZ.txt | \
    awk '{if ($2 <= 1.28) print $0}' >> CAonly

# ---- Distal elements ----
awk 'FNR==NR {x[$4];next} ($1 in x)' distal colon-H3K27ac-maxZ.txt | \
    awk '{if ($2 > 1.28) print $0}' > dELS

awk 'FNR==NR {x[$4];next} ($1 in x)' distal colon-H3K27ac-maxZ.txt | \
    awk '{if ($2 <= 1.28) print $0}' > no1

awk 'FNR==NR {x[$1];next} ($1 in x)' no1 colon-H3K4me3-maxZ.txt | \
    awk '{if ($2 > 1.28) print $0}' >> DNaseK4

awk 'FNR==NR {x[$1];next} ($1 in x)' no1 colon-H3K4me3-maxZ.txt | \
    awk '{if ($2 <= 1.28) print $0}' > no2

awk 'FNR==NR {x[$1];next} ($1 in x)' no2 colon-CTCF-maxZ.txt | \
    awk '{if ($2 > 1.28) print $0}' >> CTCFonly

awk 'FNR==NR {x[$1];next} ($1 in x)' no2 colon-CTCF-maxZ.txt | \
    awk '{if ($2 <= 1.28) print $0}' >> CAonly

# ---- Accessioning cCREs ----
echo "Accessioning cCREs..."
awk 'FNR==NR {x[$1];next} ($4 in x)' PLS      $rdhs | \
    awk '{print $1"\t"$2"\t"$3"\t"$4"\t""PLS"}'         > l.bed
awk 'FNR==NR {x[$1];next} ($4 in x)' pELS     $rdhs | \
    awk '{print $1"\t"$2"\t"$3"\t"$4"\t""pELS"}'        >> l.bed
awk 'FNR==NR {x[$1];next} ($4 in x)' dELS     $rdhs | \
    awk '{print $1"\t"$2"\t"$3"\t"$4"\t""dELS"}'        >> l.bed
awk 'FNR==NR {x[$1];next} ($4 in x)' DNaseK4  $rdhs | \
    awk '{print $1"\t"$2"\t"$3"\t"$4"\t""CA-H3K4me3"}'  >> l.bed
awk 'FNR==NR {x[$1];next} ($4 in x)' CTCFonly $rdhs | \
    awk '{print $1"\t"$2"\t"$3"\t"$4"\t""CA-CTCF"}'     >> l.bed
awk 'FNR==NR {x[$1];next} ($4 in x)' CAonly   $rdhs | \
    awk '{print $1"\t"$2"\t"$3"\t"$4"\t""CA"}'          >> l.bed


awk 'FNR==NR {x[$4];next} ($1 in x)' l.bed colon-CTCF-maxZ.txt | \
    awk '{if ($2 > 1.28) print $0}' > CTCFall

awk 'FNR==NR {x[$1];next} ($4 in x)' CTCFall l.bed | \
    awk '{print $0",CTCF-bound"}' > m.bed
awk 'FNR==NR {x[$1];next} !($4 in x)' CTCFall l.bed | \
    awk '{print $0}' >> m.bed

sort -k1,1 -k2,2n m.bed > l.bed

sort -k1,1 -k2,2n l.bed > $dataDir/colon-bulk-cCREs.bed 
echo "Done."

# ---- Summary ----
echo "Classification counts:"
awk '{print $5}' $dataDir/colon-bulk-cCREs.bed | sort | uniq -c | sort -rn
wc -l $dataDir/colon-bulk-cCREs.bed 
