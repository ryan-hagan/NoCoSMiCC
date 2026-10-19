#!/bin/bash

# Define directories
scriptDir=/Scripts
dataDir=/outputs/Bulk-colon-rCARs
peaks=/files/ENCODE-bulk-peaks.txt
Bigwigs=/files/Bigwigs
atac_signal="${Bigwigs}/BulkATAC/"
dnase_signal="${Bigwigs}/DNase/"

while IFS=$'\t' read -r dset dpeak dsig _; do

    # Create a temp directory
    tmp_dir=$(mktemp -d -t nocosmicc-XXXXXX)
    cd "$tmp_dir" || exit 1

    # Download the BigWig signal file if not local
    if [ -f "${atac_signal}${dsig}.bigWig" ]; then
        signal="${atac_signal}${dsig}.bigWig"
    elif [ -f "${dnase_signal}${dsig}.bigWig" ]; then
        signal="${dnase_signal}${dsig}.bigWig"
    else
        wget https://www.encodeproject.org/files/$dsig/@@download/$dsig.bigWig
        signal=$dsig.bigWig
    fi

    # Prepare the chromatin accessible region (CAR) file
    car=$dataDir/bulkCARs/$dpeak.CARs.bed
    cp "$car" bed
    awk '{print $1 "\t" $2 "\t" $3 "\t" "'"$dpeak"'-"NR "\t" $4}' bed | sort -k4,4 > new
    awk '{print $1 "\t" $2 "\t" $3 "\t" $4}' new > new.bed
    bigWigAverageOverBed "$signal" new.bed out.tab
    sort -k1,1 out.tab > tmp.1
    paste new tmp.1 | awk '{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $10}' > output.$dsig

    # Move the output to the data dir
    mkdir -p "$dataDir/Processed-CARs/"
    mv output.$dsig "$dataDir/Processed-CARs/"

    # Clean up the temporary dir
    rm -r "$tmp_dir"

done < "$peaks"
