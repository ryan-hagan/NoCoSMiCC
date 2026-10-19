#!/bin/bash

# Define file with sample info
peaks=/files/ENCODE-bulk-peaks.txt

# Set path
dataDir=/outputs/Bulk-colon-rCARs
cd $dataDir

# Iterate over each line in the sample info file
while IFS=$'\t' read -r exp peaks _; do

    # Check if the output file already exists
    if [ ! -f "$dataDir/bulkCARs/$peaks.CARs.bed" ]; then

        tmp_dir=$(mktemp -d -t nocosmicc-XXXXXX)
        cd "$tmp_dir" || exit 1

        wget "https://www.encodeproject.org/files/$peaks/@@download/$peaks.bed.gz"

        gunzip "$peaks.bed.gz"

        awk '{if ($3-$2 >= 50 && $3-$2 <= 1000) print $0}' "$peaks.bed" > peaks

        mkdir -p "$dataDir/bulkCARs"

        mv peaks "$dataDir/bulkCARs/$peaks.CARs.bed"

        rm -r "$tmp_dir"

    fi

done < "$peaks"
