#!/bin/bash

markers=(ACT CAL ITS LNS2 LSU PGK RPB1 RPB2 SSU TEF1 TOP1 TUB2)

for marker in "${markers[@]}"
do
    echo "======================="
    echo "Running ${marker}"
    echo "======================="

    mkdir -p output/publication_scale/no_complexes/${marker}

    python dnabarcoder.py predict \
    -i input/publication_scale_clean/${marker}.fasta \
    -c input/publication_scale_classification/${marker}.classification \
    -st 0.9 -et 1 -s 0.001 \
    -rank species \
    -higherrank genus \
    -ml 400 \
    -prefix ${marker}_no_complexes \
    -o output/publication_scale/no_complexes/${marker}

    if [ $? -ne 0 ]; then
        echo "FAILED: ${marker}" >> failed_markers.txt
        echo "Skipping..."
        continue
    fi

    echo "${marker} finished" >> completed_markers.txt

done

echo "Finished all analyses"
