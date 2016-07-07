#!/bin/bash

# Gather raw counts from the rna-seq datasets using featureCounts.

DATA_DIR=$1   # Path to data directory
GENES_DIR=$2  # Directory containing gtf gene annotations

for SPP in hg19 mm9; do
    >&2 echo

    if [ "$SPP" = "hg19" ]; then
        CELLS=("GM12878" "K562")
    else
        CELLS=("CH12" "MEL")
    fi

    for CELL in ${CELLS[@]}; do
        >&2 echo "Running featureCounts on $SPP $CELL..."
	>&2 echo

	FILES=($DATA_DIR/$CELL/*.bam)

	featureCounts -T 8 -t exon -g gene_id -a $GENES_DIR/knownGenes.$SPP.by-symbol.gtf -o $DATA_DIR/$SPP/$CELL/raw_counts.txt ${FILES[0]} ${FILES[1]}
	
        >&2 echo "Done with $SPP $CELL."
	>&2 echo
    done

done

>&2 echo "Done!"
