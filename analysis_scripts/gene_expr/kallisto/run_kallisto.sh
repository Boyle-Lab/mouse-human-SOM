#!/bin/bash

# Generate rna-seq quantifications with kallisto quant

BASEDIR=$1    # Working directory -- should contain kallisto transcript indeces
DATA_DIR=$2   # Location of fastq input files

###
# Run kallisto quant for all species/cells
for SPP in hg19 mm9; do
    >&2 echo

    if [ "$SPP" = "hg19" ]; then
        CELLS=("GM12878" "K562")
    else
        CELLS=("CH12" "MEL")
    fi


    for CELL in ${CELLS[@]}; do
        >&2 echo "Running kallisto quant on $SPP $CELL..."

	i=1
	for file in $(ls $DATA_DIR/$SPP/$CELL/*.fastq.gz); do
	    >&2 echo "    Replicate $i..."

	    OUT=$(basename $file .fasta.gz)

	    kallisto quant -i transcripts.$SPP.idx -o $BASEDIR/$SPP/$CELL/$OUT -b 100 $file --single -l 200 -s 80 -t 8
	    ((i++))
	done

        >&2 echo "Done with $SPP $CELL."
	>&2 echo
    done
    
done
>&2 echo "Done!"
