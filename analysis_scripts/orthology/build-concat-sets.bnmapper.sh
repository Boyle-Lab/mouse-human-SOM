#!/bin/bash

# Build various sets of concatenated peaks for downstream analyses.

IN_DIR=$1
OUT_DIR=$2

# By species and cell type
for SPP in hg19 mm9; do

    if [ $SPP == hg19 ]; then

	for CELL in K562 GM12878 ; do
	    grep $CELL $OUT_DIR/ALL.$SPP.native.bed > $OUT_DIR/$SPP.$CELL.native.bed
	    grep $CELL $OUT_DIR/ALL.$SPP.lifted.bed > $OUT_DIR/$SPP.$CELL.lifted.bed
	done

    else

	for CELL in CH12 MEL ; do
	    grep $CELL $OUT_DIR/ALL.$SPP.native.bed > $OUT_DIR/$SPP.$CELL.native.bed
	    grep $CELL $OUT_DIR/ALL.$SPP.lifted.bed > $OUT_DIR/$SPP.$CELL.lifted.bed
	done
    fi
done
