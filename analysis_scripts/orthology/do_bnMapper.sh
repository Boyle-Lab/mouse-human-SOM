#!/bin/bash

# Do bnMapper.py mapping between species using the UCSC liftOver chains.

IFS=$'\n'
x="$(date +"%s")"
x=$(( $x * $RANDOM ))
TMP_F="$x.tmp"

DATA_DIR=$1
FN_ROOT=$2
CHAIN_DIR=$3
OUT_DIR=$4

if [ ! -z "$5" ]; then
    RBEST=$5  # Use reciprocal best chains
else
    RBEST=0
fi

if [ ! -d $OUT_DIR/lifted ]; then
    mkdir -pv $OUT_DIR/lifted
fi

if [ ! -d $OUT_DIR/unmapped ]; then
    mkdir -pv $OUT_DIR/unmapped
fi

# bnMapper options...
THRESH=0.1   # --threshold - required overlap to call a match
GAP=-1       # --gap - reject matches with gaps larger than this (-1 to ignore)

for FILE in $DATA_DIR/$FN_ROOT.{hg19,mm9}.bed ; do
# Test Loop
#for FILE in $DATA_DIR/a_C.out.mm9.bed ; do

    >&2 echo
    >&2 echo "Working on $FILE..."
    
    SP1=$(echo $FILE | awk '/mm9/{ print "mm9" }')

    if [ -z $SP1 ]; then
	SP1=$(echo $FILE | awk '/hg19/{ print "hg19" }')
    fi

    if [ $SP1 == "mm9" ]; then
	SP2="hg19"
    else
	SP2="mm9"
    fi
#    echo $FILE $SP1  $SP2
    
    # First run outputs the actual mapping
    >&2 echo
    >&2 echo "Mapping coordinates from $SP1 to $SP2 reference..."

    if [ $RBEST = 0 ]; then
	CHAIN=$CHAIN_DIR/$SP1.$SP2.over.chain.gz
	>&2 echo "Using liftOver chains"
    else
	CHAIN=$CHAIN_DIR/$SP1.$SP2.rbest.chain.gz
	>&2 echo "Using reciprocal-best chains"
    fi

    python /home/adadiehl/james_taylor-bx-python-da37e3aa45dc/scripts/bnMapper.py $FILE $CHAIN -f BED12 -t $THRESH -g $GAP > $OUT_DIR/lifted/$(basename $FILE .bed).lifted.bed

    >&2 echo
    >&2 echo "Finding unmapped features..."
    # Find the unmapped features in the input file based on the feature name.
    # WARNING: This assumes each SOM peak has a unique name in the bed file!
    
    for line in $(cat $OUT_DIR/lifted/$(basename $FILE .bed).lifted.bed) ; do

	LNAME=$(echo $line | awk '{print $4}')
	# exact_grep is a simple, grep-like program that searches for
	# exact word matches (using c's standard strcmp() function) within a
	# line. It is available from https://github.com/adadiehl/exact_grep
	exact_grep "$LNAME" $FILE >> $TMP_F

    done

    diff --new-line-format="" --unchanged-line-format="" <(sort $FILE) <(sort $TMP_F) > $OUT_DIR/unmapped/$(basename $FILE .bed).unmapped.bed
#    rm $TMP_F

    # Save the native sequences -- these may be needed for intersectPhastCons.sh
    mv $TMP_F $OUT_DIR/lifted/$(basename $FILE .bed).native.bed
  
    # Sanity check...
    N_ORIG=$(grep -v "#" $FILE | wc -l)
    N_LIFTED=$(grep -v "#" $OUT_DIR/lifted/$(basename $FILE .bed).lifted.bed | wc -l)
    N_UNMAPPED=$(grep -v "#" $OUT_DIR/unmapped/$(basename $FILE .bed).unmapped.bed | wc -l)
    N_TOTAL=$(($N_LIFTED + $N_UNMAPPED))

    if [ $N_TOTAL != $N_ORIG ]; then
	>&2 echo "WARNING: Total number of lines between lifted and unmapped fractions does not equal original number of lines!"
    fi

done
