#!/bin/bash

# Prepares gene expression data normalized for batch effects following a
# protocol adapted from Gilad, Y., & Mizrahi-Man, O. (2015). F1000Res, 4, 121.

# Gather raw counts from bam alignments with featureCounts
./get_raw_counts.sh ../ENCODE_data ../genes/
# ENCODE_data should have structure ENCODE_data/{hg19/{K562,GM12878},mm9/{MEL,CH12}}
# and each terminal directory should contain bam alignments for two bioreplicates

./prepare_counts_matrix_complete.pl ../ENCODE_data/hg19/GM12878/raw_counts.txt ../ENCODE_data/hg19/K562/raw_counts.txt ../ENCODE_data/mm9/CH12/raw_counts.txt ../ENCODE_data/mm9/MEL/raw_counts.txt ../genes/hg19-mm9.1-1-orthologs.txt > counts_matrix_complete.txt

./build_gc_table_complete.pl counts_matrix_complete.txt ../genes/hg19-mm9.1-1-orthologs.txt ../genes/knownGenes.hg19.genes.gc ../genes/knownGenes.mm9.genes.gc > gc_table_complete.txt

# Effective gene lengths, from kallisto
for SPP in hg19 mm9 ; do awk -F $'\t' '{printf "%s\t%s\n", $1, $14}' ../kallisto/norm_genes.$SPP.total.tab | sed s/#// > eff_lens.$SPP.dat ; done

# The rest is in R...
./do_combat.Rscript # writes normalized TPM to normCountsComplete.tab
