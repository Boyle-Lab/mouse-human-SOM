#!/bin/bash

# Location of genomic fasta files
FASTA_DIR=/data/UCSC/FASTA

# Get fasta gene sequences for the UCSC knownGenes
WD=$(pwd)
for SPP in hg19 mm9; do
    cd ../genes 
    gffread -w knownGenes.$SPP.fa -g $FASTA_DIR/$SPP.fa knownGenes.$SPP.gtf
    gzip knownGenes.$SPP.fa
done
cd $WD

# UCSC knowngenes based transcript indeces
for SPP in hg19 mm9 ; do kallisto index -i transcripts.$SPP.idx ../genes/knownGenes.$SPP.fa.gz ; done

# Run kallisto
mkdir -pv {hg19/{GM12878,K562},mm9/{CH12,MEL}}
./run_kallisto.sh . ../ENCODE_data
# ENCODE_data should have structure ENCODE_data/{hg19/{K562,GM12878},mm9/{MEL,CH12}}
# with terminal directories containing fastq files for all replicates in the given
# cell type

# Run Sleuth
./do_sleuth.Rscript

# Add gene names
./rename_genes.pl norm_genes.hg19.tab ../genes/ucsc2geneId.hg19.tab --t-col 0 --r-col 4 > norm_genes.hg19.genes.tab
./rename_genes.pl norm_genes.mm9.tab ../genes/ucsc2geneId.mm9.tab --t-col 0 --r-col 4 > norm_genes.mm9.genes.tab

# Total up isoforms...
./total_isoforms.pl norm_genes.hg19.genes.tab > norm_genes.hg19.total.tab
./total_isoforms.pl norm_genes.mm9.genes.tab > norm_genes.mm9.total.tab

# Prepare composite data table
head -n 1 norm_genes.hg19.total.tab | sed 's/#//' | awk -v spp="species" '{printf "#%s\t%s\n", spp, $0}' > ALL.norm_genes.tab
for SPP in hg19 mm9 ; do grep -v "^#" norm_genes.$SPP.total.tab | awk -v spp=$SPP '{printf "%s\t%s\n", spp, $0}' >> ALL.norm_genes.tab ; done
./build_kallisto_table.pl ALL.norm_genes.tab ../genes/hg19-mm9.1-1-orthologs.txt > ALL.norm_genes.dat

# Because R doesn't like headers commented with '#'...
sed 's/^#//' ALL.norm_genes.dat > tmp

# Prepare quantile-normalized version (reads and writes data from tmp)
./do_qnorm.Rscript
head -n1 ALL.norm_genes.dat > ALL.norm_genes.qnorm.dat
tail -n $(( $(cat tmp | wc -l) - 1 )) tmp >> ALL.norm_genes.qnorm.dat
rm tmp
