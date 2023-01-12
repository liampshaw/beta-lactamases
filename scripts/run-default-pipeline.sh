#!/bin/bash
# Runs the default pipeline for a given beta-lactamase gene variant

GENE=$1
FASTADIR=$2
OUTPUTDIR=$3

FAMILY=$(echo $GENE | rev | cut -d '-' -f 2- | rev )
THRESHOLD="25"
#GENE_SEQUENCE=$(grep -A 1 $GENE"." gene-lists.fa | tail -n 1)
GENE_SEQUENCE=$(grep -A 1 $GENE" " ../data/nucleotide_fasta_protein_homolog_model.fasta  | tail -n 1)
echo $GENE
echo $FAMILY
echo $GENE_SEQUENCE

mkdir -p "$OUTPUTDIR"


python get-gene-cluster-accessions.py $FAMILY $GENE $THRESHOLD --singlehits | cut -d ',' -f 1 | sed -e 's/$/.fa/g' | sed -e 's~^~'$FASTADIR'/~g'> $OUTPUTDIR/accs.txt

#while read f;
#do
#name=$(echo $f | cut -d ',' -f 1)
#cat $FASTADIR/"$name".fa >> $OUTPUTDIR/seqs.fa
#echo $f
#done < $OUTPUTDIR/accs.txt
xargs cat < $OUTPUTDIR/accs.txt > $OUTPUTDIR/seqs.fa

echo ">gene\n" > $OUTPUTDIR/gene.fa
echo $GENE_SEQUENCE >> $OUTPUTDIR/gene.fa
scriptDir='./'

# cat $FASTADIR/N*.fa | awk '$1!=""' > all_seqs.fa # remove empty lines
#rm N*fa # remove original files

# Extract regions
echo "Extracting flanking regions..."
python "$scriptDir"/extract-region.py --gene $OUTPUTDIR/gene.fa --input $OUTPUTDIR/seqs.fa --upstream 5000 --downstream 5000 --complete --output $OUTPUTDIR/all_u5k_d5k
# unsure if --circular is appropriate here - depends on whether plasmids and chromosomes can be trusted to be circular
# The blast hits - some are not full-length. I don't check for this at the moment. It seems it could be caused by issues the with gap penalty etc parameters of blast. If one requests --complete, then there is an implicit check because the extracted regions must be equal to gene_length + flanking_regions

# Align gene sequences
echo "Aligning gene sequences with mafft..."
mafft --quiet $OUTPUTDIR/all_u5k_d5k_focal_gene.fa > $OUTPUTDIR/all_u5k_d5k_focal_gene.aln
# Deduplicate
seqkit rmdup -s < $OUTPUTDIR/all_u5k_d5k_focal_gene.aln > $OUTPUTDIR/all_u5k_d5k_focal_gene.dedup.aln -D $OUTPUTDIR/all_u5k_d5k_focal_gene.dedup.txt
# SNP distances
echo "SNP distances..."
snp-dists -q -m $OUTPUTDIR/all_u5k_d5k_focal_gene.aln   > $OUTPUTDIR/pangraph_all_u5k_d5k.gfa.gene.snps.tsv

# Pangraph
echo "Pangraph..."
pangraph build $OUTPUTDIR/all_u5k_d5k.fa > $OUTPUTDIR/pangraph_all_u5k_d5k.json # what should the minimum block size be?
pangraph export --edge-minimum-length 0 $OUTPUTDIR/pangraph_all_u5k_d5k.json -p pangraph_all_u5k_d5k  -o $OUTPUTDIR
# Optional: polishing (time-intensive) and export
#pangraph polish pangraph_all_u5k_d5k.json > pangraph_all_u5k_d5k_polished.json
# pangraph export --edge-minimum-length 0 pangraph_all_u5k_d5k_polished.json -p pangraph_all_u5k_d5k_polished  -o ./
# pangraph export pangraph_all_u5k_d5k.json -p pangraph_all_u5k_d5k_polished  -o ./ --export-panX --no-export-gfa


# Prepare gfa
echo "Preparing gfa..."
python "$scriptDir"/prepare-pangraph-gfa.py $OUTPUTDIR/pangraph_all_u5k_d5k.gfa
# change the .colours.csv output
# Should use proper lengths of blocks - but this is a minor addition for once the pipeline is done.
# this now also saves the counts of unique full-length paths and a representative of the most frequent path to
# a file

# Find the gene block in pangraph
# (could be incorporated into the prepare pangraph step?
echo "Finding gene block..."
makeblastdb -in $OUTPUTDIR/pangraph_all_u5k_d5k.fa -dbtype 'nucl'
geneBlock=$(blastn -query $OUTPUTDIR/gene.fa -db $OUTPUTDIR/pangraph_all_u5k_d5k.fa -outfmt 6 | cut -f 2)

# Use the compute-distances script
echo "Computing distances..."
python "$scriptDir"/compute-distances.py $OUTPUTDIR/pangraph_all_u5k_d5k.gfa $geneBlock

# Plot the distances
echo "Plotting distances..."
Rscript "$scriptDir"/plot-output-dists.R $OUTPUTDIR/pangraph_all_u5k_d5k.gfa.output_dists.csv $OUTPUTDIR/pangraph_all_u5k_d5k.gfa.most_frequent_path_representative.txt $OUTPUTDIR/all_u5k_d5k_focal_gene.dedup.txt
