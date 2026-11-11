#!/bin/bash
# Usage: bash Hybrid_assembly_genome_analysis.sh

ROOT=/storage/group/jzk303/default/data/tyler/Slovenian_Isolates
WORK=/storage/work/tfc5209/Slovenian_BC
TRIM=/scratch/tfc5209/4_trimmed_illumina
ASM=$WORK/6_Unicycler_assembly
MOB=$WORK/14_mob_recon_results
MASH_DIR=$WORK/15_mash_network
FASTQC=/storage/work/tfc5209/programs/FastQC/fastqc
TRIMMOMATIC=/storage/work/tfc5209/programs/Trimmomatic-0.39/trimmomatic-0.39.jar
ADAPTERS=/storage/work/tfc5209/programs/Trimmomatic-0.39/adapters/NexteraPE-PE.fa

# 1. FastQC on raw Illumina reads
mkdir -p $WORK/3_fastqc_Illumina && $FASTQC $ROOT/0_Illumina_reads/00_fastq/*.fastq.gz -o $WORK/3_fastqc_Illumina -t 4

# 2. Trimmomatic on paired-end Illumina reads
cd $ROOT/0_Illumina_reads/00_fastq && mkdir -p $TRIM && for r1 in *_R1_*.fastq.gz; do b=${r1%_R1_*}; [ -f $TRIM/${b}_R1_001.trimmedP.fastq.gz ] || java -jar $TRIMMOMATIC PE -threads 4 -phred33 ${b}_R1_001.fastq.gz ${b}_R2_001.fastq.gz $TRIM/${b}_R1_001.trimmedP.fastq.gz $TRIM/${b}_R1_001.trimmedS.fastq.gz $TRIM/${b}_R2_001.trimmedP.fastq.gz $TRIM/${b}_R2_001.trimmedS.fastq.gz ILLUMINACLIP:$ADAPTERS:2:30:10:2:True LEADING:3 TRAILING:3 MINLEN:36; done

# 3. Filtlong on raw Nanopore reads
mkdir -p $ROOT/2_trimmed_nanopore_reads && for fq in $ROOT/1_Nanopore_reads/*.fastq*; do s=$(basename $fq | sed 's/\.fastq.*//'); [ -f $ROOT/2_trimmed_nanopore_reads/${s}_filt.fastq.gz ] || filtlong --min_length 1000 --keep_percent 90 $fq | gzip > $ROOT/2_trimmed_nanopore_reads/${s}_filt.fastq.gz; done

# 4. NanoPlot on filtered Nanopore reads
for fq in $ROOT/2_trimmed_nanopore_reads/*.fastq.gz; do s=$(basename $fq .fastq.gz); [ -f $ROOT/3_nanoplot_results/$s/.done ] || NanoPlot --fastq $fq --outdir $ROOT/3_nanoplot_results/$s --prefix $s --threads 4 >/dev/null 2>&1 && touch $ROOT/3_nanoplot_results/$s/.done; done

# 5. Unicycler hybrid assembly (Illumina + Nanopore)
cd $TRIM && for r1 in *_R1_001.trimmedP.fastq.gz; do b=${r1%_R1_001.trimmedP.fastq.gz}; [ -d $ASM/${b}_unicycler ] || [ ! -f $ROOT/2_trimmed_nanopore_reads/${b}_filt.fastq.gz ] || unicycler -1 ${b}_R1_001.trimmedP.fastq.gz -2 ${b}_R2_001.trimmedP.fastq.gz -l $ROOT/2_trimmed_nanopore_reads/${b}_filt.fastq.gz -o $ASM/${b}_unicycler -t 20 --mode normal --kmers 99,127; done

# 6. Collect the assemblies into one folder
mkdir -p $ASM/assemblies && for d in $ASM/*_unicycler; do [ -f $d/assembly.fasta ] && cp $d/assembly.fasta $ASM/assemblies/$(basename $d _unicycler)_assembly.fasta; done

# 7. QUAST on all assemblies
quast.py $ASM/assemblies/*.fasta -o $WORK/7_quast -t 16

# 8. BTyper3 typing of each assembly, then one summary file
cd $ASM/assemblies && for f in *.fasta; do [ -d $WORK/8_btyper_results/${f%.fasta} ] || btyper3 -i $f --ani_geneflow True --ani_typestrains True -o $WORK/8_btyper_results/${f%.fasta}; done
mkdir -p $WORK/8_btyper_results/summary && cat $(find $WORK/8_btyper_results -path '*/btyper3_final_results/*_final_results.txt') > $WORK/8_btyper_results/summary/all_results.txt

# 9. MOB-recon plasmid typing of each assembly
for f in $ASM/assemblies/*.fasta; do s=$(basename $f .fasta); [ -f $MOB/$s/.done ] || mob recon -i $f -o $MOB/$s -n 4 >/dev/null 2>&1 && touch $MOB/$s/.done; done

# 10. Mash distances between plasmids (AE387 16S artifact excluded)
mkdir -p $MASH_DIR && find $MOB -name 'plasmid_*.fasta' ! -name 'plasmid_AE387.fasta' | sort > $MASH_DIR/plasmid_files.txt && mash sketch -k 21 -s 1000 -o $MASH_DIR/plasmids -l $MASH_DIR/plasmid_files.txt && mash dist $MASH_DIR/plasmids.msh $MASH_DIR/plasmids.msh > $MASH_DIR/mash_dist.tsv

# 11. kSNP4 for extraction of SNPs between each panC group
export PATH="/storage/work/tfc5209/programs/kSNP4_1_Linux_package/kSNP4_1_pkg/:$PATH"
MakeKSNP4infile -indir group_II
MakeKSNP4infile -indir group_III
MakeKSNP4infile -indir group_IV
MakeKSNP4infile -indir group_V

Kchooser4 -in group_II
Kchooser4 -in group_III
Kchooser4 -in group_IV
Kchooser4 -in group_V

kSNP4 -in group_II -k 21 -outdir group_II_out/ -core -ML -CPU 2 -debug
kSNP4 -in group_III -k 21 -outdir group_III_out/ -core -ML -CPU 2 -debug
kSNP4 -in group_IV -k 21 -outdir group_IV_out/ -core -ML -CPU 2 -debug
kSNP4 -in group_V -k 21 -outdir group_V_out/ -core -ML -CPU 2 -debug

#12. RaxML for phylogenetic tree building
for d in /storage/work/tfc5209/Slovenian_BC/10_kSNP4/ksnp4_out_*/; do [ -d "$d" ] || continue; g=$(basename "$d" | sed 's/^ksnp4_out_//'); o=/storage/work/tfc5209/Slovenian_BC/11_raxML/raxML_out_$g; mkdir -p "$o"; a=$(find "$d" -name "*core*SNP*.fasta" | head -1); if [ -f "$a" ]; then (cd "$o" && raxml-ng --all --msa "$a" --model GTR+G --prefix "group_$g" --seed 1234567890 --bs-trees 1000 --threads 2); else echo "ERROR: no core SNP matrix in $d"; fi; done
