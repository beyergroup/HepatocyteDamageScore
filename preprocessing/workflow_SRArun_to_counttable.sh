## Extraction with SRA toolkit

# Download SRA run in SRA format
#./prefetch SRAnumber

# or on whole list of samples that belong to one accession
# Download SRA run in SRA format
#./prefetch --option-file SraAccList.txt

# check for integrity: should report 'ok' and 'consistent' for all parameters
#./vdb-validate SRAnumber

# convert from sra to fastq
#./fasterq-dump SRAnumber -O /data/public/pungerav/liver_disease_score/data_disease_score/GSE111828

#questions: 1. what does the output on the console mean: exp.
#[pungerav@beyer-ws01 bin]$ ./fasterq-dump SRR5572772
#spots read      : 17,993,089
#reads read      : 35,986,178
#reads written   : 35,986,178
#
# any options I should consider adding?

## Quality Control with FastQC

#./fastqc path_containing_output_from_fasterq_dump.fastq

# output in html form stored in directory containing fastq file
# cd to directory all the fastqc output files
# summarize results from all samples with multiqc

#multiqc .

# optional trimming and filtering step

## alignment with STAR

#STAR command line has the following format:
#STAR --option1-name option1-value(s)--option2-name option2-value(s) ...
#If an option can accept multiple values, they are separated by spaces,
#and in a few cases - by commas

#Two step process:
#1)generate genome index files: supply reference genome fasta files
# and gtf annotations file
# For GSE99010 using GRCm39 reference genome:
# https://www.ncbi.nlm.nih.gov/assembly/GCF_000001635.27/
# download corresponding annotations (beware of chromsome names)
# ftp://ftp.ncbi.nlm.nih.gov/genomes//all/annotation_releases/10090/109/GCF_000001635.27_GRCm39/GCF_000001635.27_GRCm39_genomic.gtf.gz
#STAR --runMode genomeGenerate --genomeFastaFiles thefastfile1.fasta the fastafileN.fasta \
 #--sjdbGTFfile annotations_file.gtf --genomeDir GenomeDirectory/

#Example:
#STAR --runMode genomeGenerate --genomeFastaFiles \
#/data/public/pungerav/global_data/genomes/mus_musculus/GRCm39/ncbi-genomes-2021-12-01/GCF_000001635.27_GRCm39_genomic.fasta \
#--sjdbGTFfile /data/public/pungerav/global_data/genomes/mus_musculus/GRCm39/ncbi-genomes-2021-12-01/GCF_000001635.27_GRCm39_genomic.gtf \
#--genomeDir /data/public/pungerav/global_data/genomes/mus_musculus/GRCm39/ncbi-genomes-2021-12-01/index/ \
#--sjdbOverhang 100 \
#--runThreadN 20


#2) Mapping
PATH_INPUT="/data/public/pungerav/liver_disease_score/data_disease_score/GSE99010/temp/SRR55727"
#STAR --genomeDir /data/public/tpadvits/PROJECTS/PodocytePJ/Liver_PAULA/STAR_output --readFilesIn /data/public/pungerav/liver_disease_score/data_disease_score/GSE99019/temp/SRR557263_1.fastq,/data/public/pungerav/liver_disease_score/data_disease_score/GSE99019/temp/SRR557263_2.fastq,/data/public/pungerav/liver_disease_score/data_disease_score/GSE99019/temp/SRR557264_1.fastq,/data/public/pungerav/liver_disease_score/data_disease_score/GSE99019/temp/SRR557264_2.fastq,/data/public/pungerav/liver_disease_score/data_disease_score/GSE99019/temp/SRR557265_1.fastq,/data/public/pungerav/liver_disease_score/data_disease_score/GSE99019/temp/SRR557265_2.fastq,/data/public/pungerav/liver_disease_score/data_disease_score/GSE99019/temp/SRR557266_1.fastq,/data/public/pungerav/liver_disease_score/data_disease_score/GSE99019/temp/SRR557266_2.fastq,/data/public/pungerav/liver_disease_score/data_disease_score/GSE99019/temp/SRR557267_1.fastq,/data/public/pungerav/liver_disease_score/data_disease_score/GSE99019/temp/SRR557267_2.fastq,/data/public/pungerav/liver_disease_score/data_disease_score/GSE99019/temp/SRR557268_1.fastq,/data/public/pungerav/liver_disease_score/data_disease_score/GSE99019/temp/SRR557268_2.fastq,/data/public/pungerav/liver_disease_score/data_disease_score/GSE99019/temp/SRR557269_1.fastq,/data/public/pungerav/liver_disease_score/data_disease_score/GSE99019/temp/SRR557269_2.fastq,/data/public/pungerav/liver_disease_score/data_disease_score/GSE99019/temp/SRR557270_1.fastq,/data/public/pungerav/liver_disease_score/data_disease_score/GSE99019/temp/SRR557270_2.fastq,/data/public/pungerav/liver_disease_score/data_disease_score/GSE99019/temp/SRR557271_1.fastq,/data/public/pungerav/liver_disease_score/data_disease_score/GSE99019/temp/SRR557271_2.fastq,/data/public/pungerav/liver_disease_score/data_disease_score/GSE99019/temp/SRR557272_1.fastq,/data/public/pungerav/liver_disease_score/data_disease_score/GSE99019/temp/SRR557272_2.fastq \
#--quantMode GeneCounts --runThreadN 10 --outSAMtype BAM Unsorted

STAR --genomeDir /data/public/tpadvits/PROJECTS/PodocytePJ/Liver_PAULA/STAR_output --readFilesIn ${PATH_INPUT}63_1.fastq,${PATH_INPUT}63_2.fastq,${PATH_INPUT}64_1.fastq,${PATH_INPUT}64_2.fastq,${PATH_INPUT}65_1.fastq,${PATH_INPUT}65_2.fastq,${PATH_INPUT}66_1.fastq,${PATH_INPUT}66_2.fastq,${PATH_INPUT}67_1.fastq,${PATH_INPUT}67_2.fastq,${PATH_INPUT}68_1.fastq,${PATH_INPUT}68_2.fastq,${PATH_INPUT}69_1.fastq,${PATH_INPUT}69_2.fastq,${PATH_INPUT}70_1.fastq,${PATH_INPUT}70_2.fastq,${PATH_INPUT}71_1.fastq,${PATH_INPUT}71_2.fastq,${PATH_INPUT}72_1.fastq,${PATH_INPUT}72_2.fastq \
--quantMode GeneCounts --runThreadN 10 --outSAMtype BAM Unsorted

#alternative for pipeline: 
#filenames=$(ls /data/public/pungerav/liver_disease_score/data_disease_score/GSE99019/temp/)
#filename=($filenames)


#outTmpKeep None
#when running on server:
#option --runThreadN defines number of threads to be used for genome generation
# depends on number of available cores
