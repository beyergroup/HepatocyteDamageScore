
#GSE138419 using GRCm39 reference genome
#single end reads

PATH_INPUT="/data/public/pungerav/liver_disease_score/data_disease_score/GSE138419/SRR102305"
INDEX_DIRECTORY="/data/public/tpadvits/PROJECTS/PodocytePJ/Liver_PAULA/index_GRm39"
OUTPUT_NAME="SRR102305"

for ((i=41; i<=51; i++ ))
do

STAR --genomeDir ${INDEX_DIRECTORY} --readFilesIn ${PATH_INPUT}${i}.fastq --quantMode GeneCounts --runThreadN 10 --outSAMtype BAM Unsorted --outFileNamePrefix ${OUTPUT_NAME}${i}

done
