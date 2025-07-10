#GSE97234 using GRCm39 reference genome
#single end reads

PATH_INPUT="/data/public/pungerav/liver_disease_score/data_disease_score/GSE97234/SRR53981"
INDEX_DIRECTORY="/data/public/tpadvits/PROJECTS/PodocytePJ/Liver_PAULA/index_GRm39"
OUTPUT_NAME="SRR53981"

for ((i=59; i<=68; i++ ))
do

STAR --genomeDir ${INDEX_DIRECTORY} --readFilesIn ${PATH_INPUT}${i}.fastq --quantMode GeneCounts --runThreadN 10 --outSAMtype BAM Unsorted --outFileNamePrefix ${OUTPUT_NAME}${i}

done
