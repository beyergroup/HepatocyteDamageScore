#GSE124694 using GRCm39 reference genome
#single end

PATH_INPUT="/data/public/pungerav/liver_disease_score/data_disease_score/GSE124694/SRR83974"
INDEX_DIRECTORY="/data/public/tpadvits/PROJECTS/PodocytePJ/Liver_PAULA/index_GRm39"
OUTPUT_NAME="SRR83974"

for ((i=53; i<=59; i++ ))
do

STAR --genomeDir ${INDEX_DIRECTORY} --readFilesIn ${PATH_INPUT}${i}.fastq --quantMode GeneCounts --runThreadN 10 --outSAMtype BAM Unsorted --outFileNamePrefix ${OUTPUT_NAME}${i}

done

STAR --genomeDir ${INDEX_DIRECTORY} --readFilesIn ${PATH_INPUT}61.fastq --quantMode GeneCounts --runThreadN 10 --outSAMtype BAM Unsorted --outFileNamePrefix ${OUTPUT_NAME}61
