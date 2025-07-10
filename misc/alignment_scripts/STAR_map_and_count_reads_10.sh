#GSE97024 using GRCm39 reference genome
#single end

PATH_INPUT="/data/public/pungerav/liver_disease_score/data_disease_score/GSE97024/SRR7050"
INDEX_DIRECTORY="/data/public/tpadvits/PROJECTS/PodocytePJ/Liver_PAULA/index_GRm39"
OUTPUT_NAME="SRR7050"

for ((i=690; i<=701; i++ ))
do

STAR --genomeDir ${INDEX_DIRECTORY} --readFilesIn ${PATH_INPUT}${i}.fastq --quantMode GeneCounts --runThreadN 10 --outSAMtype BAM Unsorted --outFileNamePrefix ${OUTPUT_NAME}${i}

done
  
