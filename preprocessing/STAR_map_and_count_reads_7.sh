#GSE83240 using GRCm39 reference genome
#paired end

PATH_INPUT="/data/public/pungerav/liver_disease_score/data_disease_score/GSE83240/SRR3657"
INDEX_DIRECTORY="/data/public/tpadvits/PROJECTS/PodocytePJ/Liver_PAULA/index_GRm39"
OUTPUT_NAME="SRR3657"

for ((i=797; i<=799; i++ ))
do

STAR --genomeDir ${INDEX_DIRECTORY} --readFilesIn ${PATH_INPUT}${i}_1.fastq,${PATH_INPUT}${i}_2.fastq --quantMode GeneCounts --runThreadN 10 --outSAMtype BAM Unsorted --outFileNamePrefix ${OUTPUT_NAME}${i}

done


for ((i=803; i<=805; i++ ))
do

STAR --genomeDir ${INDEX_DIRECTORY} --readFilesIn ${PATH_INPUT}${i}_1.fastq,${PATH_INPUT}${i}_2.fastq --quantMode GeneCounts --runThreadN 10 --outSAMtype BAM Unsorted --outFileNamePrefix ${OUTPUT_NAME}${i}

done
