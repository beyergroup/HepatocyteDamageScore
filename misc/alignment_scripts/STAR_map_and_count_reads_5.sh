#GSE148849 using GRCm39 reference genome
#paired end

PATH_INPUT="/data/public/pungerav/liver_disease_score/data_disease_score/GSE148849/SRR115619"
INDEX_DIRECTORY="/data/public/tpadvits/PROJECTS/PodocytePJ/Liver_PAULA/index_GRm39"
OUTPUT_NAME="SRR115619"

for ((i=55; i<=72; i++ ))
do

STAR --genomeDir ${INDEX_DIRECTORY} --readFilesIn ${PATH_INPUT}${i}_1.fastq,${PATH_INPUT}${i}_2.fastq --quantMode GeneCounts --runThreadN 10 --outSAMtype BAM Unsorted --outFileNamePrefix ${OUTPUT_NAME}${i}

done
