
#GSE99010 using GRCm39 reference genome

#PATH_INPUT="/data/public/pungerav/liver_disease_score/data_disease_score/GSE99010/temp/SRR55727"
#INDEX_DIRECTORY="/data/public/tpadvits/PROJECTS/PodocytePJ/Liver_PAULA/index_GRm39"
#OUTPUT_NAME="SRR55727"

#for ((i=63; i<=72; i++ ))
#do

#STAR --genomeDir ${INDEX_DIRECTORY} --readFilesIn ${PATH_INPUT}${i}_1.fastq,${PATH_INPUT}${i}_2.fastq --quantMode GeneCounts --runThreadN 10 --outSAMtype BAM Unsorted --outFileNamePrefix ${OUTPUT_NAME}${i}

#done


#GSE111828 using GRCm39 reference genome

PATH_INPUT="/data/public/pungerav/liver_disease_score/data_disease_score/GSE111828/fastq/SRR68337"
INDEX_DIRECTORY="/data/public/tpadvits/PROJECTS/PodocytePJ/Liver_PAULA/index_GRm39"
OUTPUT_NAME="SRR68337"

for ((i=52; i<=75; i++ ))
do

STAR --genomeDir ${INDEX_DIRECTORY} --readFilesIn ${PATH_INPUT}${i}_1.fastq,${PATH_INPUT}${i}_2.fastq --quantMode GeneCounts --runThreadN 10 --outSAMtype BAM Unsorted --outFileNamePrefix ${OUTPUT_NAME}${i}

done

#GSE138419
# reads are single end
#PATH_INPUT="/data/public/pungerav/liver_disease_score/data_disease_score/GSE138419"
#INDEX_DIRECTORY="/data/public/tpadvits/PROJECTS/PodocytePJ/Liver_PAULA/index_GRm39"
#OUTPUT_NAME="SRR102305"

#for ((i=41; i<=51; i++ ))
#do

#STAR --genomeDir ${INDEX_DIRECTORY} --readFilesIn ${PATH_INPUT}${i}.fastq --quantMode GeneCounts --runThreadN 10 --outSAMtype BAM Unsorted --outFileNamePrefix ${OUTPUT_NAME}${i}

#done


 
