#sample input folder as first argument after script name, the fastq file 1 as second input 
#this is for analysing pair-end read. 

sample_folder=$1
sample_fastq1=$2

sample_name=${sample_fastq1%_R1_001.fastq}
sample_fastq2=${sample_name}_R2_001.fastq



#========================================================================================================
#To improve chimeric read detection, the reads are mapped by STAR using single-end alignment. That means read 2 needs to be flipped manually by seqtk before the alignment.


#flip R2
 seqtk seq -r ${sample_folder}/${sample_fastq2}   > ${sample_folder}/${sample_name}_R2_flipped.fastq  

 STAR_input1=${sample_folder}/${sample_fastq1}

 STAR_input2=${sample_folder}/${sample_name}_R2_flipped.fastq  


./STAR   --runMode genomeGenerate   --runThreadN 32 \
  --genomeDir genomeDir --genomeFastaFiles gb2018/genbank2018.fna  \
      --sjdbGTFfile gb2018/genbank2018.gff   --sjdbGTFtagExonParentTranscript Parent \
        --sjdbOverhang 74
 
# for single end
./STAR --runThreadN 32 --genomeDir genomeDir --readFilesIn  ${STAR_input1}   \
 --outFileNamePrefix STAR_output/${sample_name}_chimera_Single_R1_ --outReadsUnmapped Fastx \
 --outSAMattributes All --alignIntronMin 1 --alignIntronMax 1 --chimSegmentMin 5 --chimScoreJunctionNonGTAG 0 --chimOutType SeparateSAMold
 

 ./STAR --runThreadN 32 --genomeDir genomeDir --readFilesIn ${STAR_input2}   \
 --outFileNamePrefix STAR_output/${sample_name}_chimera_Single_R2_flipped --outReadsUnmapped Fastx \
 --outSAMattributes All --alignIntronMin 1 --alignIntronMax 1 --chimSegmentMin 5 --chimScoreJunctionNonGTAG 0 --chimOutType SeparateSAMold
 
 echo "STAR completed"
 #merge different STAR output
 python3 Script/Python_script/Combining_all_single_end_output.py $sample_name
echo "merge sam file completed"

echo "featureCounts begin"
#FeatureCounts for total alignment
subread-1.6.3-Linux-x86_64/bin/featureCounts \
-a feature_count_gtf.gtf  -o featureCounts_output/${sample_name}_chimera_Chimeric_Aligned_single_merged_.txt  \
STAR_output/${sample_name}_chimera_Chimeric_Aligned_single_merged_.out.sam  -F GTF -t exon -g gene_id -O -M -s 2 -d 1 -D 5000000 -T 32 -R CORE &&

echo "interaction count begins"
#interaction count only
Rscript Script/R_script/interaction_hyper_hyperExtract.R $sample_name

