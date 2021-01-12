### Subsample reads to generate test dataset
mkdir data

## Subset bam file
samtools view -b -s 0.1 /home/SHARED/PROJECTS/CHD_MARATO/data/WGS/BAMS/62422164_S4.bam > data/62422164_subset.bam
samtools view -b -s 0.1 /home/SHARED/PROJECTS/CHD_MARATO/data/WGS/BAMS/62422145_S6.bam > data/62422145_subset.bam

## Convert bams to fastq
samtools fastq -1 data/62422164_R1.fastq.gz -2 data/62422164_R2.fastq.gz -0 /dev/null -s /dev/null -n data/62422164_subset.bam
samtools fastq -1 data/62422145_R1.fastq.gz -2 data/62422145_R2.fastq.gz -0 /dev/null -s /dev/null -n data/62422145_subset.bam
