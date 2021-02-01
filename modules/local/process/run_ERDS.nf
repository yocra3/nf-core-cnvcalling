/*
 * CNV calling with ERDS
 */
 process RUN_ERDS {

   tag "$sampID"

   label 'process_medium'

   input:
   tuple val(sampID), file(bam), file(bai), file(vcf), file(vcf_idx)
   tuple file(fasta), file(fai), file(dict)

   output:
   tuple val(sampID), file("${sampID}/${sampID}.erds.vcf")

   script:
   """
   erds_pipeline.pl -b $bam -v $vcf -o ./$sampID -r $fasta -n $sampID --samtools /opt/conda/envs/nf-core-cnvcalling/bin/samtools
   """

 }
