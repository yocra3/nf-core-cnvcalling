/*
 * CNV calling with CNVnator
 */
 process RUN_CNVNATOR {
   tag "$sampID"

   label 'process_medium'

   input:
   tuple val(sampID), file(bam), file(bai), file(vcf), file(vcf_idx)
   file(fasta)
   val(genome)
   val(bin_size)

   output:
   tuple val(sampID), file("${sampID}.txt")

   script:
   """
   cnvnator -genome $genome -root ${sampID}.root -tree $bam
   cnvnator -genome $genome -root ${sampID}.root -his $bin_size
   cnvnator -root "${sampID}.root" -stat $bin_size
   cnvnator -root "${sampID}.root" -partition $bin_size -ngc
   cnvnator -root "${sampID}.root" -call $bin_size -ngc > ${sampID}.txt
   """

 }
