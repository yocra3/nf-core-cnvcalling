/*
 * Select indels in a VCF
 */
 process BCFTOOLS_SELECT_INDELS {
  tag "$sampID"

  input:
  tuple val(sampID), file(bam), file(bai), file(vcf), file(vcf_idx)

  output:
  tuple val(sampID), file("${sampID}.indels.vcf.gz") 

  script:
  """
  bcftools view -v indels  -f .,PASS $vcf -o ${sampID}.indels.vcf.gz -O z
  """

}
