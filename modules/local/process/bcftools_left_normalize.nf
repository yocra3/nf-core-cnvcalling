/*
 * Left normalize indels in VCF
 */
 process BCFTOOLS_LEFT_NORMALIZE {
  tag "$sampID"

  input:
  tuple val(sampID), file(vcf)
  tuple file(fasta), file(fastaidx)

  output:
  tuple val(sampID), file("${sampID}.leftNorm.vcf.gz")

  script:
  """
  bcftools norm --check-ref s -f $fasta $vcf -o ${sampID}.leftNorm.vcf.gz -O z
  """

}
