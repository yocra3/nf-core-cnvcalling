/*
 * Expand multiallelic sites from a VCF
 */
 process BCFTOOLS_EXPAND_MULTIALLELIC {
  tag "$sampID"

  input:
  tuple val(sampID), file(vcf), file(vcf_idx)

  output:
  tuple val(sampID), file("${sampID}.multi.vcf.gz")

  script:
  """
  bcftools norm -m-both $vcf -o ${sampID}.multi.vcf.gz -O z
  """

}
