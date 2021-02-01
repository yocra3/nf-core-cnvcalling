/*
* Sort and compress VCF
*/
process BCFTOOLS_SORT_COMPRESS {

  tag "$sampID"

  input:
  tuple val(sampID), file(vcf)

  output:
  tuple val(sampID), file("${vcf}.gz")

  script:
  """
  bcftools sort -o ${vcf}.gz -O z $vcf
  """
}
