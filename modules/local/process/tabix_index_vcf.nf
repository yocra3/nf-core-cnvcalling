/*
* Index vcf
*/
process TABIX_INDEX_VCF {

  tag "$sampID"

  input:
  tuple val(sampID), file(vcf)

  output:
  tuple val(sampID), file(vcf), file("${vcf}.tbi")

  script:
  """
  tabix -p vcf ${vcf}
  """
}
