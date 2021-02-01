/*
* Compress VCF with bgzip
*/
process BGZIP_COMPRESS_VCF {

  tag "$sampID"

  input:
  tuple val(sampID), file(vcf)

  output:
  tuple val(sampID), file("${vcf}.gz")

  script:
  """
  bgzip $vcf
  """
}
