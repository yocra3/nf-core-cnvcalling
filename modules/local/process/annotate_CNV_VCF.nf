/*
* Add annotation to CNVs
*/
process ANNOTATE_CNV_VCF {

  tag "$sampID"

  input:
  tuple val(sampID), file(vcf), file(vcftbi)
  file(commonCNV)
  file(clinvarCNV)
  file(gtfRef)
  file(omim)

  output:
  tuple val(sampID), file("${sampID}.ERDS_CNVnator_CNVs.annotated.vcf")

  """
  annotateCNVs.R $vcf $commonCNV $clinvarCNV $gtfRef $omim ${sampID}.ERDS_CNVnator_CNVs.annotated.vcf
  """

}
