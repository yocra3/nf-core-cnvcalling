/*
* Prioritize CNVs:
* - Remove CNVs with overlap > 20% with commonCNVs
* - Select CNVs with overlap > 80% pathogenic variants (subset1)
* - Select CNVs overlapping exons in OMIM genes (subset2)
* - Select CNVs overlapping exons in GENCODE genes (subset3)
*/
process PRIORITIZE_CNVS {

  tag "$sampID"
  publishDir "${params.outdir}/CNVs/XLSX/", mode: 'copy'


  input:
  tuple val(sampID), file(vcf), file(vcftbi)

  output:
  file("${sampID}.ERDS_CNVnator_CNVs.Prioritization.xlsx")

  script:
  """
  prioritizeCNVs.R $vcf ${sampID}.ERDS_CNVnator_CNVs.Prioritization.xlsx
  """

}
