/*
* Prioritize SNVs
*/
process PRIORITIZE_SNVS {

  tag "$sampID"
  publishDir "${params.outdir}/SNVs/XLSX/", pattern: '*.xlsx', mode: 'copy'
  publishDir "${params.outdir}/SNVs/logs/prioritization/", pattern: '*.log', mode: 'copy'
  publishDir "${params.outdir}/SNVs/TXT/filteredVariants/", pattern: '*.txt', mode: 'copy'

  input:
  tuple val(sampID), file(annovar)
  file(omim)
  file(omim_map)
  file(nonprivate)


  output:
  path("*.txt"), emit: txt
  path("*.log"), emit: log
  path("*.xlsx"), emit: xlsx

  script:
  """
  prioritizeSNVs.R $annovar $omim $omim_map $nonprivate ${sampID}.MosaicHunter.MosaicSNVs.Prioritization
  """
}
