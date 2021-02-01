/*
* Filter CNV calls
* - Select common CNVs between CNVnator and erds
* - Remove CNVs overlapping > 70% with low complexity regions
*/
process FILTER_CNV_CALLS {
  tag "$sampID"

  publishDir "${params.outdir}/CNVs/TXT/Filtered/", mode: 'copy', pattern: '*.txt'
  publishDir "${params.outdir}/CNVs/log/", mode: 'copy', pattern: '*.log'

  input:
  tuple val(sampID), file(combCall), val(sex)
  file(gapsRef)

  output:
  tuple val(sampID), path("${sampID}.ERDS_CNVnator_CNVs.filtered.txt"), emit: txt
  path("${sampID}.ERDS_CNVnator_CNVs.filtered.log"), emit: log

  script:
  """
  filterCNVs.R $combCall $sex $gapsRef ${sampID}.ERDS_CNVnator_CNVs.filtered.txt
  """
}
