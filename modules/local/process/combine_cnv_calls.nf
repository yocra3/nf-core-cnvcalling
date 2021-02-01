/*
* Combine CNVnator and ERDS calls
*/
process COMBINE_CNV_CALLS {

  tag "$sampID"
  publishDir "${params.outdir}/CNVs/TXT/Raw/", mode: 'copy'

  input:
  tuple val(sampID), file(cnvnatorCall), file(erdsCall)

  output:
  tuple val(sampID), file("${sampID}.ERDS_CNVnator_CNVs.raw.txt")

  script:
  """
  python2 /opt/TCAG-WGS-CNV-workflow/add_features.py -i $cnvnatorCall -a $erdsCall -o ${sampID}.ERDS_CNVnator_CNVs.raw.txt -s $sampID -c 0 -p reciprocal
  """
}
