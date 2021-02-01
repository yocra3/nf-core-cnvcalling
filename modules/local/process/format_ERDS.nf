/*
* Format ERDS
*/
process FORMAT_ERDS {
  tag "$sampID"

  input:
  tuple val(sampID), file("erds.vcf")

  output:
  tuple val(sampID), file("${sampID}.erds.txt")

  script:
  """
  python2 /opt/TCAG-WGS-CNV-workflow/format_erds_results.py erds.vcf ${sampID}.erds.txt
  """

}
