/*
* Format CNVnator
*/
process FORMAT_CNVNATOR {
  tag "$sampID"

  input:
  tuple val(sampID), file("${sampID}.cnvnator.txt")

  output:
  tuple val(sampID), file("out/${sampID}.cnvnator.txt")

  script:
  """
  mkdir out
  python2 /opt/TCAG-WGS-CNV-workflow/format_cnvnator_results.py ${sampID}.cnvnator.txt out/${sampID}.cnvnator.txt
  """

}
