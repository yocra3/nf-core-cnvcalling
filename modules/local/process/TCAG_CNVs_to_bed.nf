/*
* Convert TCAG CNVs to BED format
*/
process TCAG_CNVS_TO_BED {

  tag "$sampID"

  input:
  tuple val(sampID), file(txt)

  output:
  tuple val(sampID), file("${sampID}.cnvs.bed")

  script:
  """
  cut -f2-4 $txt | tail -n +2 - > ${sampID}.cnvs.bed
  """

}
