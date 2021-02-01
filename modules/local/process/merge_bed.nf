/*
* Concatenate two bed files
*/
process MERGE_BED {

  tag "$sampID"

  input:
  tuple val(sampID), file(indels), file(cnvs)

  output:
  tuple val(sampID), file("${sampID}.merged.bed") 

  """
  cp $indels ${sampID}.merged.bed
  cat $cnvs >> ${sampID}.merged.bed
  """

}
