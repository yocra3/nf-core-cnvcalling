/*
* Cluster CNVnator calls
*/
process CLUSTER_CNVNATOR {

  tag "$sampID"

  label 'process_low'

  input:
  tuple val(sampID), file(cnvnatorCall)
  file(gapsRef)

  output:
  tuple val(sampID), file("merged/${sampID}.cnvnator.txt.cluster.txt") 

  script:
  """
  mkdir calls

  echo $cnvnatorCall"\t"$sampID > ids.map
  grep -v "GL" $cnvnatorCall > calls/$cnvnatorCall

  mkdir merged
  python2 /opt/TCAG-WGS-CNV-workflow/merge_cnvnator_results.py -i ./calls/ -a ids.map -o ./merged/ -g $gapsRef
  """

}
