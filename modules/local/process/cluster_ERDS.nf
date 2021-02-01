/*
* Cluster ERDS calls
*/
process CLUSTER_ERDS {
  tag "$sampID"

  input:
  tuple val(sampID), file(erdsCall)
  file(gapsRef)

  output:
  tuple val(sampID), file("merged/${sampID}.erds.txt.cluster.txt")

  script:
  """
  mkdir calls

  echo $erdsCall"\t"$sampID > ids.map
  grep -v "GL" $erdsCall > calls/$erdsCall


  mkdir merged
  python2 /opt/TCAG-WGS-CNV-workflow/merge_erds_results.py -i ./calls/ -a ids.map -o ./merged/ -g $gapsRef
  """
}
