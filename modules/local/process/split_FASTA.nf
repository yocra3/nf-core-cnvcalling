/*
 * Pre-process:  Split fasta in chromosomes
 */
process SPLIT_FASTA {

  input:
  file(fasta)

  output:
  file("*.fa")

  script:
  """
  csplit -s -z $fasta '/>/' '{*}'
  for i in xx*
  do
    n=\$(sed 's/>// ; s/ .*// ; 1q' "\$i") ; \
    mv "\$i" "\$n.fa" ; \
  done
  """
}
