process CREATE_FASTA_INDEX {

  input:
  file("ref.fasta")

  output:
  tuple file("ref.fasta"), file("ref.fasta.fai")

  script:
  """
  samtools faidx ref.fasta
  """
}
