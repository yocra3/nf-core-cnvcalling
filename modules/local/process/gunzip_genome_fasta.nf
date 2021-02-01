// TODO nf-core: Simplificar parametros

// This complex logic is to prevent accessing the genome_fasta_gz variable if
// necessary indices for STAR, HiSAT2, Salmon already exist, or if
// params.transcript_fasta is provided as then the transcript sequences don't
// need to be extracted.
process GUNZIP_GENOME_FASTA {
  tag "$gz"

  input:
  file gz 

  output:
  file "${gz.baseName}"

  script:
  """
  gunzip -k --verbose --stdout --force ${gz} > ${gz.baseName}
  """
}
