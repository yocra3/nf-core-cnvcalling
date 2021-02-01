/*
* Add chromosome names in VCF header based in FASTA
*/
process BCFTOOLS_REHEADER_FASTA {

  tag "$sampID"

  input:
  tuple val(sampID), file(vcf)
  tuple file(fasta), file(fastaidx)

  output:
  tuple val(sampID), file("header.vcf")

  script:
  """
  bcftools reheader -f $fastaidx -o header.vcf $vcf
  """
}
