/*
* Convert VCF to annovar input
* Note: annovar is not distributed in docker. Input parameters point to annovar files.
*/
process VCF_TO_ANNOVAR {

  tag "$sampID"

  input:
  tuple val(sampID), file(vcf)
  file(annovarConv)

  output:
  tuple val(sampID), file("${sampID}.avinput")

  script:
  """
  perl $annovarConv --format vcf4 --filter pass --outfile ${sampID}.avinput $vcf
  """
}
