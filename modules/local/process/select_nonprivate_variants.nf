/*
* Select variants present in more than one individual
*/
process SELECT_NONPRIVATE_VARIANTS {

  input:
  file(vcf)

  output:
  file("nonprivate.variants.vcf.gz") 

  script:
  """
  bcftools view -i 'N_SAMPLES-N_MISSING>1'  $vcf -o nonprivate.variants.vcf.gz -O z
  """
}
