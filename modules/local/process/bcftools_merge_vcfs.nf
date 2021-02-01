/*
* Merge VCFs with bcftools
*/
process BCFTOOLS_MERGE_VCFS {

  input:
  file(vcf)

  output:
  file("merged.vcf.gz")

  script:
  """
  bcftools merge *.gz -m none -o merged.vcf.gz -O z
  """
}
