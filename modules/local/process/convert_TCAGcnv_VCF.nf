/*
* Convert CNVnator output to VCF
*/
process CONVERT_TCAGCNV_VCF {

  tag "$sampID"

  input:
  tuple val(sampID), file(merged)

  output:
  tuple val(sampID), file("${sampID}.CNVs.vcf")

  script:
  """
  ## Make Header
echo "##fileformat=VCFv4.3
##source=CNVnator
##INFO=<ID=END,Number=1,Type=Integer,Description='End position of the variant described in this record'>
##INFO=<ID=SVLEN,Number=1,Type=Integer,Description='Difference in length between REF and ALT alleles'>
##INFO=<ID=SVTYPE,Number=1,Type=String,Description='Type of structural variant'>
##FORMAT=<ID=GT,Number=1,Type=String,Description='Genotype'>
##FORMAT=<ID=NRD,Number=1,Type=String,Description='Normalized Read depth'>
##FORMAT=<ID=E,Number=1,Type=Float,Description='e-val2'>
##FORMAT=<ID=q0,Number=1,Type=Float,Description='q0'>
##FORMAT=<ID=NCNV,Number=1,Type=Integer,Description='Number of CNVs'>
##FORMAT=<ID=LCNVS,Number=1,Type=Integer,Description='Length of CNVS'>
##FORMAT=<ID=LG,Number=1,Type=Integer,Description='Length of gaps'>
##FORMAT=<ID=PG,Number=1,Type=Float,Description='Proportion of Gaps'>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t$sampID" > ${sampID}.CNVs.vcf

## Add calls
tail -n +2 $merged | awk  '{OFS = "\t"}!/^#/{print \$2, \$3, ".", "N", ".", "<"\$5">", "PASS", "SVTYPE="\$5";END="\$4";SVLEN="\$6,\
   "GT:NRD:E:q0:NCNV:LCNVS:LG:PG", "1/0:"\$7":"\$8":"\$9":"\$10":"\$11":"\$12":"\$13  }' - \
   >> ${sampID}.CNVs.vcf
   sed -i 's/:-/:\\./g' ${sampID}.CNVs.vcf ## Change empty values to .
   sed -i 's/\\x27/\\x22/g' ${sampID}.CNVs.vcf ## Change ' for "
  """
}
