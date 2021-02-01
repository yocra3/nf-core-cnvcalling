/*
* Convert MosaicHunter to VCF
*/
process CONVERT_MOSAICHUNTER_VCF {

  tag "$sampID"

  input:
  tuple val(sampID), file(tsv)

  output:
  tuple val(sampID), file("${sampID}.MosaicHunter.MosaicSNVs.vcf")

  script:
  """
  ## Make Header
echo "##fileformat=VCFv4.3
##source=MosaicHunter
##FORMAT=<ID=GT,Number=1,Type=String,Description='Genotype'>
##FORMAT=<ID=AF,Number=1,Type=Float,Description='Frequency of alternate allele'>
##FORMAT=<ID=DP,Number=1,Type=Float,Description='Read depth'>
##FORMAT=<ID=AD,Number=.,Type=Integer,Description='Allelic depths for the ref and alt alleles in the order listed.'>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t$sampID" > ${sampID}.MosaicHunter.MosaicSNVs.vcf

## Add calls
awk  '{OFS = "\t"}
{if (\$3=\$7)
      print \$1, \$2, ".", \$3, \$9, ".", "PASS", ".", "GT:AF:DP:AD", "0/1:"\$10/\$4":"\$4":"\$8","\$10;
else
    print \$1, \$2, ".", \$3, \$7, ".", "PASS", ".", "GT:AF:DP:AD", "0/1:"\$8/\$4":"\$4":"\$10","\$8;
   }'    $tsv >> ${sampID}.MosaicHunter.MosaicSNVs.vcf
   sed -i 's/:-/:\\./g' ${sampID}.MosaicHunter.MosaicSNVs.vcf ## Change empty values to .
   sed -i 's/\\x27/\\x22/g' ${sampID}.MosaicHunter.MosaicSNVs.vcf ## Change ' for "
  """
}
