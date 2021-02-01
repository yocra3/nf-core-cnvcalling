/*
* Convert indels from VCF to BED format
*/
 process VCF_INDELS_TO_BED {

   tag "$sampID"

   input:
   tuple val(sampID), file(vcf)

   output:
   tuple val(sampID), file("${sampID}.indels.bed")

   script:
   """
   ## Decompress VCF on the fly. Add +/- 5 bp as recommended in MosaicHunter
   awk  '{OFS = "\t"}!/^#/{print \$1, \$2 - 5, \$2 + (length(\$2)) + 5}' <(gzip -cd $vcf) > ${sampID}.indels.bed
   """

 }
