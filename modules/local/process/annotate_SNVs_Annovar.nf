/*
* Annotate SNVs with Annovar.
* Note: annovar is not distributed in docker. Input parameters point to annovar files.
*/
process ANNOTATE_SNVS_ANNOVAR {

  tag "$sampID"
  label 'process_long'
  label 'process_medium'
  cpus 16

  publishDir "${params.outdir}/SNVs/Annovar/", pattern: '*.txt', mode: 'copy'

  input:
  tuple val(sampID), file(avinput)
  file(annovar)
  file(annovarVar)
  file(annovarCod)
  file(annovarXref)
  file(annovarFold)

  output:
  tuple val(sampID), file("${sampID}.hg19_multianno.txt")

  script:
  """
  perl $annovar $avinput $annovarFold -buildver hg19 -out ${sampID} -remove \
  --xref $annovarXref \
  -protocol refGene,cytoBand,genomicSuperDups,gnomad211_exome,gnomad211_genome,avsnp150,kaviar_20150923,clinvar_20200316,dbnsfp41a \
  -operation gx,r,r,f,f,f,f,f,f -nastring . --otherinfo --thread ${task.cpus}
  """
}
