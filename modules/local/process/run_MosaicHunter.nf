/*
* Run MosaicHunter to detect mosaic SNVs
*/
process RUN_MOSAICHUNTER {

  tag "$sampID"
  label 'process_medium'
  publishDir "${params.outdir}/Mosaics/SNVs/TXT/unfilteredVariants/", pattern: '*.tsv', mode: 'copy', saveAs: { filename -> "${sampID}.MosaicHunter.MosaicSNVs.txt" }
  publishDir "${params.outdir}/Mosaics/SNVs/logs/runMosaicHunter/", pattern: '*.log', mode: 'copy', saveAs: { filename -> "${sampID}.MosaicHunter.MosaicSNVs.log" }

  input:
  tuple val(sampID), file(bam), file(bai), file(vcf), file(vcf_idx), val(sex), file(bed)
  file(fasta)
  file(conf)
  val(mode)
  file(gapsRef)


  output:
  tuple val(sampID), path("final.passed.tsv"), emit: tsv
  path("stdout*.log"), emit: log

  script:
  """
  java -jar /opt/MosaicHunter/build/mosaichunter.jar -C $conf \
  -P input_file=$bam \
  -P reference_file=$fasta \
  -P output_dir=./ \
  -P mosaic_filter.sex=$sex \
  -P mosaic_filter.mode=$mode \
  -P repetitive_region_filter.bed_file=$gapsRef \
  -P indel_region_filter.bed_file=$bed
  """

}
