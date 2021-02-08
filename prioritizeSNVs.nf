#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// Load input data
if (params.inputTable) {
      Channel
            .fromPath(params.inputTable)
            .splitCsv(sep: "\t")
            .map { row -> [ row[0], file(row[1], checkIfExists: true),  file(row[1] + ".tbi", checkIfExists: true) ] }
            .ifEmpty { exit 1, "params.inputTable was empty - no input files supplied" }
            .set { ch_input }

            Channel
              .fromPath(params.inputTable)
              .splitCsv(sep: "\t")
              .map { row -> [ row[0], row[2] ]}
              .set { ch_sex }

} else { exit 1, "--inputTable should be text file with the sample id, a vcf file and sample sex should be provided." }

if (params.fasta ) {

      Channel.fromPath(params.fasta, checkIfExists: true)
        .ifEmpty { exit 1, "Genome fasta file not found: ${params.fasta}" }
        .set { ch_fasta }
  }


/*
* Initialize parameters
*/
omim = file("${params.omim}")
omim_map = file("${params.omim_map}")
annovar = file("${params.annovarPath}/table_annovar.pl")
annovarVar = file("${params.annovarPath}/annotate_variation.pl")
annovarCod = file("${params.annovarPath}/coding_change.pl")
annovarConv = file("${params.annovarPath}/convert2annovar.pl")
annovarXref = file("${params.annovarPath}/example/gene_fullxref.txt")
annovarFold = file("${params.annovarFold}")

// Include modules
include { CREATE_FASTA_INDEX } from './modules/local/process/create_FASTA_index.nf'
include { CREATE_FASTA_DICT } from './modules/local/process/create_FASTA_dict.nf'
include { SPLIT_FASTA } from './modules/local/process/split_FASTA.nf'
include { BCFTOOLS_MERGE_VCFS } from './modules/local/process/bcftools_merge_vcfs.nf'
include { SELECT_NONPRIVATE_VARIANTS } from './modules/local/process/select_nonprivate_variants.nf'
include { BCFTOOLS_EXPAND_MULTIALLELIC } from './modules/local/process/bcftools_expand_multiallelic.nf'
include { BCFTOOLS_LEFT_NORMALIZE } from './modules/local/process/bcftools_left_normalize.nf'
include { VCF_TO_ANNOVAR } from './modules/local/process/VCF_to_annovar.nf'
include { ANNOTATE_SNVS_ANNOVAR } from './modules/local/process/annotate_SNVs_Annovar.nf'
include { PRIORITIZE_SNVS } from './modules/local/process/prioritize_SNVs.nf'
include { TABIX_INDEX_VCF } from './modules/local/process/tabix_index_vcf.nf'
include { BGZIP_COMPRESS_VCF } from './modules/local/process/bgzip_compress_vcf.nf'



workflow PREPROCESS {

  take:
  ch_fasta

  main:
  CREATE_FASTA_INDEX(ch_fasta)
  CREATE_FASTA_DICT(CREATE_FASTA_INDEX.out)
  SPLIT_FASTA(ch_fasta)

  emit:
  idx = CREATE_FASTA_INDEX.out
  dict = CREATE_FASTA_DICT.out
  split = SPLIT_FASTA.out
}

workflow FILTER_SNVS {

  take:
  input
  fasta
  annovar
  annovarVar
  annovarCod
  annovarXref
  annovarFold
  annovarConv
  omim
  omim_map

  main:
  BCFTOOLS_EXPAND_MULTIALLELIC(input)
  BCFTOOLS_LEFT_NORMALIZE(BCFTOOLS_EXPAND_MULTIALLELIC.out, fasta)

  TABIX_INDEX_VCF(BCFTOOLS_LEFT_NORMALIZE.out)
  BCFTOOLS_MERGE_VCFS(TABIX_INDEX_VCF.out.collect())
  SELECT_NONPRIVATE_VARIANTS(BCFTOOLS_MERGE_VCFS.out)

  VCF_TO_ANNOVAR(BCFTOOLS_LEFT_NORMALIZE.out, annovarConv)
  ANNOTATE_SNVS_ANNOVAR(VCF_TO_ANNOVAR.out, annovar, annovarVar, annovarCod, annovarXref, annovarFold)
  PRIORITIZE_SNVS(ANNOTATE_SNVS_ANNOVAR.out, omim, omim_map, SELECT_NONPRIVATE_VARIANTS.out)

}

workflow {
  PREPROCESS(ch_fasta)
  FILTER_SNVS(ch_input, PREPROCESS.out.idx.collect(), annovar, annovarVar, annovarCod, annovarXref, annovarFold, annovarConv, omim, omim_map)
}
