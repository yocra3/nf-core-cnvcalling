#!/usr/bin/env nextflow
/*
========================================================================================
                         nf-core/cnvcalling
========================================================================================
 nf-core/cnvcalling Analysis Pipeline.
 #### Homepage / Documentation
 https://github.com/nf-core/cnvcalling
----------------------------------------------------------------------------------------
*/

def helpMessage() {
    // TODO nf-core: Add to this help message with new command line parameters
    log.info nfcoreHeader()
    log.info"""

    Usage:

    The typical command for running the pipeline is as follows:

    nextflow run nf-core/cnvcalling --reads '*_R{1,2}.fastq.gz' -profile docker

    Mandatory arguments:
      --reads [file]                Path to input data (must be surrounded with quotes)
      -profile [str]                Configuration profile to use. Can use multiple (comma separated)
                                    Available: conda, docker, singularity, test, awsbatch, <institute> and more

    Options:
      --genome [str]                  Name of iGenomes reference
      --single_end [bool]             Specifies that the input is single-end reads

    References                        If not specified in the configuration file or you wish to overwrite any of the references
      --fasta [file]                  Path to fasta reference

    Other options:
      --outdir [file]                 The output directory where the results will be saved
      --email [email]                 Set this parameter to your e-mail address to get a summary e-mail with details of the run sent to you when the workflow exits
      --email_on_fail [email]         Same as --email, except only send mail if the workflow is not successful
      --max_multiqc_email_size [str]  Theshold size for MultiQC report to be attached in notification email. If file generated by pipeline exceeds the threshold, it will not be attached (Default: 25MB)
      -name [str]                     Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic

    AWSBatch options:
      --awsqueue [str]                The AWSBatch JobQueue that needs to be set when running on AWSBatch
      --awsregion [str]               The AWS Region for your AWS Batch job to run on
      --awscli [str]                  Path to the AWS CLI tool
    """.stripIndent()
}

// Show help message
if (params.help) {
    helpMessage()
    exit 0
}

/*
 * SET UP CONFIGURATION VARIABLES
 */

// Has the run name been specified by the user?
//  this has the bonus effect of catching both -name and --name
custom_runName = params.name
if (!(workflow.runName ==~ /[a-z]+_[a-z]+/)) {
    custom_runName = workflow.runName
}

if (workflow.profile.contains('awsbatch')) {
    // AWSBatch sanity checking
    if (!params.awsqueue || !params.awsregion) exit 1, "Specify correct --awsqueue and --awsregion parameters on AWSBatch!"
    // Check outdir paths to be S3 buckets if running on AWSBatch
    // related: https://github.com/nextflow-io/nextflow/issues/813
    if (!params.outdir.startsWith('s3:')) exit 1, "Outdir not on S3 - specify S3 Bucket to run on AWSBatch!"
    // Prevent trace files to be stored on S3 since S3 does not support rolling files.
    if (params.tracedir.startsWith('s3:')) exit 1, "Specify a local tracedir or run without trace! S3 cannot be used for tracefiles."
}

// Stage config files
ch_output_docs = file("$baseDir/docs/output.md", checkIfExists: true)

/*
 * Create a channel for input read files
 */

if (params.inputTable) {
      Channel
            .fromPath(params.inputTable)
            .splitCsv(sep: "\t")
            .map { row -> [ row[0], file(row[1], checkIfExists: true), file(row[1] + ".bai", checkIfExists: true) ,
             file(row[2], checkIfExists: true),  file(row[2] + ".tbi", checkIfExists: true) ] }
            .ifEmpty { exit 1, "params.inputTable was empty - no input files supplied" }
            .into { ch_input_cnvnator; ch_input_erds; ch_input_mosaichunter; ch_input_indels; ch_input_mrmosaic }

            Channel
              .fromPath(params.inputTable)
              .splitCsv(sep: "\t")
              .map { row -> [ row[0], row[3] ]}
              .into { ch_sex_CNVs;  ch_sex_mosaics}

} else { exit 1, "--inputTable should be text file with the sample id, a bam file and a vcf file should be provided." }


if (params.fasta && !params.skipAlignment) {
  if (hasExtension(params.fasta, 'gz')) {
    Channel.fromPath(params.fasta, checkIfExists: true)
        .ifEmpty { exit 1, "Genome fasta file not found: ${params.fasta}" }
        .set { genome_fasta_gz }
  } else {
      Channel.fromPath(params.fasta, checkIfExists: true)
        .ifEmpty { exit 1, "Genome fasta file not found: ${params.fasta}" }
        .into { ch_fasta_for_dict; ch_fasta_for_cnvnator; ch_fasta_for_mosaichunter }
  }
}




/*
* Initialize parameters
*/
date = java.time.LocalDate.now()
params.gapsRef = "/home/SHARED/DATA/REFERENCES/GRCh37/Repeats/RLCRs_no_Repeat_Masker.txt"
gapsRef = file("${params.gapsRef}")
params.binSize = 100
params.genome = "hg19"
commonCNV = file("${params.commonCNV}")
clinvarCNV = file("${params.clinvarCNV}")
gtfRef = file("${params.gtfRef}")
omim = file("${params.omim}")
omim_map = file("${params.omim_map}")
mosaichunter_config = file("${params.mosaichunter_config}")
annovar = file("${params.annovarPath}/table_annovar.pl")
annovarVar = file("${params.annovarPath}/annotate_variation.pl")
annovarCod = file("${params.annovarPath}/coding_change.pl")
annovarXref = file("${params.annovarPath}/example/gene_fullxref.txt")
annovarFold = file("${params.annovarFold}")

//
// // Header log info
// log.info nfcoreHeader()
// def summary = [:]
// if (workflow.revision) summary['Pipeline Release'] = workflow.revision
// summary['Run Name']         = custom_runName ?: workflow.runName
// // TODO nf-core: Report custom parameters here
// summary['Reads']            = params.reads
// summary['Fasta Ref']        = params.fasta
// summary['Data Type']        = params.single_end ? 'Single-End' : 'Paired-End'
// summary['Max Resources']    = "$params.max_memory memory, $params.max_cpus cpus, $params.max_time time per job"
// if (workflow.containerEngine) summary['Container'] = "$workflow.containerEngine - $workflow.container"
// summary['Output dir']       = params.outdir
// summary['Launch dir']       = workflow.launchDir
// summary['Working dir']      = workflow.workDir
// summary['Script dir']       = workflow.projectDir
// summary['User']             = workflow.userName
// if (workflow.profile.contains('awsbatch')) {
//     summary['AWS Region']   = params.awsregion
//     summary['AWS Queue']    = params.awsqueue
//     summary['AWS CLI']      = params.awscli
// }
// summary['Config Profile'] = workflow.profile
// if (params.config_profile_description) summary['Config Description'] = params.config_profile_description
// if (params.config_profile_contact)     summary['Config Contact']     = params.config_profile_contact
// if (params.config_profile_url)         summary['Config URL']         = params.config_profile_url
// if (params.email || params.email_on_fail) {
//     summary['E-mail Address']    = params.email
//     summary['E-mail on failure'] = params.email_on_fail
//     summary['MultiQC maxsize']   = params.max_multiqc_email_size
// }
// log.info summary.collect { k,v -> "${k.padRight(18)}: $v" }.join("\n")
// log.info "-\033[2m--------------------------------------------------\033[0m-"
//
// // Check the hostnames against configured profiles
// checkHostname()
//
// Channel.from(summary.collect{ [it.key, it.value] })
//     .map { k,v -> "<dt>$k</dt><dd><samp>${v ?: '<span style=\"color:#999999;\">N/A</a>'}</samp></dd>" }
//     .reduce { a, b -> return [a, b].join("\n            ") }
//     .map { x -> """
//     id: 'nf-core-cnvcalling-summary'
//     description: " - this information is collected when the pipeline is started."
//     section_name: 'nf-core/cnvcalling Workflow Summary'
//     section_href: 'https://github.com/nf-core/cnvcalling'
//     plot_type: 'html'
//     data: |
//         <dl class=\"dl-horizontal\">
//             $x
//         </dl>
//     """.stripIndent() }
//     .set { ch_workflow_summary }
//
compressedReference = hasExtension(params.fasta, 'gz')

if (compressedReference) {
  // TODO nf-core: Simplificar parametros

  // This complex logic is to prevent accessing the genome_fasta_gz variable if
  // necessary indices for STAR, HiSAT2, Salmon already exist, or if
  // params.transcript_fasta is provided as then the transcript sequences don't
  // need to be extracted.
  process gunzip_genome_fasta {
    tag "$gz"
    publishDir path: { params.saveReference ? "${params.outdir}/reference_genome" : params.outdir },
    saveAs: { params.saveReference ? it : null }, mode: 'copy'

    input:
    file gz from genome_fasta_gz

    output:
    file "${gz.baseName}" into ch_fasta_for_dict, ch_fasta_for_cnvnator

    script:
    """
    gunzip -k --verbose --stdout --force ${gz} > ${gz.baseName}
    """
  }
}
//
//
// /*
//  * Parse software version numbers
//  */
// process get_software_versions {
//     publishDir "${params.outdir}/pipeline_info", mode: 'copy',
//         saveAs: { filename ->
//                       if (filename.indexOf(".csv") > 0) filename
//                       else null
//                 }
//
//     output:
//     file 'software_versions_mqc.yaml' into ch_software_versions_yaml
//     file "software_versions.csv"
//
//     script:
//     // TODO nf-core: Get all tools to print their version number here
//     """
//     echo $workflow.manifest.version > v_pipeline.txt
//     echo $workflow.nextflow.version > v_nextflow.txt
//     fastqc --version > v_fastqc.txt
//     multiqc --version > v_multiqc.txt
//     scrape_software_versions.py &> software_versions_mqc.yaml
//     """
// }
//

/*
 * Pre-process: Create fasta dictionary
 */
 process createFastaDict {

   input:
   file("ref.fasta") from ch_fasta_for_dict

   output:
   set file("ref.fasta"), file("ref.fasta.fai"), file("ref.dict") into ch_fasta_for_ERDS, ch_fastadict_VCF

   """
   samtools faidx ref.fasta
   gatk CreateSequenceDictionary -R ref.fasta
   """
 }

 /*
  * Pre-process:  Split fasta in chromosomes
  */
 process splitFasta {

   input:
   file(fasta) from ch_fasta_for_cnvnator

   output:
   file("*.fa") into ch_fastaChr_cnvnator

   """
   csplit -s -z $fasta '/>/' '{*}'
   for i in xx*
   do
     n=\$(sed 's/>// ; s/ .*// ; 1q' "\$i") ; \
     mv "\$i" "\$n.fa" ; \
   done
   """
 }

 /*
  * STEP 4 - CNV calling with CNVnator
  */
  process runCNVnator {
    tag "$sampID"

    label 'process_medium'

    input:
    tuple val(sampID), file(bam), file(bai), file(vcf), file(vcf_idx) from ch_input_cnvnator
    file(fasta) from ch_fastaChr_cnvnator.collect()
    val(genome) from params.genome
    val(bin_size) from params.binSize

    output:
    set sampID, file("${sampID}.txt") into ch_cnvnator_calls

    """
    cnvnator -genome $genome -root ${sampID}.root -tree $bam
    cnvnator -genome $genome -root ${sampID}.root -his $bin_size
    cnvnator -root "${sampID}.root" -stat $bin_size
    cnvnator -root "${sampID}.root" -partition $bin_size -ngc
    cnvnator -root "${sampID}.root" -call $bin_size -ngc > ${sampID}.txt
    """

  }

  /*
  * STEP 4b - Format CNVnator
  */
  process formatCNVNnator {
    tag "$sampID"

    input:
    set sampID, file("${sampID}.cnvnator.txt") from ch_cnvnator_calls

    output:
    set sampID, file("out/${sampID}.cnvnator.txt") into ch_cnvnator_formatcalls

    """
    mkdir out
    python2 /root/TCAG-WGS-CNV-workflow/format_cnvnator_results.py ${sampID}.cnvnator.txt out/${sampID}.cnvnator.txt
    """

  }


  /*
  * STEP 4c - Cluster CNVnator calls
  */
  process clusterCNVnator {

    tag "$sampID"

    label 'process_low'

    input:
    set sampID, file(cnvnatorCall) from ch_cnvnator_formatcalls
    file(gapsRef) from gapsRef

    output:
    set sampID, file("merged/${sampID}.cnvnator.txt.cluster.txt") into ch_cnvnator_clustercalls

    """
    mkdir calls

    echo $cnvnatorCall"\t"$sampID > ids.map
    grep -v "GL" $cnvnatorCall > calls/$cnvnatorCall

    mkdir merged
    python2 /root/TCAG-WGS-CNV-workflow/merge_cnvnator_results.py -i ./calls/ -a ids.map -o ./merged/ -g $gapsRef
    """

  }

  /*
   * STEP 5 - CNV calling with ERDS
   */
   process runERDS {
     tag "$sampID"

     label 'process_medium'

     input:
     tuple val(sampID), file(bam), file(bai), file(vcf), file(vcf_idx) from ch_input_erds
     tuple file(fasta), file(fai), file(dict) from ch_fasta_for_ERDS.collect()

     output:
     set sampID, file("${sampID}/${sampID}.erds.vcf") into ch_ERDS_calls

     """
     erds_pipeline.pl -b $bam -v $vcf -o ./$sampID -r $fasta -n $sampID --samtools /opt/conda/envs/nf-core-cnvcalling/bin/samtools
     """

   }

   /*
   * STEP 5b - Format ERDS
   */
   process formatERDS {
     tag "$sampID"

     input:
     set sampID, file("erds.vcf") from ch_ERDS_calls

     output:
     set sampID, file("${sampID}.erds.txt") into ch_ERDS_formatted_calls
     """
     python2 /root/TCAG-WGS-CNV-workflow/format_erds_results.py erds.vcf ${sampID}.erds.txt
     """

   }

   /*
   * STEP 5c - Cluster ERDS calls
   */
   process clusterERDS {
     tag "$sampID"

     input:
     set sampID, file(erdsCall) from ch_ERDS_formatted_calls
     file(gapsRef) from gapsRef

     output:
     set sampID, file("merged/${sampID}.erds.txt.cluster.txt") into ch_erds_clustercalls

     """
     mkdir calls

     echo $erdsCall"\t"$sampID > ids.map
     grep -v "GL" $erdsCall > calls/$erdsCall


     mkdir merged
     python2 /root/TCAG-WGS-CNV-workflow/merge_erds_results.py -i ./calls/ -a ids.map -o ./merged/ -g $gapsRef
     """
   }

ch_combined_calls = ch_cnvnator_clustercalls.join(ch_erds_clustercalls)

/*
* STEP 5d - Combine CNVnator and ERDS calls
*/
process combineCalls {
  tag "$sampID"
  publishDir "${params.outdir}/CNVs/TXT/Raw/", mode: 'copy'

  input:
  set sampID, file(cnvnatorCall), file(erdsCall) from ch_combined_calls

  output:
  set sampID, file("${sampID}.ERDS_CNVnator_CNVs.raw.txt") into ch_final_calls

  """
  python2 /root/TCAG-WGS-CNV-workflow/add_features.py -i $cnvnatorCall -a $erdsCall -o ${sampID}.ERDS_CNVnator_CNVs.raw.txt -s $sampID -c 0 -p reciprocal
  """
}

ch_final_calls_sex = ch_final_calls.join(ch_sex_CNVs)

/*
* STEP 5e - Filter CNV calls
* - Select common CNVs between CNVnator and erds (reciprocal overlap > 50%)
* - Discard CNVs with CNVnator q0 > 0.5
* - Remove CNVs overlapping > 70% with low complexity regions
*/
process filterCalls {
  tag "$sampID"

  publishDir "${params.outdir}/CNVs/TXT/Filtered/", mode: 'copy', pattern: '*.txt'
  publishDir "${params.outdir}/CNVs/log/", mode: 'copy', pattern: '*.log'

  input:
  set sampID, file(combCall), val(sex) from ch_final_calls_sex
  file(gapsRef) from gapsRef

  output:
  set sampID, file("${sampID}.ERDS_CNVnator_CNVs.filtered.txt") into filtered_calls, ch_cnvs_mosaichunter
  file("${sampID}.ERDS_CNVnator_CNVs.filtered.log") into filtered_log

  """
  filterCNVs.R $combCall $sex $gapsRef ${sampID}.ERDS_CNVnator_CNVs.filtered.txt
  """
}

// Convert CNVnator output to VCF
process convert2VCF {

  tag "$sampID"

  input:
  set sampID, file(merged) from filtered_calls

  output:
  tuple sampID, file("${sampID}.CNVs.vcf") into vcfs

  """
  ## Make Header
echo "##fileformat=VCFv4.3
##fileDate=$date
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

// Compress, sort and index vcf
process prepareVCF {

  tag "$sampID"

  input:
  tuple sampID, file(vcf) from vcfs
  tuple file(fasta), file(fastaidx), file(dict) from ch_fastadict_VCF.collect()

  output:
  tuple sampID, file("${vcf}.gz"), file("${vcf}.gz.tbi")  into sortVCF

  """
  bcftools reheader -f $fastaidx -o header.vcf $vcf
  bcftools sort -o ${vcf}.gz -O z  header.vcf
  tabix -p vcf ${vcf}.gz
  """
}


// Add annotation to CNVs
process annotateVCF {

  tag "$sampID"
  publishDir "${params.outdir}/CNVs/VCF/", mode: 'copy'

  input:
  set sampID, file(vcf), file(vcftbi) from sortVCF
  file(commonCNV)
  file(clinvarCNV)
  file(gtfRef)
  file(omim)

  output:
  set sampID, file("${sampID}.ERDS_CNVnator_CNVs.annotated.vcf.gz"), file("${sampID}.ERDS_CNVnator_CNVs.annotated.vcf.gz.tbi")  into annotatedVCF

  """
  annotateCNVs.R $vcf $commonCNV $clinvarCNV $gtfRef $omim ${sampID}.ERDS_CNVnator_CNVs.annotated.vcf
  bgzip ${sampID}.ERDS_CNVnator_CNVs.annotated.vcf
  tabix -p vcf ${sampID}.ERDS_CNVnator_CNVs.annotated.vcf.gz
  """

}


/*
* Prioritize variants:
* - Remove CNVs with overlap > 20% with commonCNVs
* - Select CNVs with overlap > 80% pathogenic variants (subset1)
* - Select CNVs overlapping exons in OMIM genes (subset2)
* - Select CNVs overlapping exons in GENCODE genes (subset3)
*/
process prioritizeCNVs {

  tag "$sampID"
  publishDir "${params.outdir}/CNVs/XLSX/", mode: 'copy'


  input:
  set sampID, file(vcf), file(vcftbi) from annotatedVCF

  output:
  file("${sampID}.ERDS_CNVnator_CNVs.Prioritization.xlsx") into priorCNV

  script:
  """
  prioritizeCNVs.R $vcf ${sampID}.ERDS_CNVnator_CNVs.Prioritization.xlsx
  """

}

/*
 *
 */
 process prepareINDELs {
   tag "$sampID"

   input:
   tuple val(sampID), file(bam), file(bai), file(vcf), file(vcf_idx) from ch_input_indels

   output:
   tuple sampID, file("${sampID}.indels.vcf.gz") into ch_vcf_indels

   """
   bcftools view -v indels  -f .,PASS $vcf -o ${sampID}.indels.vcf.gz -O z
   """

 }

 process indelsToBed {

   tag "$sampID"

   input:
   tuple sampID, file(vcf) from ch_vcf_indels

   output:
   set sampID, file("${sampID}.indels.bed") into ch_indels_bed

   """
   ## Decompress VCF on the fly. Add +/- 5 bp as recommended in MosaicHunter
   awk  '{OFS = "\t"}!/^#/{print \$1, \$2 - 5, \$2 + (length(\$2)) + 5}' <(gzip -cd $vcf) > ${sampID}.indels.bed
   """

 }

 process cnvsToBed {

   tag "$sampID"

   input:
   tuple sampID, file(txt) from ch_cnvs_mosaichunter

   output:
   tuple sampID, file("${sampID}.cnvs.bed") into ch_cnvs_bed

   """
   cut -f2-4 $txt | tail -n +2 - > ${sampID}.cnvs.bed
   """

 }

ch_comb = ch_indels_bed.join(ch_cnvs_bed)

process mergeBeds {

  tag "$sampID"

  input:
  tuple sampID, file(indels), file(cnvs) from ch_comb

  output:
  tuple sampID, file("${sampID}.merged.bed") into ch_merged_bed

  """
  cp $indels ${sampID}.merged.bed
  cat $cnvs >> ${sampID}.merged.bed
  """

}

ch_input_mosaichunter_input = ch_input_mosaichunter.join(ch_sex_mosaics).join(ch_merged_bed)

process runMosaicHunter {

  tag "$sampID"
  label 'process_medium'
  publishDir "${params.outdir}/Mosaics/SNVs/TXT/unfilteredVariants/", pattern: '*.tsv', mode: 'copy', saveAs: { filename -> "${sampID}.MosaicHunter.MosaicSNVs.txt" }
  publishDir "${params.outdir}/Mosaics/SNVs/logs/runMosaicHunter/", pattern: '*.log', mode: 'copy', saveAs: { filename -> "${sampID}.MosaicHunter.MosaicSNVs.log" }

  input:
  tuple val(sampID), file(bam), file(bai), file(vcf), file(vcf_idx), val(sex), file(bed) from ch_input_mosaichunter_input
  file(fasta) from ch_fasta_for_mosaichunter.collect()
  file(conf) from mosaichunter_config
  val(mode) from params.mosaichunter_mode
  file(gapsRef) from gapsRef


  output:
  tuple sampID, file("final.passed.tsv") into ch_mosaichunter_out, ch_mosaichunter_out2
  file("stdout*.log")

  script:
  """
  java -jar ~/MosaicHunter/build/mosaichunter.jar -C $conf \
  -P input_file=$bam \
  -P reference_file=$fasta \
  -P output_dir=./ \
  -P mosaic_filter.sex=$sex \
  -P mosaic_filter.mode=$mode \
  -P repetitive_region_filter.bed_file=$gapsRef \
  -P indel_region_filter.bed_file=$bed
  """

}


process mosaicHunterVCF {

  tag "$sampID"

  publishDir "${params.outdir}/Mosaics/SNVs/VCF/",  mode: 'copy'


  input:
  tuple sampID, file(tsv) from ch_mosaichunter_out2

  output:
  tuple sampID, file("${sampID}.MosaicHunter.MosaicSNVs.vcf.gz") into ch_mosaichunter_vcf

  script:
  """
  ## Make Header
echo "##fileformat=VCFv4.3
##fileDate=$date
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
   bgzip ${sampID}.MosaicHunter.MosaicSNVs.vcf
  """
}

process indexVCF {

  tag "$sampID"

  input:
  tuple sampID, file(vcf) from ch_mosaichunter_vcf

  output:
  tuple file(vcf), file("${vcf}.tbi") into ch_mosaichunter_tbi

  script:
  """
  tabix -p vcf $vcf
  """
}


process mergeVCFs {

  input:
  file(vcf) from ch_mosaichunter_tbi.collect()

  output:
  file("merged.MosaicHunter.MosaicSNVs.vcf.gz") into ch_mosaichunter_vcf_merge

  script:
  """
  bcftools merge *.gz -m none -o merged.MosaicHunter.MosaicSNVs.vcf.gz -O z
  """
}

process selectNonPrivateVariants {

  input:
  file(vcf) from ch_mosaichunter_vcf_merge

  output:
  file("nonprivate.MosaicHunter.MosaicSNVs.vcf.gz") into ch_mosaichunter_vcf_nonprivate

  script:
  """
  bcftools view -i 'N_SAMPLES-N_MISSING>1'  $vcf -o nonprivate.MosaicHunter.MosaicSNVs.vcf.gz -O z
  """
}




process convertAnnovarinput {

  tag "$sampID"

  input:
  tuple sampID, file(tsv) from ch_mosaichunter_out

  output:
  tuple sampID, file("${sampID}.avinput.txt") into ch_mosaichunter_avinput

  script:
  """
awk  '{OFS = "\t"}
{if (\$3=\$7)
      print \$1, \$2, \$2, \$3, \$9, \$10/\$4, \$4, \$8"/"\$10;
else
    print \$1, \$2,  \$2, \$3, \$7,  \$8/\$4, \$4, \$10"/"\$8;
   }'    $tsv > ${sampID}.avinput.txt
  """

}

process annotateSNPs {

  tag "$sampID"

  input:
  tuple sampID, file(avinput) from ch_mosaichunter_avinput
  file(annovar)
  file(annovarVar)
  file(annovarCod)
  file(annovarXref)
  file(annovarFold)

  output:
  tuple sampID, file("${sampID}.hg19_multianno.txt") into ch_mosaichunter_annot

  script:
  """
  perl $annovar $avinput $annovarFold -buildver hg19 -out ${sampID} -remove \
  --xref $annovarXref \
  -protocol refGene,cytoBand,genomicSuperDups,gnomad211_exome,gnomad211_genome,avsnp150,kaviar_20150923,clinvar_20200316,dbnsfp41a \
  -operation gx,r,r,f,f,f,f,f,f -nastring . --otherinfo
  """
}


process prioritizeSNVs {

  tag "$sampID"
  publishDir "${params.outdir}/Mosaics/SNVs/XLSX/", pattern: '*.xlsx', mode: 'copy'
  publishDir "${params.outdir}/Mosaics/SNVs/logs/prioritization/", pattern: '*.log', mode: 'copy'
  publishDir "${params.outdir}/Mosaics/SNVs/TXT/filteredVariants/", pattern: '*.txt', mode: 'copy'

  input:
  tuple sampID, file(annovar) from ch_mosaichunter_annot
  file(omim)
  file(omim_map)
  file(nonprivate) from ch_mosaichunter_vcf_nonprivate


  output:
  file("*.txt")
  file("*.log")
  file("*.xlsx")

  script:
  """
  prioritizeSNVs.R $annovar $omim $omim_map $nonprivate ${sampID}.MosaicHunter.MosaicSNVs.Prioritization
  """
}



// /*
//  * STEP 3 - Output Description HTML
//  */
// process output_documentation {
//
//     publishDir "${params.outdir}/pipeline_info", mode: 'copy'
//
//     input:
//     file output_docs from ch_output_docs
//
//     output:
//     file "results_description.html"
//
//     script:
//     """
//     markdown_to_html.py $output_docs -o results_description.html
//     """
// }
//
// /*
//  * Completion e-mail notification
//  */
// workflow.onComplete {
//
//     // Set up the e-mail variables
//     def subject = "[nf-core/cnvcalling] Successful: $workflow.runName"
//     if (!workflow.success) {
//         subject = "[nf-core/cnvcalling] FAILED: $workflow.runName"
//     }
//     def email_fields = [:]
//     email_fields['version'] = workflow.manifest.version
//     email_fields['runName'] = custom_runName ?: workflow.runName
//     email_fields['success'] = workflow.success
//     email_fields['dateComplete'] = workflow.complete
//     email_fields['duration'] = workflow.duration
//     email_fields['exitStatus'] = workflow.exitStatus
//     email_fields['errorMessage'] = (workflow.errorMessage ?: 'None')
//     email_fields['errorReport'] = (workflow.errorReport ?: 'None')
//     email_fields['commandLine'] = workflow.commandLine
//     email_fields['projectDir'] = workflow.projectDir
//     email_fields['summary'] = summary
//     email_fields['summary']['Date Started'] = workflow.start
//     email_fields['summary']['Date Completed'] = workflow.complete
//     email_fields['summary']['Pipeline script file path'] = workflow.scriptFile
//     email_fields['summary']['Pipeline script hash ID'] = workflow.scriptId
//     if (workflow.repository) email_fields['summary']['Pipeline repository Git URL'] = workflow.repository
//     if (workflow.commitId) email_fields['summary']['Pipeline repository Git Commit'] = workflow.commitId
//     if (workflow.revision) email_fields['summary']['Pipeline Git branch/tag'] = workflow.revision
//     email_fields['summary']['Nextflow Version'] = workflow.nextflow.version
//     email_fields['summary']['Nextflow Build'] = workflow.nextflow.build
//     email_fields['summary']['Nextflow Compile Timestamp'] = workflow.nextflow.timestamp
//
//     // TODO nf-core: If not using MultiQC, strip out this code (including params.max_multiqc_email_size)
//     // On success try attach the multiqc report
//     def mqc_report = null
//     try {
//         if (workflow.success) {
//             mqc_report = ch_multiqc_report.getVal()
//             if (mqc_report.getClass() == ArrayList) {
//                 log.warn "[nf-core/cnvcalling] Found multiple reports from process 'multiqc', will use only one"
//                 mqc_report = mqc_report[0]
//             }
//         }
//     } catch (all) {
//         log.warn "[nf-core/cnvcalling] Could not attach MultiQC report to summary email"
//     }
//
//     // Check if we are only sending emails on failure
//     email_address = params.email
//     if (!params.email && params.email_on_fail && !workflow.success) {
//         email_address = params.email_on_fail
//     }
//
//     // Render the TXT template
//     def engine = new groovy.text.GStringTemplateEngine()
//     def tf = new File("$baseDir/assets/email_template.txt")
//     def txt_template = engine.createTemplate(tf).make(email_fields)
//     def email_txt = txt_template.toString()
//
//     // Render the HTML template
//     def hf = new File("$baseDir/assets/email_template.html")
//     def html_template = engine.createTemplate(hf).make(email_fields)
//     def email_html = html_template.toString()
//
//     // Render the sendmail template
//     def smail_fields = [ email: email_address, subject: subject, email_txt: email_txt, email_html: email_html, baseDir: "$baseDir", mqcFile: mqc_report, mqcMaxSize: params.max_multiqc_email_size.toBytes() ]
//     def sf = new File("$baseDir/assets/sendmail_template.txt")
//     def sendmail_template = engine.createTemplate(sf).make(smail_fields)
//     def sendmail_html = sendmail_template.toString()
//
//     // Send the HTML e-mail
//     if (email_address) {
//         try {
//             if (params.plaintext_email) { throw GroovyException('Send plaintext e-mail, not HTML') }
//             // Try to send HTML e-mail using sendmail
//             [ 'sendmail', '-t' ].execute() << sendmail_html
//             log.info "[nf-core/cnvcalling] Sent summary e-mail to $email_address (sendmail)"
//         } catch (all) {
//             // Catch failures and try with plaintext
//             [ 'mail', '-s', subject, email_address ].execute() << email_txt
//             log.info "[nf-core/cnvcalling] Sent summary e-mail to $email_address (mail)"
//         }
//     }
//
//     // Write summary e-mail HTML to a file
//     def output_d = new File("${params.outdir}/pipeline_info/")
//     if (!output_d.exists()) {
//         output_d.mkdirs()
//     }
//     def output_hf = new File(output_d, "pipeline_report.html")
//     output_hf.withWriter { w -> w << email_html }
//     def output_tf = new File(output_d, "pipeline_report.txt")
//     output_tf.withWriter { w -> w << email_txt }
//
//     c_green = params.monochrome_logs ? '' : "\033[0;32m";
//     c_purple = params.monochrome_logs ? '' : "\033[0;35m";
//     c_red = params.monochrome_logs ? '' : "\033[0;31m";
//     c_reset = params.monochrome_logs ? '' : "\033[0m";
//
//     if (workflow.stats.ignoredCount > 0 && workflow.success) {
//         log.info "-${c_purple}Warning, pipeline completed, but with errored process(es) ${c_reset}-"
//         log.info "-${c_red}Number of ignored errored process(es) : ${workflow.stats.ignoredCount} ${c_reset}-"
//         log.info "-${c_green}Number of successfully ran process(es) : ${workflow.stats.succeedCount} ${c_reset}-"
//     }
//
//     if (workflow.success) {
//         log.info "-${c_purple}[nf-core/cnvcalling]${c_green} Pipeline completed successfully${c_reset}-"
//     } else {
//         checkHostname()
//         log.info "-${c_purple}[nf-core/cnvcalling]${c_red} Pipeline completed with errors${c_reset}-"
//     }
//
// }
//
//
// def nfcoreHeader() {
//     // Log colors ANSI codes
//     c_black = params.monochrome_logs ? '' : "\033[0;30m";
//     c_blue = params.monochrome_logs ? '' : "\033[0;34m";
//     c_cyan = params.monochrome_logs ? '' : "\033[0;36m";
//     c_dim = params.monochrome_logs ? '' : "\033[2m";
//     c_green = params.monochrome_logs ? '' : "\033[0;32m";
//     c_purple = params.monochrome_logs ? '' : "\033[0;35m";
//     c_reset = params.monochrome_logs ? '' : "\033[0m";
//     c_white = params.monochrome_logs ? '' : "\033[0;37m";
//     c_yellow = params.monochrome_logs ? '' : "\033[0;33m";
//
//     return """    -${c_dim}--------------------------------------------------${c_reset}-
//                                             ${c_green},--.${c_black}/${c_green},-.${c_reset}
//     ${c_blue}        ___     __   __   __   ___     ${c_green}/,-._.--~\'${c_reset}
//     ${c_blue}  |\\ | |__  __ /  ` /  \\ |__) |__         ${c_yellow}}  {${c_reset}
//     ${c_blue}  | \\| |       \\__, \\__/ |  \\ |___     ${c_green}\\`-._,-`-,${c_reset}
//                                             ${c_green}`._,._,\'${c_reset}
//     ${c_purple}  nf-core/cnvcalling v${workflow.manifest.version}${c_reset}
//     -${c_dim}--------------------------------------------------${c_reset}-
//     """.stripIndent()
// }
//
// def checkHostname() {
//     def c_reset = params.monochrome_logs ? '' : "\033[0m"
//     def c_white = params.monochrome_logs ? '' : "\033[0;37m"
//     def c_red = params.monochrome_logs ? '' : "\033[1;91m"
//     def c_yellow_bold = params.monochrome_logs ? '' : "\033[1;93m"
//     if (params.hostnames) {
//         def hostname = "hostname".execute().text.trim()
//         params.hostnames.each { prof, hnames ->
//             hnames.each { hname ->
//                 if (hostname.contains(hname) && !workflow.profile.contains(prof)) {
//                     log.error "====================================================\n" +
//                             "  ${c_red}WARNING!${c_reset} You are running with `-profile $workflow.profile`\n" +
//                             "  but your machine hostname is ${c_white}'$hostname'${c_reset}\n" +
//                             "  ${c_yellow_bold}It's highly recommended that you use `-profile $prof${c_reset}`\n" +
//                             "============================================================"
//                 }
//             }
//         }
//     }
// }
//
// Check file extension
def hasExtension(it, extension) {
    it.toString().toLowerCase().endsWith(extension.toLowerCase())
}
