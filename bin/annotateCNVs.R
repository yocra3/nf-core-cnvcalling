#! /opt/conda/envs/nf-core-cnvcalling/bin/Rscript

#'#################################################################################
#'#################################################################################
#' Annotate CNVs detected from CNVnator with the following fields:
#' - commonCNV: any overlap with common calls from nstd186 (AF > 0.01)
#' - commonCNV20: reference CNV contains >20% sample CNV. Common calls are from nstd186 (AF > 0.01)
#' - commonCNV80: reference CNV contains >80% sample CNV. Common calls are from nstd186 (AF > 0.01)
#' - clinCNV: any overlap with ClinVar calls from nstd102
#' - beningCNV20: overlap with ClinVar benign or likely benign calls from nstd102. Reference CNV contains >20% sample CNV.
#' - beningCNV80: overlap with ClinVar benign or likely benign calls from nstd102. Reference CNV contains >80% sample CNV.
#' - pathoCNV20: overlap with ClinVar pathogenic or likely pathogenic calls from nstd102. Sample CNV contains >20% reference CNV.
#' - pathoCNV80: overlap with ClinVar pathogenic or likely pathogenic calls from nstd102. Sample CNV contains >80% reference CNV.
#' - gencodeGENES: Overlap with GENCODE genes
#' - gencodeEXONS: Overlap with exons
#' - omimGENES: Overlap with OMIM genes
#'#################################################################################
#'#################################################################################

## Capture arguments
args <- commandArgs(trailingOnly=TRUE)
cnvFile <- args[1]
commonCNVFile <- args[2]
clinvarCNVFile <- args[3]
gencodeFile <- args[4]
omimFile <- args[5]
outVCF <- args[6]


## Load libraries
library(VariantAnnotation)
library(rtracklayer)

## Load data
cnv <- readVcf(cnvFile)
common <- readVcf(commonCNVFile)
clinvar <- readVcf(clinvarCNVFile)

gencodeGR <- readGFFAsGRanges(gencodeFile)
omim <- read.table(gzfile(omimFile), sep = "\t", quote = "", as.is = TRUE)
colnames(omim) <- c("chromosome", "start", "end", "cytogenetic_loc", "cytogenetic_loc2",
                    "MIM_Number", "Symbols", "Gene_Name", "Approved_Symbol",
                    "EntrezID", "EnsemblID", "Comments", "Phenotypes", "Mouse Gene Symbol/ID")


createGRfromVCF <- function(vcf){

  cnvGR <- rowRanges(vcf)
  mcols(cnvGR) <- info(vcf)
  cnvGR$POS <- start(cnvGR)
  start(cnvGR) <- pmin(cnvGR$POS, cnvGR$END)
  end(cnvGR) <- pmax(cnvGR$POS, cnvGR$END)
  cnvGR
}

## Convert VCF files to GenomicRanges and add end position
cnvGR <- createGRfromVCF(cnv)
commonGR <- createGRfromVCF(common)
clinvarGR <- createGRfromVCF(clinvar)

comOver <- findOverlaps(cnvGR, commonGR)
ovList <- split(names(commonGR)[to(comOver)], names(cnvGR)[from(comOver)])
cnvGR$commonCNV <- ovList[names(cnvGR)]

## Define functions to compute overlaps ####
getOverlapRefSample <- function(refCNVGR, cnvVec, overlap){

  ov <- findOverlaps(cnvVec, refCNVGR, minoverlap = width(refCNVGR)*overlap)

  if (length(ov) > 0) {
    data.frame(cnv = names(cnvVec)[from(ov)], ref = rep(names(refCNVGR), length(ov)))
  }
}

getOverlapSampleRef <- function(cnvSamp, refCNVs, overlap){

  ov <- findOverlaps(cnvSamp, refCNVs, minoverlap = width(cnvSamp)*overlap)

  if (length(ov) > 0) {
    names(refCNVs)[to(ov)]
  } else {
    character(0)
  }
}

## Add overlaps with common CNVs ####
cnvGR$commonCNV20 <- lapply(seq_len(length(cnvGR)), function(i)
  getOverlapSampleRef(cnvGR[i], commonGR[cnvGR$commonCNV[[i]]], overlap = 0.2))

cnvGR$commonCNV80 <- lapply(seq_len(length(cnvGR)), function(i)
    getOverlapSampleRef(cnvGR[i], commonGR[cnvGR$commonCNV20[[i]]], overlap = 0.8))

## Add overlaps with clinvar CNVs ####
clinOver <- findOverlaps(cnvGR, clinvarGR)
cnvGR$clinCNV <- split(names(clinvarGR)[to(clinOver)], names(cnvGR)[from(clinOver)])[names(cnvGR)]

### Correct CLNSIG column
clinvarGR$CLNSIG <- unlist(clinvarGR$CLNSIG)
clinvarGR$CLNSIG <- gsub('\"', "", clinvarGR$CLNSIG)

## Benign/Likely Bening
beningGR <- clinvarGR[clinvarGR$CLNSIG %in% c("Benign", "Benign/Likely benign", "Likely benign")]
cnvGR$beningCNV20 <- lapply(seq_len(length(cnvGR)), function(i)
  getOverlapSampleRef(cnvGR[i], beningGR, overlap = 0.2))
cnvGR$beningCNV80 <- lapply(seq_len(length(cnvGR)), function(i)
    getOverlapSampleRef(cnvGR[i], beningGR[cnvGR$beningCNV20[[i]]], overlap = 0.8))

## Pathogenic/Likely pathogenic
pathogenicGR <- clinvarGR[clinvarGR$CLNSIG %in% c("Pathogenic", "Pathogenic/Likely pathogenic", "Likely pathogenic")]
cnvtmp <- cnvGR[lengths(cnvGR$clinCNV) > 0]
patOver20 <- lapply(seq_len(length(pathogenicGR)), function(i)
  getOverlapRefSample(pathogenicGR[i], cnvtmp, overlap = 0.20))
patOver20df <- Reduce(rbind, patOver20)
if(is.null(nrow(patOver20df))){
  cnvGR$pathoCNV20 <- character(1)
  cnvGR$pathoCNV80 <- character(1)
} else {
  cnvGR$pathoCNV20 <- split(patOver20df$ref, patOver20df$cnv)[names(cnvGR)]

  refCNVs <- pathogenicGR[unique(unlist(cnvGR$pathoCNV20))]
  cnvtmp <- cnvGR[lengths(cnvGR$pathoCNV20) > 0]
  patOver80 <- lapply(seq_len(length(refCNVs)), function(i)
    getOverlapRefSample(refCNVs[i], cnvtmp, overlap = 0.80))
  patOver80df <- Reduce(rbind, patOver80)

  if (is.null(nrow(patOver80df))){
    cnvGR$pathoCNV80 <- character(1)
  } else {
    cnvGR$pathoCNV80 <- split(patOver80df$ref, patOver80df$cnv)[names(cnvGR)]
  }
}


## Add overlaps with GENCODE genes ####
seqlevelsStyle(gencodeGR) <- "NCBI"
genesGR <- gencodeGR[gencodeGR$type == "gene"]
genesOver <- findOverlaps(cnvGR, genesGR)
genesList <- split(genesGR$gene_name[to(genesOver)], names(cnvGR)[from(genesOver)])
cnvGR$gencodeGENES <- genesList[names(cnvGR)]

exonsGR <- gencodeGR[gencodeGR$type == "exon"]
exonsOver <- findOverlaps(cnvGR, exonsGR)
exonsList <- split(exonsGR$gene_name[to(exonsOver)], names(cnvGR)[from(exonsOver)])
exonsList <- lapply(exonsList, unique)
cnvGR$gencodeEXONS <- exonsList[names(cnvGR)]

## Add overlaps with OMIM genes ####
### Use coordinates from gencode based on ENSEMBLIDs
omimGenes <- unique(as.character(omim$EnsemblID))
genesGR$ensemblID <- gsub("\\..*$", "", genesGR$gene_id)
omimGR <- genesGR[genesGR$ensemblID %in% omimGenes]
omimOver <- findOverlaps(cnvGR, omimGR)
omimList <- split(omimGR$gene_name[to(omimOver)], names(cnvGR)[from(omimOver)])
cnvGR$omimGENES <- omimList[names(cnvGR)]

omimList.ensem <- split(omimGR$ensemblID[to(omimOver)], names(cnvGR)[from(omimOver)])
omimList.pheno <- lapply(omimList.ensem, function(x) omim[omim$EnsemblID %in% x, "Phenotypes"])
cnvGR$omimPHENO <- omimList.pheno[names(cnvGR)]


## Output vcfs
cnvGR$POS <- NULL

### Convert fields to character
mcolsdf <- mcols(cnvGR)
newCols <- c("commonCNV", "commonCNV20", "commonCNV80", "clinCNV",
  "beningCNV20", "beningCNV80", "pathoCNV20", "pathoCNV80",
  "gencodeGENES", "gencodeEXONS", "omimGENES", "omimPHENO")
pasteCol <- function(col) vapply(col, paste, collapse = "|", character(1))
for(var in newCols) mcolsdf[[var]] <- pasteCol(mcolsdf[[var]])
mcolsdf$omimPHENO <- gsub(";", "|", mcolsdf$omimPHENO) ## Remove ; to avoid problems when parsing the VCF
info(cnv) <- mcolsdf


newFields <- DataFrame(Number = "1", Type = "String",
Description = c("any overlap with common calls from nstd186 (AF > 0.01)",
  "reference CNV contains >20% sample CNV. Common calls are from nstd186 (AF > 0.01)",
  "reference CNV contains >80% sample CNV. Common calls are from nstd186 (AF > 0.01)",
  "any overlap with ClinVar calls from nstd102",
  "overlap with ClinVar benign or likely benign calls from nstd102. Reference CNV contains >20% sample CNV",
  "overlap with ClinVar benign or likely benign calls from nstd102. Reference CNV contains >80% sample CNV",
  "overlap with ClinVar pathogenic or likely pathogenic calls from nstd102. Sample CNV contains >20% reference CNV.",
  "overlap with ClinVar pathogenic or likely pathogenic calls from nstd102. Sample CNV contains >80% reference CNV",
  "overlap with GENCODE genes",
  "overlap with exons",
  "overlap with OMIM genes",
  "phenotypes of OMIM genes"))
rownames(newFields) <- newCols
info(header(cnv)) <- rbind(info(header(cnv)), newFields)

writeVcf(cnv, filename = outVCF)
