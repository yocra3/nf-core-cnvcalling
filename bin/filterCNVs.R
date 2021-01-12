#! /opt/conda/envs/nf-core-cnvcalling/bin/Rscript
#'#################################################################################
#'#################################################################################
#' Filter CNV calls
#' - Select common CNVs between CNVnator and erds  (reciprocal overlap > 50%)
#' - Discard CNVs with CNVnator q0 > 0.5
#' - Remove CNVs overlapping > 70% with low complexity regions
#'#################################################################################
#'#################################################################################

## Capture arguments
args <- commandArgs(trailingOnly=TRUE)
cnvFile <- args[1]
LCRFile <- args[2]
outFile <- args[3]

## Load libraries
library(GenomicRanges)

## Load data
cnv <- read.delim(cnvFile, as.is = TRUE)
LCRs <- read.table(LCRFile, comment.char = "", header = FALSE, as.is = TRUE)

## Select common CNVs between ERDS and CNVnator 
cnv.com1 <- subset(cnv, erds_fraction != "-" & erds_fraction != "#" & cnvn_fraction != "-" & cnv_type_conflict == "-")
## Reciprocal overlap > 50%
cnv.com <- subset(cnv.com1, erds_fraction > 0.50 & cnvn_fraction > 0.50)

## Filter CNVs with q0 > 0.5
cnv.filt <- subset(cnv.com, !(q0 != "-" & as.numeric(q0) > 0.5))

## Remove CNVs overlapping > 70% with low complexity regions
cnvGR <- makeGRangesFromDataFrame(cnv.filt, seqnames.field = "m")
seqlevelsStyle(cnvGR) <- "UCSC"
LCRGR <- makeGRangesFromDataFrame(LCRs, seqnames.field = "V1", 
                                  start.field = "V2", end.field = "V3")
seqlevelsStyle(LCRGR) <- "UCSC"
overs <- findOverlaps(cnvGR, LCRGR)

computeTotalOverlap <- function(ov, GR1, GR2){
  
  overGR <- intersect(GR1[from(ov)], GR2[to(ov)])
  sum(width(overGR))/width(GR1[unique(from(ov))])
}
cnv.Over <- data.frame(ID = unique(from(overs)))
cnv.Over$prop <- vapply(cnv.Over$ID, 
                        function(x) computeTotalOverlap(overs[from(overs) == x], cnvGR, LCRGR),
                        numeric(1))
badCNVs <- cnv.Over$ID[cnv.Over$prop > 0.7]
cnv.sing <- cnv.filt[-badCNVs, ]

## Output file
write.table(cnv.sing, file = outFile, quote = FALSE, row.names = FALSE, sep = "\t")