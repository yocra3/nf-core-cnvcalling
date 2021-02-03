#! /opt/conda/envs/nf-core-cnvcalling/bin/Rscript
#'#################################################################################
#'#################################################################################
#' Prioritize SNVs using annovar data
#' Remove SNVs that fulfill these criteria:
#' - Common SNVs: present with AF > 0.001 in any population
#' - Segmental Duplications: Remove variants present in segmental Duplications
#' - Low read depth: read depth < 10. Important: read depth should be encoded in OtherInfo2 -> 7th column (2nd extra column) of annovar input
#' - Frequent in our cohort: More than 2 occurences in our cohort
#' Select SNVs:
#' - In exonic or splicing positions
#' - Highly deleterous (frameshift, non-sense) or non-synonimous
#' Generate Excel tables. For each source (but clinVar), two tables are generated:
#'     - Dominant model (Heterozygous variants)
#'     - Recessive model (double heterozygous and homozygous)
#' - clinVar: Variants indentified in clinVar as pathogenic or likely pathogenic
#' - OMIM: Genes present in OMIM
#' - Candidate genes: Genes with high pLI (pLI > 0.8) or high pREC (pREC > 0.8)
#' - Remaining genes: Variants in other genes
#' A txt with all variants passing selection will be generated (to be used for ORVAL)
#'#################################################################################
#'#################################################################################

## Capture arguments
args <- commandArgs(trailingOnly=TRUE)
annovarFile <- args[1]
omimGenes <- args[2]
omimMap <- args[3]
cohortVCF <- args[4]
outPrefix <- args[5]

### Parameters (move to a config file?)
AF_threshold <- 0.001
min_readDepth <- 10
max_internal_freq <- 2
min_pLI <- 0.8
min_pREC <- 0.8

# Load libraries
library(VariantAnnotation)
library(openxlsx)
library(tidyverse)

## Load SNVs file and modify it:
### Create column with SNV name matching vcfs
### Add a column with genotypes (homogeneize between mosaics and germinal SNVs)
#### Consider Homozygous for AF > 0.8
ini <- read_tsv(annovarFile, na = ".") %>%
  mutate(ID = paste0(Chr, ":", Start, "_", Ref, "/", Alt),
         genotype = ifelse(Otherinfo1 == "het" | Otherinfo1 < 0.8,
         "Heterozygous", "Homozygous")
)

## Load interal freq SNVs
comSNVs <- readVcf(cohortVCF)
snpSum <- snpSummary(comSNVs)

if (nrow(snpSum) > 1){
  selSNPs <- snpSum %>%
    filter(g00 + g01 + g11 > max_internal_freq) %>%
    rownames()
  ini.snp <- filter(ini, !ID %in% selSNPs)
} else {
  ini.snp <- ini
}

## Segmental Duplications
ini.segDup <- subset(ini.snp, is.na(genomicSuperDups))

## Low read depths
ini.depth <- subset(ini.segDup, Otherinfo2 > min_readDepth)

## Common in Population
af_cols <- colnames(ini)[grep("AF", colnames(ini))] ### Columns with AF information
af_cols <- af_cols[-grep("Bayes", af_cols)] ## Remove Bayes measures
ini.com <- ini.depth %>%
  filter_at(af_cols, all_vars(is.na(.) | . < AF_threshold))

## Select exonic or splicing Variants
ini.exon <- filter(ini.com, Func.refGene %in% c("exonic", "splicing"))

## Discard synonymous variants
ini.del <- filter(ini.exon, !( !is.na(ExonicFunc.refGene) & ExonicFunc.refGene %in% c("synonymous SNV", "unknown")))

## Create table for ORVAL
orval <- ini.del %>%
  select(Chr, Start, Ref, Alt, genotype)

write.table(orval, file = paste0(outPrefix, ".selVariants.txt"), quote = FALSE,
  row.names = FALSE)



# OMIM ####
## Create table with annovar gene ID, OMIM ID and OMIM phenotype
omim_match <- read_tsv(omimMap) %>%
  mutate(Gene.refGene = `Approved symbol`,
         OMIM_ID = `OMIM ID(supplied by OMIM)`) %>%
  dplyr::select(Gene.refGene, OMIM_ID)
omim <- read_tsv(omimGenes, skip = 3) %>%
  mutate(OMIM_ID = as.character(`MIM Number`),
         OMIM_Phenotype = Phenotypes) %>%
  dplyr::select(OMIM_ID, OMIM_Phenotype) %>%
  filter(!is.na(OMIM_ID) & !is.na(OMIM_Phenotype)) %>%
  inner_join(omim_match, by = "OMIM_ID")

## Add OMIM columns to all variants annotation and define categories
vars.annot <- left_join(ini.del, omim, by = "Gene.refGene") %>%
  mutate(pLI_flag = !is.na(pLi.refGene) & pLi.refGene > min_pLI,
         misZ_flag = !is.na(pRec.refGene) & pRec.refGene > min_pREC,
         cand_flag = pLI_flag | misZ_flag,
         clinVar_flag = !is.na(CLNSIG) & grepl("Pathogenic|Likely_pathogenic", CLNSIG),
         prior_tab = ifelse(clinVar_flag, "clinVar",
                            ifelse(!is.na(OMIM_Phenotype), "OMIM",
                                   ifelse(cand_flag, "Candidate genes", "Other genes"))))

# clinVar table ####
clinvar <- filter(vars.annot, prior_tab == "clinVar")
ini.omim <- filter(vars.annot, prior_tab == "OMIM")

## Variants per gene
getGenesMultiVariants <- function(df, n = 2){
  df %>%
    dplyr::select(Gene.refGene) %>%
    group_by(Gene.refGene) %>%
    summarize(count = n()) %>%
    filter(count >= n) %>%
    pull(Gene.refGene)
}
omim_genes <- getGenesMultiVariants(ini.omim)

omim.dom <- filter(ini.omim, genotype == "Heterozygous" & !Gene.refGene %in% omim_genes)
omim.rec <- filter(ini.omim, (genotype == "Heterozygous" & Gene.refGene %in% omim_genes) |
                     genotype == "Homozygous")

# Candidate genes ####
ini.cand <- filter(vars.annot, prior_tab == "Candidate genes")
cand_genes <- getGenesMultiVariants(ini.cand)

cand.dom <- filter(ini.cand, genotype == "Heterozygous" & !Gene.refGene %in% cand_genes)
cand.rec <- filter(ini.cand, (genotype == "Heterozygous" & Gene.refGene %in% cand_genes) |
                     genotype == "Homozygous")

# Remaining genes
ini.rest <- filter(vars.annot, prior_tab == "Other genes")
rest_genes <- getGenesMultiVariants(ini.rest)

rest.dom <- filter(ini.rest, genotype == "Heterozygous" & !Gene.refGene %in% rest_genes)
rest.rec <- filter(ini.rest, (genotype == "Heterozygous" & Gene.refGene %in% rest_genes) |
                     genotype == "Homozygous")



## Create variant selection log
sumTable <- data.frame(Description = c(
  "Initial Number SNVs",
  paste("SNVs with internal Freq <", max_internal_freq),
  "SNVs not in Segmental Duplications",
  paste("SNVs with read depth >", min_readDepth),
  "SNVs in exonic or splicing positions",
  "Splicing, Frameshift, non-sense or non-synonymous SNVs",
  "Pathogenic or likely pathogenic in clinVar",
  "SNVs in OMIM genes",
  "SNVs in intolerant genes",
  "SNVs in remaining genes"),
  Number = c(nrow(ini), nrow(ini.snp), nrow(ini.segDup), nrow(ini.depth),
             nrow(ini.exon), nrow(ini.del), nrow(clinvar),
             nrow(ini.omim), nrow(ini.cand), nrow(ini.rest)))
sumTable$Proportion <- round(sumTable$Number/nrow(ini)*100, 2)
write.table(sumTable, file = paste0(outPrefix, ".log"), quote = FALSE,
            row.names = FALSE)


write.xlsx(list(clinvar, omim.dom, omim.rec, cand.dom, cand.rec,
                rest.dom, rest.rec, sumTable),
           file = paste0(outPrefix, ".xlsx"),
           rowNames = FALSE,
           colNames = TRUE,
           sheetName = c("ClinVar Pathogenic",
                         "OMIM genes - Dominant",
                         "OMIM genes - Recessive",
                         "Intolerant genes - Dominant",
                         "Intolerant genes - Recessive",
                         "Other genes - Dominant",
                         "Other genes - Recessive",
                         "Prioritization summary"))
