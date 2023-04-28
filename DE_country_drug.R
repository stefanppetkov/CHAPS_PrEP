library(DESeq2)
library(readxl)
library(tidyverse)
library(cowplot)

# load in counts data
countData <- readxl::read_xlsx("BEA21P098_rawCounts_rerun.xlsx") 
temp <- countData[c(2, 6:60669), c(1, grep("TRUE", unlist(countData[5,])))]
cols <- gsub(".*_", "", temp[1,2:142])
temp <- temp[2:60665,]
temp <- column_to_rownames(temp, names(temp)[1])
names(temp) <- cols
countData <- temp

# metadadta from unblinding (group allocation)
unblinding <- readxl::read_xlsx("unblinding_groups.xlsx")

metaData_drug <- unblinding %>% dplyr::filter(country == "SA") %>% dplyr::select(group, treatment, CODE) %>% 
dplyr::filter(!CODE == c("EB-063")) # potential outliers


metaData_drug <- metaData_drug %>% column_to_rownames("CODE") %>% dplyr::mutate(treatment = as.factor(treatment))
countData_drug <- countData[, intersect(rownames(metaData_drug), names(countData))] # find the status of missing SA samples

metaData_drug <- metaData_drug %>% rownames_to_column("CODE") %>%  dplyr::filter(CODE %in% intersect(rownames(metaData_drug), names(countData))) %>% 
  column_to_rownames("CODE")

countData_drug <- data.matrix(countData_drug)

dds_drug <- DESeqDataSetFromMatrix(countData = countData_drug,
                                         colData = metaData_drug,
                                         design =  ~ treatment)

# keep only genes for which 90% of samples have at least 1 count
nrow(dds_drug)
keep <- rowSums(counts(dds_drug) >= 2) >= dim(metaData_drug)[1]*0.9

dds_drug <- dds_drug[keep,]
nrow(dds_drug)

#### run Deseq2
dds_drug$group <- factor(dds_drug$treatment)
design(dds_drug) <- ~ group

dds_drug <- DESeq(dds_drug)
resultsNames(dds_drug)

# results function
resFunk <- function(dds, grps){
  results(dds, contrast = c("group", grps), tidy = T) %>% 
    dplyr::arrange(padj) %>% 
    dplyr::inner_join(grch38, by = c("row" = "ensgene")) %>% 
    dplyr::filter(biotype == "protein_coding", complete.cases(.)) %>%
    dplyr::mutate(lowCount = if_else(row %in% names(filtered[filtered == 'TRUE']), 'yes', 'no')) %>% # low count, filtered genes
    dplyr::mutate(significant = if_else(padj <= 0.05, 'yes', 'no')) 
}

groups_drug <- c("FTC-TDF", "FTC-TAF", "CT")
combinations_drug <- combn(groups_drug, 2, simplify = F)
combinations_drug[[1]] <- rev(combinations_drug[[1]])

results_drug <- list()
for (i in 1:length(combinations_drug)) {
  print(paste("creating results", i, "out of", length(combinations_drug)))
  out <- resFunk(dds_drug, combinations_drug[[i]])
  results_drug[[i]] <- out
}
# rename list elements
names_drug <- combinations_drug %>% map(~ paste0(.x, collapse = "vs"))
names(results_drug) <- names_drug


results_drug_significant <- map(results_drug, ~{dplyr::filter(., padj <= 0.05) 
  # tab(lowCount, significant)
})


df_results_u_significant <- bind_rows(results_drug_significant, .id = "comparison") %>% dplyr::distinct(pvalue, .keep_all = T) %>% dplyr::select(comparison, symbol, log2FoldChange, padj, lowCount)

df_results_sa_significant <- bind_rows(results_drug_significant, .id = "comparison") %>% dplyr::distinct(pvalue, .keep_all = T) %>% dplyr::select(comparison, symbol, log2FoldChange, padj, lowCount)


################# volcano plots###############
vol_u_drug_TDF_CT <- volc(results_drug[["FTC-TDFvsCT"]],
                               t = 'FTC-TDF vs CT')

vol_u_drug_TAF_CT <- volc(results_drug[["FTC-TAFvsCT"]],
                               t = 'FTC-TAF vs CT')


vol_sa_drug_TDF_CT <- volc(results_drug[["FTC-TDFvsCT"]],
                          t = 'FTC-TDF vs CT')

vol_sa_drug_TAF_TDF <- volc(results_drug[["FTC-TAFvsFTC-TDF"]],
                          t = 'FTC-TAF vs FTC-TDF')
                          
# save figure
# ggsave("figure_2_.pdf", 
#        cowplot::plot_grid(vol_sa_drug_TDF_CT, vol_sa_drug_TAF_TDF, vol_u_drug_TDF_CT, vol_u_drug_TAF_CT),
#        device = cairo_pdf, width = 10, height = 10, scale = 1.2)
