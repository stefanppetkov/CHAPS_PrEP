### Analysis of merged sites *by drug-dose*

metaData_merge_dose <- unblinding  %>% dplyr::select(treatment_dose, CODE) %>% 
  dplyr::filter(!CODE == c("EB-063"))
  # unite(drug_country, c("country", "treatment"))

metaData_merge_dose <- metaData_merge_dose %>% column_to_rownames("CODE") %>% dplyr::mutate(dose = as.factor(treatment_dose))
countData_merge_dose <- countData[, intersect(rownames(metaData_merge_dose), names(countData))] # 

metaData_merge_dose <- metaData_merge_dose %>% rownames_to_column("CODE") %>%  dplyr::filter(CODE %in% intersect(rownames(metaData_merge_dose), names(countData))) %>% 
  column_to_rownames("CODE")

countData_merge_dose <- data.matrix(countData_merge_dose)

dds_merge_dose <- DESeqDataSetFromMatrix(countData = countData_merge_dose,
                              colData = metaData_merge_dose,
                              design =  ~ dose)

# keep only genes for which 90% of samples have at least 1 count
nrow(dds_merge_dose)
keep <- rowSums(counts(dds_merge_dose) >= 2) >= dim(metaData_merge_dose)[1]*0.9
# keep <- c(merged_degs_drug, merged_ctrl_genes_drug)   # pooled differences from individ countries + control genes
filtered <- !keep
# table(keep)
dds_merge_dose <- dds_merge_dose[keep,]
nrow(dds_merge_dose)

#### run Deseq2
dds_merge_dose$group <- factor(dds_merge_dose$dose)
design(dds_merge_dose) <- ~ group

dds_merge_dose <- DESeq(dds_merge_dose)
resultsNames(dds_merge_dose)
# results function

resFunk <- function(dds, grps){
  results(dds, contrast = c("group", grps), tidy = T) %>% 
      dplyr::arrange(padj) %>% 
  dplyr::inner_join(grch38, by = c("row" = "ensgene")) %>% 
  dplyr::filter(biotype == "protein_coding", complete.cases(.)) %>%
  dplyr::mutate(lowCount = if_else(row %in% names(filtered[filtered == 'TRUE']), 'yes', 'no')) %>% # low count, filtered genes
  dplyr::mutate(significant = if_else(padj <= 0.05, 'yes', 'no')) 
}

groups_merge_dose <- rev(c("CT", "1xTDF", "1xTAF", "2xTDF", "2xTAF"))
combinations_merge_dose <- combn(groups_merge_dose, 2, simplify = F)

results_merge_dose <- list()
for (i in 1:length(combinations_merge_dose)) {
  print(paste("creating results", i, "out of", length(combinations_merge_dose)))
    out <- resFunk(dds_merge_dose, combinations_merge_dose[[i]])
    results_merge_dose[[i]] <- out
}
# rename list elements
names_merge_dose <- combinations_merge_dose %>% map(~ paste0(.x, collapse = "vs"))
names(results_merge_dose) <- names_merge_dose


results_merge_dose_significant <- map(results_merge_dose, ~{dplyr::filter(., padj <= 0.05) 
    # tab(lowCount, significant)
  })

# data for word cloud
df_results_merge_dose_significant <- bind_rows(results_merge_dose_significant, .id = "comparison") %>% dplyr::distinct(pvalue, .keep_all = T) %>% dplyr::select(comparison, symbol, log2FoldChange, padj, lowCount, row) %>% tab(symbol)

############# volcano plots ##################
volc <- function(d, t, lfc = 1){
  EnhancedVolcano(d,
                  lab = d %>% pull(symbol),
                  selectLab = d %>% dplyr::filter(padj <= 0.05, abs(log2FoldChange)> lfc) %>% pull(symbol),
                  labSize = 3,
                  x = 'log2FoldChange', 
                  y = 'padj', 
                  pCutoff = 0.05,
                  drawConnectors = T,
                  widthConnectors = 0.2,
                  maxoverlapsConnectors = 50,
                  FCcutoff = 1, 
                  col = c('grey', 'pink','#00509d', '#bf0603'),
                  title = as.character(t),
                  subtitle = '',
                  caption = '',
                  pointSize = 1,
                  gridlines.major = F,
                  gridlines.minor = F,
                  legendPosition = 'none',
                  legendLabSize = 5,
                  raster = T,
                  xlim = c(-5, 5),
                  ylim = c(0, 8)
  )
}

vol_merged_dose_a <- volc(results_merge_dose[["1xTDFvsCT"]],
                t = 'FTC-TDF (one dose) vs CT')

vol_merged_dose_b <- volc(results_merge_dose[["2xTDFvsCT"]],
                t = 'FTC-TDF (two doses) vs CT')

vol_merged_dose_c <- volc(results_merge_dose[["1xTAFvsCT"]],
                t = 'FTC-TAF (one dose) vs CT')

vol_merged_dose_d <- volc(results_merge_dose[["2xTAFvsCT"]],
                t = 'FTC-TAF (two doses) vs CT')


# save figure
ggsave("merged_sites_by_dose_volcanos__.pdf",
       cowplot::plot_grid(vol_merged_dose_a, vol_merged_dose_b, vol_merged_dose_c, vol_merged_dose_d,  nrow = 2),
       device = cairo_pdf, 
       width = 12, height = 12, scale = 1)


# get representative GO terms from DEGs
dose_degs <- df_results_merge_dose_significant %>% dplyr::select(-log2FoldChange) %>% pivot_wider(names_from = comparison, values_from = padj)
  # writexl::write_xlsx("merged_by_dose.xlsx")

dose_deg_ids <- dose_degs %>% left_join(grch38, by = "symbol")

dose_deg_2xtaf <- dose_deg_ids %>% dplyr::select(`2xTAFvsCT`, entrez, symbol, ensgene) %>% na.omit() %>% mutate(symbol = replace(symbol, symbol == "KIAA0141", "DELE1"))
dose_deg_2xtdf <- dose_deg_ids %>% dplyr::select(`2xTDFvsCT`, entrez, symbol, ensgene) %>% na.omit() %>% mutate(symbol = replace(symbol, symbol == "KIAA0141", "DELE1"),
                                                                                                                symbol = replace(symbol, symbol == "ATP5J", "ATP5PF"))
dose_deg_1xtaf <- dose_deg_ids %>% dplyr::select(`1xTAFvsCT`, entrez, symbol, ensgene) %>% na.omit() %>% mutate(symbol = replace(symbol, symbol == "KIAA0141", "DELE1"))

library(clusterProfiler)

ggo_2xtaf <- groupGO(gene     = as.character(dose_deg_2xtaf$entrez),
               OrgDb    = org.Hs.eg.db,
               ont      = "BP",
               level    = 5,
               readable = TRUE)
ggo_2xtdf <- groupGO(gene     = as.character(dose_deg_2xtdf$entrez),
               OrgDb    = org.Hs.eg.db,
               ont      = "BP",
               level    = 5,
               readable = TRUE)
ggo_1xtaf <- groupGO(gene     = as.character(dose_deg_1xtaf$entrez),
               OrgDb    = org.Hs.eg.db,
               ont      = "BP",
               level    = 5,
               readable = TRUE)



xx <- ggo_2xtaf@result %>% dplyr::filter(Count == 1) %>% pull(ID)
xy <- ggo_2xtdf@result %>% dplyr::filter(Count == 2) %>% pull(ID)
paste(xx, collapse = ",") %>% write_clip()

tafx2 <- ggo_2xtaf@result %>% dplyr::filter(Count == 1)
tafx2_p <- tafx2 %>% left_join(dose_deg_2xtaf, by = c(geneID = "symbol")) %>% dplyr::rename(padj = `2xTAFvsCT`)

tdfx2 <- ggo_2xtdf@result %>% dplyr::filter(Count == 1)
tdfx2_p <- tdfx2 %>% left_join(dose_deg_2xtdf, by = c(geneID = "symbol")) %>% dplyr::rename(padj = `2xTDFvsCT`)

tafx1 <- ggo_1xtaf@result %>% dplyr::filter(Count == 1)
tafx1_p <- tafx1 %>% left_join(dose_deg_1xtaf, by = c(geneID = "symbol")) %>% dplyr::rename(padj = `1xTAFvsCT`)


# write to clipboard and analyze with web-based revigo (tiny)
ggo_2xtaf@result %>% dplyr::filter(Count == 1) %>% pull(ID) %>% write_clip()
ggo_2xtdf@result %>% dplyr::filter(Count == 1) %>% pull(ID) %>% write_clip()
ggo_1xtaf@result %>% dplyr::filter(Count == 1) %>% pull(ID) %>% write_clip()
