library(DESeq2)
library(readxl)
library(tidyverse)
library(FactoMineR)
library(factoextra)
library(RColorBrewer)
library(glmpca)

### DDS object with both countries for diagnostic plots

metaData_all <- unblinding  %>% dplyr::select(country, treatment, CODE) %>% 
  unite(drug_country, c("country", "treatment"))

metaData_all <- metaData_all %>% column_to_rownames("CODE") %>% dplyr::mutate(drug_country = as.factor(drug_country))
countData_all <- countData[, intersect(rownames(metaData_all), names(countData))] # 

metaData_all <- metaData_all %>% rownames_to_column("CODE") %>%  dplyr::filter(CODE %in% intersect(rownames(metaData_all), names(countData))) %>% 
  column_to_rownames("CODE")

countData_all <- data.matrix(countData_all)

dds_all <- DESeqDataSetFromMatrix(countData = countData_all,
                              colData = metaData_all,
                              design =  ~ drug_country)

# keep only genes for which 90% of samples have at least 1 count
nrow(dds_all)
keep <- rowSums(counts(dds_all) >= 2) >= dim(metaData_all)[1]*0.9
filtered <- !keep
# keep <- rowSums(fpm(dds_sa, robust = F) >= 1) >= 7
# table(keep)
dds_all <- dds_all[keep,]
nrow(dds_all)

#### PCA
vsd_all <- vst(dds_all, blind = FALSE)
pca_data_all <- plotPCA(vsd_all, intgroup = c("drug_country"), returnData = T)

all_vst_pca <- ggplot(pca_data_all %>% dplyr::filter(grepl("CT", drug_country)), aes(x = PC1, y = PC2)) +
  geom_point(aes(color = drug_country), size =2) + 
  geom_point(shape = 1, size = 2, color = "black") + 
  ggtitle("") +
  hrbrthemes::theme_ipsum() +
  xlim(-50, 50) +
  ylim(-50, 50) +
  scale_color_brewer(palette = "Set3")

# outlier identification
pca_data_all %>% mutate(zscore = (PC1 - mean(PC1))/sd(PC1),
                  outlier = if_else(abs(zscore)>3, "yes", "no"))
