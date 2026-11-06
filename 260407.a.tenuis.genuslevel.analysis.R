library(tidyr)
library(gclus)
library(gradientForest)
library(randomForest)
library(RDAforest)
library(stats)
library(pheatmap)
library(dendextend)
library(impute)
library(stringr)
library(readr)
library(janitor)
library(reshape2)
library(RColorBrewer)
library(xgboost)
library(tidyverse)
library(ggplot2)
library(vegan)
library(mgcv)
library(dplyr)
library(sf)
library(cowplot)
library(ggspatial)
library(ggrepel)
library(tidyr)
library(VennDiagram)
library(UpSetR)
library(MASS)

atenuis.genus.bracken.output.df <- read_tsv("/Users/karimprimov/Dropbox/UTAustin/2026/mcav.cryptic.lineages/mcav.flkeys.metagenomics/atenuis.genus.level.kraken.bracken.reports/260320.atenuis.genus.level.bracken.taxonomic.report.tsv")
a.tenuis.lineage.designation <- read_csv("/Users/karimprimov/Dropbox/UTAustin/2026/mcav.cryptic.lineages/mcav.flkeys.metagenomics/a.tenuis.metadata.csv")

atenuis.filtered.microbes <- read_csv("/Users/karimprimov/Dropbox/UTAustin/2026/mcav.cryptic.lineages/mcav.flkeys.metagenomics/260320.atenuis.strict_coral_microbiome.csv")

atenuis.genus.bracken.output.filtered.df <- atenuis.genus.bracken.output.df[atenuis.genus.bracken.output.df$name %in% atenuis.filtered.microbes$Genus, ]

keep_genera <- c(
  "Endozoicomonas",
  "Alteromonas",
  "Pseudoalteromonas",
  "Shewanella",
  "Vibrio",
  "Acinetobacter",
  "Stenotrophomonas",
  "Roseivirga",
  "Leeuwenhoekiella",
  "Flavobacterium",
  "Sphingomonas",
  "Hyphobacterium",
  "Paracoccus",
  "Synechococcus",
  "Photobacterium",
  "Tritonibacter",
  "Thalassotalea"
)


contaminants <- c(
  # Human / host contamination
  "Homo", "Staphylococcus", "Streptococcus", "Corynebacterium",
  "Cutibacterium", "Propionibacterium", "Micrococcus", "Dermacoccus",
  
  # Common kit contaminants (DNA extraction / PCR reagents)
  "Ralstonia", "Burkholderia", "Cupriavidus", "Sphingomonas",
  "Bradyrhizobium", "Methylobacterium", "Massilia",
  "Brevundimonas", "Caulobacter", "Novosphingobium",
  "Rhizobium", "Agrobacterium", "Ochrobactrum",
  
  # Water / lab environmental contaminants
  "Pseudomonas", "Acinetobacter", "Stenotrophomonas",
  "Delftia", "Comamonas", "Alcaligenes",
  
  # Soil-associated (unlikely in coral tissue microbiome)
  "Bacillus", "Paenibacillus", "Clostridium", "Clostridioides",
  "Enterococcus", "Lactobacillus",
  
  # Enteric / fecal indicators (strong contamination signal)
  "Escherichia", "Shigella", "Klebsiella", "Enterobacter",
  "Salmonella", "Serratia",
  
  # Photosynthetic / organelle contamination
  "Chloroplast", "Mitochondria", "Cyanobium", "Prochlorococcus",
  
  # Fungal / eukaryotic bleed-through (if present in taxonomy)
  "Candida", "Saccharomyces", "Aspergillus",
  
  # Weird freshwater / non-marine genera (context-dependent but usually safe to drop)
  "Aeromonas", "Chromobacterium",
  
  # Human / clinical / oral / pathogen-associated
  "Neisseria", "Leptospira", "Borrelia", "Bartonella",
  "Rickettsia", "Treponema", "Campylobacter",
  "Fusobacterium", "Listeria", "Yersinia",
  "Providencia", "Citrobacter", "Actinomyces",
  "Lactococcus", "Leuconostoc", "Lactiplantibacillus",
  "Weissella", "Bifidobacterium",
  
  # Gut / fecal-associated
  "Bacteroides", "Phocaeicola", "Alistipes",
  "Cloacibacterium", "Cloacibacillus",
  
  # Soil / terrestrial / unlikely marine
  "Nocardiopsis", "Mycobacterium", "Rhodococcus",
  "Microbacterium", "Nocardioides",
  "Lysobacter", "Xanthomonas",
  "Pedobacter", "Flavihumibacter", "Mucilaginibacter",
  
  # Kit / reagent / lab contaminants
  "Paraburkholderia", "Sphingobium", "Sphingopyxis",
  "Achromobacter", "Variovorax",
  
  # Anaerobic / odd environmental contaminants
  "Anaerotignum", "Gottschalkia", "Sporomusa",
  "Acetobacteroides",
  
  # Non-bacterial (remove if present)
  "Alphapapillomavirus",  # virus
  "Methanosarcina",       # archaea
  "Methanococcus",        # archaea
  "Natrinema",             # archaea
  
  # Human / clinical / host-associated
  "Moraxella", "Spiroplasma",
  
  # Plant / terrestrial / non-marine weird
  "Pectobacterium", "Sorangium", "Lentzea",
  "Conexibacter", "Rubrobacter",
  
  # Cyanobacteria / photosynthetic contamination (optional but recommended)
  "Nostoc", "Fischerella", "Microcystis",
  
  # Viral / non-bacterial
  "Tequatrovirus",
  
  # Soil / low-relevance Actinobacteria
  "Brevibacterium", "Zhihengliuella",
  "Planococcus", "Neobacillus",
  
  # Very questionable / unlikely coral-associated
  "Kovacikia", "Sinanaerobacter",
  "Cohnella", "Heyndrickxia",
  
  # Terrestrial / plant-associated
  "Pantoea",
  
  # Soil / low-relevance environmental Actinobacteria
  "Microbispora",
  
  # Freshwater / odd environmental outliers
  "Limnoglobus",
  
  # Very questionable / poorly characterized (likely noise at this abundance)
  "Prescottella"
)





write.csv(atenuis.genus.bracken.output.filtered.df["name"], "atenuis.genus.bracken.output.filtered.df.microbe.names.csv", row.names = FALSE)


atenuis.bacteria.only <- atenuis.genus.bracken.output.filtered.df[!(atenuis.genus.bracken.output.filtered.df$name %in% contaminants), ]

atenuis.microbe.df <- atenuis.bacteria.only %>%
  mutate(
    Sample = substr(as.character(sample), 1, 8),
    Taxa = name,
    new_est_reads = as.numeric(new_est_reads)) %>%
  dplyr::select(Sample, Taxa, new_est_reads)


cleaned.pivot.tenuis.df <- atenuis.microbe.df %>%
  pivot_wider(names_from = "Sample", values_from = "new_est_reads", values_fn = sum)

transposed.tenuis.df <- cleaned.pivot.tenuis.df%>%
  column_to_rownames("Taxa") %>%
  t() %>% 
  as.data.frame() %>%
  mutate(across(everything(), ~replace_na(., 0)))


transposed.tenuis.subset.df <- transposed.tenuis.df[rownames(transposed.tenuis.df) %in% a.tenuis.lineage.designation$Sample, ]

prevalence <- colSums(transposed.tenuis.subset.df > 0)
df.final <- transposed.tenuis.subset.df[, prevalence >= 2]

#metadata:
a.tenuis.metadata <- data.frame(Sample = rownames(transposed.tenuis.subset.df),
                                stringsAsFactors = FALSE) %>%
  left_join(a.tenuis.lineage.designation, by = "Sample") %>%
  column_to_rownames("Sample") %>%
  mutate(readcount = rowSums(transposed.tenuis.subset.df))







# counts_long <- transposed.tenuis.subset.df %>%
#   as.data.frame() %>%
#   rownames_to_column("Sample")
# 
# # 2. Join with your lineage metadata
# combined_data <- counts_long %>%
#   inner_join(a.tenuis.lineage.designation, by = "Sample")
# 
# # 3. Aggregate (sum) counts by Lineage (Cluster) and Site (Reef)
# aggregated_microbiome <- combined_data %>%
#   group_by(Cluster, Reef) %>%
#   summarise(across(where(is.numeric), sum, na.rm = TRUE), 
#             n_samples = n()) %>% # Keeps track of how many samples were pooled
#   ungroup()
# 
# # 4. Check your new group depths
# rowSums(aggregated_microbiome %>% select(where(is.numeric), -n_samples))
# 
# 
# 
# # Create a 'replicate' ID within each Site/Cluster group
# aggregated_df <- combined_data %>%
#   group_by(Cluster, Reef) %>%
#   mutate(rep_group = ceiling(row_number() / 3)) %>% # Groups of 3 samples
#   group_by(Cluster, Reef, rep_group) %>%
#   summarise(across(where(is.numeric), sum, na.rm = TRUE),
#             n_summed = n()) %>%
#   ungroup()
# 
# 
# atenuis_summary <- aggregated_df %>%
#   mutate(lineage_id = paste(Cluster, Reef, rep_group, sep = "_")) %>%
#   column_to_rownames("lineage_id")
# 
# # Metadata for your 'data =' argument in capscale
# atenuis_row_totals <- atenuis_summary %>%
#   dplyr::select(Cluster, Reef, n_summed) %>%
#   mutate(row.total = rowSums(atenuis_summary %>% dplyr::select(where(is.numeric), -Latitude, -Longitude, -n_summed, -rep_group))) %>%
#   dplyr::rename(lineage = Cluster, reef = Reef)
# 
# # 2. Isolate the Taxa Matrix for CLR/Distance
# a.tenuis.microbe.filt.df <- atenuis_summary %>%
#   dplyr::select(-Latitude, -Longitude, -Cluster, -Reef, -n_summed, -rep_group)


#colnames(a.tenuis.metadata): [1] "Reef", "locationRemarks", "Latitude", "Longitude","Cluster" 


a.tenuis.microbe.filt.stand.df <- decostand(transposed.tenuis.subset.df, method = "clr", pseudocount = 0.01)
dist_matrix <- as.matrix(vegdist(a.tenuis.microbe.filt.stand.df, method = "euclidean"))

# Partial db-RDA testing Site first
# a.tenuis.microbe.pca <- capscale(dist_matrix ~ Reef + Cluster + Condition(readcount), data = a.tenuis.metadata)
# plot(a.tenuis.microbe.pca)
# 
# 
# 
# # 1. Extract the scores (CAP axes)
# # Note: 'sites' in vegan refers to your aggregated row IDs
# site_scores <- as.data.frame(scores(a.tenuis.microbe.pca, display = "sites", choices = c(1:2)))
# 
# # 2. Combine with your new metadata
# # Since site_scores and atenuis_row_totals are both indexed by lineage_id, cbind works fine
# plot_data <- cbind(site_scores, a.tenuis.metadata)
# 
# # 3. Calculate % variance for axis labels
# perc <- round(100 * summary(a.tenuis.microbe.pca)$concont$importance[2, 1:2], 1)
# 
# # 4. Plot focused on Site (Color) and Lineage (Shape)
# ggplot(plot_data, aes(x = CAP1, y = CAP2, color = Cluster)) + 
#   geom_point(size = 4, alpha = 0.8) +
#   theme_minimal() +
#   labs(
#     title = "Partial db-RDA (Capscale) of Acropora tenuis Microbiome",
#     subtitle = "Order of effect: Site > Lineage | Conditioned on Read Depth",
#     x = paste0("CAP1 (", perc[1], "%)"), 
#     y = paste0("CAP2 (", perc[2], "%)"),
#     color = "Reef / Site",
#     shape = "Cryptic Lineage"
#   ) +
#   theme(legend.position = "right")
















#only cluster a and b:
# 1. Identify which rows belong to Cluster 1 or 2
# Note: Check if your metadata calls them "Cluster 1" or just "1"
keep_indices <- which(a.tenuis.metadata$Cluster %in% c("Cluster 1A", "Cluster 1B"))

# 2. Subset the metadata
meta_sub <- a.tenuis.metadata[keep_indices, ]

# 3. Subset the distance matrix 
# (vegan requires the matrix to be subsetted to match the metadata)
dist_sub <- as.matrix(dist_matrix)[keep_indices, keep_indices]
dist_sub <- as.dist(dist_sub)

# 4. Re-run the capscale on the subset
pca_sub <- capscale(dist_sub ~ Reef + Cluster + Condition(readcount), data = meta_sub)

# 5. Extract scores and Plot
site_scores_sub <- as.data.frame(scores(pca_sub, display = "sites"))
plot_data_sub <- cbind(site_scores_sub, meta_sub)

ggplot(plot_data_sub, aes(x = CAP1, y = CAP2, color = locationRemarks)) + 
  geom_point(size = 5) + 
  theme_minimal() +
  labs(title = "Microbiome Variation: Clusters 1 & 2 Only")











#overlap in microbes across samples:



#adonis and anova:
#260409:
#adonis 2:

adonis.site.results <- adonis2(dist_sub ~ readcount + Reef + Cluster, 
                               data = meta_sub,
                               by = "margin")   # marginal tests are best for independent effects

# 2. View results
print(adonis.site.results)


adonis.site.lineage.results <- adonis2(dist_sub ~ readcount + Reef + Cluster, 
                                       data = meta_sub,
                                       by = "margin")   # marginal tests are best for independent effects

# 2. View results
print(adonis.site.lineage.results)

####################################################################################
####################################################################################

db_site_model <- dbrda(dist_sub ~ Reef + Condition(readcount), 
                       data = meta_sub)

# 2. Get the PERMANOVA table (identical results to adonis2)
anova(db_site_model, by = "margin")



db_site_lineage_model <- dbrda(dist_sub ~ Reef + Cluster + Condition(readcount), 
                               data = meta_sub)

# 2. Get the PERMANOVA table (identical results to adonis2)
anova(db_site_lineage_model, by = "margin")





varpart(dist_sub, 
        ~ Reef, 
        ~ Cluster, 
        ~ readcount, 
        data = meta_sub)

interaction(Reef, Cluster)

interaction.site.cluster <- lm(dist_sub ~ Reef * Cluster, data = meta_sub)
tidy(interaction.site.cluster) %>% 
  filter(p.value < 0.05)
anova(interaction.site.cluster)

simple_model <- lm(dist.matrix ~ Reef + Cluster, data = meta_sub)
anova(simple_model)


disp_site <- betadisper(as.dist(dist_sub), meta_sub$Reef)
permutest(disp_site) 

# For ClusterID
disp_cluster <- betadisper(as.dist(dist_sub), meta_sub$Cluster)
permutest(disp_cluster)


table(pocil.metadata.ordered$Site, pocil.metadata.ordered$ClusterID)

# 2. Chi-Squared Test (Statistical test of association)
test <- chisq.test(table(meta_sub$Reef, meta_sub$Cluster))
print(test)
















#saving work in RData file
save.image(file = "/Users/karimprimov/Dropbox/UTAustin/2026/mcav.cryptic.lineages/mcav.flkeys.metagenomics/261031.a.tenuis.cryptic.lineage.microbiome.analysis.RData")





#260929:
# ============================================================
# ACROPORA TENUIS — FINAL 1A / 1B DATASET
# ============================================================

keep_indices <- which(
  a.tenuis.metadata$Cluster %in% c("Cluster 1A", "Cluster 1B")
)

ten_meta <- droplevels(a.tenuis.metadata[keep_indices, ])

ten_dist <- as.dist(
  as.matrix(dist_matrix)[keep_indices, keep_indices]
)

# Verify alignment
nrow(ten_meta)
attr(ten_dist, "Size")

identical(
  rownames(ten_meta),
  attr(ten_dist, "Labels")
)

# Sample sizes
table(ten_meta$Cluster)
table(ten_meta$Reef)

ten_meta$Site_unique <- interaction(
  ten_meta$Reef,
  ten_meta$Latitude,
  ten_meta$Longitude,
  drop = TRUE,
  sep = "_"
)

nlevels(ten_meta$Site_unique)

table(ten_meta$Site_unique)

table(
  ten_meta$Site_unique,
  ten_meta$Cluster
)


library(vegan)

set.seed(123)

ten_adonis_site <- adonis2(
  ten_dist ~ readcount + Site_unique,
  data = ten_meta,
  permutations = 9999,
  by = "margin"
)

ten_adonis_site


set.seed(123)

ten_adonis_lineage <- adonis2(
  ten_dist ~ readcount + Cluster,
  data = ten_meta,
  permutations = 9999,
  by = "margin"
)

ten_adonis_lineage


set.seed(123)

ten_adonis_both <- adonis2(
  ten_dist ~ readcount + Site_unique + Cluster,
  data = ten_meta,
  permutations = 9999,
  by = "margin"
)

ten_adonis_both

site_lineage_table <- table(
  ten_meta$Site_unique,
  ten_meta$Cluster
)

site_lineage_table

ten_chisq <- chisq.test(site_lineage_table)
ten_chisq

library(lsr)

ten_cramers <- cramersV(site_lineage_table)
ten_cramers



# Fisher's exact test with Monte Carlo simulation
set.seed(123)
fisher.test(
  site_lineage_table,
  simulate.p.value = TRUE,
  B = 99999
)


set.seed(123)

ten_disp_site <- betadisper(
  ten_dist,
  ten_meta$Site_unique
)

permutest(
  ten_disp_site,
  permutations = 9999
)

set.seed(123)

ten_disp_lineage <- betadisper(
  ten_dist,
  ten_meta$Cluster
)

permutest(
  ten_disp_lineage,
  permutations = 9999
)



ten_dbrda_full <- dbrda(
  ten_dist ~ readcount + Site_unique + Cluster,
  data = ten_meta
)

ten_dbrda_full


set.seed(123)

anova(
  ten_dbrda_full,
  by = "margin",
  permutations = 9999
)


# ============================================================
# A. tenuis — dbRDA FOR VISUALIZATION
# Condition out readcount
# ============================================================

ten_dbrda_plot <- dbrda(
  ten_dist ~ Cluster + Condition(readcount+Site_unique),
  data = ten_meta
)

ten_dbrda_plot

# Extract sample scores
ten_scores <- as.data.frame(
  scores(
    ten_dbrda_plot,
    display = "sites",
    choices = 1:2
  )
)

ten_scores$SampleID_plot <- rownames(ten_scores)
ten_meta$SampleID_plot <- rownames(ten_meta)

ten_plot_data <- merge(
  ten_scores,
  ten_meta,
  by = "SampleID_plot",
  all.x = TRUE
)

# % constrained variation
ten_eig <- eigenvals(
  ten_dbrda_plot,
  model = "constrained"
)

ten_db1 <- 100 * ten_eig[1] / sum(ten_eig)
ten_db2 <- 100 * ten_eig[2] / sum(ten_eig)

ten_db1
ten_db2

# Plot
ggplot(
  ten_plot_data,
  aes(
    x = dbRDA1,
    y = MDS1,
    color = Cluster
  )
) +
  geom_point(
    size = 4,
    alpha = 0.85
  ) +
  theme_minimal() +
  labs(
    title = expression(italic("Acropora tenuis")),
    x = paste0("dbRDA1 (", round(ten_db1, 1), "%)"),
    y = paste0("MDS1 (", round(ten_db2, 1), "%)"),
    color = "Lineage"
  )



set.seed(123)

adonis2(
  ten_dist ~ readcount + Site_unique + Cluster,
  data = ten_meta,
  permutations = permute::how(
    nperm = 9999,
    blocks = factor(ten_meta$Site_unique)
  ),
  by = "margin"
)

ten_clr <- as.matrix(a.tenuis.microbe.filt.stand.df)
ten_clr <- ten_clr[rownames(ten_meta), , drop = FALSE]

dim(ten_clr)
anyNA(ten_clr)


library(adegenet)

set.seed(123)

ten_cv <- xvalDapc(
  ten_clr,
  grp = ten_meta$Cluster,
  n.pca.max = 20,
  n.rep = 30,
  xval.plot = TRUE
)

ten_cv$`Mean Successful Assignment by Number of PCs of PCA`
ten_cv$`Number of PCs Achieving Highest Mean Success`
ten_cv$`Number of PCs Achieving Lowest MSE`



ten_clr <- as.matrix(a.tenuis.microbe.filt.stand.df)
ten_clr <- ten_clr[rownames(ten_meta), , drop = FALSE]
dfa.lineage <- MASS::lda(ten_clr, grouping = ten_meta$Cluster)
plot(dfa.lineage)

ten_dapc <- dapc(
  ten_clr,
  grp = ten_meta$Cluster,
  n.pca = 10,
  n.da = 3
)
ten_dapc
scatter(ten_dapc)
ten_dfa <- MASS::lda(
  ten_clr[, -ncol(ten_clr)],
  grouping = ten_meta$Cluster
)

ten_dfa
plot(ten_dfa)







a.tenuis.microbe.filt.stand.df <- decostand(transposed.tenuis.subset.df, method = "clr", pseudocount = 0.01)
dist_matrix <- as.matrix(vegdist(a.tenuis.microbe.filt.stand.df, method = "euclidean"))
dim(dist_matrix)
qr(as.matrix(dist_matrix))$rank

ten_clr <- as.matrix(a.tenuis.microbe.filt.stand.df)
ten_clr <- ten_clr[rownames(ten_meta), , drop = FALSE]

ten_dist <- dist(ten_clr)

dim(ten_clr)

set.seed(123)

ten_adonis_both <- vegan::adonis2(
  ten_dist ~ readcount + Site_unique + Cluster,
  data = ten_meta,
  permutations = 9999,
  by = "margin"
)

ten_adonis_both

ten_dbrda_plot <- vegan::dbrda(
  ten_dist ~ Cluster + Condition(readcount + Site_unique),
  data = ten_meta
)

ten_dbrda_plot


ten_disp_lineage <- vegan::betadisper(ten_dist, ten_meta$Cluster)
ten_disp_site <- vegan::betadisper(ten_dist, factor(ten_meta$Site_unique))

permutest(ten_disp_lineage, permutations = 9999)
permutest(ten_disp_site, permutations = 9999)

set.seed(123)

vegan::adonis2(
  ten_dist ~ readcount + Site_unique + Cluster,
  data = ten_meta,
  permutations = permute::how(
    nperm = 9999,
    blocks = factor(ten_meta$Site_unique)
  ),
  by = "margin"
)


ten_disp_lineage <- vegan::betadisper(ten_dist, ten_meta$Cluster)
ten_disp_site <- vegan::betadisper(ten_dist, factor(ten_meta$Site_unique))

permutest(ten_disp_lineage, permutations = 9999)
permutest(ten_disp_site, permutations = 9999)
#

ten_dfa <- MASS::lda(
  ten_clr[, -ncol(ten_clr)],
  grouping = ten_meta$Cluster
)

ten_dfa
plot(ten_dfa)


mean(predict(ten_dfa)$class == ten_meta$Cluster)

table(
  Observed = ten_meta$Cluster,
  Predicted = predict(ten_dfa)$class
)

set.seed(123)

ten_cv <- adegenet::xvalDapc(
  ten_clr,
  grp = ten_meta$Cluster,
  n.pca.max = 20,
  n.rep = 30,
  xval.plot = TRUE
)

ten_cv$`Number of PCs Achieving Highest Mean Success`
ten_cv$`Mean Successful Assignment by Number of PCs of PCA`

ten_dapc <- adegenet::dapc(
  ten_clr,
  grp = ten_meta$Cluster,
  n.pca = 10,
  n.da = 1
)

scatter(ten_dapc)


ten_dapc <- adegenet::dapc(
  ten_clr,
  grp = ten_meta$Cluster,
  n.pca = 16,
  n.da = 1
)

scatter(ten_dapc)


plot_w <- 6
plot_h <- 5
plot_dpi <- 300

ten_plot <- ggplot(
  ten_plot_data,
  aes(x = dbRDA1, y = MDS1, color = Cluster)
) +
  geom_point(size = 2.5, alpha = 0.8) +
  theme_minimal() +
  labs(
    title = expression(italic("Acropora tenuis")),
    x = "dbRDA1 (100% constrained)",
    y = "MDS1 (residual)",
    color = "Lineage"
  )

ggsave("tenuis_dbRDA.pdf", ten_plot,
       width = 6, height = 5)

ten_dapc <- adegenet::dapc(
  ten_clr,
  grp = ten_meta$Cluster,
  n.pca = 16,
  n.da = 1
)

pdf("tenuis_DAPC.pdf", width = 6, height = 5)
scatter(ten_dapc)
dev.off()

ten_pca <- prcomp(ten_clr)

ten_pca_data <- as.data.frame(ten_pca$x[, 1:2])
ten_pca_data$Cluster <- ten_meta$Cluster

ten_pca_plot <- ggplot(ten_pca_data,
                       aes(PC1, PC2, color = Cluster)) +
  geom_point(size = 2.5, alpha = 0.8) +
  theme_minimal() +
  labs(
    title = expression(italic("Acropora tenuis")),
    color = "Lineage"
  )

ggsave("tenuis_PCA.pdf", ten_pca_plot,
       width = 6, height = 5)



ten_scores <- as.data.frame(
  vegan::scores(ten_dbrda_plot, display = "sites", choices = 1:2)
)

ten_scores$SampleID_plot <- rownames(ten_scores)
ten_meta$SampleID_plot <- rownames(ten_meta)

ten_plot_data <- merge(
  ten_scores, ten_meta,
  by = "SampleID_plot"
)

names(ten_scores)

ten_plot <- ggplot(
  ten_plot_data,
  aes(x = dbRDA1, y = MDS1, color = Cluster)
) +
  geom_point(size = 2.5, alpha = 0.8) +
  theme_minimal() +
  labs(
    title = expression(italic("Acropora tenuis")),
    x = "dbRDA1 (100% constrained)",
    y = "MDS1 (residual)",
    color = "Lineage"
  )

ten_plot

ggsave("tenuis_dbRDA.pdf", plot = ten_plot,
       width = 6, height = 5, units = "in")

ten_dapc <- adegenet::dapc(
  ten_clr,
  grp = ten_meta$Cluster,
  n.pca = 16,
  n.da = 1
)

pdf("tenuis_DAPC.pdf", width = 6, height = 5)
adegenet::scatter(ten_dapc)
dev.off()




ten_scores <- as.data.frame(
  scores(ten_dbrda_plot, display = "sites", choices = 1:2)
)

ten_scores$SampleID_plot <- rownames(ten_scores)

ten_plot_data <- merge(
  ten_scores, ten_meta, by = "SampleID_plot"
)
#









#





#





















library(ggplot2)

# ==========================================
# ACROPORA TENUIS — PUBLICATION dbRDA
# ==========================================

p_tenuis <- ggplot(
  ten_plot_data,
  aes(
    x = dbRDA1,
    y = dbRDA2,
    color = Cluster
  )
) +
  geom_point(
    size = 8,
    alpha = 0.85
  ) +
  
  # Species name
  labs(
    title = expression(italic("Acropora tenuis")),
    x = paste0(
      "dbRDA1 (",
      round(ten_db1, 1),
      "%)"
    ),
    y = paste0(
      "dbRDA2 (",
      round(ten_db2, 1),
      "%)"
    ),
    color = "Lineage"
  ) +
  
  # Clean publication theme
  theme_classic(base_size = 12) +
  
  # Force the plotting PANEL to be square
  theme(
    aspect.ratio = 1,
    
    plot.title = element_text(
      size = 13,
      hjust = 0
    ),
    
    axis.title = element_text(
      size = 11
    ),
    
    axis.text = element_text(
      size = 10
    ),
    
    legend.title = element_text(
      size = 10
    ),
    
    legend.text = element_text(
      size = 9
    ),
    
    legend.position = "right"
  )

# View it
p_tenuis

ggsave(
  filename = "A_tenuis_dbRDA.png",
  plot = p_tenuis,
  width = 5,
  height = 5,
  units = "in",
  dpi = 600,
  bg = "white"
)

ggsave(
  filename = "A_tenuis_dbRDA.pdf",
  plot = p_tenuis,
  width = 5,
  height = 5,
  units = "in"
)




#261104:
library(ggplot2)

ten_pca_final <- ten_pca_plot +
  theme_classic(base_size = 12) +
  theme(
    aspect.ratio = 1,
    panel.grid = element_blank(),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    plot.title = element_text(size = 13, face = "italic"),
    legend.position = "none",
    plot.margin = margin(8, 8, 8, 8)
  )

print(ten_pca_final)

ggsave(
  "ten_PCA_final.pdf",
  plot = ten_pca_final,
  width = 5,
  height = 5,
  units = "in"
)




library(adegenet)

set.seed(123)

ten_cv_final <- xvalDapc(
  ten_clr,
  grp = ten_meta$Cluster,
  n.pca.max = 20,
  training.set = 0.7,
  result = "groupMean",
  center = TRUE,
  scale = FALSE,
  n.rep = 100,
  xval.plot = TRUE
)
ten_groups <- droplevels(factor(ten_meta$Cluster))
ten_cv_final$`Mean Successful Assignment by Number of PCs of PCA`

ten_cv_final$`Number of PCs Achieving Highest Mean Success`



ten_best_npc <- as.integer(
  ten_cv_final$`Number of PCs Achieving Highest Mean Success`
)

ten_dapc_final <- dapc(
  ten_clr,
  grp = ten_groups,
  n.pca = ten_best_npc,
  n.da = 1
)

scatter(ten_dapc_final, scree.da = FALSE)



# Export final A. tenuis DAPC plot

pdf(
  "ten_DAPC_final.pdf",
  width = 5,
  height = 5
)

scatter(
  ten_dapc_final,
  scree.da = FALSE
)

dev.off()




#re-doing the db-rda
library(ggplot2)

# Final standardized A. tenuis dbRDA

ten_dbrda_final <- ggplot(
  ten_plot_data,
  aes(
    x = dbRDA1,
    y = MDS1,
    color = Cluster
  )
) +
  geom_point(size = 5, alpha = 0.85) +
  labs(
    title = "Acropora tenuis",
    x = "dbRDA1 (0.34%)",
    y = "MDS1 (residual)"
  ) +
  theme_classic(base_size = 12) +
  theme(
    aspect.ratio = 1,
    panel.grid = element_blank(),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    plot.title = element_text(size = 13, face = "italic"),
    legend.position = "none",
    plot.margin = margin(8, 8, 8, 8)
  )

# View
print(ten_dbrda_final)

# Export
ggsave(
  "ten_dbRDA_final.pdf",
  plot = ten_dbrda_final,
  width = 5,
  height = 5,
  units = "in"
)








dim(ten_meta)
names(ten_meta)
head(ten_meta)

# Number of samples per site
table(ten_meta$Site_unique)

# Lineage composition by site
table(ten_meta$Site_unique, ten_meta$Cluster)


library(dplyr)
library(tidyr)
library(ggplot2)
library(sf)
library(scatterpie)
library(rnaturalearth)

# Summarize lineage counts at each reef
ten_site_pies <- ten_meta %>%
  as.data.frame() %>%
  group_by(Site_unique, Reef, Latitude, Longitude, Cluster) %>%
  summarise(n = n(), .groups = "drop") %>%
  pivot_wider(
    names_from = Cluster,
    values_from = n,
    values_fill = 0
  )

# Rename lineage columns for convenient plotting
ten_site_pies <- ten_site_pies %>%
  rename(
    Cluster1A = `Cluster 1A`,
    Cluster1B = `Cluster 1B`
  )

# Confirm 23 sites
nrow(ten_site_pies)

australia <- ne_countries(
  scale = "medium",
  country = "Australia",
  returnclass = "sf"
)




ten_map <- ggplot() +
  
  # Australian coastline
  geom_sf(
    data = australia,
    fill = "grey90",
    color = "grey40",
    linewidth = 0.3
  ) +
  
  # Lineage composition at each reef
  geom_scatterpie(
    data = ten_site_pies,
    aes(
      x = Longitude,
      y = Latitude,
      r = 0.4
    ),
    cols = c("Cluster1A", "Cluster1B"),
    color = "black",
    linewidth = 0.2
  ) +
  
  scale_fill_manual(
    values = c(
      Cluster1A = "#F8766D",
      Cluster1B = "#00BFC4"
    ),
    labels = c(
      Cluster1A = "Cluster 1A",
      Cluster1B = "Cluster 1B"
    )
  ) +
  
  coord_sf(
    xlim = c(141, 155),
    ylim = c(-25, -10),
    expand = FALSE
  ) +
  
  labs(
    title = "Acropora tenuis",
    x = "Longitude",
    y = "Latitude",
    fill = "Lineage"
  ) +
  
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(
      face = "italic",
      size = 13
    ),
    panel.grid = element_blank(),
    legend.position = "right"
  )

print(ten_map)









library(ggplot2)
library(ggrepel)
library(scatterpie)
library(dplyr)

# Keep original reef coordinates
ten_site_pies <- ten_site_pies %>%
  mutate(
    actual_lon = Longitude,
    actual_lat = Latitude
  )

# Use a repulsion layout to calculate new pie positions
set.seed(123)

# Start with pie coordinates at the real locations
ten_site_pies <- ten_site_pies %>%
  mutate(
    pie_lon = Longitude,
    pie_lat = Latitude
  )

# Manually offset the tightly clustered northern sites
# while leaving the other reef locations unchanged
ten_site_pies <- ten_site_pies %>%
  arrange(desc(Latitude)) %>%
  mutate(
    pie_lon = if_else(
      Latitude > -15,
      Longitude + rep(c(-0.5, 0.5, -0.7, 0.7), length.out = n()),
      Longitude
    )
  )



ten_map <- ggplot() +
  
  geom_sf(
    data = australia,
    fill = "grey90",
    color = "grey40",
    linewidth = 0.3
  ) +
  
  # Thin lines to actual reef locations
  geom_segment(
    data = ten_site_pies,
    aes(
      x = actual_lon,
      y = actual_lat,
      xend = pie_lon,
      yend = pie_lat
    ),
    color = "grey40",
    linewidth = 0.3
  ) +
  
  # Actual locations
  geom_point(
    data = ten_site_pies,
    aes(x = actual_lon, y = actual_lat),
    color = "black",
    size = 0.8
  ) +
  
  # Pies at adjusted positions
  geom_scatterpie(
    data = ten_site_pies,
    aes(
      x = pie_lon,
      y = actual_lat,
      r = 0.4
    ),
    cols = c("Cluster1A", "Cluster1B"),
    color = "black",
    linewidth = 0.2
  ) +
  
  scale_fill_manual(
    values = c(
      Cluster1A = "#F8766D",
      Cluster1B = "#00BFC4"
    ),
    labels = c(
      Cluster1A = "Cluster 1A",
      Cluster1B = "Cluster 1B"
    )
  ) +
  
  coord_sf(
    xlim = c(141, 155),
    ylim = c(-25, -10),
    expand = FALSE
  ) +
  
  labs(
    title = "Acropora tenuis",
    x = "Longitude",
    y = "Latitude",
    fill = "Lineage"
  ) +
  
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "italic", size = 13),
    panel.grid = element_blank()
  )

print(ten_map)







ggsave(
  "/Users/karimprimov/Dropbox/UTAustin/2026/mcav.cryptic.lineages/261104.plots/tenuis_sampling_map_final.pdf",
  plot = ten_map,
  width = 6,
  height = 7,
  units = "in"
)




















































































