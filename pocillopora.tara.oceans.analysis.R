library(tidyr)
library(gclus)
library(gradientForest)
library(randomForest)
library(RDAforest)
library(stats)
library(pheatmap)
library(dendextend)
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
library(pcutils)
library(pctax)
library(adegenet)
library(lsr)

pocillopora.tara.bracken.output.df <- read_tsv("/Users/karimprimov/Dropbox/UTAustin/2026/mcav.cryptic.lineages/mcav.flkeys.metagenomics/tara.oceans.bracken.reports/260422.pocillopora.effusa.tara.genus.level.bracken.taxonomic.report.tsv")
pocil.microbe.taxa <- unique(pocillopora.tara.bracken.output.df$name)
write.csv(pocil.microbe.taxa, "260427.pocillopora.microbe.taxa.screencheck.csv", row.names = FALSE)

keep_taxa <- c(
  
  # =========================
  # CORE CORAL SYMBIONTS (high-confidence reef core)
  # =========================
  "Endozoicomonas",
  "Ruegeria",
  "Phaeobacter",
  "Sulfitobacter",
  "Roseobacter",
  "Roseovarius",
  "Salipiger",
  "Labrenzia",
  "Paracoccus",
  "Devosia",
  "Jannaschia",
  "Gemmobacter",
  "Roseibium",
  "Leisingera",
  
  # =========================
  # ROSEOBACTER / RHODOBACTERACEAE EXPANSION (important pocillopora signal)
  # =========================
  "Loktanella",
  "Shimia",
  "Octadecabacter",
  "Celeribacter",
  "Amaricoccus",
  "Thalassococcus",
  "Dinoroseobacter",
  "Planktomarina",
  "Roseivivax",
  "Roseicitreum",
  "Phaeobacter",
  
  # =========================
  # ALPHA / MARINE SURFACE & BIOFILM (reef mucus layer)
  # =========================
  "Hyphomicrobium",
  "Caulobacter",
  "Brevundimonas",
  "Hyphomonas",
  "Methylobacterium",
  "Methylorubrum",
  "Rhodopseudomonas",
  "Thalassospira",
  "Octadecabacter",
  "Amaricoccus",
  "Sphingomonas",
  "Sphingopyxis",
  "Sphingobium",
  "Novosphingobium",
  "Erythrobacter",
  "Porphyrobacter",
  "Altererythrobacter",
  
  # =========================
  # CYANOBACTERIA (reef-realistic subset only)
  # =========================
  "Synechococcus",
  "Prochlorococcus",
  "Trichodesmium",
  "Nostoc",
  "Chroococcidiopsis",
  "Leptolyngbya",
  "Oscillatoria",
  "Phormidium",
  "Geitlerinema",
  
  # =========================
  # FLAVOBACTERIA / MUCUS DEGRADERS (core reef processing layer)
  # =========================
  "Flavobacterium",
  "Polaribacter",
  "Winogradskyella",
  "Formosa",
  "Maribacter",
  "Algibacter",
  "Aquimarina",
  "Cellulophaga",
  "Zobellia",
  "Nonlabens",
  "Dokdonia",
  "Kordia",
  "Tenacibaculum",
  "Flammeovirga",
  "Chryseobacterium",
  "Pedobacter",
  "Sphingobacterium",
  "Mucilaginibacter",
  "Hymenobacter",
  "Rufibacter",
  
  # =========================
  # GAMMAPROTEOBACTERIA (reef + seawater interface core)
  # =========================
  "Alteromonas",
  "Pseudoalteromonas",
  "Marinobacter",
  "Marinomonas",
  "Shewanella",
  "Colwellia",
  "Psychrobacter",
  "Halomonas",
  "Salinivibrio",
  "Aliivibrio",
  "Vibrio",
  "Photobacterium",
  "Thalassotalea",
  "Neptunomonas",
  "Idiomarina",
  "Catenovulum",
  "Ferrimonas",
  "Alcanivorax",
  "Marinobacterium",
  "Thalassolituus",
  "Agarivorans",
  "Paraglaciecola",
  "Thalassomonas",
  "Psychromonas",
  "Spongiibacter",
  
  # =========================
  # SULFUR / REDOX CYCLING (reef biogeochemistry core)
  # =========================
  "Thiomicrorhabdus",
  "Sulfurimonas",
  "Beggiatoa",
  "Desulfovibrio",
  "Desulfobacter",
  "Desulfuromonas",
  "Geobacter",
  "Desulfobulbus",
  "Syntrophus",
  "Thiothrix",
  "Halothiobacillus",
  "Thioalkalivibrio",
  "Sedimenticola",
  "Acidithiobacillus",
  
  # =========================
  # NITROGEN / METHYLOTROPHS (reef surface chemistry)
  # =========================
  "Nitrosomonas",
  "Nitrosococcus",
  "Nitrospira",
  "Methylophaga",
  "Methylotenera",
  "Methyloversatilis",
  "Methylococcus",
  "Methylomonas",
  "Hyphomicrobium",
  "Methylocapsa",
  "Methylocella",
  
  # =========================
  # ACTINOBACTERIA (reef + sediment + biofilm contributors)
  # =========================
  "Rhodococcus",
  "Nocardioides",
  "Microbacterium",
  "Corynebacterium",
  "Streptomyces",
  "Propionibacterium",
  
  # =========================
  # CONDITIONAL / STRESS-ASSOCIATED (keep but interpret carefully)
  # =========================
  "Pseudomonas",
  "Acinetobacter",
  "Burkholderia",
  "Cupriavidus",
  "Comamonas",
  "Hydrogenophaga",
  "Variovorax",
  "Delftia",
  "Stenotrophomonas",
  "Ralstonia"
)
pocil.metadata <- read_csv("/Users/karimprimov/Dropbox/UTAustin/2026/sp26/cryptic.lineage.microbiome/tara.oceans/pocillopora.effusa.metadata.csv")
pocil.metadata <- pocil.metadata %>%
  column_to_rownames(var = "Sample Name")


keep.bacteria.df <- pocillopora.tara.bracken.output.df[pocillopora.tara.bracken.output.df$name %in% keep_taxa,]

keep.bacteria.wide.join.df <- keep.bacteria.df %>%
  dplyr::select(Sample = sample, Taxa = name, readcount = new_est_reads) %>%
  pivot_wider(names_from = "Sample", values_from = "readcount", values_fn = sum) %>%
  column_to_rownames(var = "Taxa") %>%
  mutate(across(everything(), ~replace(., is.na(.), 1)))

t.keep.bacteria.wide.join.df <- t(keep.bacteria.wide.join.df)
rownames(t.keep.bacteria.wide.join.df) <- gsub(".bracken.genus.level.output.report", "", rownames(t.keep.bacteria.wide.join.df))
rownames(pocil.metadata) <- gsub(".bracken.genus.level.output.report", "", rownames(pocil.metadata))
pocil.metadata.ordered <- pocil.metadata[match(rownames(t.keep.bacteria.wide.join.df), rownames(pocil.metadata)), ]

pocil.metadata.ordered$readcount = rowSums(t.keep.bacteria.wide.join.df)




stand.microbes <- decostand(t.keep.bacteria.wide.join.df, method = "rclr")

dist.matrix <- as.matrix(dist(stand.microbes, method = "euclidean"))

capscale.bacteria <- capscale(dist.matrix ~ Site + ClusterID + Condition(readcount), data = pocil.metadata.ordered)
plot(capscale.bacteria)

microbe.scores <- scores(capscale.bacteria, display = "sites", choices=c(1:2)) 

merged.scores.metadata <- merge(microbe.scores, pocil.metadata.ordered, by = "row.names", all.x = TRUE)

ggplot(merged.scores.metadata, aes(x = MDS1, y = MDS2, color =Latitude, label = rownames(merged.scores.metadata))) +
  geom_point(size = 12) +
  theme_minimal() + 
  coord_fixed()

ggplot(merged.scores.metadata, aes(x = CAP1, y = CAP2, color = Latitude, label = rownames(merged.scores.metadata))) +
  geom_point(size = 12) +
  theme_minimal() + 
  coord_fixed()


#

adonis.site.results <- adonis2(dist.matrix ~ readcount + Site + ClusterID, 
                               data = pocil.metadata.ordered,
                               by = "margin")   # marginal tests are best for independent effects


varpart(dist.matrix ~ Site ~ ClusterID ~readcount, data = pocil.metadata.ordered)







#


pocillopora.db_readcount_model <- dbrda(pocillopora.bacteria.dist.matrix ~ ClusterID + Site + Condition(readcount), 
                                    data = pocil.metadata.ordered)

# 2. Get the PERMANOVA table (identical results to adonis2)
anova(pocillopora.db_readcount_model, by = "margin")

#variance partitioning
varpart(pocillopora.bacteria.dist.matrix, 
        ~ Site, 
        ~ ClusterID, 
        ~ readcount, 
        data = pocil.metadata.ordered)

interaction(Site, ClusterID)

interaction.site.cluster <- lm(dist.matrix ~ Site * ClusterID, data = pocil.metadata.ordered)
tidy(interaction.site.cluster) %>% 
  filter(p.value < 0.05)
anova(interaction.site.cluster)

simple_model <- lm(dist.matrix ~ Site + ClusterID, data = pocil.metadata.ordered)
anova(simple_model)


disp_site <- betadisper(as.dist(dist.matrix), pocil.metadata.ordered$Site)
permutest(disp_site) 

# For ClusterID
disp_cluster <- betadisper(as.dist(dist.matrix), pocil.metadata.ordered$ClusterID)
permutest(disp_cluster)


table(pocil.metadata.ordered$Site, pocil.metadata.ordered$ClusterID)

# 2. Chi-Squared Test (Statistical test of association)
test <- chisq.test(table(pocil.metadata.ordered$Site, pocil.metadata.ordered$ClusterID))
print(test)



pocillopora.microbiome.centroid.distances <- as.matrix(dist(disp_cluster$centroids))




#



















#Saving RData for the day!
#######################################################################################################################################
save.image("/Users/karimprimov/Dropbox/UTAustin/2026/mcav.cryptic.lineages/mcav.flkeys.metagenomics/261031.pocillopora.comparison.genuslevel.analysis.RData")















#260929:

rownames(dist.matrix)
colnames(dist.matrix)
rownames(pocil.metadata.ordered)
head(pocil.metadata.ordered)

adonis2(dist.matrix ~ readcount + Site, data = pocil.metadata.ordered, permutations = 9999, by = "margin")
adonis2(dist.matrix ~ readcount + ClusterID, data = pocil.metadata.ordered, permutations = 9999, by = "margin")
adonis2(dist.matrix ~ readcount + Site + ClusterID, data = pocil.metadata.ordered, permutations = 9999, by = "margin")

# See exactly what's missing
pocil.metadata.ordered[
  is.na(pocil.metadata.ordered$ClusterID) |
    pocil.metadata.ordered$ClusterID == "#N/A",
  c("Species", "Site", "Latitude", "Longitude", "ClusterID")
]

keep <- !is.na(pocil.metadata.ordered$ClusterID) &
  pocil.metadata.ordered$ClusterID != "#N/A"

poc_meta <- droplevels(pocil.metadata.ordered[keep, ])

poc_dist <- as.dist(
  as.matrix(dist.matrix)[keep, keep]
)


nrow(poc_meta)
attr(poc_dist, "Size")

identical(
  rownames(poc_meta),
  attr(poc_dist, "Labels")
)

table(poc_meta$ClusterID)

poc_meta$Site_unique <- interaction(
  poc_meta$Site,
  poc_meta$Latitude,
  poc_meta$Longitude,
  drop = TRUE,
  sep = "_"
)

nlevels(poc_meta$Site_unique)

table(poc_meta$Site_unique)

table(
  poc_meta$Site_unique,
  poc_meta$ClusterID
)


library(vegan)

set.seed(123)

# ---------------------------------
# 1. SITE-ONLY PERMANOVA
# ---------------------------------

adonis_site <- adonis2(
  poc_dist ~ readcount + Site_unique,
  data = poc_meta,
  permutations = 9999,
  by = "margin"
)

adonis_site


# ---------------------------------
# 2. LINEAGE-ONLY PERMANOVA
# ---------------------------------

set.seed(123)

adonis_lineage <- adonis2(
  poc_dist ~ readcount + ClusterID,
  data = poc_meta,
  permutations = 9999,
  by = "margin"
)

adonis_lineage


# ---------------------------------
# 3. JOINT PERMANOVA
# ---------------------------------

set.seed(123)

adonis_both <- adonis2(
  poc_dist ~ readcount + Site_unique + ClusterID,
  data = poc_meta,
  permutations = 9999,
  by = "margin"
)

adonis_both


# ---------------------------------
# 4. SITE × LINEAGE ASSOCIATION
# ---------------------------------

site_lineage_table <- table(
  poc_meta$Site_unique,
  poc_meta$ClusterID
)

chi_test <- chisq.test(site_lineage_table)
chi_test

library(lsr)
cramersV(site_lineage_table)


# ---------------------------------
# 5. PERMDISP — SITE
# ---------------------------------

set.seed(123)

disp_site <- betadisper(
  poc_dist,
  poc_meta$Site_unique
)

permutest(
  disp_site,
  permutations = 9999
)


# ---------------------------------
# 6. PERMDISP — LINEAGE
# ---------------------------------

set.seed(123)

disp_lineage <- betadisper(
  poc_dist,
  poc_meta$ClusterID
)

permutest(
  disp_lineage,
  permutations = 9999
)


# ---------------------------------
# 7. dbRDA
# ---------------------------------

poc_dbrda_full <- dbrda(
  poc_dist ~ readcount + Site_unique + ClusterID,
  data = poc_meta
)

poc_dbrda_full


# marginal tests
set.seed(123)

anova(
  poc_dbrda_full,
  by = "margin",
  permutations = 9999
)




library(vegan)
library(ggplot2)

# Full dbRDA model
poc_dbrda_plot <- dbrda(
  poc_dist ~ Site_unique + ClusterID + Condition(readcount),
  data = poc_meta
)

# Extract sample scores
poc_scores <- as.data.frame(
  scores(
    poc_dbrda_plot,
    display = "sites",
    choices = 1:2
  )
)

# Add sample IDs
poc_scores$SampleID_plot <- rownames(poc_scores)
poc_meta$SampleID_plot <- rownames(poc_meta)

# Merge scores with metadata
poc_plot_data <- merge(
  poc_scores,
  poc_meta,
  by = "SampleID_plot",
  all.x = TRUE
)

# Get % of constrained variation on first 2 axes
eig <- eigenvals(poc_dbrda_plot, model = "constrained")

db1 <- 100 * eig[1] / sum(eig)
db2 <- 100 * eig[2] / sum(eig)

# Plot
ggplot(
  poc_plot_data,
  aes(
    x = dbRDA1,
    y = dbRDA2,
    color = ClusterID
  )
) +
  geom_point(size = 12, alpha = 0.85) +
  theme_minimal() +
  coord_fixed() +
  labs(
    title = "Pocillopora effusa",
    x = paste0("dbRDA1 (", round(db1, 1), "%)"),
    y = paste0("dbRDA2 (", round(db2, 1), "%)"),
    color = "Lineage"
  )










library(vegan)

poc_meta$Lineage_broad <- factor(
  sub("I[0-9]+$", "", as.character(poc_meta$ClusterID))
)

table(poc_meta$Lineage_broad)

table(poc_meta$Site_unique, poc_meta$Lineage_broad)

stopifnot(
  identical(rownames(poc_meta), attr(poc_dist, "Labels"))
)

poc_adonis <- adonis2(
  poc_dist ~ readcount + Site_unique + Lineage_broad,
  data = poc_meta,
  by = "margin",
  permutations = 9999
)

poc_adonis





poc_disp <- betadisper(
  poc_dist,
  group = poc_meta$Lineage_broad
)

permutest(poc_disp, permutations = 9999)

# Within-site restricted permutations
poc_perm <- how(
  nperm = 9999,
  blocks = poc_meta$Site_unique
)

adonis2(
  poc_dist ~ readcount + Site_unique + Lineage_broad,
  data = poc_meta,
  by = "margin",
  permutations = poc_perm
)





#261009 pocillopora analysis:
library(vegan)
library(permute)

poc_remove <- c(
  "I01S01C006POC",
  "I03S01C011POC",
  "I03S01C016POC",
  "I03S01C018POC",
  "I03S01C019POC",
  "I03S01C020POC",
  "I05S02C010POC"
)

# Confirm sample alignment
stopifnot(
  identical(rownames(poc_meta), attr(poc_dist, "Labels"))
)

# Remove duplicate ramets
keep <- !poc_meta$`Cluster Pocil` %in% poc_remove

poc_meta_pruned <- droplevels(poc_meta[keep, , drop = FALSE])

poc_dist_pruned <- as.dist(
  as.matrix(poc_dist)[keep, keep]
)

# Verify sample count and group representation
nrow(poc_meta_pruned)
table(poc_meta_pruned$Lineage_broad)

# Marginal PERMANOVA
adonis2(
  poc_dist_pruned ~ readcount + Site_unique + Lineage_broad,
  data = poc_meta_pruned,
  by = "margin",
  permutations = 9999
)

adonis2(
  poc_dist_pruned ~ readcount + Site_unique + Lineage_broad,
  data = poc_meta_pruned,
  by = "margin",
  permutations = 9999
)

# Site-restricted permutation test
poc_perm_pruned <- how(
  nperm = 9999,
  blocks = poc_meta_pruned$Site_unique
)

adonis2(
  poc_dist_pruned ~ readcount + Site_unique + Lineage_broad,
  data = poc_meta_pruned,
  by = "margin",
  permutations = poc_perm_pruned
)

# Dispersion
poc_disp_pruned <- betadisper(
  poc_dist_pruned,
  group = poc_meta_pruned$Lineage_broad
)

permutest(poc_disp_pruned, permutations = 9999)





poc_meta_final <- droplevels(poc_meta[keep, , drop = FALSE])

poc_dist_final <- as.dist(
  as.matrix(poc_dist)[keep, keep]
)

poc_meta_final$Lineage_broad <- factor(
  sub("I[0-9]+$", "",
      as.character(poc_meta_final$ClusterID))
)





# Full constrained ordination
poc_cap <- capscale(
  poc_dist_final ~
    Site_unique +
    Lineage_broad +
    Condition(readcount),
  data = poc_meta_final
)

summary(poc_cap)

# Overall significance
anova(
  poc_cap,
  permutations = 9999
)

# Marginal effects
anova(
  poc_cap,
  by = "margin",
  permutations = 9999
)

# Adjusted R²
RsquareAdj(poc_cap)


poc_cap_lineage <- capscale(
  poc_dist_final ~
    Lineage_broad +
    Condition(Site_unique + readcount),
  data = poc_meta_final
)

summary(poc_cap_lineage)

anova(
  poc_cap_lineage,
  permutations = 9999
)

RsquareAdj(poc_cap_lineage)

site_perm <- how(
  nperm = 9999,
  blocks = poc_meta_final$Site_unique
)

anova(
  poc_cap_lineage,
  permutations = site_perm
)

# Partition variation between site and genetic group,
# conditional on sequencing read count.

poc_vp <- varpart(
  poc_dist_final,
  ~ Site_unique,
  ~ Lineage_broad,
  ~ readcount,
  data = poc_meta_final
)

poc_vp
plot(poc_vp)





site_lineage_table <- table(
  poc_meta_pruned$Site_unique,
  poc_meta_pruned$Lineage_broad
)

chi_test <- chisq.test(site_lineage_table)
chi_test

library(lsr)
cramersV(site_lineage_table)









library(ggplot2)

# Extract partial dbRDA site scores
poc_scores <- as.data.frame(
  scores(poc_cap_lineage, display = "sites",
         choices = 1:2, scaling = 1)
)

poc_scores$Lineage <- poc_meta_final$Lineage_broad
poc_scores$Site <- poc_meta_final$Site_unique

ggplot(
  poc_scores,
  aes(x = CAP1, y = CAP2, color = Lineage)
) +
  geom_point(size = 6, alpha = 0.85) +
  theme_classic(base_size = 13) +
  labs(
    title = "Pocillopora bacterial microbiome",
    subtitle = "Partial dbRDA controlling for site and read count",
    x = "CAP1",
    y = "CAP2",
    color = "Genetic group"
  )

eig <- eigenvals(poc_cap_lineage, model = "constrained")

db1 <- 100 * eig[1] / sum(eig)
db2 <- 100 * eig[2] / sum(eig)










#261031:
library(dplyr)
library(tidyr)
library(tibble)

poc_ids <- rownames(poc_meta_final)
poc_genera <- rownames(keep.bacteria.wide.join.df)

poc_counts <- pocillopora.tara.bracken.output.df %>%
  mutate(sample = sub(
    "\\.bracken\\.genus\\.level\\.output\\.report$",
    "", sample
  )) %>%
  filter(
    sample %in% poc_ids,
    name %in% poc_genera
  ) %>%
  group_by(sample, name) %>%
  summarise(readcount = sum(new_est_reads),
            .groups = "drop") %>%
  complete(
    sample = poc_ids,
    name = poc_genera,
    fill = list(readcount = 0)
  ) %>%
  pivot_wider(
    names_from = name,
    values_from = readcount
  ) %>%
  column_to_rownames("sample") %>%
  as.matrix()

poc_counts <- poc_counts[
  poc_ids, poc_genera, drop = FALSE
]


library(vegan)

# Preserve the original 58 x 132 count matrix
poc_counts_unfiltered <- poc_counts

# Candidate contaminants identified so far
poc_exclude <- c(
  "Stenotrophomonas",
  "Pseudomonas",
  "Streptomyces"
)

# Remove contaminants
poc_counts_clean <- poc_counts[
  , !colnames(poc_counts) %in% poc_exclude,
  drop = FALSE
]

# Verify
dim(poc_counts_clean)
intersect(colnames(poc_counts_clean), poc_exclude)

# Quantify excluded reads
100 * (1 - sum(poc_counts_clean) /
         sum(poc_counts_unfiltered))

# Update sequencing-depth covariate
poc_meta_final$readcount <- rowSums(poc_counts_clean)

# Ordinary CLR transformation
poc_clr_corrected <- decostand(
  poc_counts_clean,
  method = "clr",
  pseudocount = 0.01
)

# Aitchison distance
poc_dist_corrected <- dist(
  poc_clr_corrected,
  method = "euclidean"
)

# Final checks
dim(poc_clr_corrected)
anyNA(poc_clr_corrected)
attr(poc_dist_corrected, "Size")
identical(
  rownames(poc_clr_corrected),
  rownames(poc_meta_final)
)


library(vegan)

poc_meta_final$Lineage_broad <- factor(
  poc_meta_final$Lineage_broad
)

poc_meta_final$Site_unique <- factor(
  poc_meta_final$Site_unique
)

set.seed(123)

poc_adonis_corrected <- adonis2(
  poc_dist_corrected ~ readcount + Site_unique + Lineage_broad,
  data = poc_meta_final,
  permutations = 9999,
  by = "margin"
)

poc_adonis_corrected


poc_site_lineage <- table(
  poc_meta_final$Site_unique,
  poc_meta_final$Lineage_broad
)

poc_site_lineage

# Number of sites containing multiple lineages
sum(rowSums(poc_site_lineage > 0) > 1)

# Check model rank
poc_model_matrix <- model.matrix(
  ~ readcount + Site_unique + Lineage_broad,
  data = poc_meta_final
)

c(
  rank = qr(poc_model_matrix)$rank,
  columns = ncol(poc_model_matrix)
)


library(vegan)
library(permute)

set.seed(123)

poc_adonis_restricted_corrected <- adonis2(
  poc_dist_corrected ~ readcount + Site_unique + Lineage_broad,
  data = poc_meta_final,
  permutations = how(
    nperm = 9999,
    blocks = poc_meta_final$Site_unique
  ),
  by = "margin"
)

poc_adonis_restricted_corrected

# Dispersion by lineage
set.seed(123)

poc_disp_lineage_corrected <- betadisper(
  poc_dist_corrected,
  poc_meta_final$Lineage_broad
)

permutest(
  poc_disp_lineage_corrected,
  permutations = 9999
)

# Dispersion by site
set.seed(123)

poc_disp_site_corrected <- betadisper(
  poc_dist_corrected,
  poc_meta_final$Site_unique
)

permutest(
  poc_disp_site_corrected,
  permutations = 9999
)


library(lsr)

# Cramér's V
lsr::cramersV(poc_site_lineage)

# Fisher's exact test with simulated p-value
set.seed(123)

fisher.test(
  poc_site_lineage,
  simulate.p.value = TRUE,
  B = 99999
)



library(vegan)

poc_dbrda_corrected <- dbrda(
  poc_dist_corrected ~ Lineage_broad +
    Condition(Site_unique + readcount),
  data = poc_meta_final
)

poc_dbrda_corrected

# Percent total inertia explained by each constrained axis
100 * eigenvals(
  poc_dbrda_corrected,
  model = "constrained"
) / poc_dbrda_corrected$tot.chi




library(ggplot2)
library(vegan)

# Partial dbRDA
poc_scores_corrected <- as.data.frame(
  scores(poc_dbrda_corrected, display = "sites", choices = 1:2)
)

poc_scores_corrected$Lineage <- poc_meta_final$Lineage_broad

poc_db_eig <- 100 * eigenvals(
  poc_dbrda_corrected, model = "constrained"
) / poc_dbrda_corrected$tot.chi

poc_dbrda_figure <- ggplot(
  poc_scores_corrected,
  aes(dbRDA1, dbRDA2, color = Lineage)
) +
  geom_point(size = 7, alpha = 0.85) +
  coord_fixed() +
  theme_classic() +
  labs(
    title = expression(italic("Pocillopora effusa")),
    x = sprintf("dbRDA1 (%.2f%%)", poc_db_eig[1]),
    y = sprintf("dbRDA2 (%.2f%%)", poc_db_eig[2]),
    color = "Lineage"
  )

print(poc_dbrda_figure)

ggsave("poc_dbRDA_corrected.pdf",
       plot = poc_dbrda_figure, width = 6, height = 5)

# Unconstrained PCA
poc_pca_corrected <- prcomp(poc_clr_corrected)

poc_pca_scores <- as.data.frame(poc_pca_corrected$x[, 1:2])
poc_pca_scores$Lineage <- poc_meta_final$Lineage_broad

poc_pca_var <- 100 * poc_pca_corrected$sdev^2 /
  sum(poc_pca_corrected$sdev^2)

poc_pca_figure <- ggplot(
  poc_pca_scores,
  aes(PC1, PC2, color = Lineage)
) +
  geom_point(size = 7, alpha = 0.85) +
  coord_fixed() +
  theme_classic() +
  labs(
    title = expression(italic("Pocillopora effusa")),
    x = sprintf("PC1 (%.2f%%)", poc_pca_var[1]),
    y = sprintf("PC2 (%.2f%%)", poc_pca_var[2]),
    color = "Lineage"
  )

print(poc_pca_figure)

ggsave("poc_PCA_corrected.pdf",
       plot = poc_pca_figure, width = 6, height = 5)





library(MASS)
library(adegenet)

# In-sample DFA accuracy
poc_dfa_corrected <- MASS::lda(
  poc_clr_corrected,
  grouping = poc_meta_final$Lineage_broad
)

poc_dfa_accuracy <- mean(
  predict(poc_dfa_corrected)$class ==
    poc_meta_final$Lineage_broad
) * 100

poc_dfa_accuracy

# DAPC cross-validation
set.seed(123)

poc_cv_corrected <- xvalDapc(
  poc_clr_corrected,
  grp = poc_meta_final$Lineage_broad,
  n.pca.max = 20,
  training.set = 0.7,
  result = "groupMean",
  center = TRUE,
  scale = FALSE,
  n.rep = 100,
  xval.plot = TRUE
)

poc_cv_corrected$`Number of PCs Achieving Highest Mean Success`
poc_cv_corrected$`Mean Successful Assignment by Number of PCs of PCA`



# Find the sample furthest along dbRDA1
poc_outlier_id <- rownames(poc_scores_corrected)[
  which.max(abs(poc_scores_corrected$dbRDA1))
]

poc_outlier_id

# Metadata
poc_meta_final[
  poc_outlier_id,
  c("Lineage_broad", "Site_unique", "readcount"),
  drop = FALSE
]

# Compare its retained read count with other samples
summary(poc_meta_final$readcount)



library(adegenet)

poc_dapc_corrected <- dapc(
  poc_clr_corrected,
  grp = poc_meta_final$Lineage_broad,
  n.pca = 10,
  n.da = 4
)

scatter(
  poc_dapc_corrected,
  scree.da = FALSE,
  legend = TRUE,
  main = "Pocillopora effusa — DAPC"
)

pdf("poc_DAPC_corrected.pdf", width = 7, height = 6)

scatter(
  poc_dapc_corrected,
  scree.da = FALSE,
  legend = TRUE,
  main = "Pocillopora effusa — DAPC"
)

dev.off()







#sensitivity for dropped sample:
# Identify its original count before filtering
c(
  original = sum(poc_counts_unfiltered["pocil.tenoko3", ]),
  retained = sum(poc_counts_clean["pocil.tenoko3", ])
)

# Compare its retained genus richness
sort(rowSums(poc_counts_clean > 0))[1:10]

# Sensitivity analysis: exclude the outlier
poc_keep_sensitivity <- rownames(poc_meta_final) != "pocil.tenoko3"

poc_meta_sensitivity <- droplevels(
  poc_meta_final[poc_keep_sensitivity, , drop = FALSE]
)

poc_dist_sensitivity <- as.dist(
  as.matrix(poc_dist_corrected)[
    poc_keep_sensitivity,
    poc_keep_sensitivity
  ]
)

set.seed(123)

poc_adonis_sensitivity <- vegan::adonis2(
  poc_dist_sensitivity ~ readcount + Site_unique + Lineage_broad,
  data = poc_meta_sensitivity,
  permutations = 9999,
  by = "margin"
)

poc_adonis_sensitivity






#261104:






library(ggplot2)

#same formatting used for Acropora and Montastraea
coral_theme <- theme_classic(base_size = 12) +
  theme(
    aspect.ratio = 1,
    panel.grid = element_blank(),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    plot.title = element_text(size = 13, face = "italic"),
    legend.position = "none",
    plot.margin = margin(8, 8, 8, 8)
  )

#standardized figures
poc_dbrda_final <- poc_dbrda_figure + coral_theme
poc_pca_final <- poc_pca_figure + coral_theme

#print for inspection
print(poc_dbrda_final)
print(poc_pca_final)

#export PDFs
ggsave(
  "poc_dbRDA_final.pdf",
  plot = poc_dbrda_final,
  width = 5,
  height = 5,
  units = "in"
)

ggsave(
  "poc_PCA_final.pdf",
  plot = poc_pca_final,
  width = 5,
  height = 5,
  units = "in"
)

library(adegenet)

pdf(
  "poc_DAPC_final.pdf",
  width = 5,
  height = 5
)

scatter(
  poc_dapc_corrected,
  scree.da = FALSE
)

dev.off()














































































































































































































































































