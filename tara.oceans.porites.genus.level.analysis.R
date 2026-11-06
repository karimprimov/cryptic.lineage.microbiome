
#installing packages:
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
library(broom)

porites.tara.bracken.output.df <- read_tsv("/Users/karimprimov/Dropbox/UTAustin/2026/mcav.cryptic.lineages/mcav.flkeys.metagenomics/tara.oceans.bracken.reports/260422.porites.lobata.tara.genus.level.bracken.taxonomic.report.tsv")
porites.metadata <- read_csv("/Users/karimprimov/Dropbox/UTAustin/2026/sp26/cryptic.lineage.microbiome/tara.oceans/porites.lobata.metadata.csv")
porites.metadata <- porites.metadata %>%
  column_to_rownames(var = "...1")
#filtering out contaminant, human, and other unwanted microbe taxa:
porites.microbe.taxa <- unique(porites.tara.bracken.output.df$name)
write.csv(microbe.taxa, "260422.porites.lobata.tara.microbe.taxa.screencheck.csv", row.names = FALSE)


#filtering to keep these taxa: 

keep_taxa <- c(
  
  # =========================
  # CORE REEF MICROBIOME
  # =========================
  
  # Coral-associated (stable symbiosis / mucus / host interface)
  "Endozoicomonas",
  "Ruegeria",
  "Sulfitobacter",
  "Phaeobacter",
  "Roseobacter",
  "Roseovarius",
  "Salipiger",
  "Labrenzia",
  "Paracoccus",
  "Devosia",
  "Jannaschia",
  "Gemmobacter",
  "Roseibium",
  
  # Cyanobacteria (primary production)
  "Synechococcus",
  "Prochlorococcus",
  "Trichodesmium",
  "Nostoc",
  
  # Strong marine polysaccharide degraders
  "Flavobacterium",
  "Polaribacter",
  "Winogradskyella",
  "Formosa",
  "Maribacter",
  "Algibacter",
  "Aquimarina",
  "Cellulophaga",
  "Zobellia",
  "Dokdonia",
  "Nonlabens",
  "Kordia",
  "Tenacibaculum",
  "Flammeovirga",
  
  # Classic marine Gammaproteobacteria (carbon cycling core)
  "Alteromonas",
  "Pseudoalteromonas",
  "Marinobacter",
  "Marinomonas",
  "Shewanella",
  "Colwellia",
  "Psychrobacter",
  "Vibrio",
  "Photobacterium",
  "Halomonas",
  "Salinivibrio",
  "Aliivibrio",
  
  # Sulfur + redox cycling core
  "Thiomicrorhabdus",
  "Sulfurimonas",
  "Desulfovibrio",
  "Desulfobacter",
  "Desulfuromonas",
  "Geobacter",
  "Beggiatoa",
  
  # Nitrogen cycling core
  "Nitrosomonas",
  "Nitrosococcus",
  "Methylophaga",
  
  # Alphaproteobacteria marine generalists
  "Sphingomonas",
  "Novosphingobium",
  "Erythrobacter",
  "Parerythrobacter",
  "Caulobacter",
  "Hyphomonas",
  
  # =========================
  # EDGE CASE / CONTEXT-DEPENDENT TAXA
  # (IMPORTANT: include to avoid biasing toward shared-core signal)
  # =========================
  
  # Opportunistic / bloom / conditionally present reef bacteria
  "Pseudomonas",
  "Acinetobacter",
  "Burkholderia",
  "Paraburkholderia",
  "Cupriavidus",
  "Comamonas",
  "Hydrogenophaga",
  "Variovorax",
  "Delftia",
  "Stenotrophomonas",
  
  # Vibrios (ecologically real but highly state-dependent)
  "Vibrio",
  "Photobacterium",
  
  # Flavobacteria that can be host-associated OR transient
  "Chryseobacterium",
  "Pedobacter",
  "Sphingobacterium",
  "Mucilaginibacter",
  "Hymenobacter",
  
  # Rare-but-real marine heterotrophs (often lineage-sensitive)
  "Thalassotalea",
  "Neptunomonas",
  "Idiomarina",
  "Catenovulum",
  "Ferrimonas",
  
  # Sulfur/metal cycling edge taxa (often habitat-specific)
  "Desulfatibacillum",
  "Desulfobulbus",
  "Syntrophus",
  "Anaeromyxobacter",
  
  # Methane / one-carbon cycling (often patchy but meaningful)
  "Methylotenera",
  "Methyloversatilis",
  
  # Alphaproteobacteria that can swing between free-living and host-associated
  "Rhizobium",
  "Bradyrhizobium",
  "Azospirillum",
  "Pararhizobium",
  "Sinorhizobium",
  "Mesorhizobium",
  
  # Cyanobacteria that may reflect transient blooms
  "Microcystis",
  "Planktothrix"
)

keep.bacteria.df <- porites.tara.bracken.output.df[porites.tara.bracken.output.df$name %in% keep_taxa,]

keep.bacteria.wide.join.df <- keep.bacteria.df %>%
  dplyr::select(Sample = sample, Taxa = name, readcount = new_est_reads) %>%
  pivot_wider(names_from = "Sample", values_from = "readcount", values_fn = sum) %>%
  column_to_rownames(var = "Taxa") %>%
  mutate(across(everything(), ~replace(., is.na(.), 1)))

t.keep.bacteria.wide.join.df <- t(keep.bacteria.wide.join.df)
rownames(t.keep.bacteria.wide.join.df) <- sub(".bracken.genus.level.output.report", "", rownames(t.keep.bacteria.wide.join.df), fixed = TRUE)

porites.metadata_new <- porites.metadata[match(rownames(t.keep.bacteria.wide.join.df), rownames(porites.metadata)), ]
porites.metadata_new$readcount <- rowSums(t.keep.bacteria.wide.join.df, na.rm = TRUE)

stand.porites.microbes <- decostand(t.keep.bacteria.wide.join.df, method = "rclr")

porites.bacteria.dist.matrix <- as.matrix(dist(stand.porites.microbes, method = "euclidean"))

capscale.porites.bacteria <- capscale(porites.bacteria.dist.matrix ~ Site + ClusterID + Condition(readcount), data = porites.metadata_new)

plot(capscale.porites.bacteria)


porites.microbe.scores <- scores(capscale.porites.bacteria, display = "sites", choices=c(1:2)) 

merged.fam.porites.microbe.scores.metadata <- merge(porites.microbe.scores, porites.metadata_new, by = "row.names", all.x = TRUE)

ggplot(merged.fam.porites.microbe.scores.metadata, aes(x = CAP1, y = CAP2, color =as.factor(Longitude), label = rownames(merged.fam.porites.microbe.scores.metadata))) +
  geom_point(size = 12) +
  theme_minimal() + 
  coord_fixed()


adonis.site.results <- adonis2(porites.bacteria.dist.matrix ~ readcount + ClusterID + Site, 
                               data = porites.metadata_new,
                               by = "margin")   # marginal tests are best for independent effects


#260424:
#discriminant function analysis:

#1: cluster
cluster.lda_model <- lda(ClusterID ~ MDS1 + MDS2, data = merged.fam.porites.microbe.scores.metadata)
print(cluster.lda_model)
cluster.predictions <- predict(cluster.lda_model)$class
table(Actual = merged.fam.porites.microbe.scores.metadata$ClusterID, Predicted = cluster.predictions)
mean(predict(cluster.lda_model)$class == merged.fam.porites.microbe.scores.metadata$ClusterID)


#2: site
site.lda_model <- lda(Site ~ MDS1 + MDS2, data = merged.fam.porites.microbe.scores.metadata)
print(site.lda_model)
site.predictions <- predict(site.lda_model)$class
table(Actual = merged.fam.porites.microbe.scores.metadata$Site, Predicted = site.predictions)
mean(predict(site.lda_model)$class == merged.fam.porites.microbe.scores.metadata$Site)




adonis.site.results <- adonis2(porites.bacteria.dist.matrix ~ readcount + Site, 
                               data = porites.metadata_new,
                               by = "margin")   # marginal tests are best for independent effects

# 2. View results
print(adonis.site.results)


adonis.site.lineage.results <- adonis2(porites.bacteria.dist.matrix ~ readcount + Site + ClusterID, 
                                       data = porites.metadata_new,
                                       by = "margin")   # marginal tests are best for independent effects

adonis.lineage.results <- adonis2(porites.bacteria.dist.matrix ~ readcount + ClusterID, 
                                       data = porites.metadata_new,
                                       by = "margin")   # marginal tests are best for independent effects

# 2. View results
print(adonis.site.lineage.results)




######
#conditioning by readcount and site:
capscale.porites.bacteria.no.rc.no.site <- capscale(porites.bacteria.dist.matrix ~ ClusterID + Condition(readcount+Site), data = porites.metadata_new )
capscale.porites.bacteria.no.rc.no.lineage <- capscale(porites.bacteria.dist.matrix ~ Site + Condition(readcount+ClusterID), data = porites.metadata_new )

plot(capscale.porites.bacteria.no.rc.no.site)


porites.microbe.no.rc.no.site.scores <- scores(capscale.porites.bacteria.no.rc.no.site, display = "sites", choices=c(1:2)) 

merged.fam.porites.microbe.scores.no.rc.no.site.scores.metadata <- merge(porites.microbe.no.rc.no.site.scores, porites.metadata_new, by = "row.names", all.x = TRUE)

ggplot(merged.fam.porites.microbe.scores.no.rc.no.site.scores.metadata, aes(x = CAP1, y = CAP2, color =Longitude)) +
  geom_point(size = 12) +
  theme_minimal() + 
  coord_fixed()




###
env.porites.metadata_new <- porites.metadata_new[,c(2,3,4,12)]
rownames(env.porites.metadata_new) <- env.porites.metadata_new $Site


pheatmap(env.porites.metadata_new)




#260424:
#seeing if site and lineage correlate:
site.cluster.metadata <- table(meta75$Site_unique, meta75$ClusterID)
chi_test.site.cluster <- chisq.test(site.cluster.metadata)
cramersV(site.cluster.metadata) 



porites.db_readcount_model <- dbrda(porites.bacteria.dist.matrix ~ ClusterID + Site + Condition(readcount), 
                       data = porites.metadata_new)

# 2. Get the PERMANOVA table (identical results to adonis2)
anova(porites.db_readcount_model, by = "margin")

#variance partitioning
varpart(porites.bacteria.dist.matrix, 
        ~ Site, 
        ~ ClusterID, 
        ~ readcount, 
        data = porites.metadata_new)

interaction(Site, ClusterID)

interaction.site.cluster <- lm(porites.bacteria.dist.matrix ~ Site *ClusterID, data = porites.metadata_new)
tidy(interaction.site.cluster) %>% 
  filter(p.value < 0.05)
anova(interaction.site.cluster)

simple_model <- lm(porites.bacteria.dist.matrix ~ Site + ClusterID, data = porites.metadata_new)
anova(simple_model)


disp_site <- betadisper(as.dist(porites.bacteria.dist.matrix), porites.metadata_new$Site)
permutest(disp_site) 

# For ClusterID
disp_cluster <- betadisper(as.dist(porites.bacteria.dist.matrix), porites.metadata_new$ClusterID)
permutest(disp_cluster)

porites.disp.cluster.df <- data.frame(group = disp_cluster$group, distance = disp_cluster$distances, Species = "massive Porites")
plot(porites.disp.cluster.df)
write.csv(porites.disp.cluster.df, "porites.tara.beta.disp.distances.csv", row.names =TRUE)


table(porites.metadata_new$Site, porites.metadata_new$ClusterID)

# 2. Chi-Squared Test (Statistical test of association)
test <- chisq.test(table(porites.metadata_new$Site, porites.metadata_new$ClusterID))
print(test)


#genetic ditances vs. microbiome distances
porites.microbiome.centroid.distances <- as.matrix(dist(disp_cluster$centroids))
porites.genetic.distances <- read.csv("/Users/karimprimov/Dropbox/UTAustin/2026/mcav.cryptic.lineages/mcav.flkeys.metagenomics/porites.genetic.distances.csv", row.names = 1)
porites.genetic.distances.matrix <- as.matrix(porites.genetic.distances)

porites.microbiome.centroid.distances.match <- porites.microbiome.centroid.distances[rownames(porites.microbiome.centroid.distances) %in% rownames(porites.genetic.distances.matrix), colnames(porites.microbiome.centroid.distances) %in% colnames(porites.genetic.distances.matrix)]

porites.microbe.dist.ordination <- capscale(porites.microbiome.centroid.distances.match ~1)
porites.genetic.distances.ordination <- capscale(porites.genetic.distances.matrix ~1)
plot(porites.microbe.dist.ordination)
plot(porites.genetic.distances.ordination)
procrustes.microbe.genetic.distance <- protest(porites.microbe.dist.ordination, porites.genetic.distances.ordination)



# Convert matrices to vectors
gen_dist_vec <- as.vector(as.dist(porites.genetic.distances.matrix))
mic_dist_vec <- as.vector(as.dist(porites.microbiome.centroid.distances.match))

# Create a data frame for plotting
dist_df <- data.frame(
  GeneticDistance = gen_dist_vec,
  MicrobiomeDistance = mic_dist_vec
)


ggplot(dist_df, aes(x = GeneticDistance, y = MicrobiomeDistance)) +
  geom_point(color = "steelblue", size = 3, alpha = 0.7) +
  geom_smooth(method = "lm", color = "firebrick", se = TRUE) + # Add regression line
  theme_classic() +
  labs(
    title = "Isolation by Distance: Genetics vs. Microbiome",
    x = "Genetic Distance (Centroids)",
    y = "Microbiome Dissimilarity (Centroids)"
  ) +
  # Add the Mantel R and P-value manually or via annotate
  annotate("text", x = min(dist_df$GeneticDistance), y = max(dist_df$MicrobiomeDistance), 
           label = paste("Mantel r =", round(procrustes.microbe.genetic.distance$t0, 3)), 
           hjust = 0, vjust = 1, fontface = "italic")





#

#Saving RData for the day!
#######################################################################################################################################
save.image("/Users/karimprimov/Dropbox/UTAustin/2026/mcav.cryptic.lineages/mcav.flkeys.metagenomics/261031.porites.bracken.genuslevel.analysis.RData")




















#260929:
# NOTE:
# Partial capscale models return F values matching marginal adonis2
# but permutation p-values = 1 under this heavily aliased design.
# Do not report those p-values unless resolved.

porites.metadata_new$ClusterID <- factor(porites.metadata_new$ClusterID)
porites.metadata_new$Site_unique <- factor(porites.metadata_new$Site_unique)

# common subset for fair comparison
keep <- !is.na(porites.metadata_new$ClusterID)

meta75 <- droplevels(porites.metadata_new[keep, ])

dist75 <- as.dist(
  as.matrix(porites.bacteria.dist.matrix)[keep, keep]
)

set.seed(123)

adonis_site <- adonis2(
  dist75 ~ readcount + Site_unique,
  data = meta75,
  by = "margin",
  permutations = 9999
)

adonis_lineage <- adonis2(
  dist75 ~ readcount + Lineage_broad,
  data = meta75,
  by = "margin",
  permutations = 9999
)

adonis_both <- adonis2(
  dist75 ~ readcount + Site_unique + Lineage_broad,
  data = meta75,
  by = "margin",
  permutations = 9999
)

adonis_site
adonis_lineage
adonis_both


# Site dispersion
bd.site <- betadisper(
  dist75,
  meta75$Site_unique
)
permutest(bd.site, permutations = 9999)

# Lineage dispersion
bd.lineage <- betadisper(
  dist75,
  meta75$ClusterID
)
permutest(bd.lineage, permutations = 9999)




capscale.porites.bacteria.site.lineage <- capscale(
  dist75 ~ Site_unique + Lineage_broad + Condition(readcount),
  data = meta75
)

capscale.porites.bacteria.nosite <- capscale(
  dist75 ~ Lineage_broad + Condition(readcount + Site_unique),
  data = meta75
)

capscale.porites.bacteria.no.lineage <- capscale(
  dist75 ~ Site_unique + Condition(readcount + Lineage_broad),
  data = meta75
)


porites.microbe.scores.nosite <- scores(capscale.porites.bacteria.nosite, display = "sites", choices=c(1:2)) 

merged.fam.porites.porites.microbe.scores.nosite <- merge(porites.microbe.scores.nosite, meta75, by = "row.names", all.x = TRUE)

ggplot(merged.fam.porites.porites.microbe.scores.nosite, aes(x = CAP1, y = CAP2, color =ClusterID, label = rownames(merged.fam.porites.porites.microbe.scores.nosite))) +
  geom_point(size = 12) +
  theme_minimal() + 
  coord_fixed()



porites.microbe.scores.no.lineage <- scores(capscale.porites.bacteria.no.lineage, display = "sites", choices=c(1:2)) 

merged.fam.porites.microbe.scores.no.lineage <- merge(porites.microbe.scores.no.lineage, meta75, by = "row.names", all.x = TRUE)

ggplot(merged.fam.porites.microbe.scores.no.lineage, aes(x = CAP1, y = CAP2, color =Site_unique, label = rownames(merged.fam.porites.microbe.scores.no.lineage))) +
  geom_point(size = 12) +
  theme_minimal() + 
  coord_fixed()





porites.microbe.scores.site.lineage <- scores(capscale.porites.bacteria.site.lineage, display = "sites", choices=c(1:2)) 

merged.fam.porites.microbe.scores.site.lineage <- merge(porites.microbe.scores.site.lineage, meta75, by = "row.names", all.x = TRUE)

ggplot(merged.fam.porites.microbe.scores.site.lineage, aes(x = CAP1, y = CAP2, color =ClusterID, label = rownames(merged.fam.porites.microbe.scores.site.lineage))) +
  geom_point(size = 12) +
  theme_minimal() + 
  coord_fixed()


set.seed(123)

disp_site <- betadisper(
  dist75,
  meta75$Site_unique
)

permutest(disp_site, permutations = 9999)

set.seed(123)

disp_cluster <- betadisper(
  dist75,
  meta75$ClusterID
)

permutest(disp_cluster, permutations = 9999)


varpart(
  dist75,
  ~ Site_unique,
  ~ ClusterID,
  ~ readcount,
  data = meta75
)

set.seed(123)

anova(
  capscale.porites.bacteria.nosite,
  permutations = 9999
)

set.seed(123)

anova(
  capscale.porites.bacteria.no.lineage,
  permutations = 9999
)












#fixing buggy result
porites_dbrda_full <- dbrda(
  dist75 ~ readcount + Site_unique + Lineage_broad,
  data = meta75
)

set.seed(123)

anova(
  porites_dbrda_full,
  by = "margin",
  permutations = 9999
)



set.seed(123)
anova(
  capscale.porites.bacteria.nosite,
  permutations = 9999,
  model = "reduced"
)

set.seed(123)
anova(
  capscale.porites.bacteria.nosite,
  permutations = 9999,
  model = "direct"
)

set.seed(123)
anova(
  capscale.porites.bacteria.nosite,
  permutations = 9999,
  model = "full"
)




set.seed(123)
anova(
  capscale.porites.bacteria.no.lineage,
  permutations = 9999,
  model = "reduced"
)

set.seed(123)
anova(
  capscale.porites.bacteria.no.lineage,
  permutations = 9999,
  model = "direct"
)

set.seed(123)
anova(
  capscale.porites.bacteria.no.lineage,
  permutations = 9999,
  model = "full"
)



set.seed(123)

perm_lineage <- permutest(
  capscale.porites.bacteria.nosite,
  permutations = 9999,
  model = "reduced"
)

perm_lineage

summary(perm_lineage$F.perm)
range(perm_lineage$F.perm)







library(vegan)
library(ggplot2)

# Final dbRDA
porites_dbrda_plot <- dbrda(
  dist75 ~ Site_unique + ClusterID + Condition(readcount),
  data = meta75
)

# Extract site/sample scores
porites_scores <- as.data.frame(
  scores(
    porites_dbrda_plot,
    display = "sites",
    choices = 1:2
  )
)

# Preserve sample IDs
porites_scores$SampleID_plot <- rownames(porites_scores)

# Add metadata using row names
meta75$SampleID_plot <- rownames(meta75)

porites_plot_data <- merge(
  porites_scores,
  meta75,
  by = "SampleID_plot",
  all.x = TRUE
)

# Check names
colnames(porites_plot_data)

eig <- eigenvals(porites_dbrda_plot, model = "constrained")

eig

db1 <- 100 * eig[1] / sum(eig)
db2 <- 100 * eig[2] / sum(eig)

db1
db2


ggplot(
  porites_plot_data,
  aes(
    x = dbRDA1,
    y = dbRDA2,
    color = Site_unique
  )
) +
  geom_point(size = 10, alpha = 0.85) +
  theme_minimal() +
  coord_fixed() +
  labs(
    x = paste0("dbRDA1 (", round(db1, 1), "%)"),
    y = paste0("dbRDA2 (", round(db2, 1), "%)"),
    color = "Host lineage"
  )




site.cluster.metadata <- table(meta75$Site_unique, meta75$ClusterID)
chi_test.site.cluster <- chisq.test(site.cluster.metadata)
cramersV(site.cluster.metadata) 



meta75$Lineage_broad <- sub("I[0-9]+$", "", as.character(meta75$ClusterID))

table(meta75$Lineage_broad)
table(meta75$Site_unique, meta75$Lineage_broad)















keep <- meta75$Lineage_broad %in% c("K2", "K3")

meta_lobata <- droplevels(meta75[keep, ])
dist_lobata <- as.dist(as.matrix(dist75)[keep, keep])

adonis2(
  dist_lobata ~ readcount + Site_unique + Lineage_broad,
  data = meta_lobata,
  by = "margin",
  permutations = 9999
)







#261009:

library(vegan)
library(permute)

# Verify full dataset alignment
stopifnot(
  identical(
    rownames(meta75),
    attr(dist75, "Labels")
  )
)

# Check available genetic groups
table(meta75$Lineage_broad)

# Select the two cryptic P. lobata groups
keep_por <- as.character(meta75$Lineage_broad) %in%
  c("K2", "K3")

# Subset metadata and distance matrix together
por_meta <- droplevels(
  meta75[keep_por, , drop = FALSE]
)

por_dist <- as.dist(
  as.matrix(dist75)[keep_por, keep_por]
)

# Confirm alignment
stopifnot(
  identical(
    rownames(por_meta),
    attr(por_dist, "Labels")
  )
)

# Inspect the final dataset
nrow(por_meta)
table(por_meta$Lineage_broad)
length(unique(por_meta$Site_unique))



# Inspect distance metadata
attributes(por_dist)

# Check metadata columns
names(por_meta)

# Confirm predictors are present
stopifnot(
  all(c("readcount", "Site_unique",
        "Lineage_broad") %in% names(por_meta))
)

# Check for missing values
colSums(is.na(
  por_meta[, c("readcount",
               "Site_unique",
               "Lineage_broad")]
))

set.seed(123)

por_perm <- adonis2(
  por_dist ~ readcount + Site_unique + Lineage_broad,
  data = por_meta,
  by = "margin",
  permutations = 9999
)

por_perm


por_site_perm <- how(
  nperm = 9999,
  blocks = por_meta$Site_unique
)

set.seed(123)

por_perm_restricted <- adonis2(
  por_dist ~ readcount + Site_unique + Lineage_broad,
  data = por_meta,
  by = "margin",
  permutations = por_site_perm
)

por_perm_restricted

por_disp <- betadisper(
  por_dist,
  group = factor(por_meta$Lineage_broad)
)

set.seed(123)

permutest(
  por_disp,
  permutations = 9999
)

por_cap_lineage <- capscale(
  por_dist ~ Lineage_broad +
    Condition(Site_unique + readcount),
  data = por_meta
)

por_cap_lineage

set.seed(123)
anova(por_cap_lineage, permutations = por_site_perm)

RsquareAdj(por_cap_lineage)


por_disp_site <- betadisper(
  por_dist,
  group = factor(por_meta$Site_unique)
)

set.seed(123)
permutest(por_disp_site, permutations = 9999)

por_meta[, c(
  "IO_Number_Dashed",
  "Lineage_broad",
  "Site_unique"
)]

# Number of colonies in each genetic group at each site
por_site_table <- table(
  por_meta$Site_unique,
  por_meta$Lineage_broad
)

por_site_table

# Sites containing both K2 and K3
mixed_sites <- rownames(por_site_table)[
  rowSums(por_site_table > 0) == 2
]

length(mixed_sites)
mixed_sites

# Number of colonies at mixed sites
sum(por_meta$Site_unique %in% mixed_sites)

library(ggplot2)

por_scores <- as.data.frame(
  scores(
    por_cap_lineage,
    display = "sites",
    choices = 1,
    scaling = 1
  )
)

stopifnot(nrow(por_scores) == nrow(por_meta))

por_scores$Lineage <- por_meta$Lineage_broad

ggplot(
  por_scores,
  aes(x = Lineage, y = CAP1, fill = Lineage)
) +
  geom_boxplot(
    width = 0.4,
    alpha = 0.3,
    outlier.shape = NA
  ) +
  geom_jitter(
    width = 0.12,
    size = 2.5,
    alpha = 0.8
  ) +
  theme_classic(base_size = 13) +
  labs(
    title = "Porites lobata complex",
    subtitle = "Partial dbRDA controlling for site and read count",
    x = "Cryptic host species",
    y = "CAP1 (100% of constrained variation)"
  ) +
  guides(fill = "none")


library(vegan)
library(ggplot2)

# Extract site scores from the partial dbRDA
por_ord <- as.data.frame(
  scores(
    por_cap_lineage,
    display = "sites",
    choices = 1:2,
    scaling = 1
  )
)

# Inspect the available axes
names(por_ord)

# Attach metadata in the verified sample order
stopifnot(
  identical(rownames(por_ord), rownames(por_meta))
)

por_ord$Lineage <- por_meta$Lineage_broad

# Two-dimensional partial dbRDA ordination
ggplot(
  por_ord,
  aes(x = CAP1, y = MDS1, color = Lineage)
) +
  geom_point(size = 3, alpha = 0.85) +
  scale_color_manual(
    values = c("K2" = "#E69F9A",
               "K3" = "#64BFC4")
  ) +
  theme_classic(base_size = 13) +
  labs(
    title = "Porites lobata complex",
    subtitle = "Partial dbRDA controlling for site and read count",
    x = "CAP1 (100% of constrained variation)",
    y = "MDS1 (unconstrained variation)",
    color = "Cryptic species"
  ) +
  theme(
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(margin = margin(b = 12))
  )




# Keep colonies from sites containing both K2 and K3
por_mixed_keep <- por_meta$Site_unique %in% mixed_sites

por_meta_mixed <- droplevels(
  por_meta[por_mixed_keep, , drop = FALSE]
)

por_dist_mixed <- as.dist(
  as.matrix(por_dist)[por_mixed_keep, por_mixed_keep]
)

# Verify alignment
stopifnot(
  identical(
    rownames(por_meta_mixed),
    attr(por_dist_mixed, "Labels")
  )
)

# Verify sample sizes
nrow(por_meta_mixed)
table(por_meta_mixed$Lineage_broad)

# Site-restricted permutations
por_mixed_perms <- how(
  nperm = 9999,
  blocks = por_meta_mixed$Site_unique
)

set.seed(123)

por_mixed_perm <- adonis2(
  por_dist_mixed ~ readcount +
    Site_unique + Lineage_broad,
  data = por_meta_mixed,
  by = "margin",
  permutations = por_mixed_perms
)

por_mixed_perm



















#261013:
library(vegan)
library(permute)

# 1. Joint marginal PERMANOVA
set.seed(123)

por_perm <- adonis2(
  por_dist ~ readcount + Site_unique + Lineage_broad,
  data = por_meta,
  by = "margin",
  permutations = 9999
)

por_perm

# 2. Site-restricted PERMANOVA
por_site_perm <- how(
  nperm = 9999,
  blocks = por_meta$Site_unique
)

set.seed(123)

por_perm_restricted <- adonis2(
  por_dist ~ readcount + Site_unique + Lineage_broad,
  data = por_meta,
  by = "margin",
  permutations = por_site_perm
)

por_perm_restricted

# 3. Lineage dispersion
por_disp <- betadisper(
  por_dist,
  group = factor(por_meta$Lineage_broad)
)

set.seed(123)
permutest(por_disp, permutations = 9999)

# 4. Site dispersion
por_disp_site <- betadisper(
  por_dist,
  group = factor(por_meta$Site_unique)
)

set.seed(123)
permutest(por_disp_site, permutations = 9999)

# 5. Site–lineage association
site_lineage_table <- table(
  por_meta$Site_unique,
  por_meta$Lineage_broad
)

lsr::cramersV(site_lineage_table)

set.seed(123)
fisher.test(
  site_lineage_table,
  simulate.p.value = TRUE,
  B = 99999
)



por_cap_lineage <- capscale(
  por_dist ~ Lineage_broad +
    Condition(Site_unique + readcount),
  data = por_meta
)

por_cap_lineage
RsquareAdj(por_cap_lineage)

set.seed(123)
anova(por_cap_lineage, permutations = 9999)

set.seed(123)
anova(por_cap_lineage, permutations = por_site_perm)


por_cap_site <- capscale(
  por_dist ~  + Site_unique +
    Condition(Lineage_broad + readcount),
  data = por_meta
)

por_cap_site
RsquareAdj(por_cap_site)

set.seed(123)
anova(por_cap_site, permutations = 9999)


por_lin_perm <- how(
  nperm = 9999,
  blocks = por_meta$Lineage_broad
)

set.seed(123)
anova(por_cap_site, permutations = por_lin_perm)


por_cap_site_lineage <- capscale(
  por_dist ~  + Site_unique + Lineage_broad +
  Condition(readcount),
  data = por_meta
)

por_cap_site_lineage
RsquareAdj(por_cap_site_lineage)

set.seed(123)
anova(por_cap_site_lineage, permutations = 9999)



aggregate(
  por_disp$distances,
  by = list(Lineage = por_meta$Lineage_broad),
  FUN = mean
)



library(vegan)
library(permute)

# Identify sites containing both lineages
site_table <- table(
  por_meta$Site_unique,
  por_meta$Lineage_broad
)

mixed_sites <- rownames(site_table)[
  rowSums(site_table > 0) == 2
]

# Subset metadata and distance matrix
keep <- por_meta$Site_unique %in% mixed_sites

por_meta_mixed <- droplevels(
  por_meta[keep, , drop = FALSE]
)

por_dist_mixed <- as.dist(
  as.matrix(por_dist)[keep, keep]
)

# Verify sample alignment
stopifnot(identical(
  rownames(por_meta_mixed),
  attr(por_dist_mixed, "Labels")
))

# Check sample sizes
nrow(por_meta_mixed)
table(por_meta_mixed$Lineage_broad)

# Permute samples only within sites
mixed_perms <- how(
  nperm = 9999,
  blocks = por_meta_mixed$Site_unique
)

# Joint marginal PERMANOVA
set.seed(123)

por_mixed_perm <- adonis2(
  por_dist_mixed ~ readcount + Site_unique + Lineage_broad,
  data = por_meta_mixed,
  by = "margin",
  permutations = mixed_perms
)

por_mixed_perm





#261031:

library(dplyr)
library(tidyr)
library(tibble)

# Final K2/K3 sample IDs
por_ids <- rownames(por_meta)

# Genus list from your original Porites script
por_keep_taxa <- c(
  "Endozoicomonas", "Ruegeria", "Sulfitobacter",
  "Phaeobacter", "Roseobacter", "Roseovarius",
  "Salipiger", "Labrenzia", "Paracoccus",
  "Devosia", "Jannaschia", "Gemmobacter", "Roseibium",
  "Synechococcus", "Prochlorococcus",
  "Trichodesmium", "Nostoc",
  "Flavobacterium", "Polaribacter", "Winogradskyella",
  "Formosa", "Maribacter", "Algibacter",
  "Aquimarina", "Cellulophaga", "Zobellia",
  "Dokdonia", "Nonlabens", "Kordia",
  "Tenacibaculum", "Flammeovirga",
  "Alteromonas", "Pseudoalteromonas",
  "Marinobacter", "Marinomonas", "Shewanella",
  "Colwellia", "Psychrobacter", "Vibrio",
  "Photobacterium", "Halomonas", "Salinivibrio",
  "Aliivibrio", "Thiomicrorhabdus", "Sulfurimonas",
  "Desulfovibrio", "Desulfobacter", "Desulfuromonas",
  "Geobacter", "Beggiatoa", "Nitrosomonas",
  "Nitrosococcus", "Methylophaga", "Sphingomonas",
  "Novosphingobium", "Erythrobacter",
  "Parerythrobacter", "Caulobacter", "Hyphomonas",
  "Pseudomonas", "Acinetobacter", "Burkholderia",
  "Paraburkholderia", "Cupriavidus", "Comamonas",
  "Hydrogenophaga", "Variovorax", "Delftia",
  "Stenotrophomonas", "Chryseobacterium",
  "Pedobacter", "Sphingobacterium",
  "Mucilaginibacter", "Hymenobacter",
  "Thalassotalea", "Neptunomonas",
  "Idiomarina", "Catenovulum", "Ferrimonas",
  "Desulfatibacillum", "Desulfobulbus",
  "Syntrophus", "Anaeromyxobacter",
  "Methylotenera", "Methyloversatilis",
  "Rhizobium", "Bradyrhizobium", "Azospirillum",
  "Pararhizobium", "Sinorhizobium",
  "Mesorhizobium", "Microcystis", "Planktothrix"
)

# Reconstruct observed read counts, filling absent taxa with zero
por_counts <- porites.tara.bracken.output.df %>%
  mutate(sample = sub(
    ".bracken.genus.level.output.report",
    "", sample, fixed = TRUE
  )) %>%
  filter(
    sample %in% por_ids,
    name %in% por_keep_taxa
  ) %>%
  group_by(sample, name) %>%
  summarise(
    readcount = sum(new_est_reads),
    .groups = "drop"
  ) %>%
  complete(
    sample = por_ids,
    name = unique(por_keep_taxa),
    fill = list(readcount = 0)
  ) %>%
  pivot_wider(
    names_from = name,
    values_from = readcount
  ) %>%
  column_to_rownames("sample") %>%
  as.matrix()

por_counts <- por_counts[
  por_ids, unique(por_keep_taxa), drop = FALSE
]

dim(por_counts)
anyNA(por_counts)
sum(por_counts == 0)
sum(por_counts == 1)

identical(
  rownames(por_counts),
  rownames(por_meta)
)

sort(colSums(por_counts), decreasing = TRUE)[1:15]

library(vegan)

# Preserve original reconstructed counts
por_counts_unfiltered <- por_counts

# Candidate contaminant exclusions
por_exclude <- c(
  "Pseudomonas",
  "Stenotrophomonas"
)

# Remove candidate contaminants and entirely absent genera
por_counts_clean <- por_counts[
  ,
  !colnames(por_counts) %in% por_exclude &
    colSums(por_counts) > 0,
  drop = FALSE
]

# Quantify filtering
dim(por_counts_clean)

100 * (
  1 - sum(por_counts_clean) /
    sum(por_counts_unfiltered)
)

# Inspect retained sequencing depths
summary(rowSums(por_counts_clean))
sort(rowSums(por_counts_clean))[1:5]

# Ensure lineage and site are factors
por_meta$Lineage_broad <- factor(por_meta$Lineage_broad)
por_meta$Site_unique <- factor(por_meta$Site_unique)

# Read-count covariate
por_meta$readcount <- rowSums(por_counts_clean)

# Ordinary CLR with 0.01 pseudocount
por_clr_corrected <- decostand(
  por_counts_clean,
  method = "clr",
  pseudocount = 0.01
)

# Euclidean distance on CLR = Aitchison distance
por_dist_corrected <- dist(
  por_clr_corrected,
  method = "euclidean"
)

# Verify
dim(por_clr_corrected)
anyNA(por_clr_corrected)
attr(por_dist_corrected, "Size")

identical(
  rownames(por_clr_corrected),
  rownames(por_meta)
)



library(vegan)
library(lsr)

# 1. Marginal PERMANOVA
set.seed(123)

por_adonis_corrected <- adonis2(
  por_dist_corrected ~
    readcount + Site_unique + Lineage_broad,
  data = por_meta,
  permutations = 9999,
  by = "margin"
)

por_adonis_corrected

# 2. Site × lineage contingency table
por_site_lineage_corrected <- table(
  por_meta$Site_unique,
  por_meta$Lineage_broad
)

por_site_lineage_corrected

# 3. Sites containing both lineages
sum(rowSums(por_site_lineage_corrected > 0) == 2)

# 4. Cramér's V
lsr::cramersV(por_site_lineage_corrected)

# 5. Fisher's exact test
set.seed(123)

fisher.test(
  por_site_lineage_corrected,
  simulate.p.value = TRUE,
  B = 99999
)

# 6. Check model rank
por_model_matrix <- model.matrix(
  ~ readcount + Site_unique + Lineage_broad,
  data = por_meta
)

c(
  rank = qr(por_model_matrix)$rank,
  columns = ncol(por_model_matrix)
)



library(vegan)
library(permute)

# Site-restricted marginal PERMANOVA
set.seed(123)

por_adonis_restricted_corrected <- adonis2(
  por_dist_corrected ~
    readcount + Site_unique + Lineage_broad,
  data = por_meta,
  permutations = how(
    nperm = 9999,
    blocks = por_meta$Site_unique
  ),
  by = "margin"
)

por_adonis_restricted_corrected

# Lineage dispersion
set.seed(123)

por_disp_lineage_corrected <- betadisper(
  por_dist_corrected,
  por_meta$Lineage_broad
)

permutest(
  por_disp_lineage_corrected,
  permutations = 9999
)

# Site dispersion
set.seed(123)

por_disp_site_corrected <- betadisper(
  por_dist_corrected,
  por_meta$Site_unique
)

permutest(
  por_disp_site_corrected,
  permutations = 9999
)



por_dbrda_corrected <- dbrda(
  por_dist_corrected ~ Lineage_broad +
    Condition(Site_unique + readcount),
  data = por_meta
)

por_dbrda_corrected

100 * eigenvals(
  por_dbrda_corrected,
  model = "constrained"
) / por_dbrda_corrected$tot.chi




library(vegan)
library(ggplot2)

por_scores_corrected <- as.data.frame(
  scores(por_dbrda_corrected, display = "sites", choices = 1:2)
)

names(por_scores_corrected)


por_scores_corrected$Lineage <- por_meta$Lineage_broad

por_db1 <- 100 * eigenvals(
  por_dbrda_corrected, model = "constrained"
)[1] / por_dbrda_corrected$tot.chi

por_mds1 <- 100 * eigenvals(
  por_dbrda_corrected, model = "unconstrained"
)[1] / por_dbrda_corrected$tot.chi

por_dbrda_figure <- ggplot(
  por_scores_corrected,
  aes(dbRDA1, MDS1, color = Lineage)
) +
  geom_point(size = 7, alpha = 0.85) +
  coord_fixed() +
  theme_classic() +
  labs(
    title = expression(italic("Porites lobata")),
    x = sprintf("dbRDA1 (%.2f%%)", por_db1),
    y = sprintf("MDS1 (%.2f%%)", por_mds1),
    color = "Lineage"
  )

print(por_dbrda_figure)

ggsave(
  "por_dbRDA_corrected.pdf",
  plot = por_dbrda_figure,
  width = 6, height = 5
)

# Unconstrained PCA
por_pca_corrected <- prcomp(por_clr_corrected)

por_pca_scores <- as.data.frame(
  por_pca_corrected$x[, 1:2]
)
por_pca_scores$Lineage <- por_meta$Lineage_broad

por_pca_var <- 100 * por_pca_corrected$sdev^2 /
  sum(por_pca_corrected$sdev^2)

por_pca_figure <- ggplot(
  por_pca_scores,
  aes(PC1, PC2, color = Lineage)
) +
  geom_point(size = 7, alpha = 0.85) +
  coord_fixed() +
  theme_classic() +
  labs(
    title = expression(italic("Porites lobata")),
    x = sprintf("PC1 (%.2f%%)", por_pca_var[1]),
    y = sprintf("PC2 (%.2f%%)", por_pca_var[2]),
    color = "Lineage"
  )

print(por_pca_figure)

ggsave(
  "por_PCA_corrected.pdf",
  plot = por_pca_figure,
  width = 6, height = 5
)






library(MASS)
library(adegenet)

por_dfa_corrected <- MASS::lda(
  por_clr_corrected,
  grouping = por_meta$Lineage_broad
)

por_dfa_accuracy <- mean(
  predict(por_dfa_corrected)$class ==
    por_meta$Lineage_broad
) * 100

por_dfa_accuracy

set.seed(123)

por_cv_corrected <- xvalDapc(
  por_clr_corrected,
  grp = por_meta$Lineage_broad,
  n.pca.max = 20,
  training.set = 0.7,
  result = "groupMean",
  center = TRUE,
  scale = FALSE,
  n.rep = 100,
  xval.plot = TRUE
)

por_cv_corrected$`Number of PCs Achieving Highest Mean Success`
por_cv_corrected$`Mean Successful Assignment by Number of PCs of PCA`




library(adegenet)

por_dapc_corrected <- dapc(
  por_clr_corrected,
  grp = por_meta$Lineage_broad,
  n.pca = 8,
  n.da = 1
)

# Print DAPC
scatter(
  por_dapc_corrected,
  scree.da = FALSE,
  legend = TRUE,
  main = "Porites lobata — DAPC"
)

# Save PDF
pdf("por_DAPC_corrected.pdf", width = 7, height = 6)

scatter(
  por_dapc_corrected,
  scree.da = FALSE,
  legend = TRUE,
  main = "Porites lobata — DAPC"
)

dev.off()



#261104:
library(ggplot2)

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

# Apply identical formatting
por_dbrda_final <- por_dbrda_figure + coral_theme
por_pca_final <- por_pca_figure + coral_theme

# Inspect
print(por_dbrda_final)
print(por_pca_final)

# Export both PDFs
ggsave(
  "por_dbRDA_final.pdf",
  plot = por_dbrda_final,
  width = 5,
  height = 5,
  units = "in"
)

ggsave(
  "por_PCA_final.pdf",
  plot = por_pca_final,
  width = 5,
  height = 5,
  units = "in"
)


library(adegenet)

pdf("por_DAPC_final.pdf", width = 5, height = 5)

scatter(
  por_dapc_corrected,
  scree.da = FALSE
)

dev.off()





















































































