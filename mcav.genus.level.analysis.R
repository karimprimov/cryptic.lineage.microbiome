#loading libraries
install.packages("/Users/karimprimov/Dropbox/UTAustin/2023/rpackages/RDAforest_2.6.9.tar.gz")
install.packages("extendedForest", repos="http://R-Forge.R-project.org")

install.packages("gradientForest", repos="http://R-Forge.R-project.org")



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
library(pcutils)
library(pctax)
library(adegenet)


mcav.genus.bracken.output.df <- read_tsv("/Users/karimprimov/Dropbox/UTAustin/2026/mcav.cryptic.lineages/mcav.flkeys.metagenomics/260415.mcav.genus.level.bracken.taxonomic.report.tsv")
#filtering out contaminant, human, and other unwanted microbe taxa:
microbe.taxa <- unique(mcav.genus.bracken.output.df$name)
write.csv(microbe.taxa, "260415.microbe.taxa.screencheck.csv", row.names = FALSE)

mcav.genus.bracken.output.df$taxonomy_lvl[mcav.genus.bracken.output.df$taxonomy_lvl =='FALSE'] <- 'Genus'
head(mcav.genus.bracken.output.df)
mcav.genus.bracken.output.df


raw.bacteria.list <- c(
  
  # =========================
  # TIER 1 — CORE SYMBIONTS
  # =========================
  
  "Endozoicomonas",
  "Ruegeria",
  "Phaeobacter",
  "Sulfitobacter",
  "Roseovarius",
  "Leisingera",
  "Octadecabacter",
  "Labrenzia",
  "Loktanella",
  "Jannaschia",
  "Hyphomicrobium",
  "Hyphobacterium",
  "Rhodobacter",
  "Paracoccus",
  
  # =========================
  # SELECT TIER 2 (STRICTED)
  # only reef-relevant, commonly reported
  # =========================
  
  "Alteromonas",
  "Pseudoalteromonas",
  "Vibrio",
  "Photobacterium",
  "Marinomonas",
  "Neptunomonas",
  "Idiomarina",
  "Shewanella",
  "Thalassotalea",
  "Colwellia",
  
  "Marinobacter",
  "Halomonas",
  "Cobetia",
  
  "Thalassospira",
  "Alcanivorax",
  
  "Polaribacter",
  "Tenacibaculum",
  "Maribacter",
  "Aquimarina",
  "Cellulophaga",
  "Winogradskyella",
  "Kordia",
  "Formosa",
  "Algibacter",
  "Nonlabens",
  "Lacinutrix",
  "Lutibacter",
  "Zunongwangia",
  "Seonamhaeicola",
  
  "Flavobacterium",
  
  "Sphingomonas",
  "Novosphingobium",
  "Sphingobium",
  "Sphingopyxis",
  "Erythrobacter",
  
  "Paracoccus",
  "Roseomonas",
  
  "Nitrosomonas",
  "Nitrosospira",
  "Nitrospira",
  
  "Thiobacillus",
  "Thiomicrospira",
  "Thiomicrorhabdus",
  "Halothiobacillus",
  
  "Desulfovibrio",
  "Desulfobacter",
  "Desulfosarcina",
  
  "Bdellovibrio",
  "Myxococcus"
)

genus.microbe_pattern <- paste(readLines(textConnection(raw.bacteria.list)), collapse = "|")

mcav.genus.microbe.updated.df <- mcav.genus.bracken.output.df %>%
  transmute(Sample = gsub(".genus.family.output.report", "", Sample),
            Taxon = name,
            readcount = new_est_reads, 
  ) %>%
  pivot_wider(names_from = "Sample", values_from = "readcount", values_fn = sum) %>%
  filter(grepl(genus.microbe_pattern,Taxon)) %>%
  column_to_rownames("Taxon") %>%
  filter(rowSums(!is.na(.)) >= nrow(.) /43) %>% #removing rows found in <50% of samples
  mutate(across(everything(), ~replace(., is.na(.), 1)))

mcav.microbe.df.t <- t(mcav.genus.microbe.updated.df)
mcav.microbe.df.t.df <- data.frame(mcav.microbe.df.t)

mcav.genotype.data.csv <- read.csv("/Users/karimprimov/Dropbox/UTAustin/2026/mcav.cryptic.lineages/mcav.flkeys.metagenomics/260415.mcav.lineage.sample.info.csv", header = TRUE, row.names = 1)
mcav.tobesequenced <- read.csv("/Users/karimprimov/Dropbox/UTAustin/2026/mcav.cryptic.lineages/mcav.flkeys.metagenomics/To.be.Sequenced.MCav.csv", header = TRUE, row.names = 1)

common_ids <- intersect(rownames(mcav.genotype.data.csv), rownames(mcav.microbe.df.t.df))

# 2. Subset and reorder both dataframes to match that list exactly
mcav.genotype.sorted <- mcav.genotype.data.csv[common_ids, ]
mcav.microbe.sorted  <- mcav.microbe.df.t.df[common_ids, ]

# 3. Double-check they match perfectly
identical(rownames(mcav.genotype.sorted), rownames(mcav.microbe.sorted))


#adding readcounts to metadata to condition readcount and batch from the PCA
mcav.genotype.sorted$readcount <- rowSums(mcav.microbe.sorted)

stand.mcav.microbes <- decostand(mcav.microbe.sorted, method = "rclr")
mcav.microbe.matrix <- as.matrix(vegdist(stand.mcav.microbes, method = "euclidean", na.rm=TRUE))

rownames(mcav.genotype.sorted) <- sub("_.*", "", rownames(mcav.genotype.sorted))
# <- subset(mcav.genotype.sorted, select = -depth)

merged.data.metadata.readcount.df <- left_join(mcav.genotype.sorted, mcav.tobesequenced, by = "row.names")


merged.data.metadata.readcount.df <- mcav.genotype.sorted %>%
  rownames_to_column("rowid") %>%
  left_join(rownames_to_column(mcav.tobesequenced, "rowid"), by = "rowid") %>%
  column_to_rownames("rowid")
merged.data.metadata.readcount.df$depth <- as.numeric(merged.data.metadata.readcount.df$Depth.ft.)
merged.data.metadata.readcount.df$diameter <- as.numeric(merged.data.metadata.readcount.df$Diameter.cm.)
merged.data.metadata.readcount.df$height <- as.numeric(merged.data.metadata.readcount.df$Height.cm.)
merged.data.metadata.readcount.df$depth_bin <- cut(merged.data.metadata.readcount.df$depth, breaks = c(0,10,20,30,50,80))
rownames()

#condition batch and readcount!
#condition off readcount and batch and plot with respect to site
rownames(merged.data.metadata.readcount.df) <- rownames(mcav.microbe.matrix)
pca.mcav.depth.mcirobe <- capscale(mcav.microbe.matrix ~ Lineage.x + Condition(readcount), data = merged.data.metadata.readcount.df)

plot(pca.mcav.depth.mcirobe, scaling = 1)





mcav.microbe.scores <- scores(pca.mcav.depth.mcirobe, display = "sites", choices=c(1:6)) 
rownames(mcav.microbe.scores) <- rownames(mcav.microbe.matrix)

merged.fam.mcav.microbe.scores.metadata <- cbind(mcav.microbe.scores, merged.data.metadata.readcount.df, by = "row.names", all.x = TRUE)

ggplot(merged.fam.mcav.microbe.scores.metadata, aes(x = CAP1, y = CAP2, color =as.factor(Lineage.x), label = rownames(merged.fam.mcav.microbe.scores.metadata))) +
  geom_point(size = 12) +
  theme_minimal() + 
  coord_fixed()


#260506:
rownames(merged.fam.mcav.microbe.scores.metadata) <- merged.fam.mcav.microbe.scores.metadata$Row.names

merged.fam.mcav.microbe.scores.metadata.new <- merged.fam.mcav.microbe.scores.metadata[match(rownames(mcav.microbe.matrix), rownames(merged.fam.mcav.microbe.scores.metadata)), ]



pca.mcav.depth.mcirobe <- capscale(mcav.microbe.matrix ~ Lineage + Condition(readcount), data = merged.fam.mcav.microbe.scores.metadata.new)






adonis.site.results <- adonis2(mcav.microbe.matrix ~ readcount + Lineage + Site, 
                               data = merged.fam.mcav.microbe.scores.metadata.new,
                               by = "margin")   # marginal tests are best for independent effects


disp_site <- betadisper(as.dist(mcav.microbe.matrix), merged.fam.mcav.microbe.scores.metadata.new$Site)
permutest(disp_site) 

# For ClusterID
disp_cluster <- betadisper(as.dist(mcav.microbe.matrix), merged.fam.mcav.microbe.scores.metadata.new$Lineage)
permutest(disp_cluster)


table(pocil.metadata.ordered$Site, pocil.metadata.ordered$ClusterID)

# 2. Chi-Squared Test (Statistical test of association)
test <- chisq.test(table(merged.fam.mcav.microbe.scores.metadata.new$Site, merged.fam.mcav.microbe.scores.metadata.new$Lineage))
print(test)


























beta.transect <- betadisper(as.dist(mcav.microbe.matrix), merged.fam.mcav.microbe.scores.metadata$Lineage.x)
anova(beta.transect)

#is there correlation of genetic and microbiome variation
#measure lineage by lineage
#strongest genetic differentiation of host is what differentiates microbiome too
#genetic structure 
#cryptic lineages are cryptic to the microbiome 



# 2. View results
print(adonis.site.results)


adonis.site.lineage.results <- adonis2(mcav.microbe.matrix ~ readcount + Site + Lineage, 
                                       data = mcav.genotype.sorted,
                                       by = "margin")   # marginal tests are best for independent effects

# 2. View results
print(adonis.site.lineage.results)

####################################################################################
####################################################################################

db_site_model <- dbrda(mcav.microbe.matrix ~ Site + Condition(readcount), 
                       data = mcav.genotype.sorted)

# 2. Get the PERMANOVA table (identical results to adonis2)
anova(db_site_model, by = "margin")



db_site_lineage_model <- dbrda(mcav.microbe.matrix ~ Site + Lineage + Condition(readcount), 
                               data = mcav.genotype.sorted)

# 2. Get the PERMANOVA table (identical results to adonis2)
anova(db_site_lineage_model, by = "margin")


############################

adonis.depth.results <- adonis2(mcav.microbe.matrix ~ readcount + depth, 
                               data = mcav.genotype.sorted,
                               by = "margin")   # marginal tests are best for independent effects

# 2. View results
print(adonis.depth.results)


















#

























hc.microbe=hclust(as.dist(mcav.microbe.matrix),method="ave")
plot(hc.microbe)








#
mcav.microbe.df.t.df

prev <- colSums(mcav.microbe.df.t.df > 0)

n_samples <- nrow(mcav.microbe.df.t.df)
cat("Number of Mcav samples:", n_samples, "\n")
cat("Number of ASVs total:", ncol(mcav.microbe.df.t.df), "\n")

asvs_remaining <- sapply(1:n_samples, function(k) {
  sum(prev >= k)
})

curve_df <- data.frame(
  min_individuals_shared = 1:n_samples,
  ASVs_remaining = asvs_remaining
)

ggplot(curve_df,
       aes(x = min_individuals_shared,
           y = ASVs_remaining)) +
  geom_line(size = 1.2) +
  theme_minimal() +
  labs(
    x = "Minimum number of Mcav individuals sharing an ASV",
    y = "Number of ASVs remaining",
    title = "Mcav ASV prevalence (sharedness) curve"
  )



cutoff <- ceiling(0.5 * n_samples)
keep_asvs <- prev >= cutoff
mcav.microbe.df.t.df <- mcav.microbe.df.t.df[, keep_asvs]


#formula based on lineage

mcav.genotype.data.csv <- read.csv("/Users/karimprimov/Dropbox/UTAustin/2026/mcav.cryptic.lineages/mcav.flkeys.metagenomics/mcav.genotype.hybrid.data.csv", header = TRUE, row.names = 1)

mcav.genotype.data.csv2 <- mcav.genotype.data.csv[match(rownames(mcav.microbe.df.t.df), rownames(mcav.genotype.data.csv)), ]



mcav.genus.decostand.hybrid.nohybrid <- decostand(mcav.microbe.df.t.df, method = "rclr")
mcav.pca <- rda(mcav.genus.decostand ~ Lineage, data = mcav.genotype.data.csv2)
plot(mcav.pca)

mcav.genus.microbe.scores <- scores(mcav.pca, display = "sites", choices=c(1:2)) 

merged.mcav.microbe.scores.metadata <- merge(mcav.genus.microbe.scores, mcav.genotype.data.csv2, by = "row.names", all.x = TRUE)

ggplot(merged.mcav.microbe.scores.metadata, aes(x = CAP1, y = CAP2, color = Lineage)) +
  geom_point(size = 12) +
  geom_label_repel(aes(label = Row.names),
                   box.padding   = 0.35, 
                   point.padding = 0.5,
                   segment.color = 'grey50') +
  theme_minimal()




#protest
#pca.mcav.mcirobe 
#&
#mcav.pca
mcav.fam.genus.procrust <-protest(pca.mcav.mcirobe, mcav.pca, permutations = 1000)
mcav.fam.genus.procrust





#remove hybrids and make pca 
mcav.genotype.data.nohybrids <- mcav.genotype.data.csv2[mcav.genotype.data.csv2$Genotype != "Hybrid", ]

mcav.microbe.df.t.df.nohybrids <- mcav.microbe.df.t.df[rownames(mcav.microbe.df.t.df) %in% rownames(mcav.genotype.data.nohybrids), ]


mcav.genus.decostand.nohybrids <- decostand(mcav.microbe.df.t.df.nohybrids, method = "rclr")
vegdist.mcav.nohybrids.matrix <- as.matrix(vegdist(mcav.genus.decostand, method = "euclidean", na.rm = TRUE))

mcav.nohybrids.pca <- rda(mcav.genus.decostand.nohybrids ~ Lineage, data = mcav.genotype.data.nohybrids)
plot(mcav.nohybrids.pca)

mcav.genus.nohybrids.microbe.scores <- scores(mcav.nohybrids.pca, display = "sites", choices=c(1:4)) 

merged.mcav.microbe.scores.metadata.nohybrids <- merge(mcav.genus.nohybrids.microbe.scores, mcav.genotype.data.nohybrids, by = "row.names", all.x = TRUE)

ggplot(merged.mcav.microbe.scores.metadata.nohybrids, aes(x = RDA1, y = RDA2, color = Lineage)) +
  geom_point(aes(shape = as.factor(reef)), size = 12) +
  geom_label_repel(aes(label = Row.names),
                   box.padding   = 0.35, 
                   point.padding = 0.5,
                   segment.color = 'grey50') +
  theme(aspect.ratio = 1) 

#260209:
#RDA1 and RDA2 percent variation
perc_rda1 <- (mcav.nohybrids.pca$CCA$eig[1] / mcav.nohybrids.pca$tot.chi) * 100
#RDA1 
#14.58422 
perc_rda2 <- (mcav.nohybrids.pca$CCA$eig[2] / mcav.nohybrids.pca$tot.chi) * 100
#RDA2 
#8.528175 
# Call: rda(formula = mcav.genus.decostand.nohybrids ~ Lineage, data = mcav.genotype.data.nohybrids)
# 
# Inertia Proportion Rank
# Total         128.7032     1.0000     
# Constrained    39.4521     0.3065    3
# Unconstrained  89.2511     0.6935    7
# Inertia is variance 
# 
# Eigenvalues for constrained axes:
#   RDA1   RDA2   RDA3 
# 18.770 10.976  9.706 
# 
# Eigenvalues for unconstrained axes:
#   PC1    PC2    PC3    PC4    PC5    PC6    PC7 
# 20.698 16.225 14.290 12.018  9.351  8.808  7.861 

shannon.div <- diversity(mcav.microbe.df.t.df.nohybrids, index = "shannon")
simpson.div <- diversity(mcav.microbe.df.t.df.nohybrids, index = "simpson")
richness <- specnumber(mcav.microbe.df.t.df.nohybrids)

div.results <- data.frame(Sample = rownames(mcav.microbe.df.t.df.nohybrids),
                          Shannon = shannon.div,
                          Simpson = simpson.div,
                          Richness = richness)

merge.metadata.div.stats <- merge(mcav.genotype.data.nohybrids, div.results, by = "row.names", all.x = TRUE)
kruskal.test.mcav <- kruskal.test(shannon.div ~ Lineage, data = merge.metadata.div.stats)

ggplot(merge.metadata.div.stats, aes(x = Sample, y = Shannon, fill = Lineage)) +
  geom_boxplot(outlier.shape = NA) +          # Hide boxplot outliers since jitter shows them
  geom_jitter(width = 0.2, alpha = 0.5) +     # Add transparency (alpha) to see overlapping points
  facet_wrap(~Lineage, scales = "free_x") +   # 'free_x' lets each facet show only its samples
  theme_minimal() +
  labs(title = "Shannon Diversity by Sample", 
       y = "Shannon Index", 
       x = "Sample ID") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) # Tilt labels if you have many samples




plot_data <- merge.metadata.div.stats %>%
  pivot_longer(cols = c(Shannon, Simpson, Richness), 
               names_to = "Metric", 
               values_to = "Value")

plot_data$Lineage <- as.factor(plot_data$Lineage)

ggplot(plot_data, aes(x = Lineage, y = Value, fill = Lineage)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6) +
  # Use position_dodge if you add more categories later, 
  # but jitter is great for low N (n=2 or n=3)
  geom_jitter(width = 0.1, size = 2) + 
  facet_wrap(~Metric, scales = "free_y") +
  theme_bw() +
  labs(title = "Alpha Diversity across Lineages",
       x = "Lineage Assignment",
       y = "Score") +
  theme(legend.position = "none") # Lineage is already on the X-axis

top10microbes.nonhybrids <- sort(colSums(mcav.microbe.df.t.df.nohybrids), decreasing = TRUE)[1:10]
top25microbes.nonhybrids <- sort(colSums(mcav.microbe.df.t.df.nohybrids), decreasing = TRUE)[1:25]









#envfit:
mcav.genotype.data.nohybrids <- mcav.genotype.data.csv2[mcav.genotype.data.csv2$Genotype != "Hybrid", ]

mcav.microbe.df.t.df.nohybrids <- mcav.microbe.df.t.df[rownames(mcav.microbe.df.t.df) %in% rownames(mcav.genotype.data.nohybrids), ]


mcav.genus.decostand.nohybrids <- decostand(mcav.microbe.df.t.df.nohybrids, method = "rclr")
vegdist.mcav.nohybrids.matrix <- as.matrix(vegdist(mcav.genus.decostand, method = "euclidean", na.rm = TRUE))

mcav.nohybrids.pca <- rda(mcav.genus.decostand.nohybrids ~ Lineage, data = mcav.genotype.data.nohybrids)
plot(mcav.nohybrids.pca)


X.rda.microbe.environment = as.matrix(mcav.microbe.df.t.df.nohybrids)
X.hel <- decostand(X.rda.microbe.environment, method = "hellinger")
Y.rda.microbe.environment = scale(cdd.metadata.enviro.table)

rda.microbe.environment <- rda(Y.rda.microbe.environment, data = X.hel)
plot(rda.microbe.environment, scaling = 2, display=c("species","bp"))

microbe.fit <- envfit(rda.microbe.environment, X.hel, permutations = 999)
plot(microbe.fit, p.max = 0.01, col = "blue", cex = 0.6)
summary(rda.microbe.environment)








#










#












#260206
#project hybrids!!
#program predict

predicted.hybrid.ordination <- predict(mcav.nohybrids.pca, newdata = mcav.genus.decostand.hybrid.nohybrid, type = "wa", scaling = "sites")
predicted.hybrid.ordination.scores <- scores(predicted.hybrid.ordination, display = "sites", choices=c(1:2)) 
merged.mcav.microbe.scores.metadata.hybrids.and.nohybrids <- merge(predicted.hybrid.ordination.scores, mcav.genotype.data.csv2, by = "row.names", all.x = TRUE)
ggplot(merged.mcav.microbe.scores.metadata.hybrids.and.nohybrids, aes(x = RDA1, y = RDA2, color = Lineage)) +
  geom_point(size = 12) +
  geom_label_repel(aes(label = Row.names),
                   box.padding   = 0.35, 
                   point.padding = 0.5,
                   segment.color = 'grey50') +
  theme(aspect.ratio = 1) 



#adding 4.7 to pca!
mcav.genotype.data.nohybrids.4.7 <- mcav.genotype.data.csv2[mcav.genotype.data.csv2$Lineage %in% c("4", "7", "21", "30", "4.7"), ]
mcav.microbe.df.t.df.nohybrids.4.7 <- mcav.microbe.df.t.df[rownames(mcav.microbe.df.t.df) %in% rownames(mcav.genotype.data.nohybrids.4.7), ]
mcav.genus.decostand.nohybrids.4.7 <- decostand(mcav.microbe.df.t.df.nohybrids.4.7, method = "rclr")

predicted.hybrid.ordination.4.7 <- predict(mcav.nohybrids.pca, newdata = mcav.genus.decostand.nohybrids.4.7, type = "wa", scaling = "sites")
predicted.hybrid.ordination.scores.4.7 <- scores(predicted.hybrid.ordination.4.7, display = "sites", choices=c(1:2)) 
merged.mcav.microbe.scores.metadata.nohybrids.4.7 <- merge(predicted.hybrid.ordination.scores.4.7 , mcav.genotype.data.nohybrids.4.7, by = "row.names", all.x = TRUE)
ggplot(merged.mcav.microbe.scores.metadata.nohybrids.4.7, aes(x = RDA1, y = RDA2, color = Lineage)) +
  geom_point(size = 12) +
  geom_label_repel(aes(label = Row.names),
                   box.padding   = 0.35, 
                   point.padding = 0.5,
                   segment.color = 'grey50') +
  theme(aspect.ratio = 1) 



#adding 4.30 to pca!
mcav.genotype.data.nohybrids.4.7.4.3 <- mcav.genotype.data.csv2[mcav.genotype.data.csv2$Lineage %in% c("4", "7", "21", "30", "4.7", "4.3"), ]
mcav.microbe.df.t.df.nohybrids.4.7.4.3 <- mcav.microbe.df.t.df[rownames(mcav.microbe.df.t.df) %in% rownames(mcav.genotype.data.nohybrids.4.7.4.3), ]
mcav.genus.decostand.nohybrids.4.7.4.3 <- decostand(mcav.microbe.df.t.df.nohybrids.4.7.4.3, method = "rclr")

predicted.hybrid.ordination.4.7.4.3 <- predict(mcav.nohybrids.pca, newdata = mcav.genus.decostand.nohybrids.4.7.4.3, type = "wa", scaling = "sites")
predicted.hybrid.ordination.scores.4.7.4.3 <- scores(predicted.hybrid.ordination.4.7.4.3, display = "sites", choices=c(1:2)) 
merged.mcav.microbe.scores.metadata.nohybrids.4.7.4.3 <- merge(predicted.hybrid.ordination.scores.4.7.4.3, mcav.genotype.data.nohybrids.4.7.4.3, by = "row.names", all.x = TRUE)
ggplot(merged.mcav.microbe.scores.metadata.nohybrids.4.7.4.3, aes(x = RDA1, y = RDA2, color = Lineage)) +
  geom_point(size = 12) +
  geom_label_repel(aes(label = Row.names),
                   box.padding   = 0.35, 
                   point.padding = 0.5,
                   segment.color = 'grey50') +
  theme(aspect.ratio = 1) 

#adding 4.21.30 to pca!
mcav.genotype.data.nohybrids.4.7.4.3.4.21.30 <- mcav.genotype.data.csv2[mcav.genotype.data.csv2$Lineage %in% c("4", "7", "21", "30", "4.7", "4.3", "4.21.30"), ]
mcav.microbe.df.t.df.nohybrids.4.7.4.3.4.21.30 <- mcav.microbe.df.t.df[rownames(mcav.microbe.df.t.df) %in% rownames(mcav.genotype.data.nohybrids.4.7.4.3.4.21.30), ]
mcav.genus.decostand.nohybrids.4.7.4.3.4.21.30 <- decostand(mcav.microbe.df.t.df.nohybrids.4.7.4.3.4.21.30, method = "rclr")
predicted.hybrid.ordination.4.7.4.3.4.21.30 <- predict(mcav.nohybrids.pca, newdata = mcav.genus.decostand.nohybrids.4.7.4.3.4.21.30, type = "wa", scaling = "sites")
predicted.hybrid.ordination.scores.4.7.4.3.4.21.30 <- scores(predicted.hybrid.ordination.4.7.4.3.4.21.30, display = "sites", choices=c(1:2)) 
merged.mcav.microbe.scores.metadata.nohybrids.4.7.4.3.4.21.30 <- merge(predicted.hybrid.ordination.scores.4.7.4.3.4.21.30, mcav.genotype.data.nohybrids.4.7.4.3.4.21.30, by = "row.names", all.x = TRUE)
ggplot(merged.mcav.microbe.scores.metadata.nohybrids.4.7.4.3.4.21.30, aes(x = RDA1, y = RDA2, color = Lineage)) +
  geom_point(size = 12) +
  geom_label_repel(aes(label = Row.names),
                   box.padding   = 0.35, 
                   point.padding = 0.5,
                   segment.color = 'grey50') +
  theme(aspect.ratio = 1) 








#260306

#microbiome venn diagram: 
mc7.venn <- colSums(mcav_clean_df[c(1,2,11),])
mc21.venn <- colSums(mcav_clean_df[3:5,])
mc30.venn <- colSums(mcav_clean_df[c(6,8,9),])
mc4.venn <- colSums(mcav_clean_df[c(7,10),])

mc7.microbes <- names(mc7.venn[mc7.venn>0])
mc21.microbes <- names(mc21.venn[mc21.venn>0])
mc30.microbes <- names(mc30.venn[mc30.venn>0])
mc4.microbes <- names(mc4.venn[mc4.venn>0])

write.csv(mcav_clean_df, "mcav.microbe.df.t.df.nohybrids.csv", row.names = TRUE)
mcav.ven.list <- list(MC7 = mc7.microbes, MC21 = mc21.microbes, MC30 = mc30.microbes, MC4 = mc4.microbes)

mcav.venn_plot <- venn.diagram(
  x = mcav.ven.list,
  category.names = c("MC7", "MC21", "MC30", "MC4"),
  filename = NULL,
  fill = c("skyblue", "pink", "orange", "limegreen"),
  cex = 4,
  cat.cex = 4
)

# Display the plot
grid.draw(mcav.venn_plot)



upset(fromList(mcav.ven.list), 
      order.by = "freq", 
      sets.bar.color = c("skyblue", "pink", "orange", "limegreen"),
      nsets = 4,
      main.bar.color = "gray30",
      matrix.color = "black",
      text.scale = 3)


#260310:
#adonis:
rownames(merged.mcav.microbe.scores.metadata.nohybrids) <- merged.mcav.microbe.scores.metadata.nohybrids$Row.names
mcav.lineage.adonis <- adonis2(vegdist.mcav.nohybrids.matrix ~ Lineage, data=merged.mcav.microbe.scores.metadata.nohybrids)
mcav.lineage.adonis

mcav.depth.adonis <- adonis2(vegdist.mcav.nohybrids.matrix ~ depth, data=merged.mcav.microbe.scores.metadata.nohybrids)
mcav.depth.adonis

mcav.site.adonis <- adonis2(vegdist.mcav.nohybrids.matrix ~ reef, data=merged.mcav.microbe.scores.metadata.nohybrids)
mcav.site.adonis


#260324:
#sloan neutrality test:
contaminants_to_remove <- c(
  "Acinetobacter", "Ralstonia", "Burkholderia", "Cupriavidus", 
  "Stenotrophomonas", "Delftia", "Comamonas", "Variovorax", 
  "Methylobacterium", "Sphingomonas", "Novosphingobium", 
  "Sphingobium", "Sphingopyxis", "Sodalis", "Escherichia", 
  "Salmonella", "Klebsiella", "Enterobacter", "Citrobacter", 
  "Serratia", "Buchnera", "Bradyrhizobium", "Mesorhizobium", 
  "Achromobacter", "Alcaligenes","Rhizobium", "Sinorhizobium", "Paraburkholderia", "Bacteroides", 
  "Parabacteroides", "Haemophilus", "Proteus", "Kosakonia", 
  "Neptunomonas", "Bartonella", "Polynucleobacter")
mcav_clean_df <- mcav.microbe.df.t.df.nohybrids[, !(colnames(mcav.microbe.df.t.df.nohybrids) %in% contaminants_to_remove)]


t.mcav.microbe.df.t.df.nohybrids <- t(mcav_clean_df)
ncm.mcav.nohybrids <- ncm(t.mcav.microbe.df.t.df.nohybrids)
plot(ncm.mcav.nohybrids)
above_taxa <- ncm.mcav.nohybrids$ncm_data[ncm.mcav.nohybrids$ncm_data$group == "Above", ]
below_taxa <- ncm.mcav.nohybrids$ncm_data[ncm.mcav.nohybrids$ncm_data$group == "Below", ]
in_taxa <- ncm.mcav.nohybrids$ncm_data[ncm.mcav.nohybrids$ncm_data$group == "In", ]
head(above_taxa, 30)

mcav_clean_df.read.counts <- data.frame(rowSums(mcav_clean_df))
mcav.genotype.data.nohybrids$readcount <- mcav_clean_df.read.counts$rowSums.mcav_clean_df.
mcav.genotype.data.nohybrids$site <- c("Hasluns", "Hasluns", "Looe", "Looe", "Looe", "Deep", "Offshore", "Deep", "Deep", "Nearshore", "Nearshore")

mcav.genus.decostand.nohybrids <- decostand(mcav_clean_df, method = "rclr")
vegdist.mcav.nohybrids.matrix <- as.matrix(vegdist(mcav.genus.decostand, method = "euclidean", na.rm = TRUE))


mcav.nohybrids.pca <- rda(mcav.genus.decostand.nohybrids ~ Condition(readcount), data = mcav.genotype.data.nohybrids)
plot(mcav.nohybrids.pca)
mcav.nohybrids.pca


mcav.genus.microbe.scores <- scores(mcav.nohybrids.pca, display = "sites", choices=c(1:4)) 

merged.mcav.microbe.scores.metadata <- merge(mcav.genus.microbe.scores, mcav.genotype.data.nohybrids, by = "row.names", all.x = TRUE)

ggplot(merged.mcav.microbe.scores.metadata, aes(x = RDA1, y = RDA2, color = Lineage, shape = site)) +
  geom_point(size = 12) +
  geom_label_repel(aes(label = Row.names),
                   box.padding   = 0.35, 
                   point.padding = 0.5,
                   segment.color = 'grey50') +
  theme_minimal()



#260325:
#trim tara oceans data
#fasterq-dump the amillepora data
#get vcf for first Mcav run to maybe get a PCA
#neutrality test on Amillepora?!
#filter crappy microbes from the mcav and a tenuis and a hyacinthus data
#run PCAs with site AND lineage as constraining variables (readcount as condition)








#Saving RData for the day!
#######################################################################################################################################
save.image("/Users/karimprimov/Dropbox/UTAustin/2026/mcav.cryptic.lineages/mcav.flkeys.metagenomics/261031.mcav.diagnostic.nohybrids.RData")



#saving genus names to text file to remove non-bacterial taxa: 
write.table(
x = rownames(mcav.microbe.df),
file = "260209.microbe.genus.output_column.txt",
sep = "\n",          # Use newline as separator to put each value on a new line
row.names = FALSE,   # Exclude R's default row numbers
col.names = FALSE,   # Exclude the column name/header
quote = FALSE        # Prevent surrounding character data with quotation marks
)




#260402:
mcav.genus.microbe.updated.zeros.df <- mcav.genus.bracken.output.df %>%
  transmute(Sample = gsub("_S.*", "", sample),
            Taxon = name,
            readcount = new_est_reads, 
  ) %>%
  pivot_wider(names_from = "Sample", values_from = "readcount", values_fn = sum) %>%
  filter(grepl(genus.microbe_pattern,Taxon)) %>%
  column_to_rownames("Taxon") %>%
  filter(rowSums(!is.na(.)) >= nrow(.) /111) %>% #removing rows found in <50% of samples
  mutate(across(everything(), ~replace(., is.na(.), 0)))

mcav.genus.microbe.updated.df.zeros.nohybrids <- mcav.genus.microbe.updated.zeros.df[, !grepl("hybrid", names(mcav.genus.microbe.updated.df), ignore.case = TRUE)]

mcav.genus.microbe.updated.df.zeros.nohybrids.t <- t(mcav.genus.microbe.updated.df.zeros.nohybrids)
mcav.genus.microbe.updated.df.zeros.nohybrids.t.df <- data.frame(mcav.genus.microbe.updated.df.zeros.nohybrids.t)




mc7.venn <- colSums(mcav.genus.microbe.updated.df.zeros.nohybrids.t.df[c(1,3,4,7,8,10,13,27,32,33,36,40,41,42,54),])
mc21.venn <- colSums(mcav.genus.microbe.updated.df.zeros.nohybrids.t.df[c(9,16,18,24,29,34,35,44,45,46,47,48),])
mc30.venn <- colSums(mcav.genus.microbe.updated.df.zeros.nohybrids.t.df[c(5,20,21,23,37,43,49,51,52),])
mc4.venn <- colSums(mcav.genus.microbe.updated.df.zeros.nohybrids.t.df[c(2,6,11,12,15,17,25,26,28,30,31,38,39,50,53),])
mc25.venn <- colSums(mcav.genus.microbe.updated.df.zeros.nohybrids.t.df[c(14,19,22),])






mc7.microbes <- names(mc7.venn[mc7.venn>0])
mc21.microbes <- names(mc21.venn[mc21.venn>0])
mc30.microbes <- names(mc30.venn[mc30.venn>0])
mc4.microbes <- names(mc4.venn[mc4.venn>0])
mc25.microbes <- names(mc25.venn[mc25.venn>0])

write.csv(mcav_clean_df, "mcav.microbe.df.t.df.nohybrids.csv", row.names = TRUE)
mcav.ven.list <- list(MC7 = mc7.microbes, MC21 = mc21.microbes, MC30 = mc30.microbes, MC4 = mc4.microbes, MC25 = mc25.microbes)

mcav.venn_plot <- venn.diagram(
  x = mcav.ven.list,
  category.names = c("MC7", "MC21", "MC30", "MC4", "MC25"),
  filename = NULL,
  fill = c("skyblue", "pink", "orange", "limegreen", "purple"),
  cex = 2,
  cat.cex = 2
)

# Display the plot
grid.draw(mcav.venn_plot)



upset(fromList(mcav.ven.list), 
      order.by = "freq", 
      sets.bar.color = c("skyblue", "pink", "orange", "limegreen", "purple"),
      nsets = 5,
      main.bar.color = "gray30",
      matrix.color = "black",
      text.scale = 2)











#260421:
#pca, adonis, dbrda, and lda
#outputs: 1) regular capscale (conditioned by readcount), 
         #2) site capscale (conditioned by readcount), 
         #3) Lineage (conditioned by readcount), 
         #4)depth (conditioned by readcount)



#1) regular
pca.mcav.regular.mcirobe <- capscale(mcav.microbe.matrix ~ Condition(readcount), data = mcav.genotype.sorted)

plot(pca.mcav.regular.mcirobe, scaling = 1)

mcav.regular.microbe.scores <- scores(pca.mcav.regular.mcirobe, display = "sites", choices=c(1:2)) 

mcav.microbe.scores.regular.metadata <- merge(mcav.regular.microbe.scores, mcav.genotype.sorted, by = "row.names", all.x = TRUE)

ggplot(mcav.microbe.scores.regular.metadata, aes(x = MDS1, y = MDS2, color =Lineage, label = mcav.microbe.scores.regular.metadata$`Row.names`)) +
  geom_point(size = 12) +
  geom_label_repel(aes(label = Row.names),
                   box.padding   = 0.35, 
                   point.padding = 0.5,
                   segment.color = 'grey50') +
  theme_minimal() + 
  coord_fixed()


#2) site
pca.mcav.site.mcirobe <- capscale(mcav.microbe.matrix ~ Site + Condition(readcount+ Lineage + depth), data = mcav.genotype.sorted)

plot(pca.mcav.site.mcirobe, scaling = 1)

mcav.site.microbe.scores <- scores(pca.mcav.site.mcirobe, display = "sites", choices=c(1:2)) 

mcav.microbe.scores.site.metadata <- merge(mcav.site.microbe.scores, mcav.genotype.sorted, by = "row.names", all.x = TRUE) %>%
  column_to_rownames(var = "Row.names")

mcav.microbe.scores.site.metadata.new <- mcav.microbe.scores.site.metadata[match(rownames(mcav.microbe.matrix), rownames(mcav.microbe.scores.site.metadata)), ]


ggplot(mcav.microbe.scores.site.metadata.new, aes(x = CAP1, y = CAP2, color =Site, label = rownames(mcav.microbe.scores.site.metadata.new))) +
  geom_point(size = 12) +
  geom_label_repel(aes(label = rownames(mcav.microbe.scores.site.metadata.new)),
                   box.padding   = 0.35, 
                   point.padding = 0.5,
                   segment.color = 'grey50') +
  theme_minimal() + 
  coord_fixed()

#3

pca.mcav.lineage.mcirobe <- capscale(mcav.microbe.matrix ~ Lineage + Condition(readcount + Site + depth), data = mcav.genotype.sorted)

plot(pca.mcav.lineage.mcirobe, scaling = 1)

mcav.lineage.microbe.scores <- scores(pca.mcav.lineage.mcirobe, display = "sites", choices=c(1:2)) 

mcav.microbe.scores.lineage.metadata <- merge(mcav.lineage.microbe.scores, mcav.genotype.sorted, by = "row.names", all.x = TRUE)

ggplot(mcav.microbe.scores.lineage.metadata, aes(x = CAP1, y = MDS1, color =as.factor(Lineage), label = `Row.names`)) +
  geom_point(size = 12) +
  geom_label_repel(aes(label = Row.names),
                   box.padding   = 0.35, 
                   point.padding = 0.5,
                   segment.color = 'grey50') +
  theme_minimal() + 
  coord_fixed()

#4
pca.mcav.depth.mcirobe <- capscale(mcav.microbe.matrix ~ depth + Condition(readcount + Site + Lineage), data = mcav.genotype.sorted)

plot(pca.mcav.depth.mcirobe, scaling = 1)

plot(pca.mcav.depth.mcirobe$CA$eig/pca.mcav.depth.mcirobe$CCA$eig, )

mcav.depth.microbe.scores <- scores(pca.mcav.depth.mcirobe, display = "sites", choices=c(1:2)) 

mcav.microbe.scores.depth.metadata <- merge(mcav.depth.microbe.scores, mcav.genotype.sorted, by = "row.names", all.x = TRUE) %>%
  column_to_rownames(var = "Row.names")

mcav.microbe.scores.depth.metadata.new <- mcav.microbe.scores.depth.metadata[match(rownames(mcav.microbe.matrix), rownames(mcav.microbe.scores.depth.metadata)), ]

ggplot(mcav.microbe.scores.depth.metadata.new, aes(x = CAP1, y = MDS1, color =depth, label = rownames(mcav.microbe.scores.depth.metadata.new))) +
  geom_point(size = 12) +
  geom_label_repel(aes(label = rownames(mcav.microbe.scores.depth.metadata.new)),
                   box.padding   = 0.35, 
                   point.padding = 0.5,
                   segment.color = 'grey50') +
  theme_minimal() + 
  coord_fixed()



#discriminant function analysis:
#site:
dfa.model <- lda(Site + depth + Lineage + Condition(readcount) ~ mcav.microbe.matrix, data = mcav.microbe.scores.depth.metadata.new)


rda.model <- rda(mcav.microbe.sorted ~ Site + depth + Lineage + Condition(readcount), 
                 data = mcav.microbe.scores.depth.metadata.new)
anova(rda.model, by = "margin")




#DAPC:

# 1. DAPC by Site
dapc_site <- dapc(mcav.microbe.sorted, grp = mcav.microbe.scores.depth.metadata.new$Site, n.pca = 15, n.da = 2)

# 2. DAPC by Lineage
dapc_lineage <- dapc(mcav.microbe.sorted, grp = mcav.microbe.scores.depth.metadata.new$Lineage, n.pca = 15, n.da = 2)

# 3. DAPC by Depth (if categorical)
dapc_depth <- dapc(mcav.microbe.sorted, grp = as.factor(mcav.microbe.scores.depth.metadata.new$depth), n.pca = 15, n.da = 1)

# Visualize which microbes drive the Site separation
loadingplot(dapc_site$var.contr, threshold = 0.05)
loadingplot(dapc_lineage$var.contr, threshold = 0.05)
loadingplot(dapc_depth_binned$var.contr, threshold = 0.05)


# Site Plot (Clear clusters usually appear here)
scatter(dapc_site, bg="white", scree.da=FALSE, legend=TRUE, title="DAPC: Site")

# Lineage Plot
scatter(dapc_lineage, bg="white", scree.da=FALSE, legend=TRUE, title="DAPC: Lineage")
summary(dapc_lineage)
compoplot(dapc_lineage, 
          show.lab = TRUE,           # Forces the sample names to show
          csub = 0.8,                # Size of the labels
          posi = "none",         # Move the legend to the right
          main = "Membership Probability per Sample")

# Depth Plot (Note: This will look busy with 16 groups)
scatter(dapc_depth, bg="white", scree.da=FALSE, legend=TRUE, title="DAPC: Depth")



depth_cat <- cut(mcav.microbe.scores.depth.metadata.new$depth, 
                 breaks = c(0, 19, 26, Inf), 
                 labels = c("Shallow", "Mid", "Deep"),
                 right = TRUE) # This makes 20 and 25 the upper bounds of their bins


# 2. Re-run the DAPC with these categories
dapc_depth_binned <- dapc(mcav.microbe.sorted, 
                          grp = depth_cat, 
                          n.pca = 15, 
                          n.da = 2)

# 3. Now the scatter plot will work perfectly
scatter(dapc_depth_binned, 
        bg = "white", 
        scree.da = FALSE, 
        legend = TRUE, 
        title = "DAPC: Depth (Binned)")

summary(dapc_depth_binned)$assign.per.pop
# The "Compoplot" (very popular for advisors)
compoplot(dapc_depth_binned, 
            show.lab = TRUE,           # Forces the sample names to show
            csub = 0.8,                # Size of the labels
            posi = "none",         # Move the legend to the right
            main = "Membership Probability per Sample")



#
adonis2(mcav.microbe.sorted ~ depth, data = mcav.microbe.scores.depth.metadata.new)
adonis2(mcav.microbe.sorted ~ Site, data = mcav.microbe.scores.depth.metadata.new)
adonis2(mcav.microbe.sorted ~ Lineage, data = mcav.microbe.scores.depth.metadata.new)



dapc_depth_clean <- dapc(stand.mcav.microbes, grp = depth_cat, n.pca = 15, n.da = 2)
scatter(dapc_depth_clean, 
        bg = "white", 
        scree.da = FALSE, 
        legend = TRUE, 
        title = "DAPC: Depth (Binned)")

compoplot(dapc_depth_clean, 
          show.lab = TRUE,           # Forces the sample names to show
          csub = 0.8,                # Size of the labels
          posi = "none",         # Move the legend to the right
          main = "Membership Probability per Sample")





#


summary(dapc_site)
summary(dapc_lineage)
summary(dapc_depth_binned)


















#
scores_data <- as.data.frame(scores(rda.model, display = "wa"))

# 2. Add your metadata back to these coordinates
# Make sure the row order matches!
scores_data$depth <- mcav.microbe.scores.depth.metadata.new$depth
scores_data$Site  <- mcav.microbe.scores.depth.metadata.new$Site

# 3. Create the plot
ggplot(scores_data, aes(x = RDA1, y = RDA2, color = depth)) +
  geom_point(size = 7) +
  theme_minimal() +
  labs(title = "Constrained Ordination (RDA)",
       x = "RDA1", y = "RDA2") +
  scale_color_viridis_c()














#

#260925:
pca.mcav.mcirobe.cond.readcount <- capscale(mcav.microbe.matrix ~ Condition(readcount), data = merged.data.metadata.readcount.df)

plot(pca.mcav.depth.mcirobe, scaling = 1)




mcav.microbe.scores <- scores(pca.mcav.depth.mcirobe, display = "sites", choices=c(1:6)) 
rownames(mcav.microbe.scores) <- sub("_.*", "", rownames(mcav.microbe.scores))

merged.fam.mcav.microbe.scores.metadata <- cbind(mcav.microbe.scores, merged.data.metadata.readcount.df, by = "row.names", all.x = TRUE)

ggplot(merged.fam.mcav.microbe.scores.metadata, aes(x = CAP1, y = MDS1, color =depth, label = rownames(merged.fam.mcav.microbe.scores.metadata))) +
  geom_point(size = 12) +
  theme_minimal() + 
  coord_fixed()












#




















#










pca.mcav.depth.mcirobe <- capscale(mcav.microbe.matrix ~ depth + Condition(readcount), data = mcav.genotype.sorted)

plot(pca.mcav.depth.mcirobe, scaling = 1)




mcav.microbe.scores <- scores(pca.mcav.depth.mcirobe, display = "sites", choices=c(1:2)) 

merged.fam.mcav.microbe.scores.metadata <- merge(mcav.microbe.scores, mcav.genotype.sorted, by = "row.names", all.x = TRUE)

ggplot(merged.fam.mcav.microbe.scores.metadata, aes(x = CAP1, y = MDS1, color =depth, label = rownames(merged.fam.mcav.microbe.scores.metadata))) +
  geom_point(size = 12) +
  geom_label_repel(aes(label = Row.names),
                   box.padding   = 0.35, 
                   point.padding = 0.5,
                   segment.color = 'grey50') +
  theme_minimal() + 
  coord_fixed()





adonis.site.results <- adonis2(mcav.microbe.matrix ~ readcount + Site, 
                               data = mcav.genotype.sorted,
                               by = "margin")   # marginal tests are best for independent effects

# 2. View results
print(adonis.site.results)


adonis.site.lineage.results <- adonis2(mcav.microbe.matrix ~ readcount + Site + Lineage, 
                                       data = mcav.genotype.sorted,
                                       by = "margin")   # marginal tests are best for independent effects

# 2. View results
print(adonis.site.lineage.results)

####################################################################################
####################################################################################

db_site_model <- dbrda(mcav.microbe.matrix ~ Site + Condition(readcount), 
                       data = mcav.genotype.sorted)

# 2. Get the PERMANOVA table (identical results to adonis2)
anova(db_site_model, by = "margin")



db_site_lineage_model <- dbrda(mcav.microbe.matrix ~ Site + Lineage + Condition(readcount), 
                               data = mcav.genotype.sorted)

# 2. Get the PERMANOVA table (identical results to adonis2)
anova(db_site_lineage_model, by = "margin")


############################

adonis.depth.results <- adonis2(mcav.microbe.matrix ~ readcount + depth, 
                                data = mcav.genotype.sorted,
                                by = "margin")   # marginal tests are best for independent effects

# 2. View results
print(adonis.depth.results)






#260929:
mcav_meta <- mcav.genotype.data.nohybrids

mcav_dist <- as.dist(
  vegdist.mcav.nohybrids.matrix
)

# Check alignment
nrow(mcav_meta)
attr(mcav_dist, "Size")

identical(
  rownames(mcav_meta),
  attr(mcav_dist, "Labels")
)

# Sample sizes
table(mcav_meta$Lineage)
table(mcav_meta$site)

# Depth
summary(mcav_meta$depth)

# Missingness
colSums(
  is.na(
    mcav_meta[, c(
      "Lineage",
      "site",
      "depth",
      "readcount"
    )]
  )
)

# Site × lineage
table(
  mcav_meta$site,
  mcav_meta$Lineage
)

# Site × depth
aggregate(
  depth ~ site,
  data = mcav_meta,
  FUN = function(x)
    c(
      n = length(x),
      min = min(x),
      mean = mean(x),
      max = max(x)
    )
)




#####################################################################################################################################
#####################################################################################################################################
#####################################################################################################################################
library(vegan)
library(lsr)

# ==========================================
# 1. SITE-ONLY PERMANOVA
# ==========================================

set.seed(123)

mcav_adonis_site <- adonis2(
  mcav_dist ~ readcount + Site,
  data = mcav_meta,
  permutations = 9999,
  by = "margin"
)

mcav_adonis_site


# ==========================================
# 2. LINEAGE-ONLY PERMANOVA
# ==========================================

set.seed(123)

mcav_adonis_lineage <- adonis2(
  mcav_dist ~ readcount + Lineage,
  data = mcav_meta,
  permutations = 9999,
  by = "margin"
)

mcav_adonis_lineage


# ==========================================
# 3. JOINT PERMANOVA
# ==========================================

set.seed(123)

mcav_adonis_both <- adonis2(
  mcav_dist ~ readcount + Site + Lineage,
  data = mcav_meta,
  permutations = 9999,
  by = "margin"
)

mcav_adonis_both


# ==========================================
# 4. SITE × LINEAGE ASSOCIATION
# ==========================================

mcav_site_lineage_table <- table(
  mcav_meta$Site,
  mcav_meta$Lineage
)

mcav_site_lineage_table

mcav_chisq <- chisq.test(
  mcav_site_lineage_table
)

mcav_chisq

library(lsr)
mcav_cramers <- cramersV(mcav_site_lineage_table)

mcav_cramers

set.seed(123)
fisher.test(
  mcav_site_lineage_table,
  simulate.p.value = TRUE,
  B = 99999
)


# ==========================================
# 5. PERMDISP — SITE
# ==========================================

set.seed(123)

mcav_disp_site <- betadisper(
  mcav_dist,
  mcav_meta$Site
)

permutest(
  mcav_disp_site,
  permutations = 9999
)


# ==========================================
# 6. PERMDISP — LINEAGE
# ==========================================

set.seed(123)

mcav_disp_lineage <- betadisper(
  mcav_dist,
  mcav_meta$Lineage
)

permutest(
  mcav_disp_lineage,
  permutations = 9999
)


# ==========================================
# 7. FULL dbRDA
# ==========================================

mcav_dbrda_full <- dbrda(
  mcav_dist ~ readcount + Site + Lineage,
  data = mcav_meta
)

mcav_dbrda_full


# Marginal tests
set.seed(123)

anova(
  mcav_dbrda_full,
  by = "margin",
  permutations = 9999
)


library(ggplot2)

mcav_dbrda_plot <- dbrda(
  mcav_dist ~ Site + Lineage + Condition(readcount),
  data = mcav_meta
)

mcav_dbrda_plot

# Sample scores
mcav_scores <- as.data.frame(
  scores(
    mcav_dbrda_plot,
    display = "sites",
    choices = 1:2
  )
)

mcav_scores$SampleID_plot <- rownames(mcav_scores)
mcav_meta$SampleID_plot <- rownames(mcav_meta)

mcav_plot_data <- merge(
  mcav_scores,
  mcav_meta,
  by = "SampleID_plot",
  all.x = TRUE
)

# % constrained variation represented by each axis
mcav_eig <- eigenvals(
  mcav_dbrda_plot,
  model = "constrained"
)

mcav_db1 <- 100 * mcav_eig[1] / sum(mcav_eig)
mcav_db2 <- 100 * mcav_eig[2] / sum(mcav_eig)

mcav_db1
mcav_db2

# Plot by lineage
ggplot(
  mcav_plot_data,
  aes(
    x = dbRDA1,
    y = dbRDA2,
    color = Lineage
  )
) +
  geom_point(
    size = 12,
    alpha = 0.85
  ) +
  theme_minimal() +
  coord_fixed() +
  labs(
    title = expression(italic("Montastraea cavernosa")),
    x = paste0(
      "dbRDA1 (",
      round(mcav_db1, 1),
      "%)"
    ),
    y = paste0(
      "dbRDA2 (",
      round(mcav_db2, 1),
      "%)"
    ),
    color = "Lineage"
  )










library(vegan)

# ============================================================
# MCAV — WITHIN-SITE LINEAGE TESTS
# ============================================================


# ------------------------------------------------------------
# 1. CARLY'S — LINEAGE 4 vs 7
# ------------------------------------------------------------

carly_keep <- mcav_meta$Site == "Carly's"

carly_meta <- droplevels(
  mcav_meta[carly_keep, ]
)

carly_dist <- as.dist(
  as.matrix(mcav_dist)[carly_keep, carly_keep]
)

# Check design
nrow(carly_meta)
table(carly_meta$Lineage)

identical(
  rownames(carly_meta),
  attr(carly_dist, "Labels")
)


# PERMANOVA
set.seed(123)

mcav_carly_adonis <- adonis2(
  carly_dist ~ readcount + Lineage,
  data = carly_meta,
  permutations = 9999,
  by = "margin"
)

mcav_carly_adonis


# PERMDISP
set.seed(123)

mcav_carly_disp <- betadisper(
  carly_dist,
  carly_meta$Lineage
)

permutest(
  mcav_carly_disp,
  permutations = 9999
)



# ------------------------------------------------------------
# 2. HASLUNS — LINEAGE 4 vs 7
# ------------------------------------------------------------

hasluns_keep <- mcav_meta$Site == "Hasluns"

hasluns_meta <- droplevels(
  mcav_meta[hasluns_keep, ]
)

hasluns_dist <- as.dist(
  as.matrix(mcav_dist)[hasluns_keep, hasluns_keep]
)

# Check design
nrow(hasluns_meta)
table(hasluns_meta$Lineage)

identical(
  rownames(hasluns_meta),
  attr(hasluns_dist, "Labels")
)


# PERMANOVA
set.seed(123)

mcav_hasluns_adonis <- adonis2(
  hasluns_dist ~ readcount + Lineage,
  data = hasluns_meta,
  permutations = 9999,
  by = "margin"
)

mcav_hasluns_adonis


# PERMDISP
set.seed(123)

mcav_hasluns_disp <- betadisper(
  hasluns_dist,
  hasluns_meta$Lineage
)

permutest(
  mcav_hasluns_disp,
  permutations = 9999
)



# ------------------------------------------------------------
# 3. LOOE — LINEAGES 21 vs 25 vs 30
# ------------------------------------------------------------

looe_keep <- mcav_meta$Site == "Looe"

looe_meta <- droplevels(
  mcav_meta[looe_keep, ]
)

looe_dist <- as.dist(
  as.matrix(mcav_dist)[looe_keep, looe_keep]
)

# Check design
nrow(looe_meta)
table(looe_meta$Lineage)

identical(
  rownames(looe_meta),
  attr(looe_dist, "Labels")
)


# PERMANOVA
set.seed(123)

mcav_looe_adonis <- adonis2(
  looe_dist ~ readcount + Lineage,
  data = looe_meta,
  permutations = 9999,
  by = "margin"
)

mcav_looe_adonis


# PERMDISP
set.seed(123)

mcav_looe_disp <- betadisper(
  looe_dist,
  looe_meta$Lineage
)

permutest(
  mcav_looe_disp,
  permutations = 9999
)














































































































#261031:

# Align abundance data with metadata
mcav_counts <- as.matrix(
  mcav.microbe.df.t.df[rownames(mcav_meta), , drop = FALSE]
)

# Verify alignment
identical(rownames(mcav_counts), rownames(mcav_meta))

# Check counts
dim(mcav_counts)
anyNA(mcav_counts)
range(mcav_counts)




library(vegan)

# Ordinary CLR transformation
mcav_clr <- decostand(
  mcav_counts,
  method = "clr",
  pseudocount = 0.01
)

# Aitchison distances
mcav_dist_clr <- dist(mcav_clr, method = "euclidean")

# Verification
dim(mcav_clr)
anyNA(mcav_clr)
attr(mcav_dist_clr, "Size")



set.seed(123)

mcav_adonis_clr <- adonis2(
  mcav_dist_clr ~ readcount + Site + Lineage,
  data = mcav_meta,
  permutations = 9999,
  by = "margin"
)

mcav_adonis_clr


library(vegan)

set.seed(123)
mcav_disp_lineage_clr <- betadisper(
  mcav_dist_clr,
  factor(mcav_meta$Lineage)
)

permutest(mcav_disp_lineage_clr, permutations = 9999)

set.seed(123)
mcav_disp_site_clr <- betadisper(
  mcav_dist_clr,
  factor(mcav_meta$Site)
)

permutest(mcav_disp_site_clr, permutations = 9999)


library(permute)

set.seed(123)

mcav_perm_site <- how(
  nperm = 9999,
  blocks = factor(mcav_meta$Site)
)

mcav_adonis_restricted <- adonis2(
  mcav_dist_clr ~ readcount + Site + Lineage,
  data = mcav_meta,
  permutations = mcav_perm_site,
  by = "margin"
)

mcav_adonis_restricted

mcav_dbrda_clr <- dbrda(
  mcav_dist_clr ~ Lineage +
    Condition(Site + readcount),
  data = mcav_meta
)

mcav_dbrda_clr

# Percentage of total inertia
mcav_axes <- eigenvals(
  mcav_dbrda_clr,
  model = "constrained"
)

100 * mcav_axes / mcav_dbrda_clr$tot.chi

library(ggplot2)

mcav_scores_clr <- as.data.frame(
  scores(mcav_dbrda_clr, display = "sites",
         choices = 1:2)
)

mcav_scores_clr$Lineage <- factor(mcav_meta$Lineage)
mcav_scores_clr$Site <- factor(mcav_meta$Site)

ggplot(mcav_scores_clr,
       aes(x = dbRDA1, y = dbRDA2, color = Lineage)) +
  geom_point(size = 3, alpha = 0.85) +
  coord_fixed() +
  theme_classic() +
  labs(
    title = expression(italic("Montastraea cavernosa")),
    x = "dbRDA1 (2.52%)",
    y = "dbRDA2 (2.21%)"
  )


library(MASS)

mcav_dfa <- lda(
  x = mcav_clr,
  grouping = factor(mcav_meta$Lineage)
)

mcav_dfa_pred <- predict(mcav_dfa)$class

mean(mcav_dfa_pred == factor(mcav_meta$Lineage)) * 100

mcav_dfa

library(adegenet)

set.seed(123)

mcav_xval <- xvalDapc(
  mcav_clr,
  grp = factor(mcav_meta$Lineage),
  n.pca.max = 20,
  training.set = 0.7,
  result = "groupMean",
  center = TRUE,
  scale = FALSE,
  n.rep = 100,
  xval.plot = TRUE
)

mcav_xval$`Number of PCs Achieving Highest Mean Success`
mcav_xval$`Highest Mean Success`






#261031 re-start
# Match microbiome IDs to genotype IDs
mcav_ids <- rownames(mcav.microbe.df.t.df)

mcav_match <- match(
  sub("_S.*$", "", mcav_ids),
  rownames(mcav.genotype.sorted)
)

# Check matching
sum(is.na(mcav_match))
length(mcav_match)

mcav_meta <- mcav.genotype.sorted[mcav_match, , drop = FALSE]
rownames(mcav_meta) <- mcav_ids

# Verify
dim(mcav_meta)
identical(rownames(mcav_meta), mcav_ids)
table(mcav_meta$Lineage)
table(mcav_meta$Site)


mcav_counts_corrected <- mcav.genus.bracken.output.df %>%
  mutate(Sample = sub(
    "\\.genus\\.family\\.output\\.report$",
    "", Sample
  )) %>%
  filter(
    Sample %in% mcav_ids,
    name %in% colnames(mcav.microbe.df.t.df)
  ) %>%
  group_by(Sample, name) %>%
  summarise(readcount = sum(new_est_reads), .groups = "drop") %>%
  complete(
    Sample = mcav_ids,
    name = colnames(mcav.microbe.df.t.df),
    fill = list(readcount = 0)
  ) %>%
  pivot_wider(
    names_from = name,
    values_from = readcount
  ) %>%
  column_to_rownames("Sample") %>%
  as.matrix()

mcav_counts_corrected <- mcav_counts_corrected[
  mcav_ids,
  colnames(mcav.microbe.df.t.df),
  drop = FALSE
]

colSums(mcav_counts_corrected)
dim(mcav_counts_corrected)
anyNA(mcav_counts_corrected)
sum(mcav_counts_corrected == 0)
sum(mcav_counts_corrected == 1)

identical(
  rownames(mcav_counts_corrected),
  rownames(mcav_meta)
)

library(vegan)

# Update read counts
mcav_meta$readcount <- rowSums(mcav_counts_corrected)

# Ordinary CLR transformation
mcav_clr_corrected <- decostand(
  mcav_counts_corrected,
  method = "clr",
  pseudocount = 0.01
)

# Aitchison distance
mcav_dist_corrected <- dist(
  mcav_clr_corrected,
  method = "euclidean"
)

# Verify
dim(mcav_clr_corrected)
anyNA(mcav_clr_corrected)
attr(mcav_dist_corrected, "Size")

set.seed(123)

mcav_adonis_corrected <- adonis2(
  mcav_dist_corrected ~ readcount + Site + Lineage,
  data = mcav_meta,
  permutations = 9999,
  by = "margin"
)

mcav_adonis_corrected


# Ensure categorical predictors are factors
mcav_meta$Lineage <- factor(mcav_meta$Lineage)
mcav_meta$Site <- factor(mcav_meta$Site)

# Verify
str(mcav_meta[, c("Lineage", "Site", "readcount")])

set.seed(123)

mcav_adonis_corrected <- vegan::adonis2(
  mcav_dist_corrected ~ readcount + Site + Lineage,
  data = mcav_meta,
  permutations = 9999,
  by = "margin"
)

mcav_adonis_corrected


library(vegan)
library(permute)

set.seed(123)
mcav_disp_lineage_corrected <- betadisper(
  mcav_dist_corrected, mcav_meta$Lineage
)

mcav_disp_site_corrected <- betadisper(
  mcav_dist_corrected, mcav_meta$Site
)

permutest(mcav_disp_lineage_corrected, permutations = 9999)
permutest(mcav_disp_site_corrected, permutations = 9999)

# Site-restricted PERMANOVA
set.seed(123)

mcav_adonis_restricted_corrected <- adonis2(
  mcav_dist_corrected ~ readcount + Site + Lineage,
  data = mcav_meta,
  permutations = how(
    nperm = 9999,
    blocks = mcav_meta$Site
  ),
  by = "margin"
)

mcav_adonis_restricted_corrected




mcav_dbrda_corrected <- vegan::dbrda(
  mcav_dist_corrected ~ Lineage +
    Condition(Site + readcount),
  data = mcav_meta
)

mcav_dbrda_corrected

# Axis percentages relative to total inertia
100 * vegan::eigenvals(
  mcav_dbrda_corrected,
  model = "constrained"
) / mcav_dbrda_corrected$tot.chi



library(ggplot2)

mcav_scores_corrected <- as.data.frame(
  scores(mcav_dbrda_corrected,
         display = "sites",
         choices = 1:2)
)

mcav_scores_corrected$Lineage <- mcav_meta$Lineage

mcav_plot_corrected <- ggplot(
  mcav_scores_corrected,
  aes(dbRDA1, dbRDA2, color = Lineage)
) +
  geom_point(size = 7, alpha = 0.85) +
  coord_fixed() +
  theme_classic() +
  labs(
    title = expression(italic("Montastraea cavernosa")),
    x = "dbRDA1 (2.54%)",
    y = "dbRDA2 (2.29%)",
    color = "Lineage"
  )

mcav_plot_corrected

library(MASS)
library(adegenet)

# DFA training classification
mcav_dfa_corrected <- MASS::lda(
  mcav_clr_corrected,
  grouping = mcav_meta$Lineage
)

mcav_dfa_accuracy <- mean(
  predict(mcav_dfa_corrected)$class == mcav_meta$Lineage
) * 100

mcav_dfa_accuracy

# DAPC cross-validation
set.seed(123)

mcav_cv_corrected <- adegenet::xvalDapc(
  mcav_clr_corrected,
  grp = mcav_meta$Lineage,
  n.pca.max = 20,
  training.set = 0.7,
  result = "groupMean",
  center = TRUE,
  scale = FALSE,
  n.rep = 100,
  xval.plot = TRUE
)

mcav_cv_corrected$`Number of PCs Achieving Highest Mean Success`

mcav_cv_corrected$`Mean Successful Assignment by Number of PCs of PCA`


library(adegenet)

mcav_dapc_corrected <- dapc(
  mcav_clr_corrected,
  grp = mcav_meta$Lineage,
  n.pca = 10,
  n.da = 4
)

scatter(
  mcav_dapc_corrected,
  scree.da = FALSE,
  legend = TRUE,
  main = "Montastraea cavernosa — DAPC"
)


mcav_pca_corrected <- prcomp(mcav_clr_corrected)

mcav_pca_data <- as.data.frame(mcav_pca_corrected$x[, 1:2])
mcav_pca_data$Lineage <- mcav_meta$Lineage

mcav_pca_var <- 100 *
  mcav_pca_corrected$sdev^2 /
  sum(mcav_pca_corrected$sdev^2)

mcav_pca_figure <- ggplot(mcav_pca_data, aes(PC1, PC2, color = Lineage)) +
  geom_point(size = 3, alpha = 0.85) +
  coord_fixed() +
  theme_classic() +
  labs(
    title = expression(italic("Montastraea cavernosa")),
    x = paste0("PC1 (", round(mcav_pca_var[1], 2), "%)"),
    y = paste0("PC2 (", round(mcav_pca_var[2], 2), "%)")
  )

ggsave("mcav_dbRDA_corrected.pdf",
       mcav_plot_corrected, width = 6, height = 5)

ggsave("mcav_PCA_corrected.pdf",
       width = 6, height = 5)

pdf("mcav_DAPC_corrected.pdf", width = 6, height = 5)
scatter(mcav_dapc_corrected,
        scree.da = FALSE, legend = TRUE)
dev.off()



#261104:
ggsave(
  "mcav_PCA_corrected.pdf",
  plot = mcav_pca_figure,
  width = 5,
  height = 5,
  units = "in"
)



#standardized dbRDA plot for Montastraea cavernosa

library(ggplot2)
library(vegan)

mcav_scores_final <- as.data.frame(
  scores(mcav_dbrda_corrected, display = "sites", choices = 1:2)
)

#add lineage labels
mcav_scores_final$Lineage <- factor(mcav_meta$Lineage)

#calculate percentage of total variation
mcav_eig <- eigenvals(
  mcav_dbrda_corrected,
  model = "constrained"
)

mcav_db1 <- 100 * mcav_eig[1] / mcav_dbrda_corrected$tot.chi
mcav_db2 <- 100 * mcav_eig[2] / mcav_dbrda_corrected$tot.chi

#plot
mcav_dbrda_final <- ggplot(
  mcav_scores_final,
  aes(x = dbRDA1, y = dbRDA2, color = Lineage)
) +
  geom_point(size = 7, alpha = 0.85) +
  labs(
    title = "Montastraea cavernosa",
    x = paste0("dbRDA1 (", round(mcav_db1, 2), "%)"),
    y = paste0("dbRDA2 (", round(mcav_db2, 2), "%)")
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

print(mcav_dbrda_final)

#export at the same dimensions as the other species
ggsave(
  "mcav_dbRDA_final.pdf",
  plot = mcav_dbrda_final,
  width = 5,
  height = 5,
  units = "in"
)



library(dplyr)
library(tidyr)

mcav_site_pies <- mcav_meta %>%
  as.data.frame() %>%
  count(Site, Lineage) %>%
  pivot_wider(
    names_from = Lineage,
    values_from = n,
    values_fill = 0,
    names_prefix = "L"
  )

print(mcav_site_pies)




















library(dplyr)
library(tidyr)

mcav_site_pies <- mcav_meta %>%
  as.data.frame() %>%
  count(Site, Lineage) %>%
  pivot_wider(
    names_from = Lineage,
    values_from = n,
    values_fill = 0,
    names_prefix = "L"
  )

print(mcav_site_pies)


mcav_coords <- data.frame(
  Site = c("Carly's", "Hasluns", "Looe"),
  Longitude = c(-81.4294, -81.4376, -81.4141),
  Latitude = c(24.6071, 24.5532, 24.5410)
)

mcav_site_pies <- left_join(
  mcav_site_pies,
  mcav_coords,
  by = "Site"
)

mcav_site_pies



mcav_site_pies <- mcav_site_pies %>%
  mutate(
    pie_lon = c(-81.455, -81.455, -81.395),
    pie_lat = c(24.620, 24.555, 24.540)
  )

library(ggplot2)
library(scatterpie)
library(maps)

florida_map <- map_data("state") %>%
  dplyr::filter(region == "florida")


mcav_map <- ggplot() +
  
  geom_polygon(
    data = florida_map,
    aes(x = long, y = lat, group = group),
    fill = "grey90",
    color = "grey40",
    linewidth = 0.3
  ) +
  
  geom_segment(
    data = mcav_site_pies,
    aes(
      x = Longitude,
      y = Latitude,
      xend = pie_lon,
      yend = pie_lat
    ),
    color = "grey30",
    linewidth = 0.4
  ) +
  
  geom_point(
    data = mcav_site_pies,
    aes(x = Longitude, y = Latitude),
    color = "black",
    size = 2
  ) +
  
  geom_scatterpie(
    data = mcav_site_pies,
    aes(
      x = pie_lon,
      y = pie_lat,
      r = 0.010
    ),
    cols = c("L4", "L7", "L21", "L25", "L30"),
    color = "black",
    linewidth = 0.3
  ) +
  
  geom_text(
    data = mcav_site_pies,
    aes(
      x = pie_lon,
      y = pie_lat,
      label = Site
    ),
    nudge_y = 0.012,
    size = 4
  ) +
  
  coord_fixed(
    xlim = c(-81.48, -81.38),
    ylim = c(24.50, 24.66),
    expand = FALSE
  ) +
  
  labs(
    title = "Montastraea cavernosa",
    x = "Longitude",
    y = "Latitude",
    fill = "Lineage"
  ) +
  
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "italic", size = 13),
    panel.grid = element_blank()
  )

print(mcav_map)






























































































































