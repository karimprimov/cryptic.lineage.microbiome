#loading libraries
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
library(tidytext)
library(indicspecies)
library(adegenet)
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")
library(DESeq2)

#reading in family level tsv file:
a.hyacinthus.genus.bracken.output.df <- read_tsv("/Users/karimprimov/Dropbox/UTAustin/2026/mcav.cryptic.lineages/mcav.flkeys.metagenomics/acropora.hyacinthus.bracken.reports/genus.level/260302.a.hyacinthus.concatenated.genus.level.report")
a.hyacinthus.srr.names <- read_csv("/Users/karimprimov/Dropbox/UTAustin/2026/mcav.cryptic.lineages/mcav.flkeys.metagenomics/a.hyacinthus.srr.samplename.csv")
a.hyacinthus.lineage <- read_csv("/Users/karimprimov/Dropbox/UTAustin/2026/mcav.cryptic.lineages/mcav.flkeys.metagenomics/260302.ahyacinthus.cryptic.lineage.csv")
a.hyacinthus.lat.long.symb <- read_csv("/Users/karimprimov/Dropbox/UTAustin/2026/mcav.cryptic.lineages/mcav.flkeys.metagenomics/hyacinthus.lat.long.symbiont.csv")
a.hyacinthus.lat.long.symb <- a.hyacinthus.lat.long.symb %>%
  column_to_rownames(var = "Sample")

a.hyacinthus.lineage.new <- a.hyacinthus.lineage %>%
  column_to_rownames(var = "sample")

readcount.df <- data.frame(readcount = rowSums(cleaned.pivot.hyacinth.dom.microbe.df, na.rm = TRUE))

a.hyacinthus.genus.bracken.output.df <- a.hyacinthus.genus.bracken.output.df %>%
  separate_wider_delim(cols = "SRR12466364.bracken.genus.level.output.report name",
                       delim = " ",
                       names  = c("name", "Taxa"),
                       too_many = "merge")
                       
a.hyacinthus.genus.bracken.output.cleaned.df <- a.hyacinthus.genus.bracken.output.df %>%
  filter(!Taxa %in% c("name", "Homo", "Cutibacterium","Staphylococcus","Corynebacterium","Bacteroides","Prevotella","Blautia","Parabacteroides","Phocaeicola","Fusobacterium","Escherichia","Klebsiella","Enterobacter","Citrobacter","Serratia","Morganella","Providencia","Proteus","Cronobacter","Lactobacillus","Enterococcus","Haemophilus","Aggregatibacter","Actinomyces","Streptococcus","Clostridium","Propionibacterium",
                      "Mycobacterium","Mycolicibacterium","Nocardia","Bartonella","Campylobacter",
                      "Bradyrhizobium","Mesorhizobium","Rhizobium","Xanthomonas","Pantoea",
                      "Bacillus","Paenibacillus","Lysinibacillus","Neobacillus","Priestia",
                      "Streptomyces","Amycolatopsis","Micromonospora","Sorangium",
                      "Leptospira","Spiroplasma","Nostoc","Singulisphaera","Siccibacter","Filimonas",
                      "Tequatrovirus"))
a.hyacinthus.genus.bracken.output.cleaned.df$name <- substr(a.hyacinthus.genus.bracken.output.cleaned.df$name, 1, 11)

cleaned.hyacinthus.df <- a.hyacinthus.genus.bracken.output.cleaned.df %>%
  dplyr::select(name, Taxa, new_est_reads) %>%
  mutate(new_est_reads = as.numeric(new_est_reads)) %>%
  left_join(a.hyacinthus.srr.names, by = "name")

cleaned.sample.hyacinthus.df <- cleaned.hyacinthus.df %>%
  dplyr::select(sample, Taxa, new_est_reads)

cleaned.pivot.hyacinth.df <- cleaned.sample.hyacinthus.df %>%
  pivot_wider(names_from = "sample", values_from = "new_est_reads", values_fn = sum)

transposed.hyacinthus.df <- cleaned.pivot.hyacinth.df %>%
  column_to_rownames("Taxa") %>%
  t() %>%
  as.data.frame()

cleaned.pivot.hyacinth.dom.microbe.df <- transposed.hyacinthus.df %>%
  dplyr::select(where(~ is.numeric(.) && sum(!is.na(.)) >= (nrow(transposed.hyacinthus.df) * 0.05))) 


cleaned.hyacinth.df <- merge(cleaned.pivot.hyacinth.dom.microbe.df, a.hyacinthus.lineage.new, by = "row.names") %>%
  column_to_rownames("Row.names")

#removing largest sample "AH214", and 20+ other samples to keep read count similar across lineages
samplestoremove <- c("AH214", "AH2", "AH206", "AH231", "AH209", "AH264", "AH6", "AH7", "AH173", "AH1", "AH2", "AH151", "AH183", "AH150", "AH176", "AH121", "AH28", "AH68", "AH235", "AH239", "AH371", "AH73", "AH294", "AH372", "AH344")
cleaned.hyacinth.reduced.he.lineage.df <- cleaned.hyacinth.df[!rownames(cleaned.hyacinth.df) %in% samplestoremove,]

#trimming lineage df 
a.hyacinthus.lineage.new.trim <- a.hyacinthus.lineage.new[rownames(a.hyacinthus.lineage.new) %in% rownames(cleaned.hyacinth.reduced.he.lineage.df), , drop = FALSE]



a.hyacinthus.lineage.new.trim.samplecol <- a.hyacinthus.lineage.new.trim %>%
  rownames_to_column()


a.hyacinthus.lat.long.symb <- a.hyacinthus.lat.long.symb %>%
  mutate(rowname = a.hyacinthus.lat.long.symb$Sample)

a.hyacinthus.lat.long.symb <- a.hyacinthus.lat.long.symb %>%
  left_join(a.hyacinthus.lineage.new.trim.samplecol, by = "rowname")


hyacinthus.metadata.df <- a.hyacinthus.lat.long.symb %>%
  dplyr::select(Sample, Latitude, Longitude, Symbiont, Site, Lineage=lineage) %>%
  column_to_rownames("Sample")


cleaned.hyacinth.reduced.he.lineage.df <- cleaned.hyacinth.reduced.he.lineage.df %>%
  column_to_rownames("Sample")

microbe.df <- cleaned.hyacinth.reduced.he.lineage.df[,1:35]
microbe.df[is.na(microbe.df)] <- 0

hyacinthus.metadata.df$readcount <- rowSums(microbe.df, na.rm = TRUE)


a.hyacinthus.summary.stand.df <- decostand(microbe.df, method = "clr", pseudocount = 0.01)

a.hyacinthus.microbe.matrix <- as.matrix(vegdist(a.hyacinthus.summary.stand.df, method = "euclidean", na.rm=TRUE))

a.hyacinthus.microbe.pca <- capscale(a.hyacinthus.microbe.matrix ~ Site + Lineage + Condition(readcount), data = hyacinthus.metadata.df)
plot(a.hyacinthus.microbe.pca$CA$eig)

plot(a.hyacinthus.microbe.pca)



site_scores <- as.data.frame(scores(a.hyacinthus.microbe.pca, display = "sites", choices = c(1:5)))

# 3. Merge with your metadata to get your 'group' or 'lineage' columns for coloring
plot_data <- cbind(site_scores, hyacinthus.metadata.df)


ggplot(plot_data, aes(x = CAP1, y = CAP2, color = Lineage)) + 
  geom_point(size = 7, alpha = 0.8) +
  theme_minimal() +
  labs(
    title = "Partial db-RDA (Capscale)",
    subtitle = "Conditioned on Read Count",
    x = "CAP1", 
    y = "CAP2"
  ) +
  scale_color_brewer(palette = "Set1") # Change palette as needed

####################################################################################
####################################################################################

#260409:
#adonis 2:

adonis.site.results <- adonis2(a.hyacinthus.microbe.matrix ~ readcount + Site + Lineage, 
                          data = hyacinthus.metadata.df, 
                          method = "bray", # or your preferred distance metric
                          by = "terms")   # marginal tests are best for independent effects

# 2. View results
print(adonis.site.results)


adonis.site.lineage.results <- adonis2(a.hyacinthus.microbe.matrix ~ readcount + Site + Lineage, 
                                       data = hyacinthus.metadata.df, 
                                       method = "bray", # or your preferred distance metric
                                       by = "terms")   # marginal tests are best for independent effects

# 2. View results
print(adonis.site.lineage.results)

####################################################################################
####################################################################################

db_site_model <- dbrda(a.hyacinthus.microbe.matrix ~ Site + Condition(readcount), 
                  data = hyacinthus.metadata.df, 
                  distance = "bray")

# 2. Get the PERMANOVA table (identical results to adonis2)
anova(db_site_model, by = "margin")



db_site_lineage_model <- dbrda(a.hyacinthus.microbe.matrix ~ Site + Lineage + Condition(readcount), 
                  data = hyacinthus.metadata.df, 
                  distance = "bray")

# 2. Get the PERMANOVA table (identical results to adonis2)
anova(db_site_lineage_model, by = "margin")






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


disp_site <- betadisper(as.dist(a.hyacinthus.microbe.matrix), hyacinthus.metadata.df$Site)
permutest(disp_site) 

# For ClusterID
disp_cluster <- betadisper(as.dist(a.hyacinthus.microbe.matrix), hyacinthus.metadata.df$Lineage)
permutest(disp_cluster)


table(pocil.metadata.ordered$Site, pocil.metadata.ordered$ClusterID)

# 2. Chi-Squared Test (Statistical test of association)
test <- chisq.test(table(hyacinthus.metadata.df$Site, hyacinthus.metadata.df$Lineage))
print(test)

































































































#venn diagram of microbes found among and within each lineage
ha.venn <- colSums(a.hyacinthus.summary.df[1:4,])
hb.venn <- colSums(a.hyacinthus.summary.df[5:8,])
hc.venn <- colSums(a.hyacinthus.summary.df[9:12,])
he.venn <- colSums(r[13:16,])

ha.microbes <- names(ha.venn[ha.venn>0])
hb.microbes <- names(hb.venn[hb.venn>0])
hc.microbes <- names(hc.venn[hc.venn>0])
he.microbes <- names(he.venn[hd.venn>0])

ven.list <- list(HA = ha.microbes, HB = hb.microbes, HC = hc.microbes, HE = he.microbes)

venn_plot <- venn.diagram(
  x = ven.list,
  category.names = c("HA", "HB", "HC", "HE"),
  filename = NULL,
  fill = c("skyblue", "pink", "orange", "limegreen"),
  cex = 4,
  cat.cex = 4 
)

# Display the plot
grid.draw(venn_plot)

upset(fromList(ven.list), 
      order.by = "freq", 
      sets.bar.color = c("skyblue", "pink", "orange", "limegreen"),
      nsets = 4,
      main.bar.color = "gray30",
      matrix.color = "black",
      text.scale = 3)


#260310:
#adonis:
ahyacinthus.lineage.adonis <- adonis2(a.hyacinthus.microbe.matrix ~ lineage, data=readcount.summary.df)
ahyacinthus.lineage.adonis


ahyacinthus.lineage.microbe.abundances <- a.hyacinthus.summary.df %>%
  rownames_to_column(var = "Sample") %>%
  mutate(Lineage = substr(Sample, 1, 2)) %>%
  pivot_longer(cols = -c(Sample, Lineage), names_to = "category", values_to = "value")




ahyacinthus.lineage.microbe.abundances %>%
  # 1. Sum values for each category (microbe) within each Lineage
  group_by(Lineage, category) %>%
  summarise(total_abundance = sum(value, na.rm = TRUE), .groups = "drop") %>%
  
  # 2. Identify the top 10 microbes per Lineage based on that sum
  group_by(Lineage) %>%
  slice_max(order_by = total_abundance, n = 10, with_ties = FALSE) %>%
  ungroup() %>%
  
  # 3. Plotting
  ggplot(aes(x = reorder_within(category, -total_abundance, Lineage), y = total_abundance)) +
  geom_col(fill = "turquoise") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) + 
  scale_x_reordered() + 
  facet_wrap(~Lineage, scales = "free_x") + 
  theme_minimal() + 
  theme(
    axis.text.y = element_text(size=24),
    axis.text.x = element_text(size=24, angle = 45, hjust = 1, vjust = 1),
    strip.text = element_text(size = 30, face = "bold"),
    axis.title.y = element_blank(),
    axis.title.x = element_blank()
  ) 



#260312:
#stacked barchart: 
a.hyacinthus.summary.stand.hellinger.df <- decostand(a.hyacinthus.summary.df, method = "hellinger", pseudocount = 0.01)
a.hyacinthus.summary.stand.hellinger.df <- a.hyacinthus.summary.stand.hellinger.df %>%
  rownames_to_column(var = "Sample")

a.hyacinthus.summary.stand.hellinger.long.df <- a.hyacinthus.summary.stand.hellinger.df%>%
  pivot_longer(cols = -Sample,
               names_to = "Species",
               values_to = "Value")


a.hyacinthus.final_plot.df <- a.hyacinthus.summary.stand.hellinger.long.df %>%
  # 1. Replace any NAs with 0 so they don't break the 'sum' or 'fill'
  mutate(Value = replace_na(Value, 0)) %>%
  
  # 2. Group into 'Other' and top species
  mutate(Species_Grouped = fct_lump_n(Species, n = 12, w = Value, other_level = "Other")) %>%
  
  # 3. CRITICAL: Sum the values so there is ONLY ONE row per Sample/Species pair
  # This removes the "lines" inside the 'Other' block and fixes the height math
  group_by(Sample, Species_Grouped) %>%
  summarize(Value = sum(Value), .groups = "drop") %>%
  
  # 4. Fill in any missing combinations (e.g., if HA_1 doesn't have 'Vibrio')
  complete(Sample, Species_Grouped, fill = list(Value = 0)) %>%
  
  # 5. Force 'Other' to be the base of the stack
  mutate(Species_Grouped = fct_relevel(Species_Grouped, "Other"))


species_list <- levels(a.hyacinthus.final_plot.df$Species_Grouped)
n_species <- length(species_list)

# 2. Create the color vector (Paired has 12, we add others if needed)
# Using 'Set3' or combining 'Paired' with 'Dark2' for high contrast
library(RColorBrewer)
my_colors <- c(brewer.pal(12, "Paired"), brewer.pal(8, "Dark2"))[1:n_species]

# 3. Name the colors so ggplot matches them to the right species
names(my_colors) <- species_list

# 4. FIX VIBRIO: Change it to a dark color if it's currently too light
# This only works if "Vibrio" is one of the top 12 species in your levels
if("Vibrio" %in% names(my_colors)) {
  my_colors["Vibrio"] <- "firebrick4" # Dark Red
}

# 5. Plot using the CORRECT dataframe
ggplot(a.hyacinthus.final_plot.df, aes(x = Sample, y = Value, fill = Species_Grouped)) +
  # Using position_fill(reverse = TRUE) ensures they reach 1 and anchor at the bottom
  geom_col(position = position_fill(reverse = TRUE), width = 0.8) + 
  scale_fill_manual(values = my_colors) + 
  scale_y_continuous(expand = c(0,0)) + 
  theme_minimal() +
  theme(axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 12))









# 2. Plot using reorder(..., FUN = sum)
ggplot(ahyacinthus.lineage.microbe.abundances, aes(x = reorder_within(category, -value, Lineage, FUN = sum), y = value)) +
  geom_col() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) + 
  facet_wrap(~Lineage, scales = "free_x") + 
  theme(axis.text.x = element_text(angle=45, hjust = 1, vjust = 1)) + 
  labs(x = "Category (Ordered by Sum)", y = "Total Value")


ahyacinthus.lineage.microbe.abundances.new <- ahyacinthus.lineage.microbe.abundances %>%
  # 1. Filter for the top 10 per Lineage
  group_by(Lineage) %>%
  slice_max(order_by = value, n = 10, with_ties = FALSE) %>% 
  ungroup()

ggplot(ahyacinthus.lineage.microbe.abundances.new, aes(x = reorder_within(category, value, Lineage), y = value)) +
  geom_col() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) + 
  scale_x_reordered() + # Cleans up the "___Lineage" suffix from labels
  facet_wrap(~Lineage, scales = "free_x") + 
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) + 
  labs(x = "Top 10 Categories", y = "Total Value")

#








#

#indicator species analysis:
hyacinth.abundance.df <- ahyacinthus_summary[, 4:67]
lineage <- ahyacinthus_summary$lineage
a.hyacinthus.isa <- multipatt(hyacinth.abundance.df, lineage, func = "r.g", control = how(nperm = 999))
summary(a.hyacinthus.isa)


#acropora hyacinthus genera heat map
#deseq variance-stabilized transform and substract by mean of each row (use function in pheatmap which substracts mean count, to see variation across samples up and downs will be log-fold changes)


#normalize first with deseq (vst) then scale with pheatmap's scale
cleaned.hyacinth.reduced.he.lineage.t.df <- t(cleaned.hyacinth.reduced.he.lineage.df)
#remvoing lineage ID row:
cleaned.hyacinth.reduced.he.lineage.t.df <- cleaned.hyacinth.reduced.he.lineage.t.df[1:64,]
storage.mode(cleaned.hyacinth.reduced.he.lineage.t.df) <- "integer"
cleaned.hyacinth.reduced.he.lineage.t.df.pseudocount <- decostand(cleaned.hyacinth.reduced.he.lineage.t.df, method = "clr", pseudocount = 10)


vst.allsamples.df <- varianceStabilizingTransformation(cleaned.hyacinth.reduced.he.lineage.df)


# vst.transposed.hyacinthus.summary.df <- varianceStabilizingTransformation(transpose.a.hyacinthus.summary.df)
# 
# transposed.clr.df <- t(a.hyacinthus.summary.stand.df)

pheatmap(vst.transposed.hyacinthus.summary.df,
         clustering_distance_cols = "euclidean", 
         clustering_method = "ward.D2",
         scale = "row", # Z-score helps see which taxa define which lineage
         main = "Lineage Signal (VST Transformed)")

plot(hyacinthus.clr.heatmap)
#what makes them cluster




#beta dispersion
a.hyacinthus.microbe.distances <- vegdist(a.hyacinthus.summary.stand.df, method = "euclidean", na.rm=TRUE)
a.hyacinthus.betadisper <- betadisper(a.hyacinthus.microbe.distances, lineage)
a.hyacinthus.betadisper.permtest <- permutest(a.hyacinthus.betadisper)
a.hyacinthus.betadisper.permtest








# 1. Calculate row sums
row_sums <- rowSums(a.hyacinthus.summary.df)

# 2. Create the plot without the x-axis (xaxt = "n")
plot(row_sums, 
     xaxt = "n", 
     type = "b", # 'b' for both points and lines, 'p' for points only
     xlab = "Samples", 
     ylab = "Sum of Values",
     main = "Row Sums by Sample")

# 3. Add the rownames as the x-axis
axis(1, at = 1:length(row_sums), labels = rownames(a.hyacinthus.summary.df),cex.axis = 1.5)







#############

lineage_subgroups <- cleaned.pivot.hyacinth.dom.microbe.df %>% 
  rownames_to_column("sample") %>% 
  inner_join(a.hyacinthus.lineage, by = "sample") %>%
  group_by(lineage) %>%
  # Assign each sample to one of 6 groups (e.g., 1, 2, 3, 4, 5, 6, 1, 2...)
  mutate(subgroup = (row_number() %% 5) + 1) %>% 
  ungroup()

# 2. Summarise by both lineage AND subgroup
# This gives you 6 rows per lineage instead of 1
lineage_summary <- lineage_subgroups %>% 
  group_by(lineage, subgroup) %>% 
  summarise(across(where(is.numeric), ~sum(.x, na.rm = TRUE)), .groups = "drop") %>%
  # Create a unique ID for each of the 5-6 points (e.g., "LineageA_1")
  unite("lineage_id", lineage, subgroup, remove = FALSE) 

# 3. CLR Transformation
lineage_summary.stand <- lineage_summary %>% 
  column_to_rownames(var="lineage_id") %>% 
  dplyr::select(where(is.numeric)) %>%
  decostand(method = "clr", pseudocount = 0.01)

# 4. Ordination (Aitchison Distance)
pca.res <- capscale(lineage_summary.stand ~ 1, distance = "euclidean")

lineage_labels <- as.factor(lineage_summary$lineage)
group_colors <- as.numeric(lineage_labels) # Converts factor levels to 1, 2, 3...

# 2. Create the empty plot
plot(pca.res, type = "n", scaling = 1, main = "PCA by Lineage Subgroups")

# 3. Add the points
points(pca.res, display = "sites", 
       col = group_colors, 
       pch = 16, 
       cex = 1.2)

# 4. Add the legend
# We use levels() for the text and 1:length(levels) for the colors
legend("topright", 
       legend = levels(as.factor(lineage_summary$lineage)),
       col = 1:length(levels(as.factor(lineage_summary$lineage))), # Uses factor indices
       pch = 16, 
       pt.cex = 1.2,
       bty = "n", 
       cex = 0.8,
       xpd = TRUE) # This keeps the legend from being cut off

site_scores <- as.data.frame(scores(pca.res, display = "sites"))
site_scores$lineage <- lineage_summary$lineage

# 2. Plot with ggplot
ggplot(site_scores, aes(x = MDS1, y = MDS2, color = lineage)) +
  geom_point(size = 3) +
  theme_minimal() +
  labs(title = "PCA of A. hyacinthus Lineages",
       subtitle = "5-6 pooled samples per lineage") +
  theme(legend.position = "right") # Legend is now permanently attached



#check for variation in coverage
#pool first and then filter taxa 
#see if im pooling big and small ones together
#take group of samples, see how sample profiles look after pooling compared to non-pooling

















a.hyacinthus.genotype.data.csv <- read.csv("/Users/karimprimov/Dropbox/UTAustin/2026/a.hyacinthus.cryptic.lineages/a.hyacinthus.flkeys.metagenomics/a.hyacinthus.genotype.hybrid.data.csv", header = TRUE, row.names = 1)


a.hyacinthus.microbe.scores <- scores(pca.a.hyacinthus.mcirobe, display = "sites", choices=c(1:2)) 

merged.a.hyacinthus.microbe.scores.metadata <- merge(a.hyacinthus.microbe.scores, a.hyacinthus.genotype.data.csv, by = "row.names", all.x = TRUE)

ggplot(merged.a.hyacinthus.microbe.scores.metadata, aes(x = MDS1, y = MDS2, color = Lineage, label = rownames(merged.a.hyacinthus.microbe.scores.metadata))) +
  geom_point(size = 12) +
  geom_label_repel(aes(label = Row.names),
                   box.padding   = 0.35, 
                   point.padding = 0.5,
                   segment.color = 'grey50') +
  theme_minimal()

hc.microbe=hclust(as.dist(a.hyacinthus.microbe.matrix),method="ave")
plot(hc.microbe)








#
a.hyacinthus.microbe.df.t.df

prev <- colSums(a.hyacinthus.microbe.df.t.df > 0)

n_samples <- nrow(a.hyacinthus.microbe.df.t.df)
cat("Number of a.hyacinthus samples:", n_samples, "\n")
cat("Number of ASVs total:", ncol(a.hyacinthus.microbe.df.t.df), "\n")

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
    x = "Minimum number of a.hyacinthus individuals sharing an ASV",
    y = "Number of ASVs remaining",
    title = "a.hyacinthus ASV prevalence (sharedness) curve"
  )



keep_core_asvs <- names(prev[prev >= 11])

# 2. Subset your original dataframe to keep only these 951 columns
a.hyacinthus.core.15.df <- a.hyacinthus.asv.abundances.trimmed.df[, keep_core_asvs]

# Verify the count
cat("Number of ASVs remaining:", ncol(a.hyacinthus.core.15.df), "\n") # Should be 951


# 3. Create the stats dataframe for the 951 ASVs
asv_stats_core <- data.frame(
  ASV = colnames(a.hyacinthus.core.15.df),
  sum_abundance = colSums(agaricia.core.15.df),
  sd_abundance = apply(agaricia.core.15.df, 2, sd)
)

# 4. Plot
library(ggplot2)

ggplot(asv_stats_core, aes(x = reorder(ASV, sum_abundance), y = sum_abundance)) +
  # 1. The Points
  geom_point(color = "red", size = 0.8, alpha = 0.6) +
  
  # 2. The Error Bars (Calculated on raw scale, then log-transformed by the scale layer)
  geom_errorbar(aes(ymin = sum_abundance - sd_abundance, 
                    ymax = sum_abundance + sd_abundance), 
                width = 0.1, color = "darkgreen", alpha = 0.3) +
  
  # 3. Use a Log Scale for the Y-axis (this handles the skewness automatically)
  scale_y_log10() +
  
  coord_flip() + 
  theme_minimal() +
  labs(
    title = "Abundance of Core ASVs (n=951)",
    subtitle = "Shared by >= 15 individuals | Y-axis log10 transformed",
    x = "ASVs (Ranked by Total Sum)",
    y = "Total Read Sum (Log10 Scale)"
  ) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid.major.y = element_blank()
  )

























#saving work in RData file
save.image(file = "/Users/karimprimov/Dropbox/UTAustin/2026/mcav.cryptic.lineages/mcav.flkeys.metagenomics/260506.a.hyacinthus.cryptic.lineage.microbiome.analysis.RData")











































# 1. Generate PCA using ONLY the consolidated lineage data
pca_model <- prcomp(lineage_summary.stand)

# 2. Project the 112 individual samples onto that PCA space
projected_samples <- predict(pca_model, newdata = samples_clr)

# 3. Convert to dataframes for plotting
df_lineages <- as.data.frame(pca_model$x) %>% rownames_to_column("lineage")
df_samples <- as.data.frame(projected_samples) %>% rownames_to_column("sample")

# Add lineage metadata to individual samples for coloring
df_samples <- df_samples %>%
  left_join(a.hyacinthus.lineage, by = "sample")








#260311:
# -------------------------------
# 1. Subset your samples and taxa
# -------------------------------

# Keep only samples you are analyzing
samples_to_keep <- rownames(cleaned.pivot.hyacinth.dom.microbe.df) %in% rownames(a.hyacinthus.lineage.new.trim)
cleaned.hyacinth.reduced.df <- cleaned.pivot.hyacinth.dom.microbe.df[samples_to_keep, ]

# Remove rare taxa (present in less than 10 samples)
taxa_prevalence <- colSums(cleaned.hyacinth.reduced.df > 0)
cleaned.hyacinth.filtered <- cleaned.hyacinth.reduced.df[, taxa_prevalence >= 10]

# -------------------------------
# 2. VST transformation (DESeq2)
# -------------------------------

library(DESeq2)

# Prepare DESeq2 dataset
dds <- DESeqDataSetFromMatrix(
  countData = t(cleaned.hyacinth.filtered),  # DESeq2 expects taxa as rows, samples as columns
  colData = data.frame(
    lineage = a.hyacinthus.lineage.new.trim$lineage[
      match(colnames(cleaned.hyacinth.filtered), rownames(a.hyacinthus.lineage.new.trim))
    ]
  ),
  design = ~ lineage
)

# Apply variance-stabilizing transformation
vst_data <- assay(varianceStabilizingTransformation(dds))

# Transpose back to samples as rows, taxa as columns for heatmap
vst_samples.df <- t(vst_data) %>% as.data.frame()

# -------------------------------
# 3. Prepare lineage annotation for heatmap
# -------------------------------

lineage_annotation <- a.hyacinthus.lineage.new.trim$lineage[
  match(rownames(vst_samples.df), rownames(a.hyacinthus.lineage.new.trim))
]

annotation_df <- data.frame(Lineage = lineage_annotation)
rownames(annotation_df) <- rownames(vst_samples.df)

# -------------------------------
# 4. Plot heatmap
# -------------------------------

library(pheatmap)

pheatmap(
  vst_samples.df,
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "ward.D2",
  annotation_row = annotation_df,
  scale = "row",           # Z-score per taxon to see variation up/down
  main = "A. hyacinthus Microbiome Heatmap (VST)",
  show_rownames = FALSE,
  show_colnames = FALSE
)



# Euclidean distance
standardized <- decostand(a.hyacinthus.summary.df, method="total")  # relative abundances
log_standardized <- log10(standardized + 0.01)  # dampens high counts

# 2. Euclidean distance clustering
d_euc <- dist(log_standardized, method = "euclidean")
hc_euc <- hclust(d_euc, method = "ward.D2")
plot(hc_euc, main = "Euclidean Clustering of Log-Transformed Relative Abundance")

# Optional: color labels by lineage
lineage_colors <- a.hyacinthus.summary.df$lineage  # make sure this matches rownames
library(dendextend)
dend <- as.dendrogram(hc_euc)
labels_colors(dend) <- as.numeric(factor(lineage_colors))
plot(dend, main = "Euclidean Clustering with Lineage Colors")

# 3. Correlation distance clustering
d_cor <- as.dist(1 - cor(t(log_standardized)))  # transpose because cor() works on columns
hc_cor <- hclust(d_cor, method = "ward.D2")
dend_cor <- as.dendrogram(hc_cor)
labels_colors(dend_cor) <- as.numeric(factor(lineage_colors))
plot(dend_cor, main = "Correlation Clustering with Lineage Colors")




































#260312:
#seeing which microbes are most variable 
#generating new, more robust metasamples to do analysis on:

cleaned.pivot.hyacinth.dom.microbe.df <- transposed.hyacinthus.df 


cleaned.hyacinth.df <- merge(cleaned.pivot.hyacinth.dom.microbe.df, a.hyacinthus.lineage.new, by = "row.names") %>%
  column_to_rownames("Row.names")

#removing largest sample "AH214", and 20+ other samples to keep read count similar across lineages
samplestoremove <- c("AH214", "AH2", "AH206", "AH231", "AH209", "AH264", "AH6", "AH7", "AH173", "AH1", "AH2", "AH151", "AH183", "AH150", "AH176", "AH121", "AH28", "AH68", "AH235", "AH239", "AH371", "AH73", "AH294", "AH372", "AH344")
cleaned.hyacinth.reduced.he.lineage.df <- cleaned.hyacinth.df[!rownames(cleaned.hyacinth.df) %in% samplestoremove,]

#designating groups and lineages
cleaned.hyacinth.reduced.he.lineage.grouped.df <- cleaned.hyacinth.reduced.he.lineage.df %>%
  rownames_to_column("ID") %>%
  left_join(rownames_to_column(a.hyacinthus.lineage.new.trim, "ID"), by = c("ID", "lineage")) %>%
              column_to_rownames("ID")
            
cleaned.hyacinth.reduced.he.lineage.grouped.df$group <- sample(1:4, nrow(cleaned.hyacinth.reduced.he.lineage.grouped.df), replace = TRUE)


ahyacinthus.grouped.summary <- cleaned.hyacinth.reduced.he.lineage.grouped.df %>% 
  group_by(lineage, group) %>% 
  summarise(across(where(is.numeric), ~sum(.x, na.rm = TRUE)), .groups = "drop") %>%
  # Create a unique ID for each of the 5-6 points (e.g., "LineageA_1")
  unite("lineage_id", lineage, group, remove = FALSE) 


total_counts <- rowSums(ahyacinthus.grouped.summary[, 4:213], na.rm = TRUE)


hyacinthus.samplefilt.df <- ahyacinthus.grouped.summary[c(1:4,6:12,14:16),]
rowSums(hyacinthus.samplefilt.df[,4:213], na.rm = TRUE)

row_sums <- rowSums(hyacinthus.samplefilt.df[,4:213], na.rm = TRUE)

# 2. Create the plot without the x-axis (xaxt = "n")
par(mar = c(5, 7, 4, 2)) 

# 2. Create the plot (keep axes = FALSE to avoid double-drawing)
plot(row_sums, 
     axes = FALSE,      
     type = "b", 
     ylim = c(0, 10000), 
     xlab = "Samples", 
     ylab = "",
     cex.lab = 2)      # Optional: makes the "Sum of Values" title larger too

# 3. Add X-axis
axis(1, at = 1:length(row_sums), labels = rownames(standardized), cex.axis = 2)

# 4. Add Y-axis (mgp moves the labels further away from the tick marks if needed)
axis(2, cex.axis = 2, las = 1) # las = 1 makes labels horizontal/easier to read
title(ylab = "Sum of Values", line = 5, cex.lab = 2)

box()


hyacinthus.samplefilt.microbe.filter.df <- hyacinthus.samplefilt.df %>%
  dplyr::select(where(~ is.numeric(.) && sum(!is.na(.)) >= (ncol(hyacinthus.samplefilt.df) * 0.01))) 

hyacinthus.microbe.var.df <- hyacinthus.samplefilt.microbe.filter.df %>%
  summarise(across(-group, var)) %>%   # Calculate variance for all columns except 'group'
  pivot_longer(everything(), names_to = "taxon", values_to = "variance") %>%
  arrange(desc(variance)) %>%          # Sort highest to lowest
  head(50)                             # Show top 10


hyacinthus.highvarmicrobe.df <- hyacinthus.samplefilt.microbe.filter.df[, colnames(hyacinthus.samplefilt.microbe.filter.df) %in% hyacinthus.microbe.var.df$taxon]
rownames(hyacinthus.highvarmicrobe.df) <- hyacinthus.samplefilt.df$lineage_id


standardized <- decostand(hyacinthus.highvarmicrobe.df, method="total")  # relative abundances
log_standardized <- log10(standardized + 0.01)  # dampens high counts

# 2. Euclidean distance clustering
d_bray <- vegan::vegdist(standardized, method = "bray")
hc_bray <- hclust(d_bray, method = "average")
plot(hc_bray, main = "Euclidean Clustering of Log-Transformed Relative Abundance")

  
hyacinthus.highvarmicrobe.df.copy <- hyacinthus.highvarmicrobe.df 
hyacinthus.highvarmicrobe.df.copy$lineage = substr(rownames(hyacinthus.highvarmicrobe.df.copy), 1, 2)


a.hyacinthus.summary.stand.df <- decostand(hyacinthus.highvarmicrobe.df, method = "clr", pseudocount = 0.01)

a.hyacinthus.microbe.matrix <- as.matrix(vegdist(a.hyacinthus.summary.stand.df, method = "euclidean", na.rm=TRUE))

a.hyacinthus.microbe.pca <- capscale(a.hyacinthus.microbe.matrix ~ lineage, data = hyacinthus.highvarmicrobe.df.copy)
plot(a.hyacinthus.microbe.pca$CA$eig)

plot(a.hyacinthus.microbe.pca)

a.hyacinthus.microbe.pca

plot(c(a.hyacinthus.microbe.pca$CCA$eig, a.hyacinthus.microbe.pca$CA$eig))

site_scores <- as.data.frame(scores(a.hyacinthus.microbe.pca, display = "sites", choices = c(1:5)))

# 3. Merge with your metadata to get your 'group' or 'lineage' columns for coloring
plot_data <- cbind(site_scores, hyacinthus.highvarmicrobe.df.copy)


ggplot(plot_data, aes(x = CAP1, y = CAP2, color = as.factor(lineage))) + 
  geom_point(size = 12, alpha = 0.8) +
  theme_minimal() +
  labs(
    title = "Partial db-RDA (Capscale)",
    subtitle = "Conditioned on Read Count",
    x = "CAP1", 
    y = "CAP2"
  ) +
  scale_color_brewer(palette = "Set1") # Change palette as needed


pheatmap(t(a.hyacinthus.summary.stand.df),
         clustering_distance_rows = "euclidean", 
         main = "Lineage Signal (CLR Transformation)",
         angle_col = 45,        # Tilts the sample names to 45 degrees
         fontsize_col = 18,     # Adjust this number to make the font larger
         fontsize_row = 10)     # Optional: adjust the feature/row font size too

plot(hyacinthus.clr.heatmap)

a.hyacinthus.hivarmiceobes.isa <- multipatt(hyacinthus.highvarmicrobe.df, hyacinthus.highvarmicrobe.df.copy$lineage, func = "r.g", control = how(nperm = 999))
summary(a.hyacinthus.hivarmiceobes.isa)


































































#saving work in RData file
save.image(file = "/Users/karimprimov/Dropbox/UTAustin/2026/mcav.cryptic.lineages/mcav.flkeys.metagenomics/261031.a.hyacinthus.cryptic.lineage.microbiome.analysis.RData")


#260929:
# ============================================================
# ACROPORA HYACINTHUS — FINAL DATASET CHECK
# ============================================================

hy_meta <- droplevels(hyacinthus.metadata.df)

hy_dist <- as.dist(a.hyacinthus.microbe.matrix)

# Dimensions/alignment
nrow(hy_meta)
attr(hy_dist, "Size")

identical(
  rownames(hy_meta),
  attr(hy_dist, "Labels")
)

# Sample sizes
table(hy_meta$Lineage, useNA = "ifany")
table(hy_meta$Site, useNA = "ifany")

# Missing metadata
colSums(
  is.na(
    hy_meta[, c(
      "Lineage",
      "Site",
      "Latitude",
      "Longitude",
      "readcount"
    )]
  )
)


aggregate(
  cbind(Latitude, Longitude) ~ Site,
  data = hy_meta,
  FUN = function(x) length(unique(x))
)

hy_meta$Site_unique <- interaction(
  hy_meta$Site,
  hy_meta$Latitude,
  hy_meta$Longitude,
  drop = TRUE,
  sep = "_"
)

nlevels(hy_meta$Site_unique)

table(hy_meta$Site_unique)

table(
  hy_meta$Site_unique,
  hy_meta$Lineage
)


hy_meta$Site <- factor(hy_meta$Site)
hy_meta$Lineage <- factor(hy_meta$Lineage)


library(vegan)

# ------------------------------
# Site-only PERMANOVA
# ------------------------------
set.seed(123)

hy_adonis_site <- adonis2(
  hy_dist ~ readcount + Site,
  data = hy_meta,
  permutations = 9999,
  by = "margin"
)

hy_adonis_site


# ------------------------------
# Lineage-only PERMANOVA
# ------------------------------
set.seed(123)

hy_adonis_lineage <- adonis2(
  hy_dist ~ readcount + Lineage,
  data = hy_meta,
  permutations = 9999,
  by = "margin"
)

hy_adonis_lineage


# ------------------------------
# Joint PERMANOVA
# ------------------------------
set.seed(123)

hy_adonis_both <- adonis2(
  hy_dist ~ readcount + Site + Lineage,
  data = hy_meta,
  permutations = 9999,
  by = "margin"
)

hy_adonis_both


# ------------------------------
# Site × lineage association
# ------------------------------
site_lineage_table <- table(
  hy_meta$Site,
  hy_meta$Lineage
)

site_lineage_table

hy_chisq <- chisq.test(site_lineage_table)
hy_chisq

library(lsr)
hy_cramers <- cramersV(site_lineage_table)
hy_cramers


# ------------------------------
# PERMDISP — Site
# ------------------------------
set.seed(123)

hy_disp_site <- betadisper(
  hy_dist,
  hy_meta$Site
)

permutest(
  hy_disp_site,
  permutations = 9999
)


# ------------------------------
# PERMDISP — Lineage
# ------------------------------
set.seed(123)

hy_disp_lineage <- betadisper(
  hy_dist,
  hy_meta$Lineage
)

permutest(
  hy_disp_lineage,
  permutations = 9999
)


# ------------------------------
# Full dbRDA
# ------------------------------
hy_dbrda_full <- dbrda(
  hy_dist ~ readcount + Site + Lineage,
  data = hy_meta
)

hy_dbrda_full

set.seed(123)

anova(
  hy_dbrda_full,
  by = "margin",
  permutations = 9999
)

library(vegan)
library(ggplot2)

# Final dbRDA model
hy_dbrda_plot <- dbrda(
  hy_dist ~ Lineage + Condition(Site + readcount),
  data = hy_meta
)

# Extract sample scores
hy_scores <- as.data.frame(
  scores(
    hy_dbrda_plot,
    display = "sites",
    choices = 1:2
  )
)

# Add sample IDs
hy_scores$SampleID_plot <- rownames(hy_scores)
hy_meta$SampleID_plot <- rownames(hy_meta)

# Merge scores with metadata
hy_plot_data <- merge(
  hy_scores,
  hy_meta,
  by = "SampleID_plot",
  all.x = TRUE
)

# Calculate % of constrained variation on first two axes
hy_eig <- eigenvals(
  hy_dbrda_plot,
  model = "constrained"
)

hy_db1 <- 100 * hy_eig[1] / sum(hy_eig)
hy_db2 <- 100 * hy_eig[2] / sum(hy_eig)

hy_db1
hy_db2

# Plot
hy_plot <- ggplot(
  hy_plot_data,
  aes(
    x = dbRDA1,
    y = dbRDA2,
    color = Lineage
  )
) +
  geom_point(
    size = 8,
    alpha = 0.85
  ) +
  theme_minimal() +
  coord_fixed() +
  labs(
    title = expression(italic("Acropora hyacinthus")),
    x = paste0("dbRDA1 (", round(hy_db1, 1), "%)"),
    y = paste0("dbRDA2 (", round(hy_db2, 1), "%)"),
    color = "Lineage"
  )

hy_plot

site_lineage_table <- table(
  hy_meta$Site,
  hy_meta$Lineage
)

site_lineage_table

hy_chisq <- chisq.test(site_lineage_table)
hy_chisq

library(lsr)
hy_cramers <- cramersV(site_lineage_table)
hy_cramers



set.seed(123)

fisher.test(
  site_lineage_table,
  simulate.p.value = TRUE,
  B = 99999
)


adonis.site.results <- adonis2(hy_dist ~ readcount + Site + Lineage, 
                               data = hyacinthus.metadata.df, 
                               method = "bray", # or your preferred distance metric
                               by = "terms")   # marginal tests are best for independent effects

adonis.site.results

set.seed(123)

hy_adonis_site <- adonis2(
  hy_dist ~ readcount + Site,
  data = hy_meta,
  permutations = 9999,
  by = "margin"
)

hy_adonis_site


# ------------------------------
# Lineage-only PERMANOVA
# ------------------------------
set.seed(123)

hy_adonis_lineage <- adonis2(
  hy_dist ~ readcount + Lineage,
  data = hy_meta,
  permutations = 9999,
  by = "margin"
)

hy_adonis_lineage


# ------------------------------
# Joint PERMANOVA
# ------------------------------
set.seed(123)

hy_adonis_both <- adonis2(
  hy_dist ~ readcount + Site + Lineage,
  data = hy_meta,
  permutations = 9999,
  by = "margin"
)

hy_adonis_both

# a.hyacinthus.summary.stand.df <- decostand(hyacinthus.highvarmicrobe.df, method = "clr", pseudocount = 0.01)
# 
# a.hyacinthus.microbe.matrix <- as.matrix(vegdist(a.hyacinthus.summary.stand.df, method = "euclidean", na.rm=TRUE))
# 
# a.hyacinthus.microbe.pca <- capscale(a.hyacinthus.microbe.matrix ~ ineage, data = hyacinthus.highvarmicrobe.df.copy)


hy_clr <- as.matrix(a.hyacinthus.summary.stand.df)
hy_clr <- hy_clr[rownames(hy_meta), , drop = FALSE]

dfa.lineage <- MASS::lda(hy_clr, grouping = hy_meta$Lineage)
plot(dfa.lineage)

hy_dapc <- dapc(
  hy_clr,
  grp = hy_meta$Lineage,
  n.pca = 10,
  n.da = 3
)

scatter(hy_dapc)


set.seed(123)

hy_cv <- xvalDapc(
  hy_clr,
  grp = hy_meta$Lineage,
  n.pca.max = 20,
  n.rep = 30,
  xval.plot = TRUE
)

hy_cv

hy_dapc <- adegenet::dapc(hy_clr, grp = hy_meta$Lineage,
                          n.pca = 8, n.da = 3)

scatter(hy_dapc)


hy_disp_lineage <- vegan::betadisper(hy_dist, hy_meta$Lineage)
permutest(hy_disp_lineage, permutations = 9999)

hy_disp_site <- vegan::betadisper(hy_dist, factor(hy_meta$Site))
permutest(hy_disp_site, permutations = 9999)


set.seed(123)

adonis2(
  hy_dist ~ readcount + Site + Lineage,
  data = hy_meta,
  permutations = permute::how(
    nperm = 9999,
    blocks = factor(hy_meta$Site)
  ),
  by = "margin"
)

#



ggsave("hyacinthus_dbRDA.pdf", plot = hy_plot,
       width = 6, height = 5, units = "in")

pdf("hyacinthus_DAPC.pdf", width = 6, height = 5)
scatter(hy_dapc)
dev.off()






















a.hyacinthus.summary.stand.df










#261026:
#stacked barchart
hy_microbe_raw <- microbe.df[
  rownames(hy_meta_aligned),
  ,
  drop = FALSE
]

hy_microbe_rel <- decostand(
  hy_microbe_raw,
  method = "total",
  MARGIN = 1
)
hy_microbe_rel <- as.data.frame(hy_microbe_rel)

hy_microbe_rel$Lineage <- hy_meta_aligned$Lineage

mean_abund <- colMeans(
  hy_microbe_rel[, colnames(hy_microbe_rel) != "Lineage"]
)

top_taxa <- names(
  sort(mean_abund, decreasing = TRUE)
)[1:12]


hy_top <- hy_microbe_rel[, top_taxa, drop = FALSE]

hy_top <- as.data.frame(hy_top)

hy_top$Other <- 1 - rowSums(hy_top)
hy_top$SampleID <- rownames(hy_top)
hy_top$Lineage <- hy_meta_aligned$Lineage

hy_long <- hy_top %>%
  pivot_longer(
    cols = -c(SampleID, Lineage),
    names_to = "Taxon",
    values_to = "RelAbund"
  )

library(scales)

ggplot(
  hy_long,
  aes(x = SampleID, y = RelAbund, fill = Taxon)
) +
  geom_col(width = 1) +
  facet_grid(~ Lineage, scales = "free_x", space = "free_x") +
  scale_y_continuous(labels = percent) +
  labs(
    x = "Coral colonies",
    y = "Relative bacterial abundance",
    fill = "Bacterial genus"
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.spacing.x = unit(0.3, "lines")
  )



























hyac_dbrda_final <- hy_plot +
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

ggsave(
  "hyac_dbRDA_final.pdf",
  plot = hyac_dbrda_final,
  width = 5,
  height = 5,
  units = "in"
)








#261104
# standardized PCA for Acropora hyacinthus

library(ggplot2)

# run PCA on corrected CLR matrix
hy_pca <- prcomp(hy_clr)

# extract sample scores
hy_pca_scores <- as.data.frame(hy_pca$x[, 1:2])

# add lineage metadata
hy_pca_scores$Lineage <- hy_meta$Lineage

# percent variance explained
hy_pca_var <- 100 * (hy_pca$sdev^2 / sum(hy_pca$sdev^2))

# make plot
hy_pca_plot <- ggplot(
  hy_pca_scores,
  aes(x = PC1, y = PC2, color = Lineage)
) +
  geom_point(size = 7, alpha = 0.85) +
  labs(
    title = "Acropora hyacinthus",
    x = paste0("PC1 (", round(hy_pca_var[1], 2), "%)"),
    y = paste0("PC2 (", round(hy_pca_var[2], 2), "%)")
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

# view plot
print(hy_pca_plot)

# export plot
ggsave(
  "hy_PCA_final.pdf",
  plot = hy_pca_plot,
  width = 5,
  height = 5,
  units = "in"
)








library(adegenet)

set.seed(123)

hy_cv_final <- xvalDapc(
  hy_clr,
  grp = hy_meta$Lineage,
  n.pca.max = 20,
  training.set = 0.7,
  result = "groupMean",
  center = TRUE,
  scale = FALSE,
  n.rep = 100,
  xval.plot = TRUE
)

# examine the cross-validation results
hy_cv_final

# inspect the selected number of PCs
hy_cv_final$`Number of PCs Achieving Highest Mean Success`
names(hy_cv_final)


hy_best_npc <- as.integer(
  hy_cv_final$`Number of PCs Achieving Highest Mean Success`
)

hy_dapc_final <- adegenet::dapc(
  hy_clr,
  grp = hy_meta$Lineage,
  n.pca = hy_best_npc,
  n.da = 3
)

#view ordination
scatter(hy_dapc_final, scree.da = FALSE)

#export figure
pdf("hy_DAPC_final.pdf", width = 5, height = 5)
scatter(hy_dapc_final, scree.da = FALSE)
dev.off()






library(dplyr)
library(tidyr)
library(ggplot2)
library(scatterpie)

# Calculate one geographic centroid per site
hy_site_locations <- hy_map_data %>%
  group_by(Site) %>%
  summarise(
    Longitude = mean(Longitude, na.rm = TRUE),
    Latitude = mean(Latitude, na.rm = TRUE),
    .groups = "drop"
  )

# Count lineages within each site
hy_site_lineages <- hy_map_data %>%
  count(Site, Lineage) %>%
  pivot_wider(
    names_from = Lineage,
    values_from = n,
    values_fill = 0
  )

# Combine geographic coordinates with lineage counts
hy_site_pies <- left_join(
  hy_site_locations,
  hy_site_lineages,
  by = "Site"
)

print(hy_site_pies)


hy_lineages <- c("HA", "HC", "HD", "HE")

hy_map_final <- ggplot() +
  geom_scatterpie(
    data = hy_site_pies,
    aes(
      x = Longitude,
      y = Latitude,
      r = 0.0008
    ),
    cols = hy_lineages,
    color = "black",
    linewidth = 0.3
  ) +
  geom_text(
    data = hy_site_pies,
    aes(
      x = Longitude,
      y = Latitude,
      label = paste("Site", Site)
    ),
    nudge_y = 0.0013,
    size = 4
  ) +
  coord_fixed(
    xlim = c(-169.664, -169.650),
    ylim = c(-14.187, -14.175),
    expand = FALSE
  ) +
  labs(
    title = "Acropora hyacinthus",
    x = "Longitude",
    y = "Latitude",
    fill = "Lineage"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(
      face = "italic",
      size = 13
    )
  )

print(hy_map_final)










library(sf)
library(ggplot2)
library(scatterpie)

# Folder containing the American Samoa shapefile
samoa_folder <- "/Users/karimprimov/Dropbox/UTAustin/2026/mcav.cryptic.lineages/261104.plots/ahyacinthus/geoBoundaries-ASM-ADM0-all"

# Find shapefile
shp_files <- list.files(
  samoa_folder,
  pattern = "\\.shp$",
  full.names = TRUE,
  recursive = TRUE
)

print(shp_files)

# Load shapefile
stopifnot(length(shp_files) > 0)

samoa_land <- st_read(shp_files[1], quiet = TRUE)

# Convert to longitude/latitude
samoa_land <- st_transform(samoa_land, 4326)

# Preview
plot(st_geometry(samoa_land))



hy_map_final <- ggplot() +
  
  # American Samoa landmass
  geom_sf(
    data = samoa_land,
    fill = "grey90",
    color = "grey40",
    linewidth = 0.4
  ) +
  
  # Lineage composition at each site
  geom_scatterpie(
    data = hy_site_pies,
    aes(
      x = Longitude,
      y = Latitude,
      r = 0.0008
    ),
    cols = c("HA", "HC", "HD", "HE"),
    color = "black",
    linewidth = 0.3
  ) +
  
  # Site labels
  geom_text(
    data = hy_site_pies,
    aes(
      x = Longitude,
      y = Latitude,
      label = paste("Site", Site)
    ),
    nudge_y = 0.0013,
    size = 4
  ) +
  
  # Local geographic extent
  coord_sf(
    xlim = c(-169.680, -169.640),
    ylim = c(-14.200, -14.165),
    expand = FALSE
  ) +
  
  labs(
    title = "Acropora hyacinthus",
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
    panel.grid = element_blank()
  )

print(hy_map_final)








































library(ggplot2)
library(dplyr)
library(scatterpie)
library(sf)

# 1. Isolate Tutuila from the American Samoa shapefile
#    using the island closest to your sampling sites

samoa_land <- st_transform(samoa_land, 4326)

island_parts <- st_cast(
  st_union(st_geometry(samoa_land)),
  "POLYGON"
)

island_sf <- st_sf(
  geometry = island_parts
)

site_center <- st_sfc(
  st_point(c(
    mean(hy_site_pies$Longitude),
    mean(hy_site_pies$Latitude)
  )),
  crs = 4326
)

# Select the island polygon nearest to the sample locations
nearest_island <- which.min(
  as.numeric(st_distance(
    st_transform(island_sf, 32702),
    st_transform(site_center, 32702)
  ))
)

tutuila <- island_sf[nearest_island, ]

# 2. Determine the extent of the entire island
bb <- st_bbox(tutuila)

map_width <- as.numeric(bb["xmax"] - bb["xmin"])
map_height <- as.numeric(bb["ymax"] - bb["ymin"])

# 3. Position large pies below the island
hy_site_pies_plot <- hy_site_pies %>%
  mutate(
    pie_lon = case_when(
      Site == 1 ~ as.numeric(bb["xmin"]) + 0.35 * map_width,
      Site == 2 ~ as.numeric(bb["xmin"]) + 0.70 * map_width
    ),
    pie_lat = as.numeric(bb["ymin"]) - 0.18 * map_height
  )

# Increase pie radius to remain visible on the full-island map
pie_radius <- 0.065 * map_width

# 4. Create full-island map
hy_map_final <- ggplot() +
  
  geom_sf(
    data = tutuila,
    fill = "grey90",
    color = "grey40",
    linewidth = 0.4
  ) +
  
  # Lines connecting the pies to actual sampling locations
  geom_segment(
    data = hy_site_pies_plot,
    aes(
      x = Longitude,
      y = Latitude,
      xend = pie_lon,
      yend = pie_lat
    ),
    color = "grey30",
    linewidth = 0.4
  ) +
  
  # Actual geographic locations
  geom_point(
    data = hy_site_pies_plot,
    aes(x = Longitude, y = Latitude),
    size = 2,
    color = "black"
  ) +
  
  # Large lineage pie charts
  geom_scatterpie(
    data = hy_site_pies_plot,
    aes(
      x = pie_lon,
      y = pie_lat,
      r = pie_radius
    ),
    cols = c("HA", "HC", "HD", "HE"),
    color = "black",
    linewidth = 0.3
  ) +
  
  # Pie chart labels
  geom_text(
    data = hy_site_pies_plot,
    aes(
      x = pie_lon,
      y = pie_lat,
      label = paste("Site", Site)
    ),
    nudge_y = pie_radius * 1.4,
    size = 4
  ) +
  
  coord_sf(
    xlim = c(
      as.numeric(bb["xmin"]) - 0.08 * map_width,
      as.numeric(bb["xmax"]) + 0.08 * map_width
    ),
    ylim = c(
      as.numeric(bb["ymin"]) - 0.36 * map_height,
      as.numeric(bb["ymax"]) + 0.10 * map_height
    ),
    expand = FALSE
  ) +
  
  labs(
    title = "Acropora hyacinthus",
    fill = "Lineage",
    x = "Longitude",
    y = "Latitude"
  ) +
  
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(
      face = "italic",
      size = 13
    ),
    panel.grid = element_blank()
  )

print(hy_map_final)





pie_lon = case_when(
  Site == 1 ~ as.numeric(bb["xmin"]) + 0.35 * map_width,
  Site == 2 ~ as.numeric(bb["xmin"]) + 0.70 * map_width
)


ggsave(
  "/Users/karimprimov/Dropbox/UTAustin/2026/mcav.cryptic.lineages/261104.plots/ahyacinthus/hyacinthus_sampling_map_final.pdf",
  plot = hy_map_final,
  width = 6,
  height = 5,
  units = "in"
)





























