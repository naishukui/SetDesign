# this code generates manhattan plots under different settings

library(qqman)
library(dplyr)
library(STAAR)
library(STAARpipeline)
library(ggplot2)

manhattan <- qqman::manhattan
fix(manhattan)

#genes <- genes_info
genes<-get(load("/Users/nkui/Library/CloudStorage/OneDrive-InsideMDAnderson(2)/STAAR/output/combined/gene.Rdata"))
#final setting
setting1<-get(load("/Users/nkui/Library/CloudStorage/OneDrive-InsideMDAnderson(2)/STAAR/output/output6-10/coding_combined6.Rdata"))
setting2<-get(load("/Users/nkui/Library/CloudStorage/OneDrive-InsideMDAnderson(2)/STAAR/output/output6-10/coding_combined7.Rdata"))
setting3<-get(load("/Users/nkui/Library/CloudStorage/OneDrive-InsideMDAnderson(2)/STAAR/output/output6-10/coding_combined8.Rdata"))
setting4<-get(load("/Users/nkui/Library/CloudStorage/OneDrive-InsideMDAnderson(2)/STAAR/output/output6-10/coding_combined9.Rdata"))
#original STAAR pipeline setting with splitted data
setting5<-get(load("/Users/nkui/Library/CloudStorage/OneDrive-InsideMDAnderson(2)/STAAR/output/output6-10/coding_combined10.Rdata"))
#original STAAR pipeline setting with unsplitted data
setting6<-get(load("/Users/nkui/Library/CloudStorage/OneDrive-InsideMDAnderson(2)/STAAR/output/output6-10/coding_combined11.Rdata"))


colnames(genes) <- c("gene", "chr", "start_position", "end_position")

#-------------------------------------------manhattan plot for setting 1 and setting 6-----------------------------------------------#

combo <- merge(setting1, genes, by = c("chr", "gene"))
combo$chr <- as.numeric(as.character(combo$chr))  # Ensure chr is numeric
pskat <- as.numeric(combo[, 8])  # Convert p-values to numeric
gene <- as.character(combo[, 2])  # Keep gene names as character
start_position <- as.numeric(combo[, 9])  # Convert start position to numeric
allGeneP <- data.frame(chr = combo$chr, pskat = pskat, start_position = start_position, gene = gene)
bychr <- allGeneP[order(allGeneP$chr), ]  # Order by chromosome
df <- bychr %>% group_by(chr) %>% mutate(counter = row_number())
df <- as.data.frame(df)

combo2 <- merge(setting6, genes, by = c("chr", "gene"))
pskat2 <- as.numeric(unlist(combo2[, 8]))
chr2 <- as.numeric(unlist(combo2[, 1]))
gene2 <- as.character(unlist(combo2[, 2]))
start_position2 <- as.numeric(unlist(combo2[, 9]))
allGeneP2 <- data.frame(chr2 = chr2, pskat2 = pskat2, start_position2 = start_position2, gene2 = gene2)
bychr2 <- allGeneP2[order(allGeneP2$chr2), ]  # Order by chromosome
df2 <- bychr2 %>% group_by(chr2) %>% mutate(counter = row_number())
df2 <- as.data.frame(df2)

par(mfrow = c(2, 1))
custom_colors <- c("#74c2e7","#dea3e5")#"#9bbae8", "#f8c67d",)
#Manhattan plot showing associations of gene-centric analysis for pancreatic cancer versus −log10(P) of SKAT
#（n=7435）#

manhattan2(df,
          chr = "chr",
          bp = "counter",
          snp = "gene",
          p = "pskat",
          col = custom_colors,  # Use custom colors
          #suggestiveline = -log10(1e-5),  # Optional: significance threshold
          genomewideline = -log10(5e-8),  # Optional: genome-wide significance threshold
          annotatePval = 0.001,
          ylim = c(0, 6),
          main = "Setting 1")  # Optional: add title

manhattan2(df2,
          chr = "chr2",
          bp = "counter",
          snp = "gene2",
          p = "pskat2",
          col = custom_colors,
          #suggestiveline = -log10(1e-5),
          genomewideline = -log10(5e-8),
          annotatePval = 0.001,
          ylim = c(0, 6),
          main = "Setting 6")

outlier_idx <- which(-log10(combo2$SKAT1_1) > 10)
outlier_x <- which(df2$gene2=="ULK1")
outlier_y <- -log10(combo2$SKAT1_1)[outlier_idx]
outlier_snp <- "ULK1"

points(outlier_x, rep(5.5, length(outlier_x)), pch = 17, col = "red", cex = 1)
text(outlier_x, rep(5.5, length(outlier_x)), labels = paste0(outlier_snp, " (~", round(outlier_y, 2), ")"), pos = 4, col = "red", cex = 1)

#-------------------------------------------manhattan plot for setting 2 to setting 5-----------------------------------------------#

# Create a named list for settings 1 to 6
settings_list <- list(
  # "Setting 1" = setting1,
  "Setting 2" = setting2,
  "Setting 3" = setting3,
  "Setting 4" = setting4,
  "Setting 5" = setting5
  #  "Setting 6" = setting6,

)

# Custom colors for Manhattan plots
custom_colors <- c("#74c2e7", "#dea3e5")

# Set up plotting layout: 2 rows x 1 columns (adjust as needed)
par(mfrow = c(4,1))

# Loop over each setting to create a Manhattan plot
for (name in names(settings_list)) {
  current_setting <- settings_list[[name]]

  # Merge the current setting with gene data using "chr" and "gene" columns
  combo <- merge(current_setting, genes, by = c("chr", "gene"))

  # Ensure chromosome is numeric
  combo$chr <- as.numeric(as.character(combo$chr))

  # Extract values (adjust column indices if needed)
  pskat         <- as.numeric(combo[, 8])    # p-value column (assumed at col 8)
  gene_names    <- as.character(combo[, 2])  # gene names (assumed at col 2)
  start_position<- as.numeric(combo[, 9])    # start position (assumed at col 9)

  # Build a data frame for Manhattan plot
  allGeneP <- data.frame(chr = combo$chr,
                         pskat = pskat,
                         start_position = start_position,
                         gene = gene_names)

  # Order by chromosome and add a counter to simulate base-pair position
  bychr <- allGeneP[order(allGeneP$chr), ]
  df <- bychr %>%
    group_by(chr) %>%
    mutate(counter = row_number()) %>%
    as.data.frame()

  # Create Manhattan plot for this setting
  manhattan(df,
            chr = "chr",
            bp = "counter",
            snp = "gene",
            p = "pskat",
            col = custom_colors,
            genomewideline = -log10(5e-8),
            annotatePval = 0.001,
            cex.annotate = 3,
            ylim = c(0, 6),
            main = name)
}



