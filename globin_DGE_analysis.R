
#prepare gene annotations
#define what mart and dataset to use
human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
listFilters(human)
ensembl_dataset = useDataset("hsapiens_gene_ensembl", mart = human)
genome <- listDatasets(ensembl_dataset)[(listDatasets(ensembl_dataset)=="hsapiens_gene_ensembl"),]$version
attributes = listAttributes(ensembl_dataset)
ensembl = getBM(attributes = c("ensembl_gene_id","external_gene_name","entrezgene_id","transcript_length","gene_biotype","percentage_gene_gc_content","chromosome_name","start_position","end_position"),
                mart = ensembl_dataset)
ensembl = dplyr::rename(ensembl, gene_id = ensembl_gene_id, gene_name = external_gene_name, entrez = entrezgene_id, length = transcript_length, biotype = gene_biotype, GC_content = percentage_gene_gc_content, Chr = chromosome_name, GeneStart = start_position, GeneEnd = end_position)

#create edgeR object on merged count matrix (pre-CombatSeq)
dge <- DGEList(fc, group = meta$LN.Lung)
dim(dge)
#60664 288

#add gene info to object
gene_id <- rownames(dge)
my_genes <- ensembl[match(gene_id, ensembl$gene_id),]
dge$genes <- my_genes

#barplot of raw library sizes
globin_cols <- c(rep("blue",131),rep("red",65),rep("green",64),rep("red",30))
par(mfrow = c(1,1))
barplot(dge$samples$lib.size, names.arg = rownames(dge$samples), col = globin_cols, las = 2, cex.names = 0.5)
abline(h=2e07, col="blue")
title("Barplot of library sizes")

#define globin genes and remove from edgeR object
#pc_genes <- rownames(dge$counts[grep("protein_coding",dge$genes$biotype),])
pre_gd_counts <- dge$counts
pre_gd_libs <- dge$samples$lib.size
names(pre_gd_libs) <- colnames(dge)
dge <- dge[-grep("HBA1|HBA2|HBB|HBBP1|HBD|HBE1|HBG1|HBG2|HBM|HBQ1|HBZ|HBZP1",dge$genes$gene_name),, keep.lib.sizes=FALSE]
dim(dge)
#60641 288
gd_counts <- dge$counts
gd_libs <- dge$samples$lib.size
names(gd_libs) <- colnames(dge)

#barplot of library sizes post-bioinformatic globin-depletion
par(mfrow = c(1,1))
barplot(dge$samples$lib.size, names.arg = rownames(dge$samples), col = globin_cols, las = 2, cex.names = 0.5)
abline(h=2e07, col="blue")
title("Barplot of library sizes (post-globin-depletion)")

#check distribution of counts
lcpm <- cpm(dge, log = T)
par(mfrow = c(1,1))
boxplot(lcpm, xlab = "", ylab = "Log2 counts per million", las = 2, cex.axis = 0.5, outline = F)
abline(h = median(lcpm), col = "blue")
title("Boxplots of log-CPMs (unnormalised)")

#perform CombatSeq to neutralise batch effects
#combatseq_sub <- sva::ComBat_seq(counts=gd_counts, batch=meta$Globin, group=meta$LN.Lung,
                                  #shrink=FALSE, shrink.disp=FALSE)

#compare pre- vs post-CombatSeq by PCA
library(factoextra)
raw.res.pca <- prcomp(t(cpm(gd_counts, log = T)))
colors <- as.factor(meta$Globin)
factoextra::fviz_pca_ind(raw.res.pca,
             col.ind = colors,
             palette = "Set1",
             repel = TRUE,    # Avoid text overlapping
             title = "PCA plot of pre-CombatSeq data"
)
combat.res.pca <- prcomp(t(cpm(combatseq_sub, log = T)))
factoextra::fviz_pca_ind(combat.res.pca,
                         axes = c(1,2),
                         col.ind = colors,
                         palette = "Set1",
                         repel = TRUE,    # Avoid text overlapping
                         title = "PCA plot of post-CombatSeq data"
)
factoextra::fviz_pca_ind(combat.res.pca,
                         axes = c(2,3),
                         col.ind = colors,
                         palette = "Set1",
                         repel = TRUE,    # Avoid text overlapping
                         title = "PCA plot of post-CombatSeq data"
)
factoextra::fviz_pca_ind(combat.res.pca,
                         axes = c(3,4),
                         col.ind = colors,
                         palette = "Set1",
                         repel = TRUE,    # Avoid text overlapping
                         title = "PCA plot of post-CombatSeq data"
)

#re-create edgeR object using corrected data
#dge <- DGEList(combatseq_sub, group = meta$Globin)

#normalise
dge <- calcNormFactors(dge, method = "TMM")

#create summarized experiment for later
gi2gn <- function(x){
  gene_id <- rownames(x)
  genes <- ensembl[match(gene_id, ensembl$gene_id),]
  genes <- genes[!duplicated(genes$gene_name),]
  genes <- genes[!is.na(genes$gene_name),]
  x <- x[rownames(x) %in% genes$gene_id,]
  genes <- genes[genes$gene_id %in% rownames(x),]
  rownames(x) <- genes$gene_name
  x
}
library(SummarizedExperiment)
se_dge <- SummarizedExperiment(assays = list(counts = as.matrix(gi2gn(dge$counts))), colData = dge$samples)
library(TBSignatureProfiler)
se_dge <- mkAssay(se_dge, input_name = "counts", log = TRUE, counts_to_CPM = TRUE)

#check distribution of counts post-normalisation
lcpm <- cpm(dge, log = T)
par(mfrow = c(1,1))
boxplot(lcpm, xlab = "", ylab = "Log2 counts per million", las = 2, cex.axis = 0.5, outline = F)
abline(h = median(lcpm), col = "blue")
title("Boxplots of log-CPMs (TMM-normalised)")

#check clustering of samples post-normalisation
tmm <- cpm(dge)
tmm.res.pca <- prcomp(t(tmm))
colors <- as.factor(meta$Globin)
factoextra::fviz_pca_ind(tmm.res.pca,
                         axes = c(1,2),
                         col.ind = colors,
                         palette = "Set1",
                         repel = TRUE,    # Avoid text overlapping
                         title = "PCA plot of post-CombatSeq data"
)

#restrict to samples to those from participants present in old and new batches
meta_old <- meta[c(1:131),]
meta_new <- meta[c(132:nrow(meta)),]
meta_both <- meta[meta$PID %in% intersect(meta_old$PID,meta_new$PID),]
dge_both <- dge[,rownames(meta_both$PID)]

#check clustering based on participant ID
tmm_both <- cpm(dge_both)
tmm.both.pca <- prcomp(t(tmm_both))
colors <- as.factor(meta_both$PID)
c25 <- c(
  "dodgerblue2", "#E31A1C", # red
  "green4",
  "#6A3D9A", # purple
  "#FF7F00", # orange
  "black", "gold1",
  "skyblue2", "#FB9A99", # lt pink
  "palegreen2",
  "#CAB2D6", # lt purple
  "#FDBF6F", # lt orange
  "gray70", "khaki2",
  "maroon", "orchid1", "deeppink1", "blue1", "steelblue4",
  "darkturquoise", "green1", "yellow4", "yellow3",
  "darkorange4", "brown"
)
factoextra::fviz_pca_ind(tmm.both.pca,
                         axes = c(1,2),
                         col = c25,
                         col.ind = colors,
                         geom = "text",
                         repel = TRUE,    # Avoid text overlapping
                         title = "PCA plot of post-CombatSeq data"
)

#test expression across duplicates (globin depleted vs non-depleted)
dups <- c("PCM035","PCM071","PCM073","PCM074","PCM075","PCM076","PCM077","PCM078","PCM079","PCM080","PCM081","PCM083","PCM084","PCM085","PCM086","PCM087","PCM088","PCM089","PCM090","PCM091","PCM092","PCM093","PCM094","PCM095","PCM096","PCM097","PCM098","PCM099","PCM100","PCM035_dup","PCM071_dup","PCM073_dup","PCM074_dup","PCM075_dup","PCM076_dup","PCM077_dup","PCM078_dup","PCM079_dup","PCM080_dup","PCM081_dup","PCM083_dup","PCM084_dup","PCM085_dup","PCM086_dup","PCM087_dup","PCM088_dup","PCM089_dup","PCM090_dup","PCM091_dup","PCM092_dup","PCM093_dup","PCM094_dup","PCM095_dup","PCM096_dup","PCM097_dup","PCM098_dup","PCM099_dup","PCM100_dup")
rownames(meta_new) <- meta_new$RNA.seq.ID
meta_dup <- meta_new[dups,]
meta_dup$PCM <- rownames(meta_dup)
meta_dup$PCM <- gsub("_dup", "", meta_dup$PCM)
dge_dups <- dge[,rownames(meta_dup)]
meta_dup
meta_dup$PID

#compare QC stats
glob_qc_stats <- data.frame(sample = c("PCM035","PCM035","PCM071","PCM071","PCM073","PCM073","PCM074","PCM074","PCM075","PCM075","PCM076","PCM076","PCM077","PCM077","PCM078","PCM078","PCM079","PCM079","PCM080","PCM080","PCM081","PCM081","PCM083","PCM083","PCM084","PCM084","PCM085","PCM085","PCM086","PCM086","PCM087","PCM087","PCM088","PCM088","PCM089","PCM089","PCM090","PCM090","PCM091","PCM091","PCM092","PCM092","PCM093","PCM093","PCM094","PCM094","PCM095","PCM095","PCM096","PCM096","PCM097","PCM097","PCM098","PCM098","PCM099","PCM099","PCM0100","PCM0100","PCM035_dup","PCM035_dup","PCM071_dup","PCM071_dup","PCM073_dup","PCM073_dup","PCM074_dup","PCM074_dup","PCM075_dup","PCM075_dup","PCM076_dup","PCM076_dup","PCM077_dup","PCM077_dup","PCM078_dup","PCM078_dup","PCM079_dup","PCM079_dup","PCM080_dup","PCM080_dup","PCM081_dup","PCM081_dup","PCM083_dup","PCM083_dup","PCM084_dup","PCM084_dup","PCM085_dup","PCM085_dup","PCM086_dup","PCM086_dup","PCM087_dup","PCM087_dup","PCM088_dup","PCM088_dup","PCM089_dup","PCM089_dup","PCM090_dup","PCM090_dup","PCM091_dup","PCM091_dup","PCM092_dup","PCM092_dup","PCM093_dup","PCM093_dup","PCM094_dup","PCM094_dup","PCM095_dup","PCM095_dup","PCM096_dup","PCM096_dup","PCM097_dup","PCM097_dup","PCM098_dup","PCM098_dup","PCM099_dup","PCM099_dup","PCM0100_dup","PCM0100_dup"), globin = c(rep("Globin-depleted",58),rep("Non-depleted",58)), read = rep(c("F","R"),58), dups = c(62,60.6,61,59.4,57,56.1,59.6,60.5,58.7,59.1,59.1,58.1,58.1,57.4,57.9,58.3,47.4,47.7,52.9,52.9,48.5,48.3,56.4,56.5,53.1,53.7,52.7,52.7,57.8,57.6,48.8,49.1,52.4,52.7,53.6,53.9,63.1,63.4,50.9,51.2,53.3,53.2,51.5,52,52,52.6,49.6,50.7,49.9,50.2,54.1,54.6,53.6,52.5,60.8,59.4,55.4,52.4,70.4,70,83,63.1,70.5,50.7,76,57.3,83.7,63.1,67.2,46.4,67.7,66.8,80.4,78.9,67,65.2,73.5,72.1,66.8,65.4,72.4,71,79,77.4,76.2,74.5,79.4,77.1,65.6,63.6,68.8,67,72.7,70.9,72.5,71.9,76.7,74.4,74.6,72.7,63.6,62.7,71.9,69.9,70.7,68.5,70.9,69.2,72.3,70.7,71.7,71,66.5,64.6,69.5,58.4))
glob_qc_stats_melt <- reshape2::melt(glob_qc_stats)
colnames(glob_qc_stats_melt)[2] <- "Globin"

ggplot(glob_qc_stats_melt, aes(x = Globin, y = value, group = Globin)) +
  geom_boxplot(aes(fill = Globin), alpha = 0.7, outlier.shape = NA) +
  scale_fill_manual(values = c("#4DAF4A","#377EB8")) +
  geom_point(position = position_dodge(0.2)) +
  theme_bw() +
  xlab("Globin status") +
  ylab("% Duplicates") +
  stat_compare_means(aes(label = "..p.adj.."), method = "t.test", label = "p.signif", paired = T, label.y = 80) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_text(colour = "black", size = 12), axis.title.y = element_text(size = 18, face = "bold"), axis.title.x = element_text(size = 18, face = "bold"), legend.text = element_text(size = 14), legend.title = element_text(face = "bold", size = 16), strip.text.x = element_text(size = 14))

ggplot(glob_qc_stats_melt, aes(x = Globin, y = value, group = Globin)) +
  geom_boxplot(aes(fill = Globin), alpha = 0.7, outlier.shape = NA) +
  scale_fill_manual(values = c("#4DAF4A","#377EB8")) +
  geom_point(position = position_dodge(0.2)) +
  facet_wrap(~read) +
  theme_bw() +
  xlab("Globin status") +
  ylab("% Duplicates") +
  stat_compare_means(aes(label = "..p.adj.."), method = "t.test", label = "p.signif", paired = T, label.y = 80) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_text(colour = "black", size = 12), axis.title.y = element_text(size = 18, face = "bold"), axis.title.x = element_text(size = 18, face = "bold"), legend.text = element_text(size = 14), legend.title = element_text(face = "bold", size = 16), legend.position = "bottom", strip.text.x = element_text(size = 14))

nd_f_read <- glob_qc_stats_melt[glob_qc_stats_melt$Globin=="Non-depleted" & glob_qc_stats_melt$read=="F",]
mean(nd_f_read$value)
nd_r_read <- glob_qc_stats_melt[glob_qc_stats_melt$Globin=="Non-depleted" & glob_qc_stats_melt$read=="R",]
mean(nd_r_read$value)

#run TBSignatureProfiler
library(SummarizedExperiment)
library(TBSignatureProfiler)
library(ggpubr)
library(dplyr)
dup_counts <- as.data.frame(dge_dups$counts)
gene_id <- rownames(dup_counts)
my_genes <- ensembl[match(gene_id, ensembl$gene_id),]
my_genes <- my_genes[!duplicated(my_genes$gene_name),]
dup_counts <- dup_counts[rownames(dup_counts) %in% my_genes$gene_id,]
my_genes <- my_genes[my_genes$gene_id %in% rownames(dup_counts),]
dup_counts$gene_names <- my_genes$gene_name
dup_counts <- as.data.frame(dup_counts %>% group_by(gene_names) %>% summarise_if(is.numeric, sum))
dup_counts <- dup_counts[c(2:nrow(dup_counts)),]
rownames(dup_counts) <- dup_counts$gene_names
dup_counts <- dup_counts[2:nrow(dup_counts),2:ncol(dup_counts)]
se_dups <- SummarizedExperiment(assays = list(counts = as.matrix(dup_counts)), colData = meta_dup)
se_dups <- mkAssay(se_dups, input_name = "counts", log = TRUE, counts_to_CPM = TRUE)

#run scoring algorithm
names(TBsignatures)
TBsignatures_orig <- TBsignatures[c(1,2,3,4,6,5,7,16,15,17,24,25,32,34,35,36,39,41,43,44,48,49,50,51,52,54,55,57,58,59,60,61,65,67,66,68)]
names(TBsignatures_orig)
scores_GSVA <- runTBsigProfiler(se_dups, useAssay = "log_counts_cpm", signatures = TBsignatures, algorithm = "GSVA")
signatureBoxplot(inputData = scores_GSVA, 
                 name = "GSVA",
                 signatureColNames = names(TBsignatures),
                 annotationColName = "Globin", includePoints = TRUE) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

cal_z_score <- function(x){
  (x - mean(x)) / sd(x)
}

gsva_df <- as.data.frame(colData(scores_GSVA)[,c(5,109:ncol(colData(scores_GSVA)))])
gsva_glob <- as.data.frame(t(gsva_df[c(2:30),])[2:71,])
gsva_glob <- gsva_glob[,c(29,1:28)]
gsva_glob[] <- lapply(gsva_glob, as.numeric)
gsva_nonglob <- as.data.frame(t(gsva_df[c(1,31:nrow(gsva_df)),])[2:71,])
gsva_nonglob[] <- lapply(gsva_nonglob, as.numeric)

gsva_df <- as.data.frame(colData(scores_GSVA)[,c(5,109:ncol(colData(scores_GSVA)))])
gsva_glob <- as.data.frame(t(gsva_df[c(2:30),])[2:37,])
gsva_glob <- gsva_glob[,c(29,1:28)]
gsva_glob[] <- lapply(gsva_glob, as.numeric)
gsva_nonglob <- as.data.frame(t(gsva_df[c(1,31:nrow(gsva_df)),])[2:37,])
gsva_nonglob[] <- lapply(gsva_nonglob, as.numeric)

#correlate globin-depleted and non-depleted GSVA results
cor_mat <- cor(gsva_glob,gsva_nonglob)
round_df <- function(x, digits) {
  # round all numeric variables
  # x: data frame 
  # digits: number of digits to round
  numeric_columns <- sapply(x, mode) == 'numeric'
  x[numeric_columns] <-  round(x[numeric_columns], digits)
  x
}

cor_mat <- round_df(cor_mat, 2)

melted_cor_mat <- melt(cor_mat)

#get lower triangle of the correlation matrix
get_lower_tri<-function(cormat){
  cormat[upper.tri(cormat)] <- NA
  return(cormat)
}
lower_tri <- get_lower_tri(cor_mat)

melted_cor_mat <- melt(lower_tri, na.rm = TRUE)

cor_pairs <- melted_cor_mat[c(1,30,58,85,111,136,160,183,205,226,246,265,283,300,316,331,345,358,370,381,391,400,408,415,421,426,430,433,435),]

#plot correlation scores for full signature list
ggheatmap <- ggplot(cor_pairs, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0.5, limit = c(0,1), space = "Lab", 
                       name="Pearson\nCorrelation") +
  theme_minimal()+ # minimal theme
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1),
        axis.text.y = element_text(size = 12))+
  coord_fixed()
ggheatmap + 
  geom_text(aes(Var2, Var1, label = value), color = "black", size = 4) +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank(),
    legend.justification = c(1, 0),
    legend.position = c(0.6, 0.7),
    legend.direction = "horizontal")+
  guides(fill = guide_colorbar(barwidth = 7, barheight = 1,
                               title.position = "top", title.hjust = 0.5))

#cor_pairs_D0 <- cor_pairs[grep("PCM078_dup|PCM080_dup|PCM084_dup|PCM086_dup|PCM088_dup|PCM090_dup|PCM092_dup|PCM094_dup|PCM096_dup|PCM098_dup|PCM100_dup",cor_pairs$Var2),]
cor_pairs_D0 <- cor_pairs[grep("PCM074_dup|PCM076_dup|PCM078_dup|PCM080_dup|PCM084_dup|PCM086_dup|PCM088_dup|PCM090_dup|PCM092_dup|PCM094_dup|PCM096_dup|PCM098_dup|PCM100_dup",cor_pairs$Var2),]
min(cor_pairs_D0$value)
max(cor_pairs_D0$value)
median(cor_pairs_D0$value)
cor_pairs_M6 <- cor_pairs[-grep("PCM075_dup|PCM077_dup|PCM078_dup|PCM080_dup|PCM084_dup|PCM086_dup|PCM088_dup|PCM090_dup|PCM092_dup|PCM094_dup|PCM096_dup|PCM098_dup|PCM100_dup",cor_pairs$Var2),]
min(cor_pairs_M6$value)
max(cor_pairs_M6$value)
median(cor_pairs_M6$value)

#assess effect of globins on libraries
libs_df <- cbind(as.data.frame(pre_gd_libs),as.data.frame(gd_libs))
libs_df$glob_count <- (libs_df$pre_gd_libs - libs_df$gd_libs)
libs_df$prop_glob <- (libs_df$glob_count/libs_df$pre_gd_libs)*100
libs_df$Var1 <- rownames(libs_df)

#restrict to duplicates to plot globin read proportions first
libs_df_dups <- libs_df[dups,]
libs_df_dups <- libs_df_dups[c(30,2:29,1,31:58),]
libs_df_dups$sample <- dups
libs_df_dups$Globin <- c(rep("Globin-depleted",29),rep("Non-depleted",29))
rownames(libs_df_dups) <- c()
libs_df_dups_prop <- libs_df_dups[,c(6,7,4)]
libs_df_dups <- libs_df_dups[,c(6,7,1,2)]
colnames(libs_df_dups) <- c("sample","Globin","Pre-bioinformatic globin-depletion","Post-bioinformatic globin-depletion")
libs_df_dups_prop_melt <- reshape2::melt(libs_df_dups_prop)
libs_df_dups_melt <- reshape2::melt(libs_df_dups)
sum(libs_df_dups_melt[libs_df_dups_melt$Globin=="Globin-depleted" & libs_df_dups_melt$variable=="Post-bioinformatic globin-depletion",]$value)
gd_pre_gd <- libs_df_dups_melt[libs_df_dups_melt$Globin=="Globin-depleted" & libs_df_dups_melt$variable=="Pre-bioinformatic globin-depletion",]
rownames(gd_pre_gd) <- gd_pre_gd$sample
sum(gd_pre_gd[colnames(dge_dups_D0)[1:11],]$value)
sum(gd_pre_gd[colnames(dge_dups_M6)[1:11],]$value)
gd_prop <- libs_df_dups_prop_melt[libs_df_dups_prop_melt$Globin=="Globin-depleted",]
min(gd_prop$value)
max(gd_prop$value)
median(gd_prop$value)
sum(libs_df_dups_melt[libs_df_dups_melt$Globin=="Non-depleted" & libs_df_dups_melt$variable=="Post-bioinformatic globin-depletion",]$value)
nd_pre_gd <- libs_df_dups_melt[libs_df_dups_melt$Globin=="Non-depleted" & libs_df_dups_melt$variable=="Pre-bioinformatic globin-depletion",]
rownames(nd_pre_gd) <- nd_pre_gd$sample
sum(nd_pre_gd[colnames(dge_dups_D0)[1:11],]$value)
sum(nd_pre_gd[colnames(dge_dups_M6)[1:11],]$value)
nd_prop <- libs_df_dups_prop_melt[libs_df_dups_prop_melt$Globin=="Non-depleted",]
min(nd_prop$value)
max(nd_prop$value)
median(nd_prop$value)

libs_df_dups_prop_melt$Hb <- meta_dup$Hb
libs_df_dups_prop_melt$sex <- meta_dup$Sex
libs_df_dups_prop_melt$visit <- meta_dup$Visit
nd_libs_df_dups_prop_melt <- libs_df_dups_prop_melt[libs_df_dups_prop_melt$Globin=="Non-depleted",]
nd_libs_df_dups_prop_melt <- nd_libs_df_dups_prop_melt[grep("D0|M6",nd_libs_df_dups_prop_melt$visit),]

ggplot(nd_libs_df_dups_prop_melt, aes(x = visit, y = value)) +
  geom_boxplot(aes(fill = visit), alpha = 0.7, outlier.shape = NA) +
  scale_fill_manual(values = c("#4DAF4A","#377EB8")) +
  geom_point(position = position_dodge(0.2), aes(colour = sex)) +
  theme_bw() +
  xlab("Visit") +
  ylab("Reads mapped to globin genes (%)") +
  stat_compare_means(aes(label = "..p.adj.."), method = "t.test", label = "p.signif", paired = T, label.y = 80) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_text(colour = "black", size = 12), axis.title.y = element_text(size = 16, face = "bold"), axis.title.x = element_text(size = 18, face = "bold"), strip.text.x = element_text(size = 14))

nd_libs_df_dups_prop_melt2 <- nd_libs_df_dups_prop_melt[-grep("N/A",nd_libs_df_dups_prop_melt$Hb),]
nd_libs_df_dups_prop_melt2$Hb <- as.numeric(nd_libs_df_dups_prop_melt2$Hb)
ggplot(nd_libs_df_dups_prop_melt2, aes(x = visit, y = Hb)) +
  geom_boxplot(aes(fill = visit), alpha = 0.7, outlier.shape = NA) +
  scale_fill_manual(values = c("#4DAF4A","#377EB8")) +
  geom_point(position = position_dodge(0.2), aes(colour = sex)) +
  theme_bw() +
  xlab("Visit") +
  ylab("Haemoglobin count (g/dL)") +
  stat_compare_means(aes(label = "..p.adj.."), method = "t.test", label = "p.signif", paired = F) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_text(colour = "black", size = 12), axis.title.y = element_text(size = 16, face = "bold"), axis.title.x = element_text(size = 18, face = "bold"), strip.text.x = element_text(size = 14))

ggplot(libs_df_dups_prop_melt, aes(x = Globin, y = value, group = Globin)) +
  geom_boxplot(aes(fill = Globin), alpha = 0.7, outlier.shape = NA) +
  scale_fill_manual(values = c("#4DAF4A","#377EB8")) +
  geom_point(position = position_dodge(0.2)) +
  theme_bw() +
  xlab("Globin status") +
  ylab("Reads mapped to globin genes (%)") +
  stat_compare_means(aes(label = "..p.adj.."), method = "t.test", label = "p.signif", paired = T, label.y = 80) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_text(colour = "black", size = 12), axis.title.y = element_text(size = 16, face = "bold"), axis.title.x = element_text(size = 18, face = "bold"), legend.text = element_text(size = 14), legend.title = element_text(face = "bold", size = 16), legend.position = "bottom", strip.text.x = element_text(size = 14))

compare_means(value ~ variable, data = libs_df_dups_melt)

ggplot(libs_df_dups_melt, aes(x = sample, y = value, fill = Globin, group = Globin)) +
  geom_bar(stat = "identity") +
  facet_wrap(~variable, nrow = 2, ncol = 1) +
  scale_fill_manual(values = c("#4DAF4A","#377EB8")) +
  theme_bw() +
  xlab("Globin status") +
  ylab("Library size") +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_text(colour = "black", size = 12), axis.title.y = element_text(size = 18, face = "bold"), axis.title.x = element_blank(), legend.text = element_text(size = 14), legend.title = element_text(face = "bold", size = 16), legend.position = "bottom", strip.text.x = element_text(size = 14))

#assess correlation between proportion of globins and value from correlations between pairs
libs_df_dups <- libs_df[dups,]
colnames(libs_df_dups)[5] <- "Var2"
dup_cor_lib <- merge(cor_pairs, libs_df_dups, by = "Var2")
dup_cor_lib$sex <- meta_dup[dup_cor_lib$Var1,]$Sex
dup_cor_lib$sex <- factor(dup_cor_lib$sex)

meta_dup_cor <- meta_dup[levels(dup_cor_lib$Var2),]
meta_dup_cor <- meta_dup_cor[,c(2,29)]
meta_dup_cor$Var2 <- rownames(meta_dup_cor)
rownames(dup_cor_lib) <- dup_cor_lib$Var2
glob_cor <- merge(dup_cor_lib, meta_dup_cor, by = "Var2")
glob_cor <- glob_cor[c(2:9,11:18,20,23,24,26:29),]

colnames(glob_cor)[8] <- "Sex"
glob_cor$Hb <- as.numeric(glob_cor$Hb)
ggplot(glob_cor, aes(x = Sex, y = Hb, group = Sex)) +
  geom_boxplot(aes(fill = Sex), alpha = 0.5, outlier.shape = NA) +
  stat_compare_means(aes(label = "..p.adj.."), method = "t.test", label = "p.signif", ref.group = "Female", hide.ns = F) +
  geom_jitter(aes(colour = Sex), width = 0.15, size = 2) +
  ylab("Whole blood haemoglobin count (g/dL)") +
  theme_bw() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_text(colour = "black", size = 12), axis.title.y = element_text(size = 18, face = "bold"), axis.title.x = element_blank(), legend.text = element_text(size = 14), legend.title = element_text(face = "bold", size = 16), legend.position = "bottom", strip.text.x = element_text(size = 14))

glob_cor[glob_cor$Sex=="Male",]

glob_cor$Hb <- as.numeric(glob_cor$Hb)
compare_means(Hb ~ Sex, data = glob_cor)

colnames(dup_cor_lib)[8] <- "Sex"
ggplot(dup_cor_lib, aes(x = Sex, y = prop_glob, group = Sex)) +
  geom_boxplot(aes(fill = Sex), alpha = 0.5, outlier.shape = NA) +
  stat_compare_means(aes(label = "..p.adj.."), method = "t.test", label = "p.signif", ref.group = "Female", hide.ns = T) +
  geom_jitter(aes(colour = Sex), width = 0.15, size = 2) +
  ylab("Percentage of reads mapping to haemoglobin genes") +
  theme_bw() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_text(colour = "black", size = 12), axis.title.y = element_text(size = 18, face = "bold"), axis.title.x = element_blank(), legend.text = element_text(size = 14), legend.title = element_text(face = "bold", size = 16), legend.position = "bottom", strip.text.x = element_text(size = 14))

dup_cor_lib_D0 <- dup_cor_lib[grep("PCM074|PCM076|PCM078|PCM080|PCM084|PCM086|PCM088|PCM090|PCM092|PCM094|PCM096|PCM098|PCM100|PCM074_dup|PCM076_dup|PCM078_dup|PCM080_dup|PCM084_dup|PCM086_dup|PCM088_dup|PCM090_dup|PCM092_dup|PCM094_dup|PCM096_dup|PCM098_dup|PCM100_dup",rownames(dup_cor_lib)),]
ggplot(dup_cor_lib_D0, aes(x = Sex, y = prop_glob, group = Sex)) +
  geom_boxplot(aes(fill = Sex), alpha = 0.5, outlier.shape = NA) +
  stat_compare_means(aes(label = "..p.adj.."), method = "t.test", label = "p.signif", ref.group = "Female", hide.ns = T) +
  geom_jitter(aes(colour = Sex), width = 0.15, size = 2) +
  ylab("Percentage of reads mapping to haemoglobin genes") +
  theme_bw() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_text(colour = "black", size = 12), axis.title.y = element_text(size = 18, face = "bold"), axis.title.x = element_blank(), legend.text = element_text(size = 14), legend.title = element_text(face = "bold", size = 16), legend.position = "bottom", strip.text.x = element_text(size = 14))

dup_cor_lib_M6 <- dup_cor_lib[grep("PCM075|PCM077|PCM079|PCM081|PCM083|PCM085|PCM087|PCM089|PCM091|PCM093|PCM095|PCM097|PCM099|PCM075_dup|PCM077_dup|PCM079_dup|PCM081_dup|PCM083_dup|PCM085_dup|PCM087_dup|PCM089_dup|PCM091_dup|PCM093_dup|PCM095_dup|PCM097_dup|PCM099_dup",rownames(dup_cor_lib)),]
ggplot(dup_cor_lib_M6, aes(x = Sex, y = prop_glob, group = Sex)) +
  geom_boxplot(aes(fill = Sex), alpha = 0.5, outlier.shape = NA) +
  stat_compare_means(aes(label = "..p.adj.."), method = "t.test", label = "p.signif", ref.group = "Female", hide.ns = T) +
  geom_jitter(aes(colour = Sex), width = 0.15, size = 2) +
  ylab("Percentage of reads mapping to haemoglobin genes") +
  theme_bw() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_text(colour = "black", size = 12), axis.title.y = element_text(size = 18, face = "bold"), axis.title.x = element_blank(), legend.text = element_text(size = 14), legend.title = element_text(face = "bold", size = 16), legend.position = "bottom", strip.text.x = element_text(size = 14))

#plot correlation for all signatures
dup_cor_lib2 <- rbind(dup_cor_lib_D0, dup_cor_lib_M6)
ggscatter(dup_cor_lib2, x = "prop_glob", y = "value", size = 4,
          cor.coef = TRUE,
          cor.coeff.args = list(method = "pearson", label.y.npc = "bottom",  size = 6),
          xlab = "Percentage of reads corresponding to globin mRNA", ylab = "Pearson correlation with globin-depleted duplicate") +
  geom_smooth(method=lm, size = 1) +
  theme(axis.text.y = element_text(colour = "black", size = 12), axis.title.y = element_text(size = 18, face = "bold"), axis.title.x = element_text(size = 18, face = "bold"), legend.text = element_text(size = 14), legend.title = element_text(face = "bold", size = 16), strip.text.x = element_text(size = 14))

dup_cor_lib_u60 <- dup_cor_lib[dup_cor_lib$prop_glob<60,]
dup_cor_lib_o60 <- dup_cor_lib[dup_cor_lib$prop_glob>60,]

ggscatter(dup_cor_lib_u60, x = "prop_glob", y = "value",
          cor.coef = TRUE,
          cor.coeff.args = list(method = "pearson", label.y.npc = "bottom",  size = 6),
          xlab = "Percentage of reads corresponding to globin mRNA", ylab = "Pearson correlation with globin-depleted duplicate") +
  geom_smooth(method=lm) +
  theme(axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14))
ggscatter(dup_cor_lib_o60, x = "prop_glob", y = "value",
          cor.coef = TRUE,
          cor.coeff.args = list(method = "pearson", label.y.npc = "bottom",  size = 6),
          xlab = "Percentage of reads corresponding to globin mRNA", ylab = "Pearson correlation with globin-depleted duplicate") +
  geom_smooth(method=lm) +
  theme(axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14))

dup_cor_lib_D0 <- dup_cor_lib[grep("PCM074|PCM076|PCM078|PCM080|PCM084|PCM086|PCM088|PCM090|PCM092|PCM094|PCM096|PCM098|PCM100|PCM074_dup|PCM076_dup|PCM078_dup|PCM080_dup|PCM084_dup|PCM086_dup|PCM088_dup|PCM090_dup|PCM092_dup|PCM094_dup|PCM096_dup|PCM098_dup|PCM100_dup",rownames(dup_cor_lib)),]
dup_cor_lib_M6 <- dup_cor_lib[grep("PCM075|PCM077|PCM079|PCM081|PCM083|PCM085|PCM087|PCM089|PCM091|PCM093|PCM095|PCM097|PCM099|PCM075_dup|PCM077_dup|PCM079_dup|PCM081_dup|PCM083_dup|PCM085_dup|PCM087_dup|PCM089_dup|PCM091_dup|PCM093_dup|PCM095_dup|PCM097_dup|PCM099_dup",rownames(dup_cor_lib)),]

ggscatter(dup_cor_lib_D0, x = "prop_glob", y = "value",
          cor.coef = TRUE,
          cor.coeff.args = list(method = "pearson", label.y.npc = "bottom",  size = 6),
          xlab = "Percentage of reads corresponding to globin mRNA", ylab = "Pearson correlation with globin-depleted duplicate") +
  geom_smooth(method=lm) +
  ylim(c(0.4, 1.3)) +
  theme(axis.text.y = element_text(colour = "black", size = 12), axis.title.y = element_text(size = 18, face = "bold"), axis.title.x = element_text(size = 18, face = "bold"), legend.text = element_text(size = 14), legend.title = element_text(face = "bold", size = 16), strip.text.x = element_text(size = 14))

ggscatter(dup_cor_lib_M6, x = "prop_glob", y = "value",
          cor.coef = TRUE,
          cor.coeff.args = list(method = "pearson", label.y.npc = "bottom",  size = 6),
          xlab = "Percentage of reads corresponding to globin mRNA", ylab = "Pearson correlation with globin-depleted duplicate") +
  geom_smooth(method=lm) +
  ylim(c(0.4, 1.3)) +
  theme(axis.text.y = element_text(colour = "black", size = 12), axis.title.y = element_text(size = 18, face = "bold"), axis.title.x = element_text(size = 18, face = "bold"), legend.text = element_text(size = 14), legend.title = element_text(face = "bold", size = 16), strip.text.x = element_text(size = 14))

#correlation between sex and globin reads
cor.test(as.numeric(glob_cor$Hb), glob_cor$glob_count)
cor.test(log(as.numeric(glob_cor$Hb)), glob_cor$prop_glob)
cor.test(as.numeric(glob_cor[glob_cor$Sex=="Male",]$Hb), glob_cor[glob_cor$Sex=="Male",]$glob_count)
cor.test(as.numeric(glob_cor[glob_cor$Sex=="Male",]$Hb), glob_cor[glob_cor$Sex=="Male",]$prop_glob)
cor.test(as.numeric(glob_cor[glob_cor$Sex=="Female",]$Hb), glob_cor[glob_cor$Sex=="Female",]$glob)

rownames(glob_cor) <- glob_cor$Var2
cor.test(as.numeric(glob_cor[rownames(glob_cor)%in%colnames(dge_dups_D0),]$Hb), glob_cor[rownames(glob_cor)%in%colnames(dge_dups_D0),]$glob_count)
cor.test(as.numeric(glob_cor[rownames(glob_cor)%in%colnames(dge_dups_M6),]$Hb), glob_cor[rownames(glob_cor)%in%colnames(dge_dups_M6),]$glob_count)

#correlation between globin reads and globin levels
cor.test(glob_cor$Hb, glob_cor$prop_glob)
cor.test(glob_cor$Hb, glob_cor$glob_count)

#plot GSVA scores between globin-depleted and non-depleted samples
dups_GSVA <- as.data.frame(colData(scores_GSVA)[,c(5,108:ncol(colData(scores_GSVA)))])
dups_GSVA <- reshape2::melt(dups_GSVA, id.vars = c("Globin","PCM"))

c30 <- c(
  "dodgerblue2", "#E31A1C", # red
  "green4",
  "#6A3D9A", # purple
  "#FF7F00", # orange
  "black", "gold1",
  "skyblue2", "#FB9A99", # lt pink
  "palegreen2",
  "#CAB2D6", # lt purple
  "#FDBF6F", # lt orange
  "gray70", "khaki2",
  "maroon", "orchid1", "deeppink1", "blue1", "steelblue4",
  "darkturquoise", "green1", "yellow4", "yellow3",
  "darkorange4", "brown",
  "turquoise", "maroon",
  "pink", "cyan",
  "magenta"
)

ggplot(dups_GSVA, aes(x = Globin, y = value, group = Globin)) +
  geom_boxplot(aes(fill = Globin), alpha = 0.7, outlier.shape = NA) +
  scale_fill_manual(values = c("#377EB8", "#4DAF4A")) +
  geom_point(aes(group = PCM, colour = PCM), position = position_dodge(0.2)) +
  geom_line(aes(group = PCM, colour = PCM), alpha = 0.6, position = position_dodge(0.2)) +
  scale_color_manual(values = c30) +
  theme_bw() +
  xlab("Group") +
  ylab("Score") +
  facet_wrap(~ variable, ncol = 7, nrow = 10) +
  stat_compare_means(aes(label = "..p.adj.."), method = "t.test", label = "p.signif", ref.group = "Depleted", hide.ns = T, label.y = 0.8) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_text(colour = "black", size = 12), axis.title.y = element_text(size = 18, face = "bold"), axis.title.x = element_text(size = 18, face = "bold"), legend.text = element_text(size = 14), legend.title = element_text(face = "bold", size = 16), strip.text.x = element_text(size = 14))

#perform an example DGE analysis between TB controls and post-treatment timepoint
gi2gn <- function(x){
  gene_id <- rownames(x)
  genes <- ensembl[match(gene_id, ensembl$gene_id),]
  genes <- genes[!duplicated(genes$gene_name),]
  genes <- genes[!is.na(genes$gene_name),]
  x <- x[rownames(x) %in% genes$gene_id,]
  genes <- genes[genes$gene_id %in% rownames(x),]
  rownames(x) <- genes$gene_name
  x
}
dge_dups$samples$Visit <- meta_dup$Visit
dge_dups$samples$Globin <- meta_dup$Globin
dge_dups$samples$PID <- meta_dup$PID
dge_dups_TB <- dge_dups[,grep("D0|M6",dge_dups$samples$Visit)]
dge_dups_TB_glob <- dge_dups_TB[,dge_dups_TB$samples$Globin=="Depleted"]
dge_dups_TB_nonglob <- dge_dups_TB[,dge_dups_TB$samples$Globin=="Non-depleted"]

fdr.cutoff <- 0.05
lfc.cutoff <- 0.58

sig_res <- function(x){
  res <- topTags(x, n=nrow(x$table), adjust.method = "BH", sort.by = "PValue")
  res <- as.data.frame(res)
  sig.res <- res[res$FDR < fdr.cutoff,]
  sig.res <- sig.res[abs(sig.res$logFC) > lfc.cutoff,]
  gene_id <- rownames(sig.res)
  gene_id <- ensembl[match(gene_id, ensembl$gene_id),]
  gene_id <- gene_id[!duplicated(gene_id$gene_name),]
  sig.res <- sig.res[rownames(sig.res) %in% gene_id$gene_id,]
  gene_id <- gene_id[gene_id$gene_id %in% rownames(sig.res),]
  rownames(sig.res) <- gene_id$gene_name
  sig.res$gene_id <- gene_id$gene_id
  sig.res$transcript_id <- gene_id$transcript_id
  sig.res$entrez <- gene_id$entrez
  sig.res$gene_biotype <- gene_id$gene_biotype
  sig.res$chromosome <- gene_id$chromosome
  sig.res
}

des_dups_TB_glob <- model.matrix(~ dge_dups_TB_glob$samples$Visit)
dge_dups_TB_glob <- estimateDisp(dge_dups_TB_glob, design = des_dups_TB_glob)
fit_dups_TB_glob <- glmQLFit(dge_dups_TB_glob, des_dups_TB_glob)
lrt_dups_TB_glob <- glmLRT(fit_dups_TB_glob, coef = 2)
summary(decideTests(lrt_dups_TB_glob))
#       dge_dups_TB_glob$samples$VisitM6
#Down                                 46
#NotSig                            60097
#Up                                  498

tt_dups_TB_glob <- topTags(lrt_dups_TB_glob, n=nrow(lrt_dups_TB_glob$table), adjust.method = "BH", sort.by = "PValue")
tt_dups_TB_glob <- gi2gn(as.data.frame(tt_dups_TB_glob))
fg_dups_TB_glob <- rownames(tt_dups_TB_glob)[tt_dups_TB_glob$FDR < 0.05 & abs(tt_dups_TB_glob$logFC) > 1]
tt_dups_TB_glob_up <- tt_dups_TB_glob[tt_dups_TB_glob$logFC>0,]
fg_dups_TB_glob_up <- rownames(tt_dups_TB_glob)[tt_dups_TB_glob$FDR < 0.05 & tt_dups_TB_glob$logFC > 1]
tt_dups_TB_glob_down <- tt_dups_TB_glob[tt_dups_TB_glob$logFC<0,]
fg_dups_TB_glob_down <- rownames(tt_dups_TB_glob)[tt_dups_TB_glob$FDR < 0.05 & tt_dups_TB_glob$logFC < 1]

des_dups_TB_nonglob <- model.matrix(~ dge_dups_TB_nonglob$samples$Visit)
dge_dups_TB_nonglob <- estimateDisp(dge_dups_TB_nonglob, design = des_dups_TB_nonglob)
fit_dups_TB_nonglob <- glmQLFit(dge_dups_TB_nonglob, des_dups_TB_nonglob)
lrt_dups_TB_nonglob <- glmLRT(fit_dups_TB_nonglob, coef = 2)
summary(decideTests(lrt_dups_TB_nonglob))
#       dge_dups_TB_nonglob$samples$VisitM6
#Down                                    42
#NotSig                               60204
#Up                                     395

tt_dups_TB_nonglob <- topTags(lrt_dups_TB_nonglob, n=nrow(lrt_dups_TB_nonglob$table), adjust.method = "BH", sort.by = "PValue")
tt_dups_TB_nonglob <- gi2gn(as.data.frame(tt_dups_TB_nonglob))
fg_dups_TB_nonglob <- rownames(tt_dups_TB_nonglob)[tt_dups_TB_nonglob$FDR < 0.05 & abs(tt_dups_TB_nonglob$logFC) > 1]
tt_dups_TB_nonglob_up <- tt_dups_TB_nonglob[tt_dups_TB_nonglob$logFC>0,]
fg_dups_TB_nonglob_up <- rownames(tt_dups_TB_nonglob_up)[tt_dups_TB_nonglob_up$FDR < 0.05 & tt_dups_TB_nonglob_up$logFC > 1]
tt_dups_TB_nonglob_down <- tt_dups_TB_nonglob[tt_dups_TB_nonglob$logFC<0,]
fg_dups_TB_nonglob_down <- rownames(tt_dups_TB_nonglob_down)[tt_dups_TB_nonglob_down$FDR < 0.05 & tt_dups_TB_nonglob_down$logFC < 1]

glob_nonglob_M6_olap <- venn(list(fg_dups_TB_glob, fg_dups_TB_nonglob))
glob_nonglob_M6_olap_up <- venn(list(fg_dups_TB_glob_up, fg_dups_TB_nonglob_up))
glob_nonglob_M6_olap_down <- venn(list(fg_dups_TB_glob_down, fg_dups_TB_nonglob_down))

tmod_uni_glob <- tmodHGtest(fg=as.factor(attr(glob_nonglob_M6_olap,"intersections")$A), bg=as.factor(rownames(tt_dups_TB_glob)))
tmod_uni_nonglob <- tmodHGtest(fg=as.factor(attr(glob_nonglob_M6_olap,"intersections")$B), bg=as.factor(rownames(tt_dups_TB_glob)))
tmodHGtest(fg=as.factor(attr(glob_nonglob_M6_olap,"intersections")$`A:B`), bg=as.factor(rownames(tt_dups_TB_glob)))

#repeat regressing out globin levels
dge_dups_TB_glob$samples$prop_glob <- libs_df_dups[colnames(dge_dups_TB_glob),]$prop_glob
des_dups_TB_glob <- model.matrix(~ dge_dups_TB_glob$samples$Visit + dge_dups_TB_glob$samples$prop_glob)
dge_dups_TB_glob <- estimateDisp(dge_dups_TB_glob, design = des_dups_TB_glob)
fit_dups_TB_glob <- glmQLFit(dge_dups_TB_glob, des_dups_TB_glob)
lrt_dups_TB_glob <- glmLRT(fit_dups_TB_glob, coef = 2)
summary(decideTests(lrt_dups_TB_glob))
#       dge_dups_TB_glob$samples$VisitM6
#Down                                 54
#NotSig                            60079
#Up                                  508

tt_dups_TB_glob <- topTags(lrt_dups_TB_glob, n=nrow(lrt_dups_TB_glob$table), adjust.method = "BH", sort.by = "PValue")
tt_dups_TB_glob <- gi2gn(as.data.frame(tt_dups_TB_glob))
fg_dups_TB_glob <- rownames(tt_dups_TB_glob)[tt_dups_TB_glob$FDR < 0.05 & abs(tt_dups_TB_glob$logFC) > 1]
tt_dups_TB_glob_up <- tt_dups_TB_glob[tt_dups_TB_glob$logFC>0,]
fg_dups_TB_glob_up <- rownames(tt_dups_TB_glob)[tt_dups_TB_glob$FDR < 0.05 & tt_dups_TB_glob$logFC > 1]
tt_dups_TB_glob_down <- tt_dups_TB_glob[tt_dups_TB_glob$logFC<0,]
fg_dups_TB_glob_down <- rownames(tt_dups_TB_glob)[tt_dups_TB_glob$FDR < 0.05 & tt_dups_TB_glob$logFC < 1]

dge_dups_TB_nonglob$samples$prop_glob <- libs_df_dups[colnames(dge_dups_TB_nonglob),]$prop_glob
dge_dups_TB_nonglob$samples$prop_glob_thresh <- rep("Yes",22)
dge_dups_TB_nonglob$samples[colnames(dge_dups_TB_nonglob) %in% rownames(dup_cor_lib_o60),]$prop_glob_thresh <- "No"
rownames(dup_cor_lib_o60)

des_dups_TB_nonglob <- model.matrix(~ dge_dups_TB_nonglob$samples$Visit + dge_dups_TB_nonglob$samples$prop_glob)
dge_dups_TB_nonglob <- estimateDisp(dge_dups_TB_nonglob, design = des_dups_TB_nonglob)
fit_dups_TB_nonglob <- glmQLFit(dge_dups_TB_nonglob, des_dups_TB_nonglob)
lrt_dups_TB_nonglob <- glmLRT(fit_dups_TB_nonglob, coef = 2)
summary(decideTests(lrt_dups_TB_nonglob))
#       dge_dups_TB_nonglob$samples$VisitM6
#Down                                    42
#NotSig                               60204
#Up                                     386

tt_dups_TB_nonglob <- topTags(lrt_dups_TB_nonglob, n=nrow(lrt_dups_TB_nonglob$table), adjust.method = "BH", sort.by = "PValue")
tt_dups_TB_nonglob <- gi2gn(as.data.frame(tt_dups_TB_nonglob))
fg_dups_TB_nonglob <- rownames(tt_dups_TB_nonglob)[tt_dups_TB_nonglob$FDR < 0.05 & abs(tt_dups_TB_nonglob$logFC) > 1]
tt_dups_TB_nonglob_up <- tt_dups_TB_nonglob[tt_dups_TB_nonglob$logFC>0,]
fg_dups_TB_nonglob_up <- rownames(tt_dups_TB_nonglob_up)[tt_dups_TB_nonglob_up$FDR < 0.05 & tt_dups_TB_nonglob_up$logFC > 1]
tt_dups_TB_nonglob_down <- tt_dups_TB_nonglob[tt_dups_TB_nonglob$logFC<0,]
fg_dups_TB_nonglob_down <- rownames(tt_dups_TB_nonglob_down)[tt_dups_TB_nonglob_down$FDR < 0.05 & tt_dups_TB_nonglob_down$logFC < 1]

glob_nonglob_M6_olap <- venn(list(fg_dups_TB_glob, fg_dups_TB_nonglob))
glob_nonglob_M6_olap_up <- venn(list(fg_dups_TB_glob_up, fg_dups_TB_nonglob_up))
glob_nonglob_M6_olap_down <- venn(list(fg_dups_TB_glob_down, fg_dups_TB_nonglob_down))

tmod_uni_glob <- tmodHGtest(fg=as.factor(attr(glob_nonglob_M6_olap,"intersections")$A), bg=as.factor(rownames(tt_dups_TB_glob)))
tmod_uni_nonglob <- tmodHGtest(fg=as.factor(attr(glob_nonglob_M6_olap,"intersections")$B), bg=as.factor(rownames(tt_dups_TB_glob)))

gi2ez <- function(x){
  gene_id <- rownames(x)
  genes <- ensembl[match(gene_id, ensembl$gene_id),]
  genes <- genes[!duplicated(genes$entrez),]
  genes <- genes[!is.na(genes$entrez),]
  x <- x[rownames(x) %in% genes$gene_id,]
  genes <- genes[genes$gene_id %in% rownames(x),]
  rownames(x) <- genes$entrez
  x
}

res_dups_TB_nonglob_comp <- res_dups_TB_nonglob[res_dups_TB_nonglob$Symbol %in% res_dups_TB_glob$Symbol,]
res_dups_TB_glob_comp <- res_dups_TB_glob[res_dups_TB_glob$Symbol %in% res_dups_TB_nonglob$Symbol,]
res_dups_TB_glob_v_nonglob_comp <- merge(res_dups_TB_glob_comp, res_dups_TB_nonglob_comp, by = "Symbol")

rownames(res_dups_TB_glob_v_nonglob_comp) <- res_dups_TB_glob_v_nonglob_comp$Symbol
res_dups_TB_glob_v_nonglob_comp <- res_dups_TB_glob_v_nonglob_comp[,c("logFC.x","logFC.y","glob_signif","nonglob_signif")]
res_dups_TB_glob_v_nonglob_comp$comp <- ifelse(res_dups_TB_glob_v_nonglob_comp$glob_signif == TRUE & res_dups_TB_glob_v_nonglob_comp$nonglob_signif == FALSE, "globin_only", ifelse(res_dups_TB_glob_v_nonglob_comp$glob_signif == FALSE & res_dups_TB_glob_v_nonglob_comp$nonglob_signif == TRUE, "nonglobin_only", ifelse(res_dups_TB_glob_v_nonglob_comp$glob_signif == TRUE & res_dups_TB_glob_v_nonglob_comp$nonglob_signif == TRUE, "Both", "Neither")))
res_dups_TB_glob_v_nonglob_comp_genes <- rownames(res_dups_TB_glob_v_nonglob_comp)
gene_labs <- subset(res_dups_TB_glob_v_nonglob_comp, comp == "Both")

sig_res_dups_TB_glob_v_nonglob_comp <- res_dups_TB_glob_v_nonglob_comp[grep("Both|globin_only|nonglobin_only",res_dups_TB_glob_v_nonglob_comp$comp),]
cor.test(sig_res_dups_TB_glob_v_nonglob_comp$logFC.x, sig_res_dups_TB_glob_v_nonglob_comp$logFC.y, method = "spearman")
#Spearman's rank correlation rho

#data:  sig_res_dups_TB_glob_v_nonglob_comp$logFC.x and sig_res_dups_TB_glob_v_nonglob_comp$logFC.y
#S = 3430397540, p-value < 2.2e-16
#alternative hypothesis: true rho is not equal to 0
#sample estimates:
#      rho 
#0.5372118

library(ggrepel)
ggplot(res_dups_TB_glob_v_nonglob_comp) +
  geom_point(data = subset(res_dups_TB_glob_v_nonglob_comp, comp = "Both"), aes(x = logFC.x, y = logFC.y), color = "red", alpha = 0.8, size = 2) +
  geom_point(data = subset(res_dups_TB_glob_v_nonglob_comp, comp == "Neither"), aes(x = logFC.x, y = logFC.y), color = "black", alpha = 0.1, size = 2) +
  geom_point(data = subset(res_dups_TB_glob_v_nonglob_comp, comp == "globin_only"), aes(x = logFC.x, y = logFC.y), color = "lightgreen", alpha = 0.5, size = 2) +
  geom_point(data = subset(res_dups_TB_glob_v_nonglob_comp, comp == "nonglobin_only"), aes(x = logFC.x, y = logFC.y), color = "lightblue", alpha = 0.5, size = 2) +
  geom_hline(yintercept = c(-1,1), linetype = "dashed", color = "red", alpha = 0.5) +
  geom_vline(xintercept = c(-1,1), linetype = "dashed", color = "red", alpha = 0.5) +
  labs(x = expression(paste(log[2], " fold change for globin-depleted")), y = expression(paste(log[2], " fold change for non-globin-depleted"))) +
  geom_smooth(data = subset(res_dups_TB_glob_v_nonglob_comp, comp == "Both"), aes(x = logFC.x, y = logFC.y), method = "lm", linetype = "dashed", colour = "blue", alpha = 0.3, size = 0.5, se = F, fullrange = T) +
  annotate("text", x=2, y=5, size = 8, label = "rho = 0.5372118") +
  annotate("text", x=2, y=4.5, size = 8, label = "p-value < 2.2e-16") +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14))

#directly compare globin-depleted and non-depleted groups
gi2ez <- function(x){
  gene_id <- rownames(x)
  genes <- ensembl[match(gene_id, ensembl$gene_id),]
  genes <- genes[!duplicated(genes$entrez),]
  genes <- genes[!is.na(genes$entrez),]
  x <- x[rownames(x) %in% genes$gene_id,]
  genes <- genes[genes$gene_id %in% rownames(x),]
  rownames(x) <- genes$entrez
  x
}

#both timepoints combined (globin-depleted vs non-depleted)
des_dups_TB <- model.matrix(~dge_dups_TB$samples$PID + dge_dups_TB$samples$Globin)
dge_dups_TB <- estimateDisp(dge_dups_TB, des_dups_TB)
fit_dup_TB <- glmQLFit(dge_dups_TB, des_dups_TB)
qlf_dup_TB <- glmQLFTest(fit_dup_TB, coef = 12)
#GO analysis
qlf_dup_TB <- gi2ez(qlf_dup_TB)
go_dup_TB <- goana(qlf_dup_TB, gene_id = rownames(qlf_dup_TB), FDR = 0.05, species="Hs")
topGO(go_dup_TB, sort = "up")
topGO(go_dup_TB, sort = "down")

#add proportion of globin reads as a covariate
dge_dups_TB$samples$prop_glob <- libs_df_dups[rownames(dge_dups_TB$samples),]$prop_glob
des_dups_TB_pg <- model.matrix(~dge_dups_TB$samples$PID + dge_dups_TB$samples$Globin + dge_dups_TB$samples$prop_glob)
dge_dups_TB_pg <- estimateDisp(dge_dups_TB, des_dups_TB_pg)
fit_dup_TB_pg <- glmQLFit(dge_dups_TB_pg, des_dups_TB_pg)
qlf_dup_TB_pg <- glmQLFTest(fit_dup_TB_pg, coef = 13)
#GO analysis
qlf_dup_TB_pg <- gi2ez(qlf_dup_TB_pg)
go_dup_TB_pg <- goana(qlf_dup_TB_pg, gene_id = rownames(qlf_dup_TB_pg), FDR = 0.05, species="Hs")
topGO(go_dup_TB_pg, sort = "up")
topGO(go_dup_TB_pg, sort = "down")

#KEGG analysis
keg_dup_TB <- kegga(qlf_dup_TB, gene_id = rownames(qlf_dup_TB), FDR = 0.05, species="Hs")
topKEGG(keg_dup_TB, sort = "up")
topKEGG(keg_dup_TB, sort = "down")
qlf_dup_TB <- glmQLFTest(fit_dup_TB, coef = 12)
#topTags and volcano plot
tt_dup_TB <- as.data.frame(topTags(qlf_dup_TB, n = nrow(qlf_dup_TB), adjust.method = "BH", sort.by = "PValue"))
tt_dup_TB <- tt_dup_TB[-3,]
tt_dup_TB <- gi2gn(tt_dup_TB)
library(EnhancedVolcano)
EnhancedVolcano(tt_dup_TB,
                lab = rownames(tt_dup_TB),
                x = 'logFC',
                y = 'PValue')

#point of diagnosis (D0, globin-depleted vs non-depleted)
dge_dups_D0 <- dge_dups_TB[,grep("D0",dge_dups_TB$samples$Visit)]
des_dups_D0 <- model.matrix(~dge_dups_D0$samples$PID + dge_dups_D0$samples$Globin)
dge_dups_D0 <- estimateDisp(dge_dups_D0, des_dups_D0)
fit_dup_D0 <- glmQLFit(dge_dups_D0, des_dups_D0)
qlf_dup_D0 <- glmQLFTest(fit_dup_D0, coef = 12)
#GO analysis
qlf_dup_D0 <- gi2ez(qlf_dup_D0)
go_dup_D0 <- goana(qlf_dup_D0, gene_id = rownames(qlf_dup_D0), FDR = 0.05, species="Hs")
topGO(go_dup_D0, sort = "up")
topGO(go_dup_D0, sort = "down")
#KEGG analysis
keg_dup_D0 <- kegga(qlf_dup_D0, gene_id = rownames(qlf_dup_D0), FDR = 0.05, species="Hs")
topKEGG(keg_dup_D0, sort = "up")
topKEGG(keg_dup_D0, sort = "down")
qlf_dup_D0 <- glmQLFTest(fit_dup_D0, coef = 12)
#topTags and volcano plot
tt_dup_D0 <- as.data.frame(topTags(qlf_dup_D0, n = nrow(qlf_dup_D0), adjust.method = "BH", sort.by = "PValue"))
tt_dup_D0 <- tt_dup_D0[-3,]
tt_dup_D0 <- gi2gn(tt_dup_D0)
EnhancedVolcano(tt_dup_D0,
                lab = rownames(tt_dup_D0),
                x = 'logFC',
                y = 'PValue')

#6 months post-treatment (M6, globin-depleted vs non-depleted)
dge_dups_M6 <- dge_dups_TB[,grep("M6",dge_dups_TB$samples$Visit)]
des_dups_M6 <- model.matrix(~dge_dups_M6$samples$PID + dge_dups_M6$samples$Globin)
dge_dups_M6 <- estimateDisp(dge_dups_M6, des_dups_M6)
fit_dup_M6 <- glmQLFit(dge_dups_M6, des_dups_M6)
qlf_dup_M6 <- glmQLFTest(fit_dup_M6, coef = 12)
#GO analysis
qlf_dup_M6 <- gi2ez(qlf_dup_M6)
go_dup_M6 <- goana(qlf_dup_M6, gene_id = rownames(qlf_dup_M6), FDR = 0.05, species="Hs")
topGO(go_dup_M6, sort = "down")
#KEGG analysis
keg_dup_M6 <- kegga(qlf_dup_M6, gene_id = rownames(qlf_dup_M6), FDR = 0.05, species="Hs")
topKEGG(keg_dup_M6, sort = "down")
qlf_dup_M6 <- glmQLFTest(fit_dup_M6, coef = 12)
#topTags and volcano plot
tt_dup_M6 <- as.data.frame(topTags(qlf_dup_M6, n = nrow(qlf_dup_M6), adjust.method = "BH", sort.by = "PValue"))
tt_dup_M6 <- tt_dup_M6[-3,]
tt_dup_M6 <- gi2gn(tt_dup_M6)
EnhancedVolcano(tt_dup_M6,
                lab = rownames(tt_dup_M6),
                x = 'logFC',
                y = 'PValue')

#D0 and M6 combined and separated (from analysis)
glob_genes <- list(D0 = rownames(tt_dup_D0)[1:200], D0_down = rownames(tt_dup_D0[tt_dup_D0$logFC>0,])[1:200], D0_up = rownames(tt_dup_D0[tt_dup_D0$logFC<0,])[1:200], M6 = rownames(tt_dup_M6)[1:200], M6_down = rownames(tt_dup_M6[tt_dup_M6$logFC<0,])[1:200], M6_up = rownames(tt_dup_D0[tt_dup_D0$logFC>0,])[1:200], Both = rownames(tt_dup_TB)[1:200], Both_down = rownames(tt_dup_TB[tt_dup_D0$logFC<0,])[1:200], Both_up = rownames(tt_dup_TB[tt_dup_D0$logFC>0,])[1:200])
se_dups_TB <- SummarizedExperiment(assays = list(counts = as.matrix(gi2gn(dge_dups_TB$counts))), colData = dge_dups_TB$samples)
se_dups_TB <- mkAssay(se_dups_TB, input_name = "counts", log = TRUE, counts_to_CPM = TRUE)
scores_GSVA <- runTBsigProfiler(se_dups_TB, useAssay = "log_counts_cpm", signatures = glob_genes, algorithm = "GSVA")
signatureBoxplot(inputData = scores_GSVA, 
                 name = "GSVA",
                 signatureColNames = names(glob_genes),
                 annotationColName = "Globin", includePoints = TRUE) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

#D0 and M6 separated only (post-analysis)
se_dups_TB <- SummarizedExperiment(assays = list(counts = as.matrix(gi2gn(dge_dups_TB$counts))), colData = dge_dups_TB$samples)
se_dups_TB <- mkAssay(se_dups_TB, input_name = "counts", log = TRUE, counts_to_CPM = TRUE)
glob_genes <- list(D0 = rownames(tt_dup_D0)[1:200], D0_down = rownames(tt_dup_D0[tt_dup_D0$logFC<0,])[1:200], D0_up = rownames(tt_dup_D0[tt_dup_D0$logFC>0,])[1:200], M6 = rownames(tt_dup_M6)[1:200], M6_down = rownames(tt_dup_M6[tt_dup_M6$logFC<0,])[1:200], M6_up = rownames(tt_dup_D0[tt_dup_D0$logFC>0,])[1:200])
library(data.table)
lapply(glob_genes, function(x) write.table( data.frame(x), 'glob_genes.csv'  , append= T, sep=',' ))
D0_M6_GSVA <- runTBsigProfiler(se_dups_TB, useAssay = "log_counts_cpm", signatures = glob_genes, algorithm = "GSVA")
D0_M6_GSVA_1 <- as.data.frame(colData(D0_M6_GSVA)[,c(5,7:ncol(colData(D0_M6_GSVA)))])
D0_M6_GSVA_1 <- reshape2::melt(D0_M6_GSVA_1, id.vars = c("Globin"))
D0_M6_GSVA_2 <- as.data.frame(colData(D0_M6_GSVA)[,c(5,8,9,11,12)])
D0_M6_GSVA_2 <- reshape2::melt(D0_M6_GSVA_2, id.vars = c("Globin"))
D0_M6_GSVA_3 <- as.data.frame(colData(D0_M6_GSVA)[,c(1,2,5,8,9,11,12)])
rownames(libs_df_dups_prop) <- libs_df_dups_prop$sample
D0_M6_GSVA_3$prop_glob <- libs_df_dups_prop[rownames(D0_M6_GSVA_3),]$prop_glob
D0_M6_GSVA_3 <- reshape2::melt(D0_M6_GSVA_3, id.vars = c("Globin","lib.size","group","prop_glob"))

ggplot(D0_M6_GSVA_2, aes(x = Globin, y = value, group = Globin)) +
  geom_boxplot(aes(fill = Globin), alpha = 0.7, outlier.shape = NA) +
  scale_fill_manual(values = c("#4DAF4A","#377EB8")) +
  geom_point(position = position_dodge(0.2)) +
  geom_line(alpha = 0.6, position = position_dodge(0.2)) +
  theme_bw() +
  xlab("Group") +
  ylab("Score") +
  facet_wrap(~ variable, ncol = 2, nrow = 2) +
  stat_compare_means(aes(label = "..p.adj.."), method = "t.test", label = "p.signif", ref.group = "Depleted", hide.ns = T, label.y = 0.8) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_text(colour = "black", size = 12), axis.title.y = element_text(size = 18, face = "bold"), axis.title.x = element_blank(), legend.text = element_text(size = 14), legend.title = element_text(face = "bold", size = 16), strip.text.x = element_text(size = 14))

mid <- median(log10(D0_M6_GSVA_3[D0_M6_GSVA_3$group == "Non-depleted",]$lib.size))
mid <- mean(D0_M6_GSVA_3[D0_M6_GSVA_3$group == "Non-depleted",]$prop_glob)
ggplot(D0_M6_GSVA_3, aes(x = Globin, y = value, group = Globin)) +
  geom_boxplot(aes(fill = Globin), alpha = 0.7, outlier.shape = NA) +
  scale_fill_manual(values = c("#4DAF4A","#377EB8")) +
  geom_point(position = position_dodge(0.2), aes(colour = prop_glob, size = 12)) +
  scale_colour_gradient2(midpoint = mid, low = "red",
                        mid = "white",
                        high = "blue") +
  geom_line(alpha = 0.6, position = position_dodge(0.2)) +
  theme_bw() +
  xlab("Group") +
  ylab("Score") +
  facet_wrap(~ variable, ncol = 2, nrow = 2) +
  stat_compare_means(aes(label = "..p.adj.."), method = "t.test", label = "p.signif", ref.group = "Depleted", hide.ns = T, label.y = 0.8) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_text(colour = "black", size = 12), axis.title.y = element_text(size = 18, face = "bold"), axis.title.x = element_blank(), legend.text = element_text(size = 14), legend.title = element_text(face = "bold", size = 16), strip.text.x = element_text(size = 14))

#validate on full dataset
se_dge$Globin <- meta$Globin
scores_GSVA <- runTBsigProfiler(se_dge, useAssay = "log_cpm", signatures = glob_genes, algorithm = "GSVA")
signatureBoxplot(inputData = scores_GSVA, 
                 name = "GSVA",
                 signatureColNames = names(glob_genes),
                 annotationColName = "Globin", includePoints = TRUE) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

#validate on published dataset
glob_val <- read.table("glob_valid.txt", header = T, row.names = 1)
colnames(glob_val) <- gsub("X.stornext.Home.data.users.allstaff.sheerin.d.globin_sra.BAM.", "", colnames(glob_val))
colnames(glob_val) <- gsub("_1.fastq.subread.BAM", "", colnames(glob_val))
colnames(glob_val) <- gsub(".sra.fastq.subread.BAM", "", colnames(glob_val))
glob_val_meta <- read.csv("glob_val_meta.csv", header = T, row.names = 1)
glob_val_meta <- glob_val_meta[c(45,44,43,42,41,40,39,38,37,36,35,34,33,32,31,30,29,28,27,26,25,24,23,22,21,20,19,18,17,16,15,14,13,12,11,10,9,8,7,6,5,4,3,2,1),]
all(rownames(glob_val_meta) %in% colnames(glob_val))
#TRUE
all(rownames(glob_val_meta) == colnames(glob_val))
#TRUE

dge_gv <- DGEList(glob_val, group = glob_val_meta$group)
dge_gv$samples$input <- glob_val_meta$input
dge_gv$samples$kit <- glob_val_meta$kit

gene_id <- rownames(dge_gv)
my_genes <- ensembl[match(gene_id, ensembl$gene_id),]
dge_gv$genes <- my_genes

#define globin genes and remove from edgeR object
pc_genes <- rownames(dge_gv$counts)
#pc_genes <- rownames(dge_gv$counts[grep("protein_coding",dge_gv$genes$biotype),])
pre_gd_counts_gv <- dge_gv$counts
pre_gd_libs_gv <- dge_gv$samples$lib.size
names(pre_gd_libs_gv) <- colnames(dge_gv)
dge_gv <- dge_gv[-grep("HBA1|HBA2|HBB|HBBP1|HBD|HBE1|HBG1|HBG2|HBM|HBQ1|HBZ|HBZP1",dge_gv$genes$gene_name),, keep.lib.sizes=FALSE]
dim(dge_gv)
#60641 45
gd_counts_gv <- dge_gv$counts
gd_libs_gv <- dge_gv$samples$lib.size
names(gd_libs_gv) <- colnames(dge_gv)

#normalise
dge_gv <- calcNormFactors(dge_gv, method = "TMM")

#create summarized experiment for later
se_dge_gv <- SummarizedExperiment(assays = list(counts = as.matrix(gi2gn(dge_gv$counts))), colData = dge_gv$samples)
se_dge_gv <- mkAssay(se_dge_gv, input_name = "counts", log = TRUE, counts_to_CPM = TRUE)

se_dge_gv$kit <- glob_val_meta$kit
scores_GSVA_gv <- runTBsigProfiler(se_dge_gv, useAssay = "log_cpm", signatures = glob_genes, algorithm = "GSVA")
signatureBoxplot(inputData = scores_GSVA_gv, 
                 name = "GSVA",
                 signatureColNames = names(glob_genes),
                 annotationColName = "kit", includePoints = TRUE) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

harrington_GSVA <- as.data.frame(colData(scores_GSVA_gv)[,c(5,7,8,10,11)])
harrington_GSVA <- reshape2::melt(harrington_GSVA, id.vars = c("kit"))
harrington_GSVA$kit[harrington_GSVA$kit == "G"] <- "Globin-Zero"
harrington_GSVA$kit[harrington_GSVA$kit == "RZG"] <- "Ribo-Zero Gold"

ggplot(harrington_GSVA, aes(x = kit, y = value, group = kit)) +
  geom_boxplot(aes(fill = kit), alpha = 0.7, outlier.shape = NA) +
  geom_point(position = position_dodge(0.2)) +
  geom_line(alpha = 0.6, position = position_dodge(0.2)) +
  scale_fill_manual(values = c("#4DAF4A","#377EB8")) +
  theme_bw() +
  xlab("Group") +
  ylab("Score") +
  facet_wrap(~ variable, ncol = 2, nrow = 2) +
  stat_compare_means(aes(label = "..p.adj.."), method = "t.test", label = "p.signif", ref.group = "Globin-Zero", hide.ns = T, label.y = 0.8) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_text(colour = "black", size = 12), axis.title.y = element_text(size = 18, face = "bold"), axis.title.x = element_blank(), legend.text = element_text(size = 14), legend.title = element_text(face = "bold", size = 16), legend.position = "bottom", strip.text.x = element_text(size = 14))

#validate recovery signatures on D0 M6 subset
pre_gd_combatseq_sub <- sva::ComBat_seq(counts=pre_gd_counts, batch=meta$Globin, group=meta$L..Lung, shrink=FALSE, shrink.disp=FALSE)

dge_pre_gd <- DGEList(pre_gd_combatseq_sub, group = meta$LN.Lung)

#normalise
dge_pre_gd <- calcNormFactors(dge_pre_gd, method = "TMM")

dge_pre_gd_dups <- dge_pre_gd[,rownames(meta_dup)]

dge_pre_gd_dups$samples$Visit <- meta_dup$Visit
dge_pre_gd_dups$samples$Globin <- meta_dup$Globin
dge_pre_gd_dups$samples$PID <- meta_dup$PID

dge_pre_gd_dups_TB <- dge_pre_gd_dups[,grep("D0|M6",dge_pre_gd_dups$samples$Visit)]

se_pre_gd_dups_TB <- SummarizedExperiment(assays = list(counts = as.matrix(gi2gn(dge_pre_gd_dups_TB$counts))), colData = dge_pre_gd_dups_TB$samples)
se_pre_gd_dups_TB <- mkAssay(se_pre_gd_dups_TB, input_name = "counts", log = TRUE, counts_to_CPM = TRUE)

treat_sigs <- TBsignatures[c(8,9,61)]
treat_sigs$Tabone_TREAT_TB27 <- c("BAK1","CASP5","GRAMD1C","SCARF1","TIFA","GBP6","ICAM1","FCGR1A","SOCS3","LIMK1","GBP5","PRTN3","SEPT4","C1QC","PDCD1LG2","FAM20A","SDC3","FCGR1B","CD274","CDCP1","GRIN3A","NRN1","C1QB","GBP1P1","IGF2BP3","P2RY14","FCGR1CP")
treat_sigs$Tabone_EarlyRESP_TB25 <- c("IFITM1","IRF1","CD274","UBE2L6","APOL1","FCER1G","EFCAB2","ANXA3","TIFA","DUSP3","POLB","PML","VAMP5","SORT1","GPR84","SEPT4","TIMM10","PSME2","RP11-468E2.5","RARRES3","RP11-476D10.1","PDCD1LG2","CASP5","APOL6","SOCS1")

treat_GSVA <- runTBsigProfiler(se_dups, useAssay = "log_cpm", signatures = treat_sigs, algorithm = "GSVA")
treat_GSVA <- as.data.frame(colData(treat_GSVA))
treat_GSVA <- treat_GSVA[c(colnames(dge_dups_D0),colnames(dge_dups_M6)),]
treat_GSVA_glob <- treat_GSVA[grep("Depleted",treat_GSVA$Globin),]
treat_GSVA_nonglob <- treat_GSVA[grep("Non-depleted",treat_GSVA$Globin),]
treat_GSVA_glob <- treat_GSVA_glob[,c(3,109:ncol(treat_GSVA_glob))]
treat_GSVA_glob <- reshape2::melt(treat_GSVA_glob, id.vars = c("Visit"))
treat_GSVA_nonglob <- treat_GSVA_nonglob[,c(3,109:ncol(treat_GSVA_nonglob))]
treat_GSVA_nonglob <- reshape2::melt(treat_GSVA_nonglob, id.vars = c("Visit"))

treat_GSVA_glob <- treat_GSVA_glob[grep("Bloom_RES_558|Tabone_EarlyRESP_TB25",treat_GSVA_glob$variable),]
treat_GSVA_nonglob <- treat_GSVA_nonglob[grep("Bloom_RES_558|Tabone_EarlyRESP_TB25",treat_GSVA_nonglob$variable),]

gg_glob <- ggplot(treat_GSVA_glob, aes(x = Visit, y = value, group = Visit)) +
  geom_boxplot(aes(fill = Visit), alpha = 0.7, outlier.shape = NA) +
  geom_point(position = position_dodge(0.2), colour = "#4DAF4A") +
  geom_line(alpha = 0.6, position = position_dodge(0.2)) +
  theme_bw() +
  xlab("Group") +
  ylab("Score") +
  facet_wrap(~ variable, ncol = 1) +
  stat_compare_means(aes(label = "..p.adj.."), method = "t.test", label = "p.signif", ref.group = "D0", hide.ns = T, label.y = 0.8) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_text(colour = "black", size = 12), axis.title.y = element_blank(), axis.title.x = element_text(size = 18, face = "bold"), legend.position = "none", strip.text.x = element_text(size = 14))

gg_non_glob <- ggplot(treat_GSVA_nonglob, aes(x = Visit, y = value, group = Visit)) +
  geom_boxplot(aes(fill = Visit), alpha = 0.7, outlier.shape = NA) +
  geom_point(position = position_dodge(0.2), colour = "#377EB8") +
  geom_line(alpha = 0.6, position = position_dodge(0.2)) +
  theme_bw() +
  xlab("Group") +
  ylab("Score") +
  facet_wrap(~ variable, ncol = 1) +
  stat_compare_means(aes(label = "..p.adj.."), method = "t.test", label = "p.signif", ref.group = "D0", hide.ns = T, label.y = 0.8) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_text(colour = "black", size = 12), axis.title.y = element_blank(), axis.title.x = element_text(size = 18, face = "bold"), legend.position = "none", strip.text.x = element_text(size = 14))

gg_non_glob_leg <- ggplot(treat_GSVA_nonglob, aes(x = Visit, y = value, group = Visit)) +
  geom_boxplot(aes(fill = Visit), alpha = 0.7, outlier.shape = NA) +
  geom_point(position = position_dodge(0.2), colour = "#377EB8") +
  geom_line(alpha = 0.6, position = position_dodge(0.2)) +
  theme_bw() +
  xlab("Group") +
  ylab("Score") +
  facet_wrap(~ variable, ncol = 1) +
  stat_compare_means(aes(label = "..p.adj.."), method = "t.test", label = "p.signif", ref.group = "D0", hide.ns = T, label.y = 0.8) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_text(colour = "black", size = 12), axis.title.y = element_text(size = 18, face = "bold"), axis.title.x = element_text(size = 18, face = "bold"), legend.text = element_text(size = 14), legend.title = element_text(face = "bold", size = 16), legend.position = "bottom", strip.text.x = element_text(size = 14))

extract_legend <- function(my_ggp) {
  step1 <- ggplot_gtable(ggplot_build(my_ggp))
  step2 <- which(sapply(step1$grobs, function(x) x$name) == "guide-box")
  step3 <- step1$grobs[[step2]]
  return(step3)
}

shared_legend <- extract_legend(gg_non_glob_leg)

gg_grid <- grid.arrange(arrangeGrob(gg_glob, gg_non_glob, ncol = 2),
             shared_legend, nrow = 2, heights = c(10, 0.5))
ggsave(gg_grid, filename = "Fig 4D.pdf", device = cairo_pdf, 
       width = 12, height = 8, units = "in")


dup_counts <- as.data.frame(dge_dups$counts)
gene_id <- rownames(dup_counts)
my_genes <- ensembl[match(gene_id, ensembl$gene_id),]
my_genes <- my_genes[!duplicated(my_genes$gene_name),]
dup_counts <- dup_counts[rownames(dup_counts) %in% my_genes$gene_id,]
my_genes <- my_genes[my_genes$gene_id %in% rownames(dup_counts),]
dup_counts$gene_names <- my_genes$gene_name
dup_counts <- as.data.frame(dup_counts %>% group_by(gene_names) %>% summarise_if(is.numeric, sum))
dup_counts <- dup_counts[c(2:nrow(dup_counts)),]
rownames(dup_counts) <- dup_counts$gene_names
dup_counts <- dup_counts[2:nrow(dup_counts),2:ncol(dup_counts)]

glob_genes_all <- c(glob_genes$D0_up, glob_genes$D0_down, glob_genes$M6_up, glob_genes$M6_down, glob_genes$Both_up, glob_genes$Both_down, glob_genes$D0, glob_genes$M6, glob_genes$Both)
'%notin%' <- Negate('%in%')
dup_counts <- dup_counts[rownames(dup_counts) %notin% glob_genes_all,]

se_dups <- SummarizedExperiment(assays = list(counts = as.matrix(dup_counts)), colData = meta_dup)
se_dups <- mkAssay(se_dups, input_name = "counts", log = TRUE, counts_to_CPM = TRUE)

treat_GSVA <- runTBsigProfiler(se_dups, useAssay = "log_cpm", signatures = treat_sigs, algorithm = "GSVA")
treat_GSVA <- as.data.frame(colData(treat_GSVA))
treat_GSVA <- treat_GSVA[c(colnames(dge_dups_D0),colnames(dge_dups_M6)),]
treat_GSVA_glob <- treat_GSVA[grep("Depleted",treat_GSVA$Globin),]
treat_GSVA_nonglob <- treat_GSVA[grep("Non-depleted",treat_GSVA$Globin),]
treat_GSVA_glob <- treat_GSVA_glob[,c(3,109:ncol(treat_GSVA_glob))]
treat_GSVA_glob <- reshape2::melt(treat_GSVA_glob, id.vars = c("Visit"))
treat_GSVA_nonglob <- treat_GSVA_nonglob[,c(3,109:ncol(treat_GSVA_nonglob))]
treat_GSVA_nonglob <- reshape2::melt(treat_GSVA_nonglob, id.vars = c("Visit"))

treat_GSVA_glob <- treat_GSVA_glob[grep("Bloom_RES_558|Tabone_EarlyRESP_TB25",treat_GSVA_glob$variable),]
treat_GSVA_nonglob <- treat_GSVA_nonglob[grep("Bloom_RES_558|Tabone_EarlyRESP_TB25",treat_GSVA_nonglob$variable),]

ggplot(treat_GSVA_glob, aes(x = Visit, y = value, group = Visit)) +
  geom_boxplot(aes(fill = Visit), alpha = 0.7, outlier.shape = NA) +
  geom_point(position = position_dodge(0.2)) +
  geom_line(alpha = 0.6, position = position_dodge(0.2)) +
  theme_bw() +
  xlab("Group") +
  ylab("Score") +
  facet_wrap(~ variable) +
  stat_compare_means(aes(label = "..p.adj.."), method = "t.test", label = "p.signif", ref.group = "D0", hide.ns = T, label.y = 0.8) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_text(colour = "black", size = 12), axis.title.y = element_text(size = 18, face = "bold"), axis.title.x = element_text(size = 18, face = "bold"), legend.text = element_text(size = 14), legend.title = element_text(face = "bold", size = 16), strip.text.x = element_text(size = 14))

ggplot(treat_GSVA_nonglob, aes(x = Visit, y = value, group = Visit)) +
  geom_boxplot(aes(fill = Visit), alpha = 0.7, outlier.shape = NA) +
  geom_point(position = position_dodge(0.2)) +
  geom_line(alpha = 0.6, position = position_dodge(0.2)) +
  theme_bw() +
  xlab("Group") +
  ylab("Score") +
  facet_wrap(~ variable) +
  stat_compare_means(aes(label = "..p.adj.."), method = "t.test", label = "p.signif", ref.group = "D0", hide.ns = T, label.y = 0.8) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_text(colour = "black", size = 12), axis.title.y = element_text(size = 18, face = "bold"), axis.title.x = element_text(size = 18, face = "bold"), legend.text = element_text(size = 14), legend.title = element_text(face = "bold", size = 16), strip.text.x = element_text(size = 14))

D0_M6_GSVA_df <- as.data.frame(colData(D0_M6_GSVA))
"Bloom_RES_558|Tabone_EarlyRESP_TB25"

treat_GSVA_nd <- treat_GSVA[treat_GSVA$Globin == "Non-depleted",]
treat_GSVA$D0_down <- D0_M6_GSVA_df[rownames(treat_GSVA),]$D0_down
treat_GSVA$D0_up <- D0_M6_GSVA_df[rownames(treat_GSVA),]$D0_up
treat_GSVA$M6_down <- D0_M6_GSVA_df[rownames(treat_GSVA),]$M6_down
treat_GSVA$M6_up <- D0_M6_GSVA_df[rownames(treat_GSVA),]$M6_up

treat_GSVA_df <- treat_GSVA[,c(110,113)]
treat_GSVA_df <- reshape2::melt(treat_GSVA_df)
treat_GSVA_df$sex <- c(rep(treat_GSVA$Sex,2))
treat_GSVA_df$globin <- c(rep(treat_GSVA$Globin,2))
treat_GSVA_df$visit <- c(rep(treat_GSVA$Visit,2))
treat_GSVA_df$D0_down <- c(rep(treat_GSVA$D0_down,2))
treat_GSVA_df$D0_up <- c(rep(treat_GSVA$D0_up,2))
treat_GSVA_df$M6_down <- c(rep(treat_GSVA$M6_down,2))
treat_GSVA_df$M6_up <- c(rep(treat_GSVA$M6_up,2))

treat_GSVA_df_nd <- treat_GSVA_df[treat_GSVA_df$globin == "Non-depleted",]
treat_GSVA_df_gd <- treat_GSVA_df[treat_GSVA_df$globin == "Globin-depleted",]

ggplot(treat_GSVA_df_nd, aes(x = visit, y = value)) +
  geom_boxplot(aes(fill = visit), alpha = 0.7, outlier.shape = NA) +
  geom_point(aes(colour = sex)) +
  #scale_colour_continuous(type = "viridis") +
  theme_bw() +
  xlab("Group") +
  ylab("Score") +
  facet_wrap(~ variable) +
  stat_compare_means(aes(label = "..p.adj.."), method = "t.test", label = "p.signif", ref.group = "D0", hide.ns = T, label.y = 0.8) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_text(colour = "black", size = 12), axis.title.y = element_text(size = 18, face = "bold"), axis.title.x = element_text(size = 18, face = "bold"), legend.text = element_text(size = 14), legend.title = element_text(face = "bold", size = 16), strip.text.x = element_text(size = 14))

treat_GSVA_df <- as.data.frame(t(treat_GSVA[,c(109:113)]))
treat_GSVA_df_nd <- treat_GSVA_df[,rownames(treat_GSVA[treat_GSVA$Globin == "Non-depleted",])]
treat_GSVA_df_nd[] <- lapply(treat_GSVA_df_nd, as.numeric)
treat_GSVA_df_gd <- treat_GSVA_df[,rownames(treat_GSVA[treat_GSVA$Globin == "Depleted",])]
treat_GSVA_df_gd[] <- lapply(treat_GSVA_df_gd, as.numeric)

cor_mat2 <- cor(treat_GSVA_df_nd[],treat_GSVA_df_gd[])
round_df <- function(x, digits) {
  # round all numeric variables
  # x: data frame 
  # digits: number of digits to round
  numeric_columns <- sapply(x, mode) == 'numeric'
  x[numeric_columns] <-  round(x[numeric_columns], digits)
  x
}

cor_mat2 <- round_df(cor_mat2, 2)

melted_cor_mat2 <- melt(cor_mat2)

#get lower triangle of the correlation matrix
get_lower_tri<-function(cormat){
  cormat[upper.tri(cormat)] <- NA
  return(cormat)
}
lower_tri2 <- get_lower_tri(cor_mat2)

melted_cor_mat2 <- melt(lower_tri2, na.rm = TRUE)
rownames(melted_cor_mat2) <- c()

cor_pairs2 <- melted_cor_mat2[c(1,24,47,70,93,116,139,162,185,208,231,254,277,300,323,346,369,392,415,438,461,484),]

rownames(libs_df_dups) <- libs_df_dups$sample 
libs_df_dups_D0M6 <- libs_df_dups[grep("D0|M6",meta_dup$Visit),]
colnames(libs_df_dups_D0M6)[5] <- "Var1"
dup_cor_lib2 <- merge(cor_pairs2, libs_df_dups_D0M6, by = "Var1")
dup_cor_lib2$sex <- meta_dup[dup_cor_lib2$Var1,]$Sex
dup_cor_lib2$sex <- factor(dup_cor_lib2$sex)

meta_dup_cor2 <- meta_dup[levels(dup_cor_lib2$Var2),]
meta_dup_cor2 <- meta_dup_cor2[,c(2,29)]
meta_dup_cor2$Var2 <- rownames(meta_dup_cor2)
rownames(dup_cor_lib2) <- dup_cor_lib2$Var2
glob_cor2 <- merge(dup_cor_lib2, meta_dup_cor2, by = "Var2")
glob_cor2 <- glob_cor2[c(2:9,11:18,20),]
D0_M6_GSVA_df$group <- rownames(D0_M6_GSVA_df)
colnames(D0_M6_GSVA_df)[1] <- "Var1"
glob_cor2 <- merge(glob_cor2, D0_M6_GSVA_df[,c(1,8:13)], by = "Var1")

glob_cor2[glob_cor2$D0_down<(-0.5)&glob_cor2$D0_up>(0.25),]$D0 <- "fail"
glob_cor2[glob_cor2$D0!="fail",]$D0 <- "pass"
glob_cor2[glob_cor2$M6_down<(-0.25)&glob_cor2$M6_up>(0.25),]$M6 <- "fail"
glob_cor2[glob_cor2$M6!="fail",]$M6 <- "pass"

dim(glob_cor2[glob_cor2$D0=="fail"|glob_cor2$M6=="fail",])
glob_cor2[glob_cor2$D0=="fail"&glob_cor2$M6=="fail",]$sample <- "fail"
glob_cor2[glob_cor2$sample!="fail",]$sample <- "pass"

ggplot(glob_cor2) +
  geom_jitter(aes(x = sample, y = value))

ggscatter(glob_cor2, x = "gd_libs", y = "value", color = "value", shape = "D0",
          cor.coef = TRUE,
          cor.coeff.args = list(method = "pearson", label.y.npc = "bottom",  size = 6),
          xlab = "Percentage of reads corresponding to globin mRNA", ylab = "Pearson correlation with globin-depleted duplicate") +
  geom_smooth(method=lm) +
  theme(axis.text.y = element_text(colour = "black", size = 12), axis.title.y = element_text(size = 18, face = "bold"), axis.title.x = element_text(size = 18, face = "bold"), legend.text = element_text(size = 14), legend.title = element_text(face = "bold", size = 16), strip.text.x = element_text(size = 14))

#correlate with WBC globin levels
meta_D0 <- meta_dup[rownames(dge_dups_D0$samples),]
meta_D0_glob <- meta_D0[,c(10,29)]

#NOISeq analysis
pre_gd_counts_dups <- pre_gd_counts[,dups]
#pre_gd_counts_dups <- pre_gd_counts_dups[!rownames(pre_gd_counts_dups)%in%pc_genes,]
pre_gd_gene_id <- rownames(pre_gd_counts_dups)
pre_gd_genes <- ensembl[match(pre_gd_gene_id, ensembl$gene_id),]
rownames(pre_gd_genes) <- rownames(pre_gd_counts_dups)
pre_gd_genes <- pre_gd_genes[,c(2:9)]

pre_gd_length <- pre_gd_genes$length
names(pre_gd_length) <- rownames(pre_gd_genes)
pre_gd_gc <- pre_gd_genes$GC_content
names(pre_gd_gc) <- rownames(pre_gd_genes)
pre_gd_biotypes <- pre_gd_genes$biotype
names(pre_gd_biotypes) <- rownames(pre_gd_genes)
pre_gd_chroms <- pre_gd_genes[,c(6:8)]
my_factors <- as.data.frame(meta_dup$Globin)
colnames(my_factors) <- c("Globin")

pre_gd_data <- readData(data = pre_gd_counts_dups, length = pre_gd_length, gc = pre_gd_gc, biotype = pre_gd_biotypes, chromosome = pre_gd_chroms, factors = my_factors)
pre_gd_explodata <- dat(pre_gd_data, type = "biodetection")
explo.plot(pre_gd_explodata, plottype = c("comparison"), samples = c(1, 31))
my_nicedata <- dat2save(my_explodata)
pre_gdcountsbio = dat(pre_gd_data, factor = NULL, type = "countsbio")
explo.plot(pre_gdcountsbio, toplot = 1, samples = NULL, plottype = "barplot")
explo.plot(pre_gdcountsbio, toplot = 1, samples = NULL, plottype = "boxplot")

pre_gd_biotypes <- as.data.frame(pre_gd_explodata@dat$biotables$PCM035_dup[2,])
colnames(pre_gd_biotypes) <- c("PCM035_dup")
pre_gd_biotypes$PCM035 <- pre_gd_explodata@dat$biotables$PCM035[2,]
pre_gd_biotypes$PCM036_dup <- pre_gd_explodata@dat$biotables$PCM036_dup[2,]
pre_gd_biotypes$PCM036 <- pre_gd_explodata@dat$biotables$PCM036[2,]
pre_gd_biotypes$PCM071 <- pre_gd_explodata@dat$biotables$PCM071[2,]
pre_gd_biotypes$PCM071_dup <- pre_gd_explodata@dat$biotables$PCM071_dup[2,]
pre_gd_biotypes$PCM072 <- pre_gd_explodata@dat$biotables$PCM072[2,]
pre_gd_biotypes$PCM072_dup <- pre_gd_explodata@dat$biotables$PCM072_dup[2,]
pre_gd_biotypes$PCM073 <- pre_gd_explodata@dat$biotables$PCM073[2,]
pre_gd_biotypes$PCM073_dup <- pre_gd_explodata@dat$biotables$PCM073_dup[2,]
pre_gd_biotypes$PCM074 <- pre_gd_explodata@dat$biotables$PCM074[2,]
pre_gd_biotypes$PCM074_dup <- pre_gd_explodata@dat$biotables$PCM074_dup[2,]
pre_gd_biotypes$PCM075 <- pre_gd_explodata@dat$biotables$PCM075[2,]
pre_gd_biotypes$PCM075_dup <- pre_gd_explodata@dat$biotables$PCM075_dup[2,]
pre_gd_biotypes$PCM076 <- pre_gd_explodata@dat$biotables$PCM076[2,]
pre_gd_biotypes$PCM076_dup <- pre_gd_explodata@dat$biotables$PCM076_dup[2,]
pre_gd_biotypes$PCM077 <- pre_gd_explodata@dat$biotables$PCM077[2,]
pre_gd_biotypes$PCM077_dup <- pre_gd_explodata@dat$biotables$PCM077_dup[2,]
pre_gd_biotypes$PCM078 <- pre_gd_explodata@dat$biotables$PCM078[2,]
pre_gd_biotypes$PCM078_dup <- pre_gd_explodata@dat$biotables$PCM078_dup[2,]
pre_gd_biotypes$PCM079 <- pre_gd_explodata@dat$biotables$PCM079[2,]
pre_gd_biotypes$PCM079_dup <- pre_gd_explodata@dat$biotables$PCM079_dup[2,]
pre_gd_biotypes$PCM080 <- pre_gd_explodata@dat$biotables$PCM080[2,]
pre_gd_biotypes$PCM080_dup <- pre_gd_explodata@dat$biotables$PCM080_dup[2,]
pre_gd_biotypes$PCM081 <- pre_gd_explodata@dat$biotables$PCM081[2,]
pre_gd_biotypes$PCM081_dup <- pre_gd_explodata@dat$biotables$PCM081_dup[2,]
pre_gd_biotypes$PCM082 <- pre_gd_explodata@dat$biotables$PCM082[2,]
pre_gd_biotypes$PCM082_dup <- pre_gd_explodata@dat$biotables$PCM082_dup[2,]
pre_gd_biotypes$PCM083 <- pre_gd_explodata@dat$biotables$PCM083[2,]
pre_gd_biotypes$PCM083_dup <- pre_gd_explodata@dat$biotables$PCM083_dup[2,]
pre_gd_biotypes$PCM084 <- pre_gd_explodata@dat$biotables$PCM084[2,]
pre_gd_biotypes$PCM084_dup <- pre_gd_explodata@dat$biotables$PCM084_dup[2,]
pre_gd_biotypes$PCM085 <- pre_gd_explodata@dat$biotables$PCM085[2,]
pre_gd_biotypes$PCM085_dup <- pre_gd_explodata@dat$biotables$PCM085_dup[2,]
pre_gd_biotypes$PCM086 <- pre_gd_explodata@dat$biotables$PCM086[2,]
pre_gd_biotypes$PCM086_dup <- pre_gd_explodata@dat$biotables$PCM086_dup[2,]
pre_gd_biotypes$PCM087 <- pre_gd_explodata@dat$biotables$PCM087[2,]
pre_gd_biotypes$PCM087_dup <- pre_gd_explodata@dat$biotables$PCM087_dup[2,]
pre_gd_biotypes$PCM088 <- pre_gd_explodata@dat$biotables$PCM088[2,]
pre_gd_biotypes$PCM088_dup <- pre_gd_explodata@dat$biotables$PCM088_dup[2,]
pre_gd_biotypes$PCM089 <- pre_gd_explodata@dat$biotables$PCM089[2,]
pre_gd_biotypes$PCM089_dup <- pre_gd_explodata@dat$biotables$PCM089_dup[2,]
pre_gd_biotypes$PCM090 <- pre_gd_explodata@dat$biotables$PCM090[2,]
pre_gd_biotypes$PCM090_dup <- pre_gd_explodata@dat$biotables$PCM090_dup[2,]
pre_gd_biotypes$PCM091 <- pre_gd_explodata@dat$biotables$PCM091[2,]
pre_gd_biotypes$PCM091_dup <- pre_gd_explodata@dat$biotables$PCM091_dup[2,]
pre_gd_biotypes$PCM092 <- pre_gd_explodata@dat$biotables$PCM092[2,]
pre_gd_biotypes$PCM092_dup <- pre_gd_explodata@dat$biotables$PCM092_dup[2,]
pre_gd_biotypes$PCM093 <- pre_gd_explodata@dat$biotables$PCM093[2,]
pre_gd_biotypes$PCM093_dup <- pre_gd_explodata@dat$biotables$PCM093_dup[2,]
pre_gd_biotypes$PCM094 <- pre_gd_explodata@dat$biotables$PCM094[2,]
pre_gd_biotypes$PCM094_dup <- pre_gd_explodata@dat$biotables$PCM094_dup[2,]
pre_gd_biotypes$PCM095 <- pre_gd_explodata@dat$biotables$PCM095[2,]
pre_gd_biotypes$PCM095_dup <- pre_gd_explodata@dat$biotables$PCM095_dup[2,]
pre_gd_biotypes$PCM096 <- pre_gd_explodata@dat$biotables$PCM096[2,]
pre_gd_biotypes$PCM096_dup <- pre_gd_explodata@dat$biotables$PCM096_dup[2,]
pre_gd_biotypes$PCM097 <- pre_gd_explodata@dat$biotables$PCM097[2,]
pre_gd_biotypes$PCM097_dup <- pre_gd_explodata@dat$biotables$PCM097_dup[2,]
pre_gd_biotypes$PCM098 <- pre_gd_explodata@dat$biotables$PCM098[2,]
pre_gd_biotypes$PCM098_dup <- pre_gd_explodata@dat$biotables$PCM098_dup[2,]
pre_gd_biotypes$PCM099 <- pre_gd_explodata@dat$biotables$PCM099[2,]
pre_gd_biotypes$PCM099_dup <- pre_gd_explodata@dat$biotables$PCM099_dup[2,]
pre_gd_biotypes$PCM100 <- pre_gd_explodata@dat$biotables$PCM100[2,]
pre_gd_biotypes$PCM100_dup <- pre_gd_explodata@dat$biotables$PCM100_dup[2,]
pre_gd_biotypes <- pre_gd_biotypes[,c(1,3,5,7,9,11,13,15,17,19,21,23,25,27,29,31,33,35,37,39,41,43,45,47,49,51,53,55,57,2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40,42,44,46,48,50,52,54,56,58)]
pre_gd_biotypes <- pre_gd_biotypes[c(1:23,25:34,38,39),]
pre_gd_biotypes$biotypes <- rownames(pre_gd_biotypes)
pre_gd_bio_melt <- melt(pre_gd_biotypes, id.vars = "biotypes")
pre_gd_bio_melt$globin <- c(rep("Globin-depleted",1015),rep("Non-depleted",1015))
pre_gd_bio_melt$biotypes <- factor(pre_gd_bio_melt$biotypes, levels = c("protein_coding","lncRNA","processed_pseudogene","unprocessed_pseudogene","misc_RNA","transcribed_unprocessed_pseudogene","TEC","snRNA","miRNA","transcribed_processed_pseudogene","snoRNA","IG_V_gene","TR_V_gene","IG_V_pseudogene","transcribed_unitary_pseudogene","rRNA_pseudogene","TR_J_gene","unitary_pseudogene","IG_D_gene","polymorphic_pseudogene","TR_V_pseudogene","rRNA","IG_J_gene","IG_C_gene","scaRNA","pseudogene","IG_C_pseudogene","TR_C_gene","TR_D_gene","IG_J_pseudogene","TR_J_pseudogene","ribozyme","sRNA","IG_pseudogene","scRNA"))

ggplot(data = pre_gd_bio_melt, 
       aes(x = globin, y = value, fill=biotypes))+
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75), alpha = 0.5)+
  geom_jitter(aes(color = biotypes, alpha = 0.5))+
  facet_wrap(~biotypes, scales = "free_y", ncol = 5)+
  xlab("Globin status") +
  ylab("Percentage of reads corresponding to gene biotype") +
  stat_compare_means(aes(group = globin, label = "..p.adj..", size = 12), method = "t.test", paired = T, label = "p.signif", hide.ns = T) +
  theme_minimal() +
  theme(legend.position = "none",
        strip.text = element_text(size = 12),
        axis.text.x = element_text(angle = 25, vjust = 1, 
                                   size = 12, hjust = 1),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14))

