###############
# Date of creation 16 Jan 2026
# Yunjiao Lu
###############
library(VennDiagram)
library(ggvenn)
library(data.table)
library(matrixStats)
dir_path <- "/Users/ylu/Documents/MOZART"
########################### Read Raw Counts ###############################
files_samples <- list.files(paste0(dir_path, "/data/RNAseq/merged_feature_counts/"))
length(grep("genes", files_samples))
files_names_idx <- sapply(files_samples, function(x) {length(grep("genes", x))})
files_names <- files_samples[files_names_idx==1]

file_name <- files_names[1]
genes_read <- read.table(file = paste0(dir_path,  "/data/RNAseq/merged_feature_counts/20170081.genes.results"),
	header = TRUE)
iso_read <- read.table(file = paste0(dir_path,  "/data/RNAseq/merged_feature_counts/20170081.isoforms.results"),
	header = TRUE)

# sanitary check - gene count = sum of the counts of all isoforms
gene_id <- "ENSG00000000003.16"
genes_read[genes_read$gene_id == gene_id, ]
iso_read[iso_read$gene_id == gene_id, ]
mean(iso_read[iso_read$gene_id == gene_id, "effective_length"]) # the length or expected length of the gene is not exactly the mean or the median of that of the isoforms
sum(iso_read[iso_read$gene_id == gene_id, "expected_count"]) # sum of the expected counts of isoforms equals to that of the gene
sum(iso_read[iso_read$gene_id == gene_id, "TPM"])# sum of the TPM of isoforms equals to that of the gene
sum(iso_read[iso_read$gene_id == gene_id, "FPKM"])# sum of the FPKM of isoforms equals to that of the gene


# Calculate TPM for one sample
calcul_TPM <- function(count_v, length_v){
	RPK <- count_v / length_v
	TPM <- RPK/sum(RPK[!is.na(RPK) & !is.infinite(RPK)])*1000000
	TPM[is.na(RPK) | is.infinite(RPK)] <- 0
	return(TPM)
}
count_v <- genes_read$expected_count
length_v <- genes_read$effective_length # use effective_length instead of length
TPM_sample_i <- calcul_TPM(count_v, length_v)
test<- data.frame(TPM_sample_i, genes_read$TPM)
genes_read[test[, 1]-test[,2]>0.5,]
test[test[, 1]-test[,2]>0.5,]


# merge TPM all patients
ref_genes <- read.table(paste0(dir_path, "/data/gencode_gene_names.tsv"),
	sep = ";",
	header = FALSE
)
ref_genes2 <- data.frame(matrix(NA, nrow = dim(ref_genes)[1], ncol = dim(ref_genes)[2]))
colnames(ref_genes2) <- c("gene_id", "gene_name")
gene_symbol_extract <- function(str){
	str <- trimws(str)
	startp <- unlist(gregexpr(" ", str))
	new_name <- substr(str, startp+1 , nchar(str))
	return(new_name)
}
ref_genes2$gene_id <- sapply(ref_genes$V1, gene_symbol_extract)
ref_genes2$gene_name <- sapply(ref_genes$V2, gene_symbol_extract)

rename_sample2 <- function(str){
	startp <- unlist(gregexpr(".genes.results", str))
	new_name <- substr(str, 1, startp - 1)
	return(new_name)
}
# exp_rsem_tpm <- data.frame(matrix(NA, nrow = dim(ref_genes2)[1], ncol = length(files_names)))
# colnames(exp_rsem_tpm) <- sapply(files_names, rename_sample2)
exp_rsem_tpm <- ref_genes2
rename_genes <- function(str){
	startp <- unlist(gregexpr("\\.", str))
	new_name <- substr(str, 1, startp - 1)
	return(new_name)
}
exp_rsem_tpm$gene_id <- sapply(exp_rsem_tpm$gene_id, rename_genes)
for (file in files_names){
	samp <- rename_sample2(file)
	genes_count_samp_i <- read.table(file = paste0(dir_path,  "/data/RNAseq/merged_feature_counts/", file),
		header = TRUE)
	genes_count_samp_i$gene_id <- sapply(genes_count_samp_i$gene_id, rename_genes)
	exp_rsem_tpm <- merge(exp_rsem_tpm, genes_count_samp_i[, c("gene_id", "TPM")], by="gene_id", all.x=TRUE)
	colnames(exp_rsem_tpm)[dim(exp_rsem_tpm)[2]] <- samp
}


########## Compare Feature counts et RSEM counts ##########
counts_feature <- read.table(paste0(dir_path, "/data/RNAseq/merged_feature_counts/read_counts/rawCounts.txt"),
	header = TRUE,
	row.names = 1)
counts_rsem <- read.table("/Users/ylu/Documents/MOZART/data/RNAseq/merged_feature_counts/read_counts/rsem.counts.matrix",
	header = TRUE,
	row.names = 1)


rename_sample <- function(str){
	startp <- unlist(gregexpr("_rawCounts.txt", str))
	new_name <- substr(str, 1, startp - 1)
	return(new_name)
}
rename_sample2 <- function(str){
	startp <- unlist(gregexpr(".genes.results", str))
	new_name <- substr(str, 1, startp - 1)
	return(new_name)
}
new_sample_name_fc <- sapply(colnames(counts_feature), rename_sample) 
new_sample_name_rsem <- sapply(colnames(counts_rsem), rename_sample2) 
colnames(counts_feature) <- new_sample_name_fc
colnames(counts_rsem) <- new_sample_name_rsem

rename_genes <- function(str){
	startp <- unlist(gregexpr("\\.", str))
	new_name <- substr(str, 1, startp - 1)
	return(new_name)
}
gene_symbol_fc <- sapply(rownames(counts_feature), rename_genes)
gene_symbol_rsem <- sapply(rownames(counts_rsem), rename_genes)

dl_venn <- list(
	FC = gene_symbol_fc,
	RSEM = gene_symbol_rsem
)
pdf(file = "/Users/ylu/Documents/MOZART/graphs/venn_gene_symbols.pdf")
ggvenn(
  dl_venn, 
  fill_color = c("#0073C2FF", "#EFC000FF"),
  stroke_size = 0.5, set_name_size = 4
  )
dev.off()

counts_feature$gene_symbol <- gene_symbol_fc
counts_rsem$gene_symbol <- gene_symbol_rsem
common_genes <- gene_symbol_rsem[is.element(gene_symbol_rsem, gene_symbol_fc)]

cor_counts <- function(sample){
	exp_fc <- counts_feature[is.element(counts_feature$gene_symbol, common_genes), c("gene_symbol", sample)]
	exp_rsem <- counts_rsem[is.element(counts_rsem$gene_symbol, common_genes),  c("gene_symbol", sample)]
	exp_sample_i <- merge(exp_fc, exp_rsem, by = "gene_symbol")
	colnames(exp_sample_i) <- c("gene_symbol", "exp_fc", "exp_rsem")
	coef_cor <- cor(exp_sample_i$exp_fc, exp_sample_i$exp_rsem, method = "spearman")
	return(coef_cor)
}
coef_cor2 <- cor(exp_sample_i$exp_fc, exp_sample_i$exp_rsem, method = "pearson")
coef_cor_sp <- sapply(new_sample_name_fc, cor_counts)
# > summary(coef_cor_sp)
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#  0.8718  0.8814  0.8836  0.8835  0.8857  0.8926 
pdf("/Users/ylu/Documents/MOZART/graphs/cor_counts_twomethods.pdf")
	plot(exp_sample_i$exp_fc, exp_sample_i$exp_rsem, xlim = c(0, 20000), ylim = c(0, 20000))
	text(5000, 18000, paste0("Spearman = ", round(coef_cor,4), "\n", "Pearnson = ", round(coef_cor2,4)))
dev.off()

pdf("/Users/ylu/Documents/MOZART/graphs/cor_counts_twomethods_full.pdf")
	plot(exp_sample_i$exp_fc, exp_sample_i$exp_rsem,)
	text(300000, 390000, paste0("Spearman = ", round(coef_cor,4), "\n", "Pearnson = ", round(coef_cor2,4)))
dev.off()

#####################################
# PCA   #
#####################################

# Extract gene lengths
library(rtracklayer)
library(data.table)

# gtf <- import("/Users/ylu/Documents/MOZART/data/gencode.v49lift37.basic.annotation.gtf")
# genes <- gtf[gtf$type == "gene"]
# dt <- as.data.table(genes)
# gene_length <- dt[, .(
#   gene_id = gene_id,
#   length = end - start + 1
# )]
# gene_length <- unique(gene_length, by = "gene_id")
# write.table(gene_length, "/Users/ylu/Documents/MOZART/data/gene_length.txt", row.names = FALSE)
# 
# df_gene_length <- read.table(
# 	"/Users/ylu/Documents/MOZART/data/gene_length.txt",
# 	header = TRUE
# )
# df_gene_length$gene_id2 <- sapply(df_gene_length$gene_id, rename_genes)
# 
# counts_rsem$gene_symbol[is.element(counts_rsem$gene_symbol, df_gene_length$gene_id2)]
# 
# # Calculate TPM from raw count
# df_gene_length$length_kb <- df_gene_length$length / 1000
# # Ensure order matches
# counts_rsem <- merge(counts_rsem, df_gene_length[, c("length_kb", "gene_id2")], by.x = "gene_symbol", by.y = "gene_id2", all = FALSE)
# 
# gene_length_kb <- counts_rsem$length_kb
# # TPM calculation
# tpm_rsem <- apply(subset(counts_rsem, select = -c(gene_symbol, length_kb)), 2, function(x) {
#   rate <- x / gene_length_kb
#   rate / sum(rate) * 1e6
# })
# 
# # tpm_rsem <- data.frame(tpm_rsem, counts_rsem$gene_symbol)



# PCA for subgroups discovary
log_tpm_rsem <- log2(subset(exp_rsem_tpm, select = -c(gene_id, gene_name)) + 0.01)
gene_var <- rowVars(as.matrix(log_tpm_rsem))
df_log_tpm_rsem <- data.frame(log_tpm_rsem, exp_rsem_tpm$gene_name)
# Use top 500–2000 most variable genes
top_genes <- order(gene_var, decreasing = TRUE)[1:20000]


log_tpm_rsem_top <- log_tpm_rsem[top_genes, ]



pca <- prcomp(t(log_tpm_rsem_top), 
	center = TRUE, 
	scale. = TRUE)

pca_df <- data.frame(
  Sample = colnames(log_tpm_rsem_top),
  PC1 = pca$x[, 1],
  PC2 = pca$x[, 2]
)

library(ggplot2)


pdf("/Users/ylu/Documents/MOZART/graphs/pca.pdf")
ggplot(pca_df, aes(PC1, PC2)) +
  geom_point(size = 4) +
  theme_classic() +
  labs(
    title = "PCA of RNA-seq Samples (log2 TPM)",
    x = paste0(
      "PC1 (",
      round(summary(pca)$importance[2, 1] * 100, 1),
      "%)"
    ),
    y = paste0(
      "PC2 (",
      round(summary(pca)$importance[2, 2] * 100, 1),
      "%)"
    )
  )
dev.off()