###############
# Date of creation 16 Jan 2026
# Yunjiao Lu
###############
# library(VennDiagram)
# library(ggvenn)
library(data.table)
library(matrixStats)
library(xtable)
library(ggplot2)
dir_path <- "/Users/ylu/Documents/MOZART"
###########################    Functions     ###############################
source("script_mozart/functions.R")
# source("/Users/ylu/RPackages/GSVA-devel/R/gsva.R") 
# source("/Users/ylu/RPackages/GSVA-devel/R/ssgsea.R")
# source("/Users/ylu/RPackages/GSVA-devel/R/ssgseaParam.R")
rename_sample2 <- function(str){
	startp <- unlist(gregexpr(".genes.results", str))
	new_name <- substr(str, 1, startp - 1)
	return(new_name)
}

rename_sample <- function(str){
	startp <- unlist(gregexpr("_rawCounts.txt", str))
	new_name <- substr(str, 1, startp - 1)
	return(new_name)
}
rename_genes <- function(str){
	startp <- unlist(gregexpr("_", trimws(str)))
	new_name <- substr(str, 1, startp[1] - 1)
	return(new_name)
}
tidy_exp_data <- function(exp_data){ 
	# Combine rows where one gene name corresponding to two gene ID;
	# Return expression matrix with rownames as gene names
	exp_data2 <- subset(exp_data, select = -c(gene_id))
	exp_data2 <- as.data.table(exp_data2)
	exp_data2 <- exp_data2[, lapply(.SD, sum), by = gene_name]
	# These genes have two IDs
	#  [1] "AKAP17A" "ASMT"    "ASMTL"   "C4orf36" "CD99"    "CRLF2"   "CSF2RA" 
	#  [8] "DHRSX"   "DUSP13B" "FAM174C" "GTPBP6"  "HERC3"   "IL3RA"   "IL9R"   
	# [15] "KYAT1"   "P2RY8"   "PDE4C"   "PDE8B"   "PINX1"   "PLCXD1"  "POLR2J3"
	# [22] "PPP2R3B" "SHOX"    "SLC25A6" "VAMP7"   "WASH6P"  "ZBED1" 
	
	# > exp_data$gene_id[exp_data$gene_name == "AKAP17A"]
	# [1] "ENSG00000197976.13" "ENSG00000292343.2" 
	
	exp_data2 <- as.data.frame(exp_data2)
	rownames(exp_data2) <- exp_data2$gene_name
	exp_data2 <- subset(exp_data2, select = -gene_name)
	return(exp_data2)
}
########################### Read Raw Counts ###############################
files_samples <- list.files(paste0(dir_path, "/data/RNAseq/merged_feature_counts/"))
length(grep("genes", files_samples))
files_names_idx <- sapply(files_samples, function(x) {length(grep("genes", x))})
files_names <- files_samples[files_names_idx==1]

file_name <- files_names[1]
genes_read <- read.table(file = paste0(dir_path,  "/data/RNAseq/merged_feature_counts/20170081.genes.results"),
	header = TRUE)
# sum(genes_read$TPM) = 10^6 it is normalized over all genes (including lncRNA)
iso_read <- read.table(file = paste0(dir_path,  "/data/RNAseq/merged_feature_counts/20170081.isoforms.results"),
	header = TRUE)

# Read protein coding gene reference
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

# sanity check - gene count = sum of the counts of all isoforms
gene_id <- "ENSG00000000003.17_TSPAN6"
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


# merge TPM of all patients
exp_rsem_tpm <- ref_genes2[!is.element(ref_genes2$gene_name, c("LIMCH1", "A1BG", "ENSG00000255330")),]
# It takes 3-4 minutes to run the following block of code
start_time <- Sys.time()
for (file in files_names){
	samp <- rename_sample2(file)
	genes_count_samp_i <- read.table(file = paste0(dir_path,  "/data/RNAseq/merged_feature_counts/", file),
		header = TRUE)
	genes_count_samp_i$gene_id <- sapply(genes_count_samp_i$gene_id, rename_genes)
	exp_rsem_tpm <- merge(exp_rsem_tpm, genes_count_samp_i[, c("gene_id", "TPM")], by="gene_id", all.x=TRUE)
	colnames(exp_rsem_tpm)[dim(exp_rsem_tpm)[2]] <- samp
}
end_time <- Sys.time()
print(end_time - start_time)
write.table(exp_rsem_tpm, 
	"data/exp_rsem_tpm.txt",
	row.names = FALSE)

# merge FPKM of all patients
exp_rsem_fpkm <- ref_genes2[!is.element(ref_genes2$gene_name, c("LIMCH1", "A1BG", "ENSG00000255330")),]
# It takes 3-4 minutes to run the following block of code
start_time <- Sys.time()
for (file in files_names){
	samp <- rename_sample2(file)
	genes_count_samp_i <- read.table(file = paste0(dir_path,  "/data/RNAseq/merged_feature_counts/", file),
		header = TRUE)
	genes_count_samp_i$gene_id <- sapply(genes_count_samp_i$gene_id, rename_genes)
	exp_rsem_fpkm <- merge(exp_rsem_fpkm, genes_count_samp_i[, c("gene_id", "FPKM")], by="gene_id", all.x=TRUE)
	colnames(exp_rsem_fpkm)[dim(exp_rsem_fpkm)[2]] <- samp
}
end_time <- Sys.time()
print(end_time - start_time)
write.table(exp_rsem_fpkm, 
	"data/exp_rsem_fpkm.txt",
	row.names = FALSE)

# merge raw count of all patients
exp_count <- ref_genes2[!is.element(ref_genes2$gene_name, c("LIMCH1", "A1BG", "ENSG00000255330")),]
# It takes 3-4 minutes to run the following block of code
start_time <- Sys.time()
for (file in files_names){
	samp <- rename_sample2(file)
	genes_count_samp_i <- read.table(file = paste0(dir_path,  "/data/RNAseq/merged_feature_counts/", file),
		header = TRUE)
	genes_count_samp_i$gene_id <- sapply(genes_count_samp_i$gene_id, rename_genes)
	exp_count <- merge(exp_count, genes_count_samp_i[, c("gene_id", "expected_count")], by="gene_id", all.x=TRUE)
	colnames(exp_count)[dim(exp_count)[2]] <- samp
}
end_time <- Sys.time()
print(end_time - start_time)
write.table(exp_count, 
	"data/exp_count.txt",
	row.names = FALSE)
########## Compare Feature counts and RSEM counts ##########
counts_feature <- read.table(paste0(dir_path, "/data/RNAseq/merged_feature_counts/read_counts/rawCounts.txt"),
	header = TRUE,
	row.names = 1)
mrna_featureC <-  counts_feature[is.element(rownames(counts_feature), ref_genes2$gene_id),]

counts_rsem <- read.table("/Users/ylu/Documents/MOZART/data/RNAseq/merged_feature_counts/read_counts/rsem.counts.matrix",
	header = TRUE,
	row.names = 1)
rownames(counts_rsem) <- sapply(rownames(counts_rsem), rename_genes)

mrna_rsem <- counts_rsem[is.element(rownames(counts_rsem), ref_genes2$gene_id),]

new_sample_name_fc <- sapply(colnames(mrna_featureC), rename_sample) 
new_sample_name_rsem <- sapply(colnames(mrna_rsem), rename_sample2) 
colnames(mrna_featureC) <- new_sample_name_fc
colnames(mrna_rsem) <- new_sample_name_rsem

gene_symbol_fc <- rownames(mrna_featureC)
gene_symbol_rsem <- rownames(mrna_rsem)

mrna_featureC$gene_symbol <- gene_symbol_fc
mrna_rsem$gene_symbol <- gene_symbol_rsem
common_genes <- gene_symbol_rsem[is.element(gene_symbol_rsem, gene_symbol_fc)]

cor_counts <- function(sample){
	exp_fc <- mrna_featureC[, c(sample, "gene_symbol")]
	exp_rsem <- mrna_rsem[, c(sample, "gene_symbol")]
	exp_sample_i <- merge(exp_fc, exp_rsem, by = "gene_symbol")
	colnames(exp_sample_i) <- c("gene_symbol", "exp_fc", "exp_rsem")
	coef_cor <- cor(exp_sample_i$exp_fc, exp_sample_i$exp_rsem, method = "spearman")
	return(coef_cor)
}
coef_cor_sp <- sapply(new_sample_name_fc, cor_counts)
coef_cor_ps <- sapply(new_sample_name_fc, cor_counts)
coef_cor <- as.data.frame(rbind(summary(coef_cor_sp), summary(coef_cor_ps)))
rownames(coef_cor) <- c("Spearman", "Pearson")
print(xtable(coef_cor), type = "latex")
# > summary(coef_cor_sp)
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#  0.9480  0.9579  0.9592  0.9588  0.9599  0.9616 
# > summary(coef_cor_ps)
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#  0.8248  0.9580  0.9646  0.9616  0.9694  0.9967 
# > 
# sample = "SLS.040"
sample = "CB0915M"
exp_fc <- mrna_featureC[, c(sample, "gene_symbol")]
exp_rsem <- mrna_rsem[, c(sample, "gene_symbol")]
exp_sample_i <- merge(exp_fc, exp_rsem, by = "gene_symbol")
colnames(exp_sample_i) <- c("gene_symbol", "exp_fc", "exp_rsem")
coef_cor2 <- cor(exp_sample_i$exp_fc, exp_sample_i$exp_rsem, method = "pearson")

pdf(paste0("/Users/ylu/Documents/MOZART/graphs/cor_counts_",sample,".pdf"), width = 12, height = 6.5)
par(mfrow = c(1,2))
plot(exp_sample_i$exp_fc, exp_sample_i$exp_rsem, xlim = c(0, 20000), ylim = c(0, 20000), asp=1, xlab = "featureCount", ylab = "RSEM")
	text(5000, 18000, paste0("Spearman = ", round(coef_cor,4), "\n", "Pearnson = ", round(coef_cor2,4)))

plot(exp_sample_i$exp_fc, exp_sample_i$exp_rsem, asp=1, xlab = "featureCount", ylab = "RSEM")
	text(300000, 390000, paste0("Spearman = ", round(coef_cor,4), "\n", "Pearnson = ", round(coef_cor2,4)))
dev.off()

#####################################
# PCA for subgroups discovery
#####################################

# Extract gene lengths
library(rtracklayer)
library(data.table)

# log2 transformation to stablize variance
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
##########################################################################
#
# Comparison of combinations of normalization and transformation methods #
#
##########################################################################
library(edgeR)
library(DESeq2)
meth_within_sample_norm <- c("none", "tpm", "fpkm")
meth_bet_sample_norm <- c("none", "tmm", "rle")
meth_transf <- c("none", "asinh")
for (i in meth_within_sample_norm){
	if (i == "none") {
		exp_data <- read.table(
			"data/exp_count.txt",
			header = TRUE
		)
		exp_mat <- tidy_exp_data(exp_data)
		assign(paste0("exp_", i), exp_mat)	
	}else if(i == "tpm"){
		exp_data <- read.table(
			"data/exp_rsem_tpm.txt",
			header = TRUE
		)
		exp_mat <- tidy_exp_data(exp_data)
		assign(paste0("exp_", i), exp_mat)	
	}else if(i == "fpkm"){
		exp_data <- read.table(
			"data/exp_rsem_fpkm.txt",
			header = TRUE
		)
		exp_mat <- tidy_exp_data(exp_data)
		assign(paste0("exp_", i), exp_mat)	
	}
	for (j in meth_bet_sample_norm){
		if (j == "none") {
			assign(paste0("exp_", i, "_", j), get(paste0("exp_", i)))	
		}else if(j == "tmm"){
			# TMM edgeR
			exp_tmm <- calcNormFactors(get(paste0("exp_", i)), method = c("TMM"))
			assign(paste0("exp_", i, "_", j), exp_tmm)
		}else if(j == "rle"){
			# RLE edgeR
			exp_rle <- calcNormFactors(get(paste0("exp_", i)), method = c("RLE"))
			assign(paste0("exp_", i, "_", j), exp_rle)
		}
		for (k in meth_transf){
			if (k == "none") {
				assign(paste0("exp_", i, "_", j, "_", k), get(paste0("exp_", i, "_", j)))	
			}else if(k == "vst"){
				exp_vst <- varianceStabilizingTransformation(
					as.matrix(get(paste0("exp_", i, "_", j))), 
					blind = TRUE, 
					fitType = "parametric"
				)
				assign(paste0("exp_", i, "_", j, "_", k), exp_vst)
			}else if(k == "rlog"){
				exp_rlog <- rlogTransformation(
					as.matrix(get(paste0("exp_", i, "_", j))),
					blind = TRUE,
					intercept,
					betaPriorVar,
					fitType = "parametric"
				)
				assign(paste0("exp_", i, "_", j, "_", k), exp_rlog)
			}else if(k == "asinh"){
				exp_asinh <- log(get(paste0("exp_", i, "_", j)) + sqrt(get(paste0("exp_", i, "_", j))^2+1))
				assign(paste0("exp_", i, "_", j, "_", k), exp_asinh)
			}
			write.table(get(paste0("exp_", i, "_", j, "_", k)),
				paste0("data/comp_norm_trans/exp_", i, "_", j, "_", k, ".txt")
			)
		}
	}
}
     dds <- makeExampleDESeqDataSet(m=6,betaSD=1)
     rld <- rlog(dds)
     dists <- dist(t(assay(rld)))
     # plot(hclust(dists))
     
exp_data <- read.table(
	"data/exp_count.txt",
	header = TRUE
)
exp_mat1 <- tidy_exp_data(exp_data)
exp_mat <- apply(exp_mat1, 2, as.integer)
rownames(exp_mat) <- rownames(exp_mat1)

meth_transf <- c("vst")
for (k in meth_transf){
	if(k == "vst"){
		exp_vst <- varianceStabilizingTransformation(
			exp_mat, 
			blind = TRUE, 
			fitType = "parametric"
		)
		assign(paste0("exp_none_none_", k), exp_vst)
	}else if(k == "rlog"){
		exp_rlog <- rlogTransformation(
			exp_mat,
			blind = TRUE,
			intercept,
			betaPriorVar,
			fitType = "parametric"
		)
		assign(paste0("exp_none_none_", k), exp_rlog)
	}
	write.table(get(paste0("exp_none_none_", k)),
		paste0("data/comp_norm_trans/exp_none_none_", k, ".txt")
	)
}

########### Transformation ##############
# hyperbolic arcsine (asinh)
exp_mat <- subset(exp_rsem_tpm, select = -c(gene_id, gene_name))
exp_tpm_asinh <- log(exp_mat + sqrt(exp_mat^2+1))

# VST
exp_vst <- varianceStabilizingTransformation(
	object, 
	blind = TRUE, 
	fitType = "parametric"
)
# rlog
exp_rlog <- rlogTransformation(
	object,
	blind = TRUE,
	intercept,
	betaPriorVar,
	fitType = "parametric"
)


########### ssGSEA ##############
library(matrixStats)
library(circlize)
library(ComplexHeatmap)
library(data.table)
library(clusterProfiler)
library(ssGSEA2)
library(GSVA)
exp_rsem_tpm <- read.table("data/exp_rsem_tpm.txt",
	header = TRUE,
	stringsAsFactors = FALSE)



log_tpm_rsem <- as.matrix(log2(exp_rsem_tpm2 + 0.01))

gmt_file <- "data/ReactomePathways.gmt" #2830 pathways
gene_sets <- read.gmt(gmt_file)
sel_paths <- grep("B cell|T cell|JAK|FLT3|NOTCH|MAPK|PI3K|apoptosis|cell cycle", 
     unique(gene_sets$term), 
     value = TRUE)
repeted_paths <- c("Signaling by NOTCH1 PEST Domain Mutants in Cancer",
"Signaling by NOTCH1 HD+PEST Domain Mutants in Cancer", 
"Constitutive Signaling by NOTCH1 HD+PEST Domain Mutants",
"Constitutive Signaling by NOTCH1 PEST Domain Mutants",
"crenolanib-resistant FLT3 mutants", "gilteritinib-resistant FLT3 mutants",            
"lestaurtinib-resistant FLT3 mutants", "linifanib-resistant FLT3 mutants",               
"midostaurin-resistant FLT3 mutants", "pexidartinib-resistant FLT3 mutants",            
"ponatinib-resistant FLT3 mutants", "quizartinib-resistant FLT3 mutants",             
"semaxanib-resistant FLT3 mutants", "sorafenib-resistant FLT3 mutants",               
"sunitinib-resistant FLT3 mutants", "tamatinib-resistant FLT3 mutants", 
"tandutinib-resistant FLT3 mutants", "Drug resistance of FLT3 mutants",
"KW2449-resistant FLT3 mutants")
sel_paths <- sel_paths[!is.element(sel_paths, repeted_paths)]
list_genesset <- lapply(sel_paths, function(x) {return(gene_sets$gene[gene_sets$term == x])})
names(list_genesset) <- sel_paths

########## Calculate the ssGSEA score #############
ssgsea_param <- ssgseaParam(log_tpm_rsem, list_genesset, verbose=FALSE, alpha = 0.25)
ES_ssgsea <- GSVA::gsva(ssgsea_param)

ES_ssgsea <- GSVA::gsva(log_tpm_rsem, list_genesset, verbose=FALSE, alpha = 0.25)

ES_ssgsea_z = (ES_ssgsea - rowMeans(ES_ssgsea))/(rowSds(as.matrix(ES_ssgsea)))[row(ES_ssgsea)]
pdf("graphs/ssgsea.pdf", width = 20, height = 10)
Heatmap(ES_ssgsea_z, col = colorRamp2(c(-2,0,2), c("orangered", "white", "purple")),
	 row_names_gp = gpar(fontsize = 6),
	 column_names_gp = gpar(fontsize = 10))
dev.off()
########## Calculate the GSVA score #############
gsva_param <- gsvaParam(log_tpm_rsem, list_genesset, verbose=FALSE, maxDiff = TRUE, kcdf = "Gaussian")
ES_GSVA_max <- GSVA::gsva(gsva_param)
ES_GSVA_z <- (ES_GSVA_max - rowMeans(ES_GSVA_max))/(rowSds(as.matrix(ES_GSVA_max)))[row(ES_GSVA_max)]
pdf("graphs/gsva.pdf", width = 20, height = 10)
Heatmap(ES_GSVA_z, col = colorRamp2(c(-2,0,2), c("orangered", "white", "purple")),
	 row_names_gp = gpar(fontsize = 6),
	 column_names_gp = gpar(fontsize = 10))
dev.off()

# Visualize distributions in R (by column or by row??) violon plot
# Use matrix norm or distance of singular values to measure the distance between matrices

