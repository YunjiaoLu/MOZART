###############
# Date of creation 29 Jan 2026
# Yunjiao Lu
###############
library(MOFA2)
library(MOFAdata)
library(data.table)
library(ggplot2)
library(tidyverse)
dir_path <- "/Users/ylu/Documents/MOZART"
########################### Read build MOFA model ###############################
utils::data("CLL_data")       
lapply(CLL_data,dim)
CLL_metadata <- fread("ftp://ftp.ebi.ac.uk/pub/databases/mofa/cll_vignette/sample_metadata.txt")
# Create the MOFA object and train the model
MOFAobject <- create_mofa(CLL_data)

pdf(file = paste0(dir_path, "/graphs/data_overview.pdf"))
plot_data_overview(MOFAobject)
dev.off()

# Define MOFA data options
data_opts <- get_default_data_options(MOFAobject)
data_opts
# Define MOFA model options
model_opts <- get_default_model_options(MOFAobject)
model_opts$num_factors <- 15
model_opts$likelihoods["Mutations"] = "bernoulli"
# Training options
train_opts <- get_default_training_options(MOFAobject)

# Train the MOFA model

MOFAobject <- prepare_mofa(MOFAobject,
  data_options = data_opts,
  model_options = model_opts,
  training_options = train_opts
)

start_time <- Sys.time()
MOFAobject <- run_mofa(MOFAobject, 
	use_basilisk = TRUE,
	outfile=paste0(dir_path, "/MOFA2_CLL.hdf5"))
end_time <- Sys.time()

end_time - start_time # 15 seconds
################### Downstream analysis ################
# Sanity check
# stopifnot(all(sort(CLL_metadata$sample)==sort(unlist(samples_names(MOFAobject)))))

# Add sample metadata to the model
samples_metadata(MOFAobject) <- CLL_metadata
# Dimensionality of the factor matrix: 200 samples, 15 factors
dim(MOFAobject@expectations$Z$group1)
# Dimensionality of the mRNA Weight matrix: 5000 features, 15 factors
dim(MOFAobject@expectations$W$mRNA)

# Correlation between factors
pdf(file = paste0(dir_path, "/graphs/cor_factors.pdf"))
plot_factor_cor(MOFAobject)
dev.off()

# Plot variance decomposition
pdf(file = paste0(dir_path, "/graphs/var_explained.pdf"))
plot_variance_explained(MOFAobject, max_r2=15)
dev.off()

pdf(file = paste0(dir_path, "/graphs/var_total_explained.pdf"))
plot_variance_explained(MOFAobject, plot_total = T)[[2]]
dev.off()

########## Characterization of factor 1
# association with sample metadata
pdf(file = paste0(dir_path, "/graphs/cov_metadata.pdf"))
correlate_factors_with_covariates(MOFAobject, 
  covariates = c("Gender","died","age"), 
  plot="log_pval"
)
dev.off()