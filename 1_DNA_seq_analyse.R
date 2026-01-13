###################################
# Author: Yunjiao Lu
# Date: 28 Oct 2025
###################################
library(readxl)
library(data.table)
library(mixOmics)
# Create patient mutation dataset
df_dnaseq <- data.frame()
filenames <- list.files(path = "data/CPX-IPC-analyseOK 3/CPX-IPC-analyseOK")

for (file in filenames) {
  data_dnaseq <- as.data.table(read_xlsx(paste0("data/CPX-IPC-analyseOK 3/CPX-IPC-analyseOK/", file)))
  patient_id <- unique(data_dnaseq$PATIENT_ID)
  data_dnaseq_sub <- data_dnaseq[Classification_mutation %in% c("A", "B", "C")]
  max_vaf <- data_dnaseq_sub[, .(max_VAF = max(Mean_VAF)), by = Gene.refGene]
  nb_muta <- table(data_dnaseq_sub$Classification_mutation)
  
  data_patient_i <- c(nb_muta, max_vaf$max_VAF )
  names(data_patient_i) <- c(paste0("nb_muta", names(nb_muta)), paste0("VAF_", (max_vaf$Gene.refGene)))
  df_patient_i <- data.frame(patient_id, as.data.frame(t(as.matrix(data_patient_i))))
  
  colnames_df <- colnames(df_dnaseq)
  new_vars <- colnames(df_patient_i)[!is.element(colnames(df_patient_i), colnames_df)]
  df_dnaseq_comp <- as.data.frame(matrix(0, nrow = dim(df_dnaseq)[1], ncol = length(new_vars)))
  colnames(df_dnaseq_comp) <- new_vars
  df_dnaseq <- data.frame(df_dnaseq, df_dnaseq_comp)
  
  old_vars <- colnames_df[!is.element(colnames_df, colnames(df_patient_i))]
  df_patient_comp <- as.data.frame(matrix(0, nrow = 1, ncol = length(old_vars)))
  colnames(df_patient_comp) <- old_vars
  df_patient_i <- data.frame(df_patient_i, df_patient_comp)
  
  df_dnaseq <- rbind(df_dnaseq, df_patient_i)
}

# modeling mutation data with sPLS
