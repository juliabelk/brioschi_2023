
root <- paste0("./gex_outputs/")


input_files <- data.frame(
  files=c(
"data/KO_1_RNA_filtered_feature_bc_matrix.h5",
"data/KO_2_RNA_filtered_feature_bc_matrix.h5",
"data/KO_3_RNA_filtered_feature_bc_matrix.h5",
"data/KO_4_RNA_filtered_feature_bc_matrix.h5",
"data/WT_1_RNA_filtered_feature_bc_matrix.h5",
"data/WT_2_RNA_filtered_feature_bc_matrix.h5",
"data/WT_3_RNA_filtered_feature_bc_matrix.h5",
"data/WT_4_RNA_filtered_feature_bc_matrix.h5"
  ),
  names=c(
    "KO_1_male-KO_1",
    "KO_2_female-KO_2",
    "KO_3_male-KO_3",
    "KO_4_female-KO_4",
    "WT_1_male-WT_1",
    "WT_2_female-WT_2",
    "WT_3_male-WT_3",
    "WT_4_female-WT_4"
  ),
  stringsAsFactors=FALSE
)

#clusters <- "RNA_snn_res.0.2"

valid_bc <- read.delim(paste0("barcode_info/seurat_metadata_trim.csv"),sep=",")[[1]]
resolution <- 0.2



