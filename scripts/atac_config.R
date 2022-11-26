
root <- paste0("./atac_outputs/")

genome_id <- "mm10"


input_files <- data.frame(
  files=c(
"../../data/KO_1_atac_fragments.tsv.gz",
"../../data/KO_2_atac_fragments.tsv.gz",
"../../data/KO_3_atac_fragments.tsv.gz",
"../../data/KO_4_atac_fragments.tsv.gz",
"../../data/WT_1_atac_fragments.tsv.gz",
"../../data/WT_2_atac_fragments.tsv.gz",
"../../data/WT_3_atac_fragments.tsv.gz",
"../../data/WT_4_atac_fragments.tsv.gz"
  ),  
  names=c(
"M_KO_1",
"F_KO_2",
"M_KO_3",
"F_KO_4",
"M_WT_1",
"F_WT_2",
"M_WT_3",
"F_WT_4"
  ),
  stringsAsFactors=FALSE
)


main_dr <- "All_15_LSI" 
main_clusters <- "All_15_200_LSI"

sub_prefix <- "Subset"
sub_dr <- "Subset_8_LSI" 
sub_clusters <- "Subset_8_400_LSI" 


valid_bc <- read.delim(paste0("barcode_info/archr_cellColData_trim.csv"),sep=",")[[1]]
sub_sel <- function(x) {
  return( rownames(x) %in% valid_bc )
}




