suppressMessages({
  library(Seurat)
  library(dplyr)
})

proc_so <- function(so,res) {
    so <- NormalizeData(so, normalization.method = "LogNormalize", scale.factor = 10000)
    so <- FindVariableFeatures(so, selection.method = "vst", nfeatures = 2000)

    so <- ScaleData(so, features = rownames(so))
    so <- RunPCA(so, features = VariableFeatures(object = so))
    so <- FindNeighbors(so, dims = 1:10)
    so <- FindClusters(so, resolution = res)
    so <- RunUMAP(so, dims = 1:10)
    return(so)
}

merge_objects <- function(root,inp) {
  oup <- paste0(root,"/r_objects/raw.rds")
  if (file.exists(oup)) {
    message("loading raw")
    merged <- readRDS(oup)
  } else { 
 

    data <- lapply(seq_along(inp[["names"]]), function (i) { 
      proj <- inp[i,"names"]
      pth <- inp[i,"files"]
      print(i)

      samplefile <- paste0(root,"/r_objects/samples/",proj,".rds")

      if (file.exists(samplefile)) {
        return(readRDS(samplefile))
      }

      print(proj)
      print(pth)

      if (endsWith(pth,"h5")) {
        mat <- Read10X_h5(pth)
      } else {
        print(list.files(pth))
        if ("filtered_feature_bc_matrix" %in% list.files(pth)) {
          matpth <- paste0(pth,"/filtered_feature_bc_matrix")
        } else {
          matpth <- pth
        }
        mat <- Read10X(data.dir=matpth)
      }

      if ("Gene Expression" %in% names(mat)) {
        ## multi-modal data
        so <- CreateSeuratObject(
          counts = mat[["Gene Expression"]],#mat,
          project = proj
        )
      } else {
        so <- CreateSeuratObject(counts=mat, project=proj)
      }

      so[["percent.mt"]] <- PercentageFeatureSet(so, pattern = "^mt-")
      if (length(unique(so$percent.mt)) == 1) { # human
        so[["percent.mt"]] <- PercentageFeatureSet(so, pattern = "^MT-")
      }

      print(so)

      print(so)
      saveRDS(so,samplefile)
      return(so)
    })

    if (length(inp[["names"]]) == 1) {
      merged <- data[[1]]
    } else {
      merged <- merge(data[[1]], y=data[2:length(data)], add.cell.ids = inp[,"names"])
    }
    saveRDS(merged,oup)
  }
  return(merged)
}

recluster_obj <- function(root,merged,valid_bc,resolution) {

  if (file.exists(paste0(root,"/r_objects/recluster.rds"))) {
    recluster1 <- readRDS(paste0(root,"/r_objects/recluster.rds"))
  } else {
    print(length(valid_bc))
    print(length(which(rownames(merged@meta.data) %in% valid_bc)))
    merged@meta.data$bc <- rownames(merged@meta.data)
    recluster1 <- subset(merged, subset = bc %in% valid_bc) 
    recluster1 <- proc_so(recluster1,resolution)
    saveRDS(recluster1,paste0(root,"/r_objects/recluster.rds"))
  }

}

args = commandArgs(trailingOnly=TRUE)

source(paste0("scripts/",args[1],".R"))

if (!dir.exists(root)) { dir.create(root) }
if (!dir.exists(paste0(root,"/r_objects"))) { dir.create(paste0(root,"/r_objects")) }
if (!dir.exists(paste0(root,"/r_objects/samples"))) { dir.create(paste0(root,"/r_objects/samples")) }

merged <- merge_objects(root,input_files)

recluster_obj(root,merged,valid_bc,resolution)









