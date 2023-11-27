getSeuratObj <- function(raw.data = NULL, project_name = "lab1474", min.cells = 5, max.mito = 0.1, min.nGene = 500, max.nGene = 2500, min.nUMI = 1000, max.nUMI = 15000, treatment = "PBS", rm.genes = NULL, rm.cells = NULL){
  if(!is.null(rm.genes)){
    percent.rm.genes = colSums(raw.data[rownames(raw.data) %in% rm.genes,])/colSums(raw.data) 
    raw.data = raw.data[!(rownames(raw.data) %in% rm.genes),]
  }
  if(!is.null(rm.cells)){
	  for(gene in rm.cells) raw.data = raw.data[,!(raw.data[rownames(raw.data) %in% gene,] > 0)]
  }
  nUMI = colSums(raw.data)
  names(nUMI) = colnames(raw.data)
  percent.gfp = raw.data[rownames(raw.data) %in% "GFP",]/nUMI 
  percent.mito = colSums(raw.data[grep("^mt-",rownames(raw.data)),])/nUMI 
  my.obj <- CreateSeuratObject(raw.data = raw.data, project = "lab1474", min.cells = min.cells)
	my.obj <- AddMetaData(object = my.obj, metadata = percent.mito, col.name = "percent.mito")
	my.obj <- AddMetaData(object = my.obj, metadata = percent.gfp, col.name = "percent.gfp")
	my.obj <- AddMetaData(object = my.obj, metadata = percent.rm.genes, col.name = "percent.rm.genes")
	my.obj@meta.data$stim <- treatment
	my.obj <- FilterCells(object = my.obj, subset.names = c("nGene","nUMI","percent.mito"), low.thresholds = c(min.nGene,min.nUMI,-Inf), high.thresholds = c(max.nGene,max.nUMI,max.mito))
	my.obj <- NormalizeData(my.obj)
	my.obj <- ScaleData(my.obj)
	return(my.obj)
}
