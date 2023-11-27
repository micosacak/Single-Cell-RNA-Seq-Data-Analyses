# updated on 24.01.2019
{
t1 = Sys.time() # to calculate average running time!
getColors = function(){
   all.colors.to.use=c("yellow3","red","skyblue","cyan3","violet","purple",
   "gray80","orange","darkgreen","pink","gold","cadetblue1","green","darkkhaki",
   "chocolate","darkred","aquamarine","darkslategray4","blue","cyan1")
   names(all.colors.to.use) = c("Im","NN_0","NN_1","NN_2","NN_3","NN_4","NN_5",
   "NN_6","NN_7","OPCOD_1","OPCOD_2","PC_0","PC_1","PC_2","PC_3","PC_4","PC_5","PC_6","PC_7","PC_8")
   
   five.clusters.colors = c("yellow3","green","pink","gold","blue")
   names(five.clusters.colors) = c("Im","NN","OPCOD_1","OPCOD_2","PC")
   return(list("all.colors.to.use" = all.colors.to.use, "five.clusters.colors" = five.clusters.colors))
}
# get inputs
{
   # these colors has been determined after one round of analysis
   my_colors = getColors()
   all.colors.to.use = my_colors$all.colors.to.use
   five.clusters.colors = my_colors$five.clusters.colors
   cells.to.rem = c("IL4_Cells_121","IL4_Cells_175","IL4_Cells_493","IL4_Cells_510","AB42_Cells_534")
   
   # load required libraries
   library("Matrix")
   library("Seurat")
   library("dplyr")
   library("openxlsx")
   library("gplots")
   library("igraph")
   
   # some custom functions
   scs_fnxs = "scs_functions/"
   fnxs = dir(scs_fnxs)
   fnxs = fnxs[grep(".R$",fnxs)]
   for(fnx in fnxs) source(paste0(scs_fnxs,fnx))
   
   # save raw file as rda
   dir.create("rdaFiles", showWarnings = F)
   if(!("raw.data.rda" %in% dir("rdaFiles/"))){
      pbs.data <- Read10X(data.dir =  "10X_Data_AutomatedCellsCounts/PBS_Cells/")
      ab42.data <- Read10X(data.dir = "10X_Data_AutomatedCellsCounts/AB42_Cells/")
      il4.data <- Read10X(data.dir =  "10X_Data_AutomatedCellsCounts/IL4_Cells/")
      colnames(pbs.data) = paste("PBS_Cells_",seq(1,ncol(pbs.data),1), sep = "")
      colnames(ab42.data) = paste("AB42_Cells_",seq(1,ncol(ab42.data),1), sep = "")
      colnames(il4.data) = paste("IL4_Cells_",seq(1,ncol(il4.data),1), sep = "")
      raw.data = list("PBS" = pbs.data, "AB42" = ab42.data, "IL4" = il4.data)
      
      #for(trt in names(raw.data)) raw.data[[trt]] = raw.data[[trt]][,!(colnames(raw.data[[trt]]) %in% cells.to.rem)]
      
      save(raw.data, file = "rdaFiles/raw.data.rda")
   }else{
      load("rdaFiles/raw.data.rda")
   }
   
  
   
   # generate new folder to save results for each run.
   mainFolder = paste0("ranOn_",gsub("[: -]","",as.character(Sys.time())))
   dir.create(mainFolder, showWarnings = FALSE)
   dir.create(paste0(mainFolder,"/rdaFiles"),showWarnings = FALSE)
   dir.create(paste0(mainFolder,"/plots"),showWarnings = FALSE)
   setwd(mainFolder)
   
   # ribosomal protein genes to be removed!
   {
      rp.genes = c("rpl18a","rplp2l","rpl34","rpl13","rplp0","rpl36a","rpl12","rpl7a","rpl19","rpl3","rpl27","rpl23","rpl5b",
                   "rplp2","rpl5a","rpl7","rpl37","rpl24","rpl9","rpl8","rpl31","rpl18","rpl28","rpl7l1","rpl6","rpl10a",
                   "rpl13a","rpl39","rpl26","rpl4","rpl35a","rpl38","rplp1","rpl30","rpl11","rpl14","rpl10","rpl37.1","rpl35",
                   "rpl17","rpl23a","rpl29","rpl22","rpl21","rpl22l1","rpl36","rpl32","rps16","rps13","rps4x","rps17","rps6ka3a",
                   "rps6ka4","rps2","rps15a","rps11","rps19bp1","rps27a","rpsa","rps26l","rps10","rps28","rps8a","rps3a","rps6",
                   "rps27.2","rps19","rps9","rps6kc1","rps7","rps29","rps8b","rps6ka1","rps6ka5","rps6kl1","rps6kal","rps6ka2","rps24",
                   "rps3","rps27.1","rps18","rps6kb1b","rps5","rps21","rps26","rps12","rps14","rps6kb1a","rps25","rps15","rps23",
                   "rps6ka3b","RPS17","RPL41")
   }
}

# first step: this step is used to determine major cell types: Radial Glia (RG), Neuron, Oligodendrocytes (OD) and Immune Cells (Im)
{
   # some preliminary analysis to determine thresholds
   for(treatment in names(raw.data)) getHistos(raw.data[[treatment]], folderName = treatment, moreGenes = c("GFP"))
   
   # filter out cells;
   # nUMI < 1000 and nUMI>10000
   # percent.mito > 0.06 (6%)
   # nGene < 500 and nGene > 2500
   # also if a gene found in less than 5 cells were removed
   
   num.ccs = 30     # number of common correlations!
   num.dims = 30
   max.genes = 1000 # number of top dipsersed genes to used!
   min.cells = 5    # min number of cells per gene!
   
   # generate and save objects!
   # generate objects
   my.main.objs = list("PBS" = "", "AB42" = "", "IL4" = "")
   for(trt in names(my.main.objs)) my.main.objs[[trt]] =  getSeuratObj(raw.data = as.matrix(raw.data[[trt]]), project_name = "lab1474", treatment = trt, rm.genes = rp.genes, min.cells = min.cells, max.mito = 0.06)
   
   # # cells below might be mixed of two cells types! Remove them from down stream analysis!
   for(trt in names(my.main.objs)) my.main.objs[[trt]] = SubsetData(my.main.objs[[trt]], cells.use = rownames(my.main.objs[[trt]]@meta.data)[!(rownames(my.main.objs[[trt]]@meta.data) %in% cells.to.rem)])
   # 
   # generate gene plots, these gene plots can be generated before filtering as well!
   for(trt in names(my.main.objs)){
      pdf(paste0("plots/",trt,"_GenePlots.pdf"))
      par(mfrow = c(1, 3))
      GenePlot(object = my.main.objs[[trt]], gene1 = "nUMI", gene2 = "percent.mito", col.use = "blue", pch.use = 20, cex.use = 0.3)
      GenePlot(object = my.main.objs[[trt]], gene1 = "nUMI", gene2 = "nGene", col.use = "blue", pch.use = 20, cex.use = 0.3)
      GenePlot(object = my.main.objs[[trt]], gene1 = "nUMI", gene2 = "percent.gfp", col.use = "blue", pch.use = 20, cex.use = 0.3)
      dev.off()
   }
   
   # violin plots for nGene, mitochondrial genes and GFP percentages
   for(trt in names(my.main.objs)){
      pdf(paste0("plots/",trt,"_VlnPlot_nGene_nUMI_mitPercent.pdf"), width = 14)
      print(VlnPlot(object = my.main.objs[[trt]], features.plot = c("nGene", "nUMI", "percent.mito","percent.gfp"), nCol = 4, point.size.use = 0.75, cols.use = "gray"))
      dev.off()
   }
   
   # select genes
   for(trt in names(my.main.objs)){
      pdf(paste0("plots/",trt,".data.varGenes.pdf"),width = 35, height = 35)
      my.main.objs[[trt]] <- FindVariableGenes(object = my.main.objs[[trt]], mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.125, x.high.cutoff = 10, y.cutoff = 0.5)
      dev.off()
   }
   
   # genes to be used!
   genes.use = c()
   for(trt in names(my.main.objs)) genes.use = union(genes.use, rownames(my.main.objs[[trt]]@hvg.info)[1:max.genes])
   length(genes.use)
   
   for(trt in names(my.main.objs)) genes.use = intersect(genes.use, rownames(my.main.objs[[trt]]@scale.data))
   length(genes.use)
   # remove mitochondrial genes!
   genes.use = genes.use[-(grep("^mt-",genes.use))]
   length(genes.use) # 1362 genes.
   
   # save objects for later use!
   rawData.main.objs = list("pbs.obj" = my.main.objs[["PBS"]],"il4.obj" = my.main.objs[["IL4"]], "ab42.obj" = my.main.objs[["AB42"]], "genes.use" = genes.use)
   save(rawData.main.objs, file = paste0("rdaFiles/rawData.main.objs.rda"))
   save(my.main.objs, file ="rdaFiles/my.main.objs.rda")
   # Run rare non-overlapping filtering
   all.combined.main = RunMultiCCA(rawData.main.objs[1:(length(rawData.main.objs)-1)], genes.use = rawData.main.objs$genes.use, num.ccs = num.ccs)
   
   # CC Selection
   metageneBicorP.main <- MetageneBicorPlot(all.combined.main, grouping.var = "stim", dims.eval = 1:num.dims, display.progress = T)
   
   # save the data!
   pdf(paste0("plots/metageneBicorPlot.main.pdf"))
   print(metageneBicorP.main)
   dev.off()
   
   save(metageneBicorP.main, file = paste0("rdaFiles/metageneBicorP.main.rda"))
   
   num.dims.use = 10
   all.combined.main <- AlignSubspace(all.combined.main, reduction.type = "cca", grouping.var = "stim", dims.align = 1:num.dims.use)
   
   # run TSNE
   all.combined.main <- RunTSNE(all.combined.main, reduction.use = "cca.aligned", dims.use = 1:num.dims.use)
   
   resoL = 1
   all.combined.main <- FindClusters(all.combined.main, reduction.type = "cca.aligned",dims.use = 1:num.dims.use, save.SNN = T, resolution = resoL, force.recalc = T)
   
   # Visualization
   TSNEPlot(all.combined.main, do.label = T, do.return = T)
   
   p1 = TSNEPlot(all.combined.main, do.return = T, group.by = "stim")
   p2 = TSNEPlot(all.combined.main, do.label = T, do.return = T) 
   
   
   pdf(paste0("plots/all.main.tSNE_samples_",resoL,".pdf"))
   print(p1)
   dev.off()
   
   pdf(paste0("plots/all.main.tSNE_clusters_",resoL,".pdf"))
   print(p2)
   dev.off()
   
   pdf(paste0("plots/all.main.tSNE_clusters_and_samples_",resoL,".pdf"), width = 14)
   print(plot_grid(p1,p2, nrow = 1))
   dev.off()
   
   ################ add new indent for 5 major cell clusters!
   new.idents = as.character(all.combined.main@ident)
   names(new.idents) = names(all.combined.main@ident)
   #TSNEPlot(all.combined.main, do.label = T)
   new.idents[new.idents %in% 10] = "Im"
   new.idents[new.idents %in% 11] = "OPCOD_1"
   new.idents[new.idents %in% 12] = "OPCOD_2"
   
   new.idents[new.idents %in% c(1,3,4,8,9)] = "NN"
   new.idents[new.idents %in% c(0,2,5,6,7)] = "PC"
   table(new.idents)
   all.combined.main = StashIdent(all.combined.main, save.name = "res.x.idents")
   all.combined.main = SetIdent(all.combined.main, ident.use = new.idents)
   
   p3 = TSNEPlot(all.combined.main, do.label = T, do.return = T)
   pdf(paste0("plots/all.main.tSNE_newColors_clusters_and_samples_",resoL,".pdf"), width = 21)
   print(plot_grid(p1,p2,p3, nrow = 1))
   dev.off()
   
   pdf(paste0("plots/all.main.tSNE_newColors_clusters_",resoL,".pdf"), width = 7)
   print(plot_grid(p3))
   dev.off()
   
   
   # stash idents
   all.combined.main = StashIdent(all.combined.main, save.name = "main.cell.clusters")
   save(all.combined.main, file = "rdaFiles/all.combined.main.rda")
   
   main.cell.markers.list <- FindAllMarkers(object = all.combined.main, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
   main.cell.markers.list = subset(main.cell.markers.list, p_val_adj < 0.1)
   nrow(main.cell.markers.list)
   
   for(topN in c(5,10,20,30,40,50,100)){
      topNgenes <- main.cell.markers.list %>% group_by(cluster) %>% top_n(topN, avg_logFC)
      pdf(paste0("plots/MAIN_cell_clusters_heatmap_top",topN,".pdf"), width = nrow(all.combined.main@meta.data)*0.01, height = 1 +length(unique(topNgenes$gene))*0.14)
      print(DoHeatmap(object = all.combined.main, genes.use = unique(topNgenes$gene), slim.col.label = TRUE, remove.key = TRUE, group.label.rot = T))
      dev.off()
   }
   save(main.cell.markers.list, file = "rdaFiles/main.cell.markers.list.rda")    
}

# second step_1: this step is used to classify RG cells 
{
{
   # the neuron and glia cluster can not be separated clearly! As a result, re-run tSNE and clusterings
   load("../rdaFiles/raw.data.rda")
   for(trt in names(raw.data)) raw.data[[trt]] = raw.data[[trt]][, colnames(raw.data[[trt]]) %in% names(all.combined.main@ident)[all.combined.main@ident %in% "PC"]]
   
   # we filter out cells;
   # nUMI < 1000 and nUMI>10000
   # percent.mito > 0.06 (6%)
   # nGene < 500 and nGene > 2500
   # also if a gene found in less than 5 cells were removed
   
   num.ccs = 30     # number of common correlations!
   max.genes = 1000 # number of top dipsersed genes to used!
   min.cells = 5    # min number of cells per gene!
   
   # generate and save objects!
   # generate objects
   my.glia.objs = list("PBS" = "", "AB42" = "", "IL4" = "")
   for(trt in names(my.glia.objs)) my.glia.objs[[trt]] =  getSeuratObj(raw.data = as.matrix(raw.data[[trt]]), project_name = "lab1474", treatment = trt, rm.genes = rp.genes, min.cells = min.cells, max.mito = 0.06)
   
   # generate gene plots, these gene plots can be generated before filtering as well!
   for(trt in names(my.glia.objs)){
      pdf(paste0("plots/glia_",trt,"_GenePlots.pdf"))
      par(mfrow = c(1, 3))
      GenePlot(object = my.glia.objs[[trt]], gene1 = "nUMI", gene2 = "percent.mito", col.use = "blue", pch.use = 20, cex.use = 0.3)
      GenePlot(object = my.glia.objs[[trt]], gene1 = "nUMI", gene2 = "nGene", col.use = "blue", pch.use = 20, cex.use = 0.3)
      GenePlot(object = my.glia.objs[[trt]], gene1 = "nUMI", gene2 = "percent.gfp", col.use = "blue", pch.use = 20, cex.use = 0.3)
      dev.off()
   }
   
   # violin plots for nGene, mitochondrial genes and GFP percentages
   for(trt in names(my.glia.objs)){
      pdf(paste0("plots/glia_",trt,"_VlnPlot_nGene_nUMI_mitPercent.pdf"), width = 14)
      print(VlnPlot(object = my.glia.objs[[trt]], features.plot = c("nGene", "nUMI", "percent.mito","percent.gfp"), nCol = 4, point.size.use = 0.75, cols.use = "gray"))
      dev.off()
   }
   
   # select genes
   for(trt in names(my.glia.objs)){
      pdf(paste0("plots/glia_",trt,".data.varGenes.pdf"),width = 35, height = 35)
      my.glia.objs[[trt]] <- FindVariableGenes(object = my.glia.objs[[trt]], mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.125, x.high.cutoff = 10, y.cutoff = 0.5)
      dev.off()
   }
   
   # genes to be used!
   genes.use = c()
   for(trt in names(my.glia.objs)) genes.use = union(genes.use, rownames(my.glia.objs[[trt]]@hvg.info)[1:max.genes])
   length(genes.use)
   for(trt in names(my.glia.objs)) genes.use = intersect(genes.use, rownames(my.glia.objs[[trt]]@scale.data))
   length(genes.use)
   
   # remove mitochondrial genes!
   genes.use = genes.use[-(grep("^mt-",genes.use))]
   length(genes.use)
   
   # save objects for later use!
   rawData.glia.objs = list("pbs.obj" = my.glia.objs[["PBS"]],"il4.obj" = my.glia.objs[["IL4"]], "ab42.obj" = my.glia.objs[["AB42"]], "genes.use" = genes.use)
   save(rawData.glia.objs, file = paste0("rdaFiles/rawData.glia.objs.rda"))
   save(my.glia.objs, file = paste0("rdaFiles/my.glia.objs.rda"))
   
   # Run rare non-overlapping filtering
   all.combined.glia = RunMultiCCA(rawData.glia.objs[1:(length(rawData.glia.objs)-1)], genes.use = rawData.glia.objs$genes.use, num.ccs = num.ccs)
   
   # CC Selection
   metageneBicorP.glia <- MetageneBicorPlot(all.combined.glia, grouping.var = "stim", dims.eval = 1:num.dims, display.progress = T)
   
   # save the data!
   pdf(paste0("plots/glia_metageneBicorPlot.glia.pdf"))
   print(metageneBicorP.glia)
   dev.off()
   save(metageneBicorP.glia, file = paste0("rdaFiles/metageneBicorP.glia.rda"))
   num.dims.use = 20
   
   all.combined.glia <- AlignSubspace(all.combined.glia, reduction.type = "cca", grouping.var = "stim", dims.align = 1:num.dims.use)
   
   # run TSNE
   all.combined.glia <- RunTSNE(all.combined.glia, reduction.use = "cca.aligned", dims.use = 1:num.dims.use)
   
   resoL = 1
   all.combined.glia <- FindClusters(all.combined.glia, reduction.type = "cca.aligned",
                                     dims.use = 1:num.dims.use, save.SNN = T, resolution = resoL, force.recalc = T)
   # Visualization
   
    TSNEPlot(all.combined.glia, do.label = T, do.return = T)
 
   p4 = TSNEPlot(all.combined.glia, do.return = T, group.by = "stim")
   p5 = TSNEPlot(all.combined.glia, do.label = T, do.return = T) 
    
   pdf(paste0("plots/all.glia.tSNE_samples_",resoL,".pdf"))
   print(p4)
   dev.off()
   
   pdf(paste0("plots/all.glia.tSNE_clusters_",resoL,".pdf"))
   print(p5)
   dev.off()
   
   pdf(paste0("plots/all.glia.tSNE_clusters_and_samples_",resoL,".pdf"), width = 14)
   print(plot_grid(p4,p5,nrow = 1))
   dev.off()
   
   all.combined.glia = SetIdent(all.combined.glia, ident.use = paste("PC_", all.combined.glia@ident, sep = ""))
   
   p6 = TSNEPlot(all.combined.glia, do.label = T, do.return = T, colors.use = all.colors.to.use[grep("PC_",names(all.colors.to.use))]) 
   pdf(paste0("plots/all.glia.tSNE_newColors_clusters_and_samples_",resoL,".pdf"), width = 21)
   print(plot_grid(p4,p5,p6,nrow = 1))
   dev.off()
   
   pdf(paste0("plots/all.glia.tSNE_newColors_clusters_",resoL,".pdf"))
   print(plot_grid(p6))
   dev.off()
   
   save(all.combined.glia,file = "rdaFiles/all.combined.glia.rda")
   
   main.glia.cell.markers.list <- FindAllMarkers(object = all.combined.glia, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
   main.glia.cell.markers.list = subset(main.glia.cell.markers.list, p_val_adj < 0.1)
   nrow(main.glia.cell.markers.list)
   
   for(topN in c(5,10,20,30,40,50,100)){
      topNgenes.glia <- main.glia.cell.markers.list %>% group_by(cluster) %>% top_n(topN, avg_logFC)
      pdf(paste0("plots/main_glia_clusters_heatmap_top",topN,".pdf"), width = nrow(all.combined.glia@meta.data)*0.01, height = 1+length(unique(topNgenes.glia$gene))*0.14)
      print(DoHeatmap(object = all.combined.glia, genes.use = unique(topNgenes.glia$gene), slim.col.label = TRUE, remove.key = TRUE, group.label.rot = T))
      dev.off()
   }
   save(main.glia.cell.markers.list, file = "rdaFiles/main.glia.cell.markers.list.rda")
   
}

# second step_2: this step is used to classify Neurons
{
   load("../rdaFiles/raw.data.rda")
   
   for(trt in names(raw.data)) raw.data[[trt]] = raw.data[[trt]][, colnames(raw.data[[trt]]) %in% names(all.combined.main@ident)[all.combined.main@ident %in% "NN"]]
   # some preliminary analysis to determine thresholds
   
   # we filter out cells;
   # nUMI < 1000 and nUMI>10000
   # percent.mito > 0.06 (6%)
   # nGene < 500 and nGene > 2500
   # also if a gene found in less than 5 cells were removed
   
   num.ccs = 30     # number of common correlations!
   max.genes = 1000 # number of top dipsersed genes to used!
   min.cells = 5    # min number of cells per gene!
   
   # generate and save objects!
   # generate objects
   my.neuron.objs = list("PBS" = "", "AB42" = "", "IL4" = "")
   for(trt in names(my.neuron.objs)) my.neuron.objs[[trt]] =  getSeuratObj(raw.data = as.matrix(raw.data[[trt]]), project_name = "lab1474", treatment = trt, rm.genes = rp.genes, min.cells = min.cells, max.mito = 0.06)
   
   # generate gene plots, these gene plots can be generated before filtering as well!
   for(trt in names(my.neuron.objs)){
      pdf(paste0("plots/neuron_",trt,"_GenePlots.pdf"))
      par(mfrow = c(1, 3))
      GenePlot(object = my.neuron.objs[[trt]], gene1 = "nUMI", gene2 = "percent.mito", col.use = "blue", pch.use = 20, cex.use = 0.3)
      GenePlot(object = my.neuron.objs[[trt]], gene1 = "nUMI", gene2 = "nGene", col.use = "blue", pch.use = 20, cex.use = 0.3)
      GenePlot(object = my.neuron.objs[[trt]], gene1 = "nUMI", gene2 = "percent.gfp", col.use = "blue", pch.use = 20, cex.use = 0.3)
      dev.off()
   }
   
   # violin plots for nGene, mitochondrial genes and GFP percentages
   for(trt in names(my.neuron.objs)){
      pdf(paste0("plots/neuron_",trt,"_VlnPlot_nGene_nUMI_mitPercent.pdf"), width = 14)
      print(VlnPlot(object = my.neuron.objs[[trt]], features.plot = c("nGene", "nUMI", "percent.mito","percent.gfp"), nCol = 4, point.size.use = 0.75, cols.use = "gray"))
      dev.off()
   }
   
   # select genes
   for(trt in names(my.neuron.objs)){
      pdf(paste0("plots/neuron_",trt,".data.varGenes.pdf"),width = 35, height = 35)
      my.neuron.objs[[trt]] <- FindVariableGenes(object = my.neuron.objs[[trt]], mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.125, x.high.cutoff = 10, y.cutoff = 0.5)
      dev.off()
   }
   
   # genes to be used!
   genes.use = c()
   for(trt in names(my.neuron.objs)) genes.use = union(genes.use, rownames(my.neuron.objs[[trt]]@hvg.info)[1:max.genes])
   length(genes.use)
   for(trt in names(my.neuron.objs)) genes.use = intersect(genes.use, rownames(my.neuron.objs[[trt]]@scale.data))
   length(genes.use)
   
   # remove mitochondrial genes!
   genes.use = genes.use[-(grep("^mt-",genes.use))]
   length(genes.use)
   
   # save objects for later use!
   rawData.neuron.objs = list("pbs.obj" = my.neuron.objs[["PBS"]],"il4.obj" = my.neuron.objs[["IL4"]], "ab42.obj" = my.neuron.objs[["AB42"]], "genes.use" = genes.use)
   save(rawData.neuron.objs, file = paste0("rdaFiles/rawData.neuron.objs.rda"))
   save(my.neuron.objs, file = paste0("rdaFiles/my.neuron.objs.rda"))
   
   # Run rare non-overlapping filtering
   all.combined.neuron = RunMultiCCA(rawData.neuron.objs[1:(length(rawData.neuron.objs)-1)], genes.use = rawData.neuron.objs$genes.use, num.ccs = num.ccs)
   
   # CC Selection
   metageneBicorP.neuron <- MetageneBicorPlot(all.combined.neuron, grouping.var = "stim", dims.eval = 1:num.dims, display.progress = T)
   
   # save the data!
   pdf(paste0("plots/neuron_metageneBicorPlot.neuron.pdf"))
   print(metageneBicorP.neuron)
   dev.off()
   save(metageneBicorP.neuron, file = paste0("rdaFiles/metageneBicorP.neuron.rda"))
   
   num.dims.use = 20
   
   all.combined.neuron <- AlignSubspace(all.combined.neuron, reduction.type = "cca", grouping.var = "stim", dims.align = 1:num.dims.use)
   
   # run TSNE
   all.combined.neuron <- RunTSNE(all.combined.neuron, reduction.use = "cca.aligned", dims.use = 1:num.dims.use)
   
   resoL = 1
   all.combined.neuron <- FindClusters(all.combined.neuron, reduction.type = "cca.aligned", dims.use = 1:num.dims.use, save.SNN = T, resolution = resoL, force.recalc = T)
   
   # Visualization
   
   TSNEPlot(all.combined.neuron, do.label = T, do.return = T)
   
   p7 = TSNEPlot(all.combined.neuron, do.return = T, group.by = "stim")
   p8 = TSNEPlot(all.combined.neuron, do.label = T, do.return = T) 
   
   
   pdf(paste0("plots/all.neuron.tSNE_samples_",resoL,".pdf"))
   print(p7)
   dev.off()
   
   pdf(paste0("plots/all.neuron.tSNE_clusters_",resoL,".pdf"))
   print(p8)
   dev.off()
   
   pdf(paste0("plots/all.neuron.tSNE_clusters_and_samples_",resoL,".pdf"), width = 14)
   print(plot_grid(p7,p8,nrow = 1))
   dev.off()
   
   all.combined.neuron = SetIdent(all.combined.neuron, ident.use = paste("NN_", all.combined.neuron@ident, sep = ""))
   
   p9 = TSNEPlot(all.combined.neuron, do.label = T, do.return = T, colors.use = all.colors.to.use[grep("NN_",names(all.colors.to.use))]) 
   
   pdf(paste0("plots/all.NewColors.neuron.tSNE_clusters_and_samples_",resoL,".pdf"), width = 21)
   print(plot_grid(p7,p8,p9,nrow = 1))
   dev.off()
   
   pdf(paste0("plots/all.NewColors.neuron.tSNE_clusters_",resoL,".pdf"))
   print(plot_grid(p9))
   dev.off()
   
   main.neuron.cell.markers.list <- FindAllMarkers(object = all.combined.neuron, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
   nrow(main.neuron.cell.markers.list)
   
   main.neuron.cell.markers.list = subset(main.neuron.cell.markers.list, p_val_adj < 0.1)
   nrow(main.neuron.cell.markers.list)
   for(topN in c(5,10,20,30,40,50,100)){
      topNgenes.neuron <- main.neuron.cell.markers.list %>% group_by(cluster) %>% top_n(topN, avg_logFC)
      pdf(paste0("plots/main.neuron.clusters_heatmap_top",topN,".pdf"), width = nrow(all.combined.neuron@meta.data)*0.01, height = 1+length(unique(topNgenes.neuron $gene))*0.14)
      print(DoHeatmap(object = all.combined.neuron, genes.use = unique(topNgenes.neuron$gene), slim.col.label = TRUE, remove.key = TRUE, group.label.rot = T))
      dev.off()
   }
   
   save(all.combined.neuron,file = "rdaFiles/all.combined.neuron.rda")
}

# third step: Combine all cell types in one Seurat Object
{
   # now combine all cell types!
   final.idents = as.character(all.combined.main@meta.data$res.x.idents)
   names(final.idents) = rownames(all.combined.main@meta.data)
   
   load("rdaFiles/all.combined.glia.rda")
   glia.idents = as.character(all.combined.glia@ident)
   names(glia.idents) = names(all.combined.glia@ident)
   
   load("rdaFiles/all.combined.neuron.rda")
   neuron.idents = as.character(all.combined.neuron@ident)
   names(neuron.idents) = names(all.combined.neuron@ident)
   
   other.idents = as.character(all.combined.main@meta.data$main.cell.clusters)
   names(other.idents) = names(all.combined.main@ident)
   other.idents = other.idents[other.idents %in% c("Im","OPCOD_1","OPCOD_2")]
   
   comb_all_idents = c(other.idents, glia.idents, neuron.idents)
   comb_all_idents = comb_all_idents[match(names(final.idents),names(comb_all_idents))]
   
   all.combined.main = SetIdent(all.combined.main, ident.use = comb_all_idents)
   
   all.combined.main = StashIdent(all.combined.main, save.name = "final.cells")
   
   TSNEPlot(all.combined.main, do.label = T, do.return = T)

   p10 = TSNEPlot(all.combined.main, do.return = T, group.by = "stim")
   p11 = TSNEPlot(all.combined.main, do.label = T, do.return = T, colors.use = all.colors.to.use) 
   
   
   pdf(paste0("plots/all.main.final.cells.tSNE_samples_",resoL,".pdf"))
   print(p10)
   dev.off()
   
   pdf(paste0("plots/all.main.final.cells.tSNE_clusters_",resoL,".pdf"))
   print(p11)
   dev.off()
   
   pdf(paste0("plots/all.main.final.cells.tSNE_clusters_and_samples_",resoL,".pdf"), width = 14)
   print(plot_grid(p10,p11,nrow = 1))
   dev.off()
   
   

   save(all.combined.main, file = "rdaFiles/all.combined.main.rda")
}
}

# add mor metda data
{
   my.final.colors = as.character(all.combined.main@meta.data$final.cells)
   table(my.final.colors)
   for(i in 1:length(all.colors.to.use)){
      my.final.colors[grep(names(all.colors.to.use)[i], my.final.colors)] = all.colors.to.use[i]
   }
   names(my.final.colors) = rownames(all.combined.main@meta.data)
   my.final.colors = as.data.frame(my.final.colors)
   all.combined.main = AddMetaData(all.combined.main, metadata = my.final.colors)
   colnames(all.combined.main@meta.data)
   
   my.main.colors = as.character(all.combined.main@meta.data$main.cell.clusters)
   table(my.main.colors)
   for(i in 1:length(five.clusters.colors)){
      my.main.colors[my.main.colors %in% names(five.clusters.colors)[i]] = five.clusters.colors[i]
   }
   names(my.main.colors) = rownames(all.combined.main@meta.data)
   my.main.colors = as.data.frame(my.main.colors)
   all.combined.main = AddMetaData(all.combined.main, metadata = my.main.colors)
   colnames(all.combined.main@meta.data)
   
   sz.fc = as.data.frame(rep(1.0,length(all.combined.main@ident)), row.names = names(all.combined.main@ident))
   colnames(sz.fc) = "Size_Factor"
   head(sz.fc)
   all.combined.main = AddMetaData(all.combined.main, metadata = sz.fc, col.name = "Size_Factor")
   
   save(all.combined.main, file = "rdaFiles/all.combined.main.rda")
}

# run monocle analysis
{
   dir.create("monocle_analysis", showWarnings = F)
   dir.create("monocle_analysis/rdaFiles", showWarnings = F)
   dir.create("monocle_analysis/plots", showWarnings = F)
   
   library("Seurat")
   library("monocle")
   
   ph.dt = all.combined.main@meta.data
   ft.dt = as.data.frame(rownames(all.combined.main@scale.data), row.names = rownames(all.combined.main@scale.data))
   head(ft.dt)
   colnames(ft.dt) = "gene_short_name"
   pd = new(ph.dt, Class ="AnnotatedDataFrame")
   ft = new(ft.dt, Class ="AnnotatedDataFrame")
   HSMM = newCellDataSet(cellData = all.combined.main@scale.data, phenoData = pd, featureData = ft)
   HSMM@phenoData@data$Size_Factor = 1.0
   load("rdaFiles/rawData.main.objs.rda")
   ordering_genes = rawData.main.objs$genes.use
   
   HSMM <- setOrderingFilter(HSMM, ordering_genes)
   HSMM <- reduceDimension(HSMM, max_components = 10, pseudo_expr = 0, norm_method = "none")
   HSMM <- orderCells(HSMM, reverse = F)
   
   
   # custom plots
   xX = 1
   yY = 2
   lib_info_with_pseudo <- pData(HSMM)
   msDfr <- reducedDimS(HSMM)
   msDfr <- data.frame(t(msDfr[c(xX, yY), ]))
   par(mfrow = c(2,1))
   pdf("monocle_analysis/plots/MonoclePseudotime_FinalCells.pdf", width = 7, height = 7)
   plot(msDfr, pch = 20, cex = 1, col = as.character(HSMM@phenoData@data$my.final.colors))
   plot.new()
   legend("center", legend = names(all.colors.to.use), col = all.colors.to.use, pch = 20)
   dev.off()
   
   par(mfrow = c(4,1))
   pdf("monocle_analysis/plots/MonoclePseudotime_AllGroupingCells.pdf", width = 7, height = 7)
   
   print(plot_cell_trajectory(HSMM,show_branch_points = T, color_by="State", cell_size = 1, show_tree = F))
   print(plot_cell_trajectory(HSMM,show_branch_points = F, color_by = "final.cells", cell_size = 1, show_tree = F))
   print(plot_cell_trajectory(HSMM,show_branch_points = F, show_backbone = T,backbone_color = "black", color_by="res.1", cell_size = 1, show_tree = F))
   print(plot_cell_trajectory(HSMM,show_branch_points = F, color_by="orig.ident", cell_size = 1, show_tree = F))
   
   dev.off()
   
   new.idents = as.data.frame(HSMM@phenoData@data$State, row.names = rownames(HSMM@phenoData@data))
   colnames(new.idents) = "State"
   all.combined.main = AddMetaData(all.combined.main, metadata = new.idents)
   
   colnames(all.combined.main@meta.data)
   TSNEPlot(all.combined.main, group.by = "State", do.label = T)
   save(all.combined.main, file = "rdaFiles/all.combined.main.rda")
   save(HSMM, file = "rdaFiles/HSMM.rda")
}

# Now we have all the required information in all.combined.main/all.combined.glia/all.combined.neuron
{
   # check the idents
   library("Seurat")
   library("dplyr")
   load("rdaFiles/all.combined.main.rda")
   # find marker genes:
   final.cells.markers = FindAllMarkers(all.combined.main, min.pct = 0.1)
   save(final.cells.markers,file = "rdaFiles/final.cells.markers.rda")
   openxlsx::write.xlsx(final.cells.markers, file = "rdaFiles/final.cells.markers.xlsx", row.names = T)
   for(topN in c(5,10,20,30,40,50,100)){
      topNgenes <- final.cells.markers %>% group_by(cluster) %>% top_n(topN, avg_logFC)
      pdf(paste0("plots/final_cell_markers_heatmap_top",topN,".pdf"), width = nrow(all.combined.main@meta.data)*0.01, height = 1 +length(unique(topNgenes$gene))*0.14)
      print(DoHeatmap(object = all.combined.main, genes.use = unique(topNgenes$gene), slim.col.label = TRUE, remove.key = TRUE, group.label.rot = T))
      dev.off()
   }
   
   all.combined.main = SetIdent(all.combined.main, ident.use = all.combined.main@meta.data$main.cell.clusters)
   main.cell.clusters.markers = FindAllMarkers(all.combined.main, min.pct = 0.1)
   save(main.cell.clusters.markers,file = "rdaFiles/main.cell.clusters.markers.rda")
   openxlsx::write.xlsx(main.cell.clusters.markers, file = "rdaFiles/main.cell.clusters.markers.xlsx", row.names = T)
   for(topN in c(5,10,20,30,40,50,100)){
      topNgenes <- main.cell.clusters.markers %>% group_by(cluster) %>% top_n(topN, avg_logFC)
      pdf(paste0("plots/main_cell_clusters_markers_heatmap_top",topN,".pdf"), width = nrow(all.combined.main@meta.data)*0.01, height = 1 +length(unique(topNgenes$gene))*0.14)
      print(DoHeatmap(object = all.combined.main, genes.use = unique(topNgenes$gene), slim.col.label = TRUE, remove.key = TRUE, group.label.rot = T))
      dev.off()
   }
   
   # neurons
   load("rdaFiles/all.combined.neuron.rda")
   neuron.cells.markers = FindAllMarkers(all.combined.neuron, min.pct = 0.1)
   save(neuron.cells.markers,file = "rdaFiles/neuron.cells.markers.rda")
   openxlsx::write.xlsx(neuron.cells.markers, file = "rdaFiles/neuron.cells.markers.xlsx", row.names = T)
   for(topN in c(5,10,20,30,40,50,100)){
      topNgenes <- neuron.cells.markers %>% group_by(cluster) %>% top_n(topN, avg_logFC)
      pdf(paste0("plots/neuron_cell_markers_heatmap_top",topN,".pdf"), width = nrow(all.combined.neuron@meta.data)*0.01, height = 1 +length(unique(topNgenes$gene))*0.14)
      print(DoHeatmap(object = all.combined.neuron, genes.use = unique(topNgenes$gene), slim.col.label = TRUE, remove.key = TRUE, group.label.rot = T))
      dev.off()
   }
   
   # glia cells
   load("rdaFiles/all.combined.glia.rda")
   glia.cells.markers = FindAllMarkers(all.combined.glia, min.pct = 0.1)
   save(glia.cells.markers,file = "rdaFiles/glia.cells.markers.rda")
   openxlsx::write.xlsx(glia.cells.markers, file = "rdaFiles/glia.cells.markers.xlsx", row.names = T)
   for(topN in c(5,10,20,30,40,50,100)){
      topNgenes <- glia.cells.markers %>% group_by(cluster) %>% top_n(topN, avg_logFC)
      pdf(paste0("plots/glia_cell_markers_heatmap_top",topN,".pdf"), width = nrow(all.combined.glia@meta.data)*0.01, height = 1 +length(unique(topNgenes$gene))*0.14)
      print(DoHeatmap(object = all.combined.glia, genes.use = unique(topNgenes$gene), slim.col.label = TRUE, remove.key = TRUE, group.label.rot = T))
      dev.off()
   }
}

# Find DEGs after each treatment ...
{   
   all.combined.main = SetIdent(all.combined.main, ident.use = paste(all.combined.main@meta.data$stim, all.combined.main@meta.data$final.cells, sep = "_"))
   the_cells = names(table(all.combined.main@meta.data$final.cells))
   ctrl_cells = paste("PBS_", the_cells, sep = "")
   trt1_cells = paste("AB42_", the_cells, sep = "")
   trt2_cells = paste("IL4_", the_cells, sep = "")
   comps = c(paste(trt1_cells, ctrl_cells, sep = "_vs_"),paste(trt2_cells, ctrl_cells, sep = "_vs_"))
   ident_counts = table(all.combined.main@ident)
   ident_counts = ident_counts[ident_counts < 3]
   for(ident_name in names(ident_counts)) comps = comps[-grep(ident_name,comps)]
   all.main.DEGs = vector("list", length = length(comps))
   names(all.main.DEGs) = comps
   for(comp in comps){
      str_split = unlist(strsplit(comp, split = "_vs_"))
      print(comp)
      all.main.DEGs[[comp]] = FindMarkers(all.combined.main, ident.1 = str_split[1], ident.2 = str_split[2])
   }
   save(all.main.DEGs, file = "rdaFiles/all.main.DEGs.rda")
   openxlsx::write.xlsx(all.main.DEGs, file = "rdaFiles/all.main.DEGs.xlsx", row.names = T)
   
   all.combined.glia = StashIdent(all.combined.glia, save.name = "PC_Cells")
   all.combined.glia = SetIdent(all.combined.glia, ident.use = paste(all.combined.glia@meta.data$stim, all.combined.glia@ident,sep = "_"))
   the_cells = names(table(all.combined.glia@meta.data$PC_Cells))
   ctrl_cells = paste("PBS_", the_cells, sep = "")
   trt1_cells = paste("AB42_", the_cells, sep = "")
   trt2_cells = paste("IL4_", the_cells, sep = "")
   comps = c(paste(trt1_cells, ctrl_cells, sep = "_vs_"),paste(trt2_cells, ctrl_cells, sep = "_vs_"))
   ident_counts = table(all.combined.glia@ident)
   ident_counts = ident_counts[ident_counts < 3]
   for(ident_name in names(ident_counts)) comps = comps[-grep(ident_name,comps)]
   glia.DEGs = vector("list", length = length(comps))
   names(glia.DEGs) = comps
   for(comp in comps){
      str_split = unlist(strsplit(comp, split = "_vs_"))
      print(comp)
      glia.DEGs[[comp]] = FindMarkers(all.combined.glia, ident.1 = str_split[1], ident.2 = str_split[2])
   }
   save(glia.DEGs, file = "rdaFiles/glia.DEGs.rda")
   openxlsx::write.xlsx(glia.DEGs, file = "rdaFiles/glia.DEGs.xlsx", row.names = T)
   
   all.combined.neuron = StashIdent(all.combined.neuron, save.name = "Neuron_Cells")
   all.combined.neuron = SetIdent(all.combined.neuron, ident.use = paste(all.combined.neuron@meta.data$stim, all.combined.neuron@ident,sep = "_"))
   the_cells = names(table(all.combined.neuron@meta.data$Neuron_Cells))
   ctrl_cells = paste("PBS_", the_cells, sep = "")
   trt1_cells = paste("AB42_", the_cells, sep = "")
   trt2_cells = paste("IL4_", the_cells, sep = "")
   comps = c(paste(trt1_cells, ctrl_cells, sep = "_vs_"),paste(trt2_cells, ctrl_cells, sep = "_vs_"))
   ident_counts = table(all.combined.neuron@ident)
   ident_counts = ident_counts[ident_counts < 3]
   for(ident_name in names(ident_counts)) comps = comps[-grep(ident_name,comps)]
   neuron.DEGs = vector("list", length = length(comps))
   names(neuron.DEGs) = comps
   for(comp in comps){
      str_split = unlist(strsplit(comp, split = "_vs_"))
      print(comp)
      neuron.DEGs[[comp]] = FindMarkers(all.combined.neuron, ident.1 = str_split[1], ident.2 = str_split[2])
   }
   save(neuron.DEGs, file = "rdaFiles/neuron.DEGs.rda")
   openxlsx::write.xlsx(neuron.DEGs, file = "rdaFiles/neuron.DEGs.xlsx", row.names = T)
}

# interaction maps

# first calculate the percentage of genes in each cell type!
library("Seurat")
#load("rdaFiles/all.combined.main.rda")
dir.create("InteractionMaps", showWarnings = F)
setwd("InteractionMaps/")
all.combined = all.combined.main # just keep the orignal one!
all.combined = SetIdent(all.combined, ident.use = paste(all.combined@meta.data$stim, all.combined@meta.data$final.cells, sep = "_"))
dir.create("rdaFiles", showWarnings = F)
save(all.combined, file = "rdaFiles/all.combined.rda")
if(!("cell.perc.rda" %in% dir("rdaFiles/"))){
   ident.levels = levels(all.combined@ident)
   ident.numbers = table(all.combined@ident)
   rw.data = all.combined@raw.data
   rw.data = rw.data[,colnames(rw.data) %in% names(all.combined@ident)] # use the cells that have been used!
   dim(rw.data)
   cell.perc = as.data.frame(matrix(0, nrow = nrow(rw.data), ncol = (length(ident.levels))))
   colnames(cell.perc) = ident.levels
   rownames(cell.perc) = rownames(rw.data)
   dim(cell.perc)                                 
   ## be careful about the idents; here the idents is cell.names + treatment name
   idents = as.character(all.combined@ident)
   names(idents) = names(all.combined@ident)
   rw.data = rw.data[,match(names(idents), colnames(rw.data))]
   ncoLs = ncol(rw.data)
   print(Sys.time())
   for(i in 1:nrow(rw.data)){
      if(sum(rw.data[i,] > 0) != 0){
         tmp.result = table(idents[rw.data[i,] > 0])
         for(j in 1:length(tmp.result)){
            tmp.result[j] = (tmp.result[j]/ident.numbers[which(names(ident.numbers) == names(tmp.result)[j])])*100.0
         }
         for(j in 1:length(tmp.result)){
            cell.perc[i,which(colnames(cell.perc) == names(tmp.result)[j])] = tmp.result[j]
         }
      }
      if(i %% 1000 == 0) print(Sys.time())
   }
   print(Sys.time())
   save(cell.perc, file = "rdaFiles/cell.perc.rda")
   head(cell.perc)
}else{
   load("rdaFiles/cell.perc.rda")
}

# find ligand receptor pairs for zebrafishs!
if(!("zebrafish.ligand.receptor.pairs.rda" %in% dir("rdaFiles/"))){
   zforthos = read.delim("../../LPpairs/zf_orthos.tab", header = T, stringsAsFactors = F)
   zforthos = zforthos[zforthos$Zebrafish.gene.name != "",]
   LRpairs = read.delim("../../LPpairs/ligand_receptor_hsa.tab", header = T, stringsAsFactors = F)
   zforthos = zforthos[zforthos$Gene.name %in% c(LRpairs$Ligand,LRpairs$Receptor),]
   LRpairs = LRpairs[LRpairs$Ligand %in% zforthos$Gene.name,]
   zfpairs = c()
   uni.ligands = unique(LRpairs$Ligand)
   for(lg in uni.ligands){
      if(lg %in% zforthos$Gene.name){
         lgL = zforthos$Zebrafish.gene.name[zforthos$Gene.name %in% lg]
         lgL = unique(lgL)
         lgRs = LRpairs$Receptor[LRpairs$Ligand %in% lg]
         zlgRs = zforthos$Zebrafish.gene.name[zforthos$Gene.name %in% lgRs]
         if(!identical(zlgRs,character(0))){
            zfpairs = c(zfpairs, paste(lgL,zlgRs, sep = "@"))
         }
      }  
   }
   zfpairs = unique(zfpairs)
   zf.dfr = as.data.frame(matrix(",",length(zfpairs),2),stringsAsFactors = F)
   colnames(zf.dfr) = c("Ligand","Receptor")
   for(i in 1:length(zfpairs)){
      lgsp = unlist(strsplit(zfpairs[i],split = "@"))
      zf.dfr[i,1] = lgsp[1]
      zf.dfr[i,2] = lgsp[2]
   }
   head(zf.dfr)
   write.table(zf.dfr, file = "zebrafish.ligand.receptor.pairs.tab", sep = "\t")
   save(zf.dfr, file = "rdaFiles/zf.dfr.rda")
}else{
   load("rdaFiles/zf.dfr.rda")
}
# find cell pairs with receptors
zf.dfr = zf.dfr[zf.dfr$Ligand %in% rownames(cell.perc),]
zf.dfr = zf.dfr[zf.dfr$Receptor %in% rownames(cell.perc),]

#### get colors used for TSNEplot
new_colors = getColors()$all.colors.to.use
########### Now find all connections; for PBS
tmp.perc.cell = cell.perc > 20.0
LgRcPairs = paste(zf.dfr$Ligand, zf.dfr$Receptor, sep = "@")
treatments = c("PBS","AB42","IL4")
my_connections = vector("list",length = 3)
names(my_connections) = treatments
for(treatment in treatments){
   my_connections[[treatment]] = vector("list", length = length(LgRcPairs))   
   names(my_connections[[treatment]]) = LgRcPairs
}
names(my_connections$PBS)
rw.names = rownames(tmp.perc.cell)
for(treatment in treatments){
   cell.perc.tmp = tmp.perc.cell[,grep(treatment, colnames(tmp.perc.cell))]
   colnames(cell.perc.tmp) = gsub(paste0(treatment,"_"),"",colnames(cell.perc.tmp))
   for(LgRcPair in LgRcPairs){
      str_split = unlist(strsplit(LgRcPair, split = "@"))
      the_ligand = str_split[1]
      the_receptor = str_split[2]
      if((the_ligand %in% rw.names) & (the_receptor %in% rw.names)){
         Lg.Cells = names(cell.perc.tmp[(rw.names %in% the_ligand),])[cell.perc.tmp[rw.names %in% the_ligand,]]
         Rc.Cells = names(cell.perc.tmp[(rw.names %in% the_receptor),])[cell.perc.tmp[rw.names %in% the_receptor,]]
         Lg.Cells
         Rc.Cells
         if(!is.null(Lg.Cells) & !is.null(Rc.Cells)){
            for(Lg.Cell in Lg.Cells){
               for(Rc.Cell in Rc.Cells){
                  cell2cellcnx = paste(Lg.Cell,Rc.Cell, sep = ";")
                  my_connections[[treatment]][[LgRcPair]] = c(my_connections[[treatment]][[LgRcPair]],cell2cellcnx)
               }
               #pritn(LgRcPair)
            }
         }
      }      
   }
}

save(my_connections, file = "rdaFiles/my_connections.rda")

############ Now plot for the following pathways; STEP1
library("igraph")
dir.create("step1", showWarnings = F)
my_paths = c(paste0("notch#",2), paste(c("wnt","fgf","igf","cxc"),1, sep = "#"),paste(c("gnai", "agrn", "appa", "ctgfa", "edil3a", "hbegf", "penkb", "ptn", "serpine1", "serpine2"),1, sep = "#"))
for(treatment in names(my_connections)){
	the_cnxs = my_connections[[treatment]]
	for(my_path in my_paths){
		str_split = unlist(strsplit(my_path, split = "#"))
		LgRcPairs = zf.dfr[grep(paste0("^",str_split[1]),zf.dfr[,as.integer(str_split[2])]),]
		LgRcPairs = paste(LgRcPairs$Ligand, LgRcPairs$Receptor, sep = "@")
		LgRcPairs
		cnx = c()
		for(LgRcPair in LgRcPairs){
			if(LgRcPair %in% names(the_cnxs)){
				cnx = c(cnx,the_cnxs[[LgRcPair]])
			}
		}
		if(!is.null(cnx)){
			my.cx = unlist(strsplit(cnx, split =";"))
			all.intrx = graph(my.cx)
            my_layout = layout_(all.intrx, nicely())
			my.cx.names= names(edges(all.intrx)[[1]][1]) # the correct order of the name is important for coloring.
			my.cx.colors = new_colors[match(my.cx.names, names(new_colors))]
			my.cx.colors 
			pdf(paste("step1/",treatment,"_",str_split[1],".pdf", sep = ""))
			plot(all.intrx, layout =my_layout,edge.color = "black",edge.arrow.size=0.4, vertex.color=my.cx.colors, vertex.size=20,
			 vertex.frame.color="black", vertex.label.color="black",
			 vertex.label.cex=0.75, vertex.label.dist=0, edge.curved=0.2)
			dev.off()
		}
	}
}
############ Now for DEGs genes; STEP2
dir.create("step2", showWarnings = F)
for(my_path in my_paths){
   print(my_path)
   str_split = unlist(strsplit(my_path, split = "#"))
   LgRcPairs = zf.dfr[grep(paste0("^",str_split[1]),zf.dfr[,as.integer(str_split[2])]),]
   LgRcPairs = paste(LgRcPairs$Ligand, LgRcPairs$Receptor, sep = "@")
   LgRcPairs
   pbs_cnx = c()
   pbs_nms = c()
   ab42_cnx = c()
   ab42_nms = c()
   il4_cnx = c()
   il4_nms = c()
   for(LgRcPair in LgRcPairs){
      if(LgRcPair %in% names(my_connections[["PBS"]])){
         pbs_cnx = c(pbs_cnx,my_connections[["PBS"]][[LgRcPair]])
         pbs_nms = c(pbs_nms, rep(LgRcPair, length(my_connections[["PBS"]][[LgRcPair]])))
      }
      
      if(LgRcPair %in% names(my_connections[["AB42"]])){
         ab42_cnx = c(ab42_cnx,my_connections[["AB42"]][[LgRcPair]])
         ab42_nms = c(ab42_nms, rep(LgRcPair, length(my_connections[["AB42"]][[LgRcPair]])))
      }
      
      if(LgRcPair %in% names(my_connections[["IL4"]])){
         il4_cnx = c(il4_cnx,my_connections[["IL4"]][[LgRcPair]])
         il4_nms = c(il4_nms, rep(LgRcPair, length(my_connections[["IL4"]][[LgRcPair]])))
      }
   }
   il4_cnx
   ab42_cnx
   pbs_cnx
   
   pbs_mrg = paste(pbs_nms,pbs_cnx, sep = "#")
   ab42_mrg = paste(ab42_nms,ab42_cnx, sep = "#")
   il4_mrg = paste(il4_nms,il4_cnx, sep = "#")
   
   table(pbs_mrg %in% ab42_mrg)
   table(pbs_mrg %in% il4_mrg)
   
   pbs_only = pbs_mrg[!(pbs_mrg %in% ab42_mrg)]
   ab42_only = ab42_mrg[!(ab42_mrg %in% pbs_mrg)]
   pbs_ab42 = pbs_mrg[pbs_mrg %in% ab42_mrg]
   
   my_cnxs = c(pbs_only, pbs_ab42, ab42_only)
   
   if(!is.null(my_cnxs)){
      my_cnx_cls = rep("black",length(my_cnxs))
      my_cnx_cls[(my_cnxs %in% pbs_only)] = "cyan"
      my_cnx_cls[(my_cnxs %in% ab42_only)] = "purple"
      my_cnx = c()
      for(the_cnx in my_cnxs){
         my_split = unlist(strsplit(the_cnx, split = "#"))
         my_cnx = c(my_cnx, my_split[2])
      }
      my_cnx_sp = unlist(strsplit(my_cnx, split = ";"))
      
      all.intrx = graph(my_cnx_sp, directed = T)
      my_layout = layout_(all.intrx, nicely())
	  my.cx.names = names(edges(all.intrx)[[1]][1]) # the correct order of the name is important for coloring.
	  my.cx.colors = new_colors[match(my.cx.names, names(new_colors))]
	  pdf(paste0("step2/Ab42_vs_PBS_",str_split[1],".pdf"))
	  plot(all.intrx, layout = my_layout, edge.color = my_cnx_cls, edge.arrow.size=0.4, vertex.color=my.cx.colors, vertex.size=20,
           vertex.frame.color="black", vertex.label.color="black",
           vertex.label.cex=0.75, vertex.label.dist=0, curve_multiple(all.intrx, .2))
      dev.off()
   }
   
   pbs_only = pbs_mrg[!(pbs_mrg %in% il4_mrg)]
   il4_only = il4_mrg[!(il4_mrg %in% pbs_mrg)]
   pbs_il4 = pbs_mrg[pbs_mrg %in%il4_mrg]
   
   my_cnxs = c(pbs_only, pbs_il4, il4_only)
   
   if(!is.null(my_cnxs)){
      my_cnx_cls = rep("black",length(my_cnxs))
      my_cnx_cls[(my_cnxs %in% pbs_only)] = "cyan"
      my_cnx_cls[(my_cnxs %in% il4_only)] = "purple"
      my_cnx_cls[(my_cnxs %in% pbs_il4)] = "black"
      my_cnx = c()
      for(the_cnx in my_cnxs){
         my_split = unlist(strsplit(the_cnx, split = "#"))
         my_cnx = c(my_cnx, my_split[2])
      }
      my_cnx_sp = unlist(strsplit(my_cnx, split = ";"))
      
      all.intrx = graph(my_cnx_sp, directed = T)
      my_layout = layout_(all.intrx, nicely())
      my.cx.names= names(edges(all.intrx)[[1]][1]) # the correct order of the name is important for coloring.
      my.cx.colors = new_colors[match(my.cx.names, names(new_colors))]
      pdf(paste0("step2/IL4_vs_PBS_",str_split[1],".pdf"))
      plot(all.intrx, layout = my_layout, edge.color = my_cnx_cls,edge.arrow.size=0.4, vertex.color=my.cx.colors, vertex.size=20,
           vertex.frame.color="black", vertex.label.color="black",
           vertex.label.cex=0.75, vertex.label.dist=0, edge.curved = curve_multiple(all.intrx, .2))
      dev.off()
   }
   
}

setwd("..")

}
{
runMe = TRUE
# below is additional analysis
if(runMe){
   #setwd("ranOn_20190131153247/")
   #load(".RData")
   dir.create("monocle_analysis.glia", showWarnings = F)
   dir.create("monocle_analysis.glia/rdaFiles", showWarnings = F)
   dir.create("monocle_analysis.glia/plots", showWarnings = F)
   
   library("Seurat")
   library("monocle")
   
   ph.dt = all.combined.glia@meta.data
   ft.dt = as.data.frame(rownames(all.combined.glia@scale.data), row.names = rownames(all.combined.glia@scale.data))
   head(ft.dt)
   colnames(ft.dt) = "gene_short_name"
   pd = new(ph.dt, Class ="AnnotatedDataFrame")
   ft = new(ft.dt, Class ="AnnotatedDataFrame")
   HSMM = newCellDataSet(cellData = all.combined.glia@scale.data, phenoData = pd, featureData = ft)
   
   HSMM@phenoData@data$Size_Factor = 1.0
   load("rdaFiles/rawData.glia.objs.rda")
   ordering_genes = rawData.glia.objs$genes.use
   
   HSMM <- setOrderingFilter(HSMM, ordering_genes)
   HSMM <- reduceDimension(HSMM, max_components = 4, pseudo_expr = 0, norm_method = "none")
   HSMM <- orderCells(HSMM, reverse = F)
   print(plot_cell_trajectory(HSMM,show_branch_points = T, color_by="State", cell_size = 1, show_tree = F))
   
   my_rg_colors = paste("PC_",all.combined.glia@meta.data$res.1, sep = "")
   my_colors = getColors()$all.colors.to.use
   my_rg_colors = my_colors[match(my_rg_colors, names(my_colors))]
   
   # custom plots
   xX = 1
   yY = 2
   lib_info_with_pseudo <- pData(HSMM)
   msDfr <- reducedDimS(HSMM)
   msDfr <- data.frame(t(msDfr[c(xX, yY), ]))
   par(mfrow = c(2,1))
   pdf("monocle_analysis.glia/plots/MonoclePseudotime_FinalCells.pdf", width = 7, height = 7)
   plot(msDfr, pch = 20, cex = 1, col = my_rg_colors, xlab = "Component1",ylab = "Component2")
   plot.new()
   legend("center", legend = names(all.colors.to.use), col = all.colors.to.use, pch = 20)
   dev.off()
   
   
   
   
   HSMM@phenoData$final.cells.cols = my_rg_colors
   HSMM@phenoData$final.cells = names(my_rg_colors)
   
   
   par(mfrow = c(3,1))
   pdf("monocle_analysis.glia/plots/MonoclePseudotime_AllGroupingCells.pdf", width = 7, height = 7)
   p20 = print(plot_cell_trajectory(HSMM,show_branch_points = T, color_by="State", cell_size = 0.5, show_tree = F))
   p21 = print(plot_cell_trajectory(HSMM,show_branch_points = F, show_backbone = T,backbone_color = "black", color_by="res.1", cell_size = 0.5, show_tree = F, use_color_gradient = F))
   p22 = print(plot_cell_trajectory(HSMM,show_branch_points = F, show_backbone = T,backbone_color = "black", color_by="final.cells", cell_size = 1, show_tree = F, col = my_rg_colors))
   p23 = print(plot_cell_trajectory(HSMM,show_branch_points = F, color_by="orig.ident", cell_size = 0.5, show_tree = F))
   dev.off()
   pdf("monocle_analysis.glia/plots/MonoclePseudotime_AllGroupingCellsPltGrid.pdf", width = 14, height = 14)
   plot_grid(p20,p21,p22,p23, nrow = 2)
   dev.off()
   
   new.idents = as.data.frame(HSMM@phenoData@data$State, row.names = rownames(HSMM@phenoData@data))
   colnames(new.idents) = "State"
   all.combined.glia = AddMetaData(all.combined.glia, metadata = new.idents)
   
   colnames(all.combined.glia@meta.data)
   TSNEPlot(all.combined.glia, group.by = "State", do.label = T)
   save(all.combined.glia, file = "rdaFiles/all.combined.glia.rda")
   save(HSMM, file = "rdaFiles/HSMM.glia.rda")
}
{
### new add
dir.create("all_interaction_maps", showWarnings = F)
setwd("all_interaction_maps/")
dir.create("rdaFiles", showWarnings = F)

# first calculate the percentage of genes in each cell type!
library("Seurat")
load("../rdaFiles/all.combined.main.rda")
all.combined = all.combined.main # just keep the orignal one!
all.combined = SetIdent(all.combined, ident.use = paste(all.combined@meta.data$stim, all.combined@meta.data$final.cells, sep = "_"))
save(all.combined, file = "rdaFiles/all.combined.rda")
if(!("cell.perc.rda" %in% dir("rdaFiles/"))){
   ident.levels = levels(all.combined@ident)
   ident.numbers = table(all.combined@ident)
   rw.data = all.combined@raw.data
   rw.data = rw.data[,colnames(rw.data) %in% names(all.combined@ident)] # use the cells that have been used!
   dim(rw.data)
   cell.perc = as.data.frame(matrix(0, nrow = nrow(rw.data), ncol = (length(ident.levels))))
   colnames(cell.perc) = ident.levels
   rownames(cell.perc) = rownames(rw.data)
   dim(cell.perc)                                 
   ## be careful about the idents; here the idents is cell.names + treatment name
   idents = as.character(all.combined@ident)
   names(idents) = names(all.combined@ident)
   rw.data = rw.data[,match(names(idents), colnames(rw.data))]
   ncoLs = ncol(rw.data)
   print(Sys.time())
   for(i in 1:nrow(rw.data)){
      if(sum(rw.data[i,] > 0) != 0){
         tmp.result = table(idents[rw.data[i,] > 0])
         for(j in 1:length(tmp.result)){
            tmp.result[j] = (tmp.result[j]/ident.numbers[which(names(ident.numbers) == names(tmp.result)[j])])*100.0
         }
         for(j in 1:length(tmp.result)){
            cell.perc[i,which(colnames(cell.perc) == names(tmp.result)[j])] = tmp.result[j]
         }
      }
      if(i %% 1000 == 0) print(Sys.time())
   }
   print(Sys.time())
   save(cell.perc, file = "rdaFiles/cell.perc.rda")
   head(cell.perc)
}else{
   load("rdaFiles/cell.perc.rda")
}

# find ligand receptor pairs for zebrafishs!
if(!("zebrafish.ligand.receptor.pairs.rda" %in% dir("rdaFiles/"))){
   zforthos = read.delim("../../LPpairs/zf_orthos.tab", header = T, stringsAsFactors = F)
   zforthos = zforthos[zforthos$Zebrafish.gene.name != "",]
   LRpairs = read.delim("../../LPpairs/ligand_receptor_hsa.tab", header = T, stringsAsFactors = F)
   zforthos = zforthos[zforthos$Gene.name %in% c(LRpairs$Ligand,LRpairs$Receptor),]
   LRpairs = LRpairs[LRpairs$Ligand %in% zforthos$Gene.name,]
   zfpairs = c()
   uni.ligands = unique(LRpairs$Ligand)
   for(lg in uni.ligands){
      if(lg %in% zforthos$Gene.name){
         lgL = zforthos$Zebrafish.gene.name[zforthos$Gene.name %in% lg]
         lgL = unique(lgL)
         lgRs = LRpairs$Receptor[LRpairs$Ligand %in% lg]
         zlgRs = zforthos$Zebrafish.gene.name[zforthos$Gene.name %in% lgRs]
         if(!identical(zlgRs,character(0))){
            zfpairs = c(zfpairs, paste(lgL,zlgRs, sep = "@"))
         }
      }  
   }
   zfpairs = unique(zfpairs)
   zf.dfr = as.data.frame(matrix(",",length(zfpairs),2),stringsAsFactors = F)
   colnames(zf.dfr) = c("Ligand","Receptor")
   for(i in 1:length(zfpairs)){
      lgsp = unlist(strsplit(zfpairs[i],split = "@"))
      zf.dfr[i,1] = lgsp[1]
      zf.dfr[i,2] = lgsp[2]
   }
   head(zf.dfr)
   write.table(zf.dfr, file = "zebrafish.ligand.receptor.pairs.tab", sep = "\t")
   save(zf.dfr, file = "rdaFiles/zf.dfr.rda")
}else{
   load("rdaFiles/zf.dfr.rda")
}
# find cell pairs with receptors
zf.dfr = zf.dfr[zf.dfr$Ligand %in% rownames(cell.perc),]
zf.dfr = zf.dfr[zf.dfr$Receptor %in% rownames(cell.perc),]

#### get colors used for TSNEplot
new_colors = getColors()$all.colors.to.use
########### Now find all connections; for PBS
tmp.perc.cell = cell.perc > 20.0
LgRcPairs = paste(zf.dfr$Ligand, zf.dfr$Receptor, sep = "@")
treatments = c("PBS","AB42","IL4")
my_connections = vector("list",length = 3)
names(my_connections) = treatments
for(treatment in treatments){
   my_connections[[treatment]] = vector("list", length = length(LgRcPairs))   
   names(my_connections[[treatment]]) = LgRcPairs
}
names(my_connections$PBS)
rw.names = rownames(tmp.perc.cell)
for(treatment in treatments){
   cell.perc.tmp = tmp.perc.cell[,grep(treatment, colnames(tmp.perc.cell))]
   colnames(cell.perc.tmp) = gsub(paste0(treatment,"_"),"",colnames(cell.perc.tmp))
   for(LgRcPair in LgRcPairs){
      str_split = unlist(strsplit(LgRcPair, split = "@"))
      the_ligand = str_split[1]
      the_receptor = str_split[2]
      if((the_ligand %in% rw.names) & (the_receptor %in% rw.names)){
         Lg.Cells = names(cell.perc.tmp[(rw.names %in% the_ligand),])[cell.perc.tmp[rw.names %in% the_ligand,]]
         Rc.Cells = names(cell.perc.tmp[(rw.names %in% the_receptor),])[cell.perc.tmp[rw.names %in% the_receptor,]]
         Lg.Cells
         Rc.Cells
         if(!is.null(Lg.Cells) & !is.null(Rc.Cells)){
            for(Lg.Cell in Lg.Cells){
               for(Rc.Cell in Rc.Cells){
                  cell2cellcnx = paste(Lg.Cell,Rc.Cell, sep = ";")
                  my_connections[[treatment]][[LgRcPair]] = c(my_connections[[treatment]][[LgRcPair]],cell2cellcnx)
               }
               #pritn(LgRcPair)
            }
         }
      }      
   }
}

save(my_connections, file = "rdaFiles/my_connections.rda")

############ Now plot for the following pathways; STEP1
library("igraph")
dir.create("step1", showWarnings = F)
#my_paths = c(paste0("notch#",2), paste(c("wnt","fgf","igf","cxc"),1, sep = "#"),paste(c("gnai", "agrn", "appa", "ctgfa", "edil3a", "hbegf", "penkb", "ptn", "serpine1", "serpine2"),1, sep = "#"))
my_paths_ligands = paste(unique(zf.dfr$Ligand),1, sep = "#")
my_paths_receptors = paste(unique(zf.dfr$Receptor),2, sep = "#")
my_paths = c(my_paths_ligands,my_paths_receptors)

for(treatment in names(my_connections)){
   the_cnxs = my_connections[[treatment]]
   for(my_path in my_paths){
        print(my_path)
      str_split = unlist(strsplit(my_path, split = "#"))
      LgRcPairs = zf.dfr[grep(paste0("^",str_split[1]),zf.dfr[,as.integer(str_split[2])]),]
      LgRcPairs = paste(LgRcPairs$Ligand, LgRcPairs$Receptor, sep = "@")
      LgRcPairs
      cnx = c()
      for(LgRcPair in LgRcPairs){
         if(LgRcPair %in% names(the_cnxs)){
            cnx = c(cnx,the_cnxs[[LgRcPair]])
         }
      }
      if(!is.null(cnx)){
         my.cx = unlist(strsplit(cnx, split =";"))
         all.intrx = graph(my.cx)
         my_layout = layout_(all.intrx, nicely())
         my.cx.names= names(edges(all.intrx)[[1]][1]) # the correct order of the name is important for coloring.
         my.cx.colors = new_colors[match(my.cx.names, names(new_colors))]
         my.cx.colors 
         pdf(paste("step1/",treatment,"_",gsub(":","",str_split[1]),".pdf", sep = ""))
         plot(all.intrx, layout =my_layout,edge.color = "black",edge.arrow.size=0.4, vertex.color=my.cx.colors, vertex.size=20,
              vertex.frame.color="black", vertex.label.color="black",
              vertex.label.cex=0.75, vertex.label.dist=0, edge.curved=0.2)
         dev.off()
      }
   }
}


############ Now for DEGs genes; STEP2
dir.create("step2", showWarnings = F)
for(my_path in my_paths){
   print(my_path)
   str_split = unlist(strsplit(my_path, split = "#"))
   LgRcPairs = zf.dfr[grep(paste0("^",str_split[1]),zf.dfr[,as.integer(str_split[2])]),]
   LgRcPairs = paste(LgRcPairs$Ligand, LgRcPairs$Receptor, sep = "@")
   LgRcPairs
   pbs_cnx = c()
   pbs_nms = c()
   ab42_cnx = c()
   ab42_nms = c()
   il4_cnx = c()
   il4_nms = c()
   for(LgRcPair in LgRcPairs){
      if(LgRcPair %in% names(my_connections[["PBS"]])){
         pbs_cnx = c(pbs_cnx,my_connections[["PBS"]][[LgRcPair]])
         pbs_nms = c(pbs_nms, rep(LgRcPair, length(my_connections[["PBS"]][[LgRcPair]])))
      }
      
      if(LgRcPair %in% names(my_connections[["AB42"]])){
         ab42_cnx = c(ab42_cnx,my_connections[["AB42"]][[LgRcPair]])
         ab42_nms = c(ab42_nms, rep(LgRcPair, length(my_connections[["AB42"]][[LgRcPair]])))
      }
      
      if(LgRcPair %in% names(my_connections[["IL4"]])){
         il4_cnx = c(il4_cnx,my_connections[["IL4"]][[LgRcPair]])
         il4_nms = c(il4_nms, rep(LgRcPair, length(my_connections[["IL4"]][[LgRcPair]])))
      }
   }
   il4_cnx
   ab42_cnx
   pbs_cnx
   
   pbs_mrg = paste(pbs_nms,pbs_cnx, sep = "#")
   ab42_mrg = paste(ab42_nms,ab42_cnx, sep = "#")
   il4_mrg = paste(il4_nms,il4_cnx, sep = "#")
   
   table(pbs_mrg %in% ab42_mrg)
   table(pbs_mrg %in% il4_mrg)
   
   pbs_only = pbs_mrg[!(pbs_mrg %in% ab42_mrg)]
   ab42_only = ab42_mrg[!(ab42_mrg %in% pbs_mrg)]
   pbs_ab42 = pbs_mrg[pbs_mrg %in% ab42_mrg]
   
   my_cnxs = c(pbs_only, pbs_ab42, ab42_only)
   
   if(!is.null(my_cnxs) & length(my_cnxs) != 0){
      my_cnx_cls = rep("black",length(my_cnxs))
      my_cnx_cls[(my_cnxs %in% pbs_only)] = "cyan"
      my_cnx_cls[(my_cnxs %in% ab42_only)] = "purple"
      my_cnx = c()
      for(the_cnx in my_cnxs){
         my_split = unlist(strsplit(the_cnx, split = "#"))
         my_cnx = c(my_cnx, my_split[2])
      }
      my_cnx_sp = unlist(strsplit(my_cnx, split = ";"))
      
      all.intrx = graph(my_cnx_sp, directed = T)
      my_layout = layout_(all.intrx, nicely())
      my.cx.names = names(edges(all.intrx)[[1]][1]) # the correct order of the name is important for coloring.
      my.cx.colors = new_colors[match(my.cx.names, names(new_colors))]
      pdf(paste0("step2/Ab42_vs_PBS_",gsub(":","",str_split[1]),".pdf"))
      plot(all.intrx, layout = my_layout, edge.color = my_cnx_cls, edge.arrow.size=0.4, vertex.color=my.cx.colors, vertex.size=20,
           vertex.frame.color="black", vertex.label.color="black",
           vertex.label.cex=0.75, vertex.label.dist=0, curve_multiple(all.intrx, .2))
      dev.off()
   }
   
   pbs_only = pbs_mrg[!(pbs_mrg %in% il4_mrg)]
   il4_only = il4_mrg[!(il4_mrg %in% pbs_mrg)]
   pbs_il4 = pbs_mrg[pbs_mrg %in%il4_mrg]
   
   my_cnxs = c(pbs_only, pbs_il4, il4_only)
   
   if(!is.null(my_cnxs)  & length(my_cnxs) != 0){
      my_cnx_cls = rep("black",length(my_cnxs))
      my_cnx_cls[(my_cnxs %in% pbs_only)] = "cyan"
      my_cnx_cls[(my_cnxs %in% il4_only)] = "purple"
      my_cnx_cls[(my_cnxs %in% pbs_il4)] = "black"
      my_cnx = c()
      for(the_cnx in my_cnxs){
         my_split = unlist(strsplit(the_cnx, split = "#"))
         my_cnx = c(my_cnx, my_split[2])
      }
      my_cnx_sp = unlist(strsplit(my_cnx, split = ";"))
      
      all.intrx = graph(my_cnx_sp, directed = T)
      my_layout = layout_(all.intrx, nicely())
      my.cx.names= names(edges(all.intrx)[[1]][1]) # the correct order of the name is important for coloring.
      my.cx.colors = new_colors[match(my.cx.names, names(new_colors))]
      pdf(paste0("step2/IL4_vs_PBS_",gsub(":","",str_split[1]),".pdf"))
      plot(all.intrx, layout = my_layout, edge.color = my_cnx_cls,edge.arrow.size=0.4, vertex.color=my.cx.colors, vertex.size=20,
           vertex.frame.color="black", vertex.label.color="black",
           vertex.label.cex=0.75, vertex.label.dist=0, edge.curved = curve_multiple(all.intrx, .2))
      dev.off()
   }
   
}
}

   setwd("..")
   
   sink("sessionInfo.txt")
   sessionInfo()
   sink()
   t2 = Sys.time()
}
runMe = TRUE
if(runMe){
   # perform GO analysis
   # some custom functions
   scs_fnxs = "../scs_functions/"
   fnxs = dir(scs_fnxs)
   fnxs = fnxs[grep(".R$",fnxs)]
   for(fnx in fnxs) source(paste0(scs_fnxs,fnx))
   dir.create("DEGGoTermAnalysis", showWarnings = F)
   setwd("DEGGoTermAnalysis/")
   library("GOstats")
   library("GSEABase")
   library("KEGG.db")
   library("AnnotationDbi")
   library("openxlsx")
   gene.info = read.delim("../../10X_Data_AutomatedCellsCounts/PBS_Cells/genes.tsv", header = F, stringsAsFactors = F)
   uniGenes = gene.info$V1
   load("../rdaFiles/all.main.DEGs.rda")
   tmp.degs = all.main.DEGs
   names(tmp.degs)
   deg.gokegg = tmp.degs
   for(cmp in names(tmp.degs)){
      tmp.dfr = tmp.degs[[cmp]]
      deg.genes = rownames(tmp.dfr)[tmp.dfr$p_val < 0.05 & abs(tmp.dfr$avg_logFC) >= 0.5]
      selGenes = gene.info$V1[gene.info$V2 %in% deg.genes]
      deg.gokegg[[cmp]] = runGOStats(selGenes = selGenes, uniGenes = uniGenes, orgId = "Dr", geneIDtype = "ENSEMBL")
      write.xlsx(deg.gokegg[[cmp]], file = paste0(cmp,"_DEG_GOKEGG.xlsx"), row.names = T)
   }
   save(deg.gokegg, file = "deg.gokegg.rda")

   load("../rdaFiles/final.cells.markers.rda")
   top.marker.genes = final.cells.markers
   head(top.marker.genes)
   the_cells = levels(top.marker.genes$cluster) 
   final.gokegg = vector("list", length = length(the_cells))
   names(final.gokegg) = the_cells
   for(the_cell in the_cells){
      deg.genes = top.marker.genes[(top.marker.genes$cluster %in% the_cell),]
      deg.genes = subset(deg.genes, p_val_adj < 0.1 & avg_logFC >= 1)
      deg.genes = subset(deg.genes, pct.1 >= 0.2)
      print(c(the_cell, ":",dim(deg.genes)))
      if(nrow(deg.genes) >= 10){
         selGenes = gene.info$V1[gene.info$V2 %in% deg.genes$gene]
         final.gokegg[[the_cell]] = runGOStats(selGenes = selGenes, uniGenes = uniGenes, orgId = "Dr", geneIDtype = "ENSEMBL")
         write.xlsx(final.gokegg[[the_cell]], file = paste0(the_cell,"_finalmarkers_GOKEGG.xlsx"), row.names = T)
      }
   }
   save(final.gokegg, file = "final.gokegg.rda")


   load("../rdaFiles/glia.cells.markers.rda")
   top.marker.genes = glia.cells.markers
   head(top.marker.genes)
   the_cells = levels(top.marker.genes$cluster) 
   glia.gokegg = vector("list", length = length(the_cells))
   names(glia.gokegg) = the_cells
   for(the_cell in the_cells){
      deg.genes = top.marker.genes[(top.marker.genes$cluster %in% the_cell),]
      deg.genes = subset(deg.genes, p_val_adj < 0.1 & avg_logFC >= 0.25)
      #deg.genes = subset(deg.genes, pct.1 >= 0.2)
      print(c(the_cell, ":",dim(deg.genes)))
      if(nrow(deg.genes) >= 10){
         selGenes = gene.info$V1[gene.info$V2 %in% deg.genes$gene]
         glia.gokegg[[the_cell]] = runGOStats(selGenes = selGenes, uniGenes = uniGenes, orgId = "Dr", geneIDtype = "ENSEMBL")
         write.xlsx(glia.gokegg[[the_cell]], file = paste0(the_cell,"_gliemarkers_GOKEGG.xlsx"), row.names = T)
      }
   }
   save(glia.gokegg, file = "glia.gokegg.rda")

   load("../rdaFiles/neuron.cells.markers.rda")
   top.marker.genes = neuron.cells.markers
   head(top.marker.genes)
   the_cells = levels(top.marker.genes$cluster) 
   neuron.gokegg = vector("list", length = length(the_cells))
   names(neuron.gokegg) = the_cells
   for(the_cell in the_cells){
      deg.genes = top.marker.genes[(top.marker.genes$cluster %in% the_cell),]
      deg.genes = subset(deg.genes, p_val_adj < 0.1 & avg_logFC >= 0.25)
      #deg.genes = subset(deg.genes, pct.1 >= 0.2)
      print(c(the_cell, ":",dim(deg.genes)))
      if(nrow(deg.genes) >= 10){
         selGenes = gene.info$V1[gene.info$V2 %in% deg.genes$gene]
         neuron.gokegg[[the_cell]] = runGOStats(selGenes = selGenes, uniGenes = uniGenes, orgId = "Dr", geneIDtype = "ENSEMBL")
         write.xlsx(neuron.gokegg[[the_cell]], file = paste0(the_cell,"_gliemarkers_GOKEGG.xlsx"), row.names = T)
      }
   }
   save(neuron.gokegg, file = "neuron.gokegg.rda")
   
   load("deg.gokegg.rda")
   dir.create("my_deg_pdfs", showWarnings = F)
   q1 = deg.gokegg
   for(nm in names(q1)){
    q2 = q1[[nm]]
    for (nn in names(q2)){
       q3 = q2[[nn]]   
       q3 = subset(q3, Pvalue < 0.05)
       drawPie(Df = q3, outputName = paste(nm,"_", nn,".pdf", sep = ""), output_folder = "my_deg_pdfs", pvalCol = 2, countCol = 5, totalCol = 6, maxRes = 50, pvalCutoff = 0.05, labCol = 7)
    }
   }

   dir.create("forWEB", showWarnings = F)
   for(deg_comp in names(deg.gokegg)){
    str_split = unlist(strsplit(deg_comp, split = "_"))
    if(length(str_split) == 5){
       CellID = str_split[5]
       comp = paste(str_split[c(1,3,4)],sep = "", collapse = "_")
    }else{
       CellID = paste(str_split[c(6,7)],sep = "", collapse = "_")
       comp = paste(str_split[c(1,4,5)],sep = "", collapse = "_")
    }
    CellID
    comp
    dir.create(paste0("forWEB/",comp), showWarnings = F)
    for(gg in names(deg.gokegg[[deg_comp]])){
       openxlsx::write.xlsx(deg.gokegg[[deg_comp]][[gg]], file = paste0("forWEB/",comp,"/",CellID,"_",gg,".xlsx"), row.names = T)
    }
   }
   
   }