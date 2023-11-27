getHistos <- function(my.data, folderName = "out", moreGenes = NULL,mt_prefix = "^mt-"){
  # remove genes that are not in the rownames
  moreGenes = moreGenes[moreGenes %in% rownames(my.data)]
  if(identical(moreGenes,character(0))) moreGenes = NULL
  print(moreGenes)
  dir.create(paste0("plots/",folderName), showWarnings = F)
  my.data.colSums = colSums(my.data)
  my.data.genes = colSums(my.data > 0)
  my.data.mito.genes = colSums(my.data[grep(mt_prefix,rownames(my.data)),])/colSums(my.data)
  if(is.null(moreGenes)){
    pdf(paste0("plots/",folderName,"/histos.pdf"),width = 9, height = 4.5)
    par(mfrow = c(1,3))
    hist(my.data.colSums, main = folderName, xlab = "nUMI")
    hist(my.data.genes, main = folderName, xlab = "nGene")
    hist(my.data.mito.genes, main = folderName, xlab = "percent.mito")
  }else{
    pdf(paste0("plots/",folderName,"/histos.pdf"), width = (9 + length(moreGenes)*3), height = 4.5)
    par(mfrow = c(1,(3+length(moreGenes))))
    hist(my.data.colSums, main = folderName, xlab = "nUMI")
    hist(my.data.genes, main = folderName, xlab = "nGene")
    hist(my.data.mito.genes, main = folderName, xlab = "percent.mito")
    if(length(moreGenes) == 1){
      gene = moreGenes
      my.data.more.genes = my.data[grep(paste0("^",gene,"$"),rownames(my.data)),]/colSums(my.data)
      hist(my.data.more.genes, main = folderName, xlab = paste0("percent.",gene))
    }else{
      for(gene in moreGenes){
        my.data.more.genes = my.data[grep(paste0("^",gene,"$"),rownames(my.data)),]/colSums(my.data)
        hist(my.data.more.genes, main = folderName, xlab = paste0("percent.",gene))
      }
    }
  }
  dev.off()
}