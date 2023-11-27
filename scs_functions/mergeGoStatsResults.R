mergeGoStatsResults <- function(rsList){
  list_length = lapply(rsList,length)
  rsList = rsList[list_length != 0]
  length(rsList)
  if(length(rsList) <= 1){
    return(rsList[[1]])
  }else{
    cellID = names(rsList)[[1]]
    merged_dfr = rsList[[1]][,c(1:2,5:6)]
    rownames(merged_dfr) = merged_dfr[,1]
    merged_dfr[,1] = NULL
    colnames(merged_dfr) = paste(cellID,c("Pvalue","Count","Size"), sep = ":")
    for(i in 2:length(rsList)){
      cellID = names(rsList)[[i]]
      tmp_dfr = rsList[[i]][,c(1:2,5:6)]
      rownames(tmp_dfr) = tmp_dfr[,1]
      tmp_dfr[,1] = NULL
      colnames(tmp_dfr) = paste(cellID,c("Pvalue","Count","Size"), sep = ":")
      merged_dfr = merge(merged_dfr,tmp_dfr, by = "row.names", all = T)
      rownames(merged_dfr) = merged_dfr$Row.names
      merged_dfr$Row.names = NULL
    }
    return(merged_dfr)
  }
}
