# updated on 16.11.2018
runGOStats <- function(orgId = "Dr", uniGenes = NULL, selGenes = NULL, pvalCutoff = 0.05, geneIDtype = "ENTREZ", organism = "Danio rerio",mainFolder = "../../",...){	
   files = dir(paste0(mainFolder,"rdaFilesGOKEGG/"))
   if(!(orgId %in% c("Dm","Dr","Mm","Hs","Rn"))){
      stop("orgId must be 'Dr', 'Mm' or 'Hs'")
   }
   if(!(orgId %in% c("Dm","Dr","Mm","Hs","Rn"))){
      stop("orgId must be 'Dr', 'Mm' or 'Hs'")
   }
   if(geneIDtype == "ENTREZ"){
      efileName <- paste("gsc",orgId,"entrezids.rda", sep = "")
      kfileName <- paste("kgsc",orgId,"entrezids.rda", sep = "")
   }else if(geneIDtype == "ENSEMBL"){
      efileName <- paste("gsc",orgId,"ensemblids.rda", sep = "")
      kfileName <- paste("kgsc",orgId,"ensemblids.rda", sep = "")
   }else{
      stop("gene id type must be ENSEMBL or ENTREZ")
   }
   # find a way to decrease time here, when there are multiple analyses
   library(paste("org",orgId,"eg.db", sep = "."), character.only = T)
   egENS <- keys(eval(parse(text = paste("org",orgId,"eg.db", sep = "."))))
   egENS <- AnnotationDbi::select(eval(parse(text = paste("org",orgId,"eg.db", sep = "."))), keys = egENS, columns = "ENSEMBL", keytype = "ENTREZID")
  
   if(!(efileName %in% files)){
      dframe <- getGOframe(geneIDtype, eval(parse(text = paste("org",orgId,"egGO", sep = "."))), eval(parse(text = paste("org",orgId,"eg.db", sep = "."))))
      goframeData = data.frame(dframe$go_id, dframe$Evidence, dframe$gene_id)
      goFrame = GOFrame(goframeData, organism = organism)
      goAllFrame = GOAllFrame(goFrame)
      gscOrg <- GeneSetCollection(goAllFrame, setType = GOCollection())
      save(gscOrg, file = paste(paste0(mainFolder,"rdaFilesGOKEGG/"),efileName, sep = ""))
   }else{
      load(file = paste(paste0(mainFolder,"rdaFilesGOKEGG/"),efileName, sep = ""))
   }
   
  if(!(kfileName %in% files)){
    kframe = toTable(eval(parse(text = paste("org",orgId,"egPATH", sep = "."))))
    keggFrameData = data.frame(kframe$path_id, kframe$gene_id)
    keggFrame=KEGGFrame(keggFrameData,organism = organism)
    kgscOrg <- GeneSetCollection(keggFrame, setType = KEGGCollection())
    save(kgscOrg, file = paste(paste0(mainFolder,"rdaFilesGOKEGG/"),kfileName, sep = ""))
  }else{
    load(file = paste(paste0(mainFolder,"rdaFilesGOKEGG/"),kfileName, sep = ""))
  }
  goTerms <- c("BP","CC","MF")
  result <- list("BP_Up" = "NA","CC_Up" = "NA","MF_Up" = "NA", "kegg_Up" = "NA")
  GOstatssm <- GOstats::summary # summary function must be from GOstats!!!
 for(goTerm in goTerms){
   hypTest <- paste(goTerm,"_Up", sep = "")
   print(hypTest)
   params <- GSEAGOHyperGParams(name="Custom en2GO GSEABase",
                                geneSetCollection = gscOrg,
                                geneIds = selGenes,
                                universeGeneIds = uniGenes,
                                ontology = goTerm,
                                pvalueCutoff = pvalCutoff,
                                conditional = FALSE,
                                testDirection = "over")
   Over <- hyperGTest(params)
   result[[hypTest]] = GOstatssm(Over)
   }

   # perform kegg pathway analysis.
  if(geneIDtype == "ENSEMBL"){
    uniGenes <- egENS$ENTREZID[egENS$ENSEMBL %in% uniGenes]
    selGenes <- egENS$ENTREZID[egENS$ENSEMBL %in% selGenes]
  }
  print("kegg_Up")
  kparams <- GSEAKEGGHyperGParams(name = "Custom en2kegg GSEABase",
                                  geneSetCollection= kgscOrg,
                                  geneIds = selGenes,
                                  universeGeneIds = uniGenes,
                                  pvalueCutoff = pvalCutoff,
                                  testDirection = "over")
  kOver <- hyperGTest(kparams)	
  result[["kegg_Up"]] = GOstatssm(kOver)
  return(result)
}
