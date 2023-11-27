runGOStats <- function(do_GO_analysis = FALSE, orgId = "Dr", uniGenes = NULL, selGenes = NULL, output_folder = NULL, pvalCutoff = 0.05, geneIDtype = "ENTREZ", organism = "Danio rerio",...){	
  files = dir(paste(paste0(mainFolder,"/rdaFiles/"), sep = ""))
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
  egENS <- keys(eval(parse(text = paste("org",orgId,"eg.db", sep = "."))))
  egENS <- AnnotationDbi::select(eval(parse(text = paste("org",orgId,"eg.db", sep = "."))), keys = egENS, columns = "ENSEMBL", keytype = "ENTREZID")
  #
  print("step 1 ...")
  if(!(efileName %in% files)){
    dframe <- getGOframe(geneIDtype, eval(parse(text = paste("org",orgId,"egGO", sep = "."))), eval(parse(text = paste("org",orgId,"eg.db", sep = "."))))
    goframeData = data.frame(dframe$go_id, dframe$Evidence, dframe$gene_id)
    goFrame = GOFrame(goframeData, organism = organism)
    goAllFrame = GOAllFrame(goFrame)
    gscOrg <- GeneSetCollection(goAllFrame, setType = GOCollection())
    save(gscOrg, file = paste(paste0("rdaFiles/"),efileName, sep = ""))
  }else{
    load(file = paste(paste0("rdaFiles/"),efileName, sep = ""))
  }
  
  print("step 2 ...")
  if(!(kfileName %in% files)){
    kframe = toTable(eval(parse(text = paste("org",orgId,"egPATH", sep = "."))))
    keggFrameData = data.frame(kframe$path_id, kframe$gene_id)
    keggFrame=KEGGFrame(keggFrameData,organism = organism)
    kgscOrg <- GeneSetCollection(keggFrame, setType = KEGGCollection())
    save(kgscOrg, file = paste(paste0("rdaFiles/"),kfileName, sep = ""))
  }else{
    load(file = paste(paste0("rdaFiles/"),kfileName, sep = ""))
  }
  
  print("step 1 ...")
  goTerms <- c("BP","CC","MF")
  #result <- list("MF_Up" = NA,"MF_Down" = NA,"CC_Up" = NA,"CC_Down" = NA,"BP_Up" = NA,"BP_Down" = NA, "kegg_Up" = NA, "kegg_Down" = NA)
  result <- list("BP_Up" = "NA","CC_Up" = "NA","MF_Up" = "NA", "kegg_Up" = "NA")
  GOstatssm <- GOstats::summary # summary function must be from GOstats!!!
  if(do_GO_analysis){
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
  }
  # perform kegg pathway analysis.
  if(geneIDtype == "ENSEMBL"){
    uniGenes <- egENS$ENTREZID[egENS$ENSEMBL %in% uniGenes]
    selGenes <- egENS$ENTREZID[egENS$ENSEMBL %in% selGenes]
  }
  print("kegg_Up")
  kparams <- GSEAKEGGHyperGParams(name="Custom en2kegg GSEABase",
                                  geneSetCollection= kgscOrg,
                                  geneIds = selGenes,
                                  universeGeneIds = uniGenes,
                                  pvalueCutoff = pvalCutoff,
                                  testDirection = "over")
  kOver <- hyperGTest(kparams)	
  result[["kegg_Up"]] = GOstatssm(kOver)
  return(result)
}
