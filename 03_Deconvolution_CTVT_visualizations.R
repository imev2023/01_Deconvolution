library(SingleCellExperiment)
library(Seurat)
library(scater)
library(LoomExperiment)
# devtools::install_github('immunogenomics/presto')
library(presto)
# install.packages("UpSetR")
library(UpSetR)

fprefix <- "03_Deconvolution_CTVT_geneCorrelations"
wDir <- paste0("~/gitClones/01_Deconvolution/03_Deconvolution_CTVT/")
dataDir <- paste0(wDir, "/031_Data/")
resDir <- paste0(wDir, "/033_Results/", fprefix, "/")
plotDir <- paste0(wDir, "/034_Plots/", fprefix, "/")

obj <- readRDS(paste0(dataDir, "SCE_concat_CTVT_includingAllGenes_logCorrecDone_vstCorrecDone_ClustDone2_ManAnot.rds"))
obj$manual.log.anot.fine

pdf(paste0("~/gitClones/01_Deconvolution/03_Deconvolution_CTVT/034_Plots/violin/03_Deconvolution_CTVT_IGKC_expression.pdf"), 
    height=4, 
    width=6)
plotExpression(obj, 
               "IGKC", 
               x = "manual.log.anot.fine", 
               color_by = "manual.log.anot.fine")
dev.off()

goi <- "UBE2S"

pdf(paste0("~/gitClones/01_Deconvolution/03_Deconvolution_CTVT/034_Plots/violin/03_Deconvolution_CTVT_", goi,"_expression.pdf"), 
    height=4, 
    width=6)
plotExpression(obj, 
               goi, 
               x = "manual.log.anot.fine", 
               color_by = "manual.log.anot.fine")
dev.off()

plotReducedDim(obj, dimred = "TSNE_log_corrected",
                by_exprs_values = "logcounts",
               color_by = goi,
                order_by=goi) +
  scale_color_gradient(low = "grey90", high = "purple")


gois <- hostSpecDfStrict$Symbol
length(set)
intersecting <- intersect(gois, rownames(obj)) # 17 not in rownames

for(g in intersecting){
  
  vio <- plotExpression(obj, 
                        g, 
                        x = "manual.log.anot.fine", 
                        color_by = "manual.log.anot.fine") + 
    ggtitle(paste0("hostSpecDfStrict - ", g))
  
  pdf(paste0("~/gitClones/01_Deconvolution/03_Deconvolution_CTVT/034_Plots/", fprefix, "/violin/03_Deconvolution_CTVT_hostSpecDfStrict_", g,"_expression.pdf"), 
      height=4, 
      width=6)
  print(vio)
  
  dev.off()
  
  feat <- plotReducedDim(obj, dimred = "TSNE_log_corrected",
                 by_exprs_values = "logcounts",
                 color_by = g,
                 order_by=g) +
    scale_color_gradient(low = "grey90", high = "purple") + 
    ggtitle(paste0("hostSpecDfStrict - ", g))
  
  pdf(paste0("~/gitClones/01_Deconvolution/03_Deconvolution_CTVT/034_Plots/", fprefix, "/feature/03_Deconvolution_CTVT_hostSpecDfStrict_", g,"_expression.pdf"), 
      height=4, 
      width=4)
  
  print(feat)
  
  
  dev.off()
  
}


gois <- tumorSpecDfStrict$Symbol
length(setdiff(gois, rownames(obj)))
intersecting <- intersect(gois, rownames(obj)) # 17 not in rownames

for(g in intersecting){
  
  vio <- plotExpression(obj, 
                        g, 
                        x = "manual.log.anot.fine", 
                        color_by = "manual.log.anot.fine") + 
    ggtitle(paste0("tumorSpecDfStrict - ", g))
  
  pdf(paste0("~/gitClones/01_Deconvolution/03_Deconvolution_CTVT/034_Plots/", fprefix, "/violin/03_Deconvolution_CTVT_tumorSpecDfStrict_", g,"_expression.pdf"), 
      height=4, 
      width=6)
  print(vio)
  
  dev.off()
  
  feat <- plotReducedDim(obj, dimred = "TSNE_log_corrected",
                         by_exprs_values = "logcounts",
                         color_by = g,
                         order_by=g) +
    scale_color_gradient(low = "grey90", high = "purple") + 
    ggtitle(paste0("tumorSpecDfStrict - ", g))
  
  pdf(paste0("~/gitClones/01_Deconvolution/03_Deconvolution_CTVT/034_Plots/", fprefix, "/feature/03_Deconvolution_CTVT_tumorSpecDfStrict_", g,"_expression.pdf"), 
      height=4, 
      width=4)
  
  print(feat)
  
  
  dev.off()
  
}

#### 2. Investigating cell type markers ####
##### 2.1 Seurat
###### 2.1.1 Convert
assay(obj, "counts") <- as.matrix(assay(obj, "counts"))
assay(obj, "logcounts") <- as.matrix(assay(obj, "logcounts"))
seuratObj <- as.Seurat(obj, counts = "counts", data = "logcounts")
seuratObj
Idents(seuratObj) <- "manual.log.anot.fine"
seuratObjNoBad <- subset(seuratObj, idents = "bad", invert=TRUE)
table(Idents(seuratObjNoBad))
DimPlot(seuratObj, reduction = "TSNE_log_corrected")
DimPlot(seuratObjNoBad, reduction = "TSNE_log_corrected")

rm(obj, seuratObj)
gc()

###### 2.1.2 Find all markers
famWilcox <- FindAllMarkers(seuratObjNoBad, 
                      group.by = "manual.log.anot.fine", 
                      min.pct = 0.05)

write.xlsx(famWilcox, 
           paste0(resDir, "findAllMarkers/", fprefix, "_seuratObjNoBad_famWilcox.xlsx"))

gc()

famMast <- FindAllMarkers(seuratObjNoBad, 
                          test.use = "MAST",
                          group.by = "manual.log.anot.fine", 
                          min.pct = 0.05)
write.xlsx(famMast, 
           paste0(resDir, "findAllMarkers/", fprefix, "_seuratObjNoBad_famMast.xlsx"))

# famMastLv <- FindAllMarkers(seuratObjNoBad, 
#                           test.use = "MAST",
#                           latent.vars = "sample",
#                           group.by = "manual.log.anot.fine", 
#                           min.pct = 0.05)
# write.xlsx(famMast, 
#            paste0(resDir, "findAllMarkers/", fprefix, "_seuratObjNoBad_famMast.xlsx"))



###### 2.1.3 Filter CTVT FAM DEG's first. --> check against genetic deconv list. Test positive genes first. 
famWilcox <- read.xlsx(paste0(resDir, "findAllMarkers/", fprefix, "_seuratObjNoBad_famWilcox.xlsx"))
famWilcoxDegs <- famWilcox %>%
  filter(p_val_adj < 0.01 & avg_log2FC > 2.5 & cluster == "CTVT") 

intersecting <- intersect(famWilcoxDegs$gene, tumorSpecSampleGeneDf$Symbol) ### 33 genes intersecting
intersectingStrict <- intersect(famWilcoxDegs$gene, tumorSpecDfStrict$Symbol) ### None overlap. 

for(g in intersecting){

  feat <- FeaturePlot(seuratObjNoBad,
                      reduction = "TSNE_log_corrected",
                      features=g) +

    ggtitle(paste0("tumorSpecDf + Wilcox DEG - ", g))
  
  featSamp <- FeaturePlot(seuratObjNoBad, 
                          reduction = "TSNE_log_corrected",
                          features=g, 
                          split.by = "Sample") 
  
  pdf(paste0("~/gitClones/01_Deconvolution/03_Deconvolution_CTVT/034_Plots/", fprefix, "/featureSeurat/tumorSpecificIntersect/03_Deconvolution_CTVT_tumorSpecDf_famWilcoxDegs_", g,"_expression.pdf"),
      height=4,
      width=4)

  print(feat)

  dev.off()
  
  
  
  pdf(paste0("~/gitClones/01_Deconvolution/03_Deconvolution_CTVT/034_Plots/", fprefix,  "/featureSeurat/tumorSpecificIntersect/03_Deconvolution_CTVT_tumorSpecDf_famWilcoxDegs_", g,"_expressionBySample.pdf"), 
      height=4, 
      width=14)
  
  print(featSamp) 
  
  dev.off()
  
}

##### Is this preseved accross 5 samples? and is this preserved in the bulk?

###### 2.1.4 Loop through NON CTVT FAM DEG's --> check against genetic deconv list. Test positive genes first. 
for(ctype in unique(famWilcox$cluster)){
  if(ctype=="CTVT"){
    next
  }
  
  tempFamWilcoxDegs <- famWilcox %>%
    filter(p_val_adj < 0.01 & avg_log2FC > 2.5 & cluster == ctype) 
  
  tempIntersecting <- intersect(tempFamWilcoxDegs$gene, hostSpecSampleGeneDf$Symbol)
  
  print(paste0("Cell type: ", ctype, " - Intersecting: ", length(tempIntersecting)))
  
  tempPlotDir <- paste0("hostSpecificIntersect", ctype)
  
  for(tg in tempIntersecting){
    
    feat <- FeaturePlot(seuratObjNoBad,
                        reduction = "TSNE_log_corrected",
                        features=tg) +
      
      ggtitle(paste0("hostSpecDf + ", ctype, " Wilcox DEG - ", tg))
    
    featSamp <- FeaturePlot(seuratObjNoBad, 
                            reduction = "TSNE_log_corrected",
                            features=tg, 
                            split.by = "Sample") 
    
    pdf(paste0("~/gitClones/01_Deconvolution/03_Deconvolution_CTVT/034_Plots/", fprefix, "/featureSeurat/", tempPlotDir, "/03_Deconvolution_CTVT_tumorSpecDf_famWilcoxDegs_", tg,"_expression.pdf"),
        height=4,
        width=4)
    
    print(feat)
    
    dev.off()
    
    
    
    pdf(paste0("~/gitClones/01_Deconvolution/03_Deconvolution_CTVT/034_Plots/", fprefix, "/featureSeurat/", tempPlotDir, "/03_Deconvolution_CTVT_tumorSpecDf_famWilcoxDegs_", tg,"_expressionBySample.pdf"), 
        height=4, 
        width=14)
    
    print(featSamp) 
    
    dev.off()
    
  }
  
}

#### Filter top Degs + top purity correlation + top eT/eH

mergedFamWilcoxDegsCTVT <- famWilcox %>%
  filter(p_val_adj < 0.01 & avg_log2FC > 2.5 & cluster=="CTVT") %>%
  filter(gene %in% tumorSpecSampleGeneDf$Symbol) %>%
  rename("Symbol"="gene") %>%
  left_join(tumorSpecSampleGeneDf) %>%
  mutate("Annotation"=paste0("tumorSpec_", cluster, "Deg")) 


write.xlsx(mergedFamWilcoxDegsCTVT, 
           paste0(resDir, fprefix, "_mergedFamWilcoxDegs_adrianTumorGenes.xlsx"))

mergedFamWilcoxDegsAdrian <- mergedFamWilcoxDegsCTVT

ctvtSet <- mergedFamWilcoxDegsCTVT $Symbol

# # 5 genes are positively correlated with tumor purity and a CTVT DEG
# candTumor <- mergedFamWilcoxDegsCTVT %>% filter(purity_correlation > 0.5) %>% select(Symbol, medianEt, purity_correlation)


mergedFamWilcoxDegsAdrian <- data.frame()


for(ctype in unique(famWilcox$cluster)){
  if(ctype=="CTVT"){
    next
  }
  
  temp <- famWilcox %>%
    filter(p_val_adj < 0.01 & avg_log2FC > 2.5& cluster==ctype) %>%
    filter(gene %in% hostSpecSampleGeneDf$Symbol) %>%
    rename("Symbol"="gene") %>%
    left_join(hostSpecSampleGeneDf) %>%
    mutate("Annotation"=paste0("hostSpec_", ctype, "Deg"))
  
  mergedFamWilcoxDegsAdrian <- rbind(mergedFamWilcoxDegsAdrian, temp)
  
}

write.xlsx(mergedFamWilcoxDegsAdrian, 
           paste0(resDir, fprefix, "_mergedFamWilcoxDegs_adrianHostGenes.xlsx"))


# 17 genes are negative correlated with tumor purity and a Bcell DEG
print(nrow(mergedFamWilcoxDegsAdrian %>% filter(cluster == "Bcells" & purity_correlation < -0.5)))
candBcell <- mergedFamWilcoxDegsAdrian %>% 
  filter(cluster == "Bcells" & purity_correlation < -0.5) %>% 
  select(Symbol, medianEh, purity_correlation, avg_log2FC)

bcellSet <- candBcell$Symbol

# 36 genes are negatively correlated with tumor purity and a Myeloid DEG
print(nrow(mergedFamWilcoxDegsAdrian %>% filter(cluster == "Myeloid" & purity_correlation < -0.5)))
candMyeloid <- mergedFamWilcoxDegsAdrian %>% 
  filter(cluster == "Myeloid" & purity_correlation < -0.5) %>% 
  select(Symbol, medianEh, purity_correlation, avg_log2FC)

myeloidSet <- candMyeloid$Symbol

# 32 genes are negatively correlated with tumor purity and a Myeloid DEG
print(nrow(mergedFamWilcoxDegsAdrian %>% filter(cluster == "CAF.Endo" & purity_correlation < -0.5)))
candCAFEndo <- mergedFamWilcoxDegsAdrian %>% 
  filter(cluster == "CAF.Endo" & purity_correlation < -0.5) %>% 
  select(Symbol, medianEh, purity_correlation, avg_log2FC)

cafendoSet <- candCAFEndo$Symbol


# 21 genes are negatively correlated with tumor purity and a Plasma DEG
print(nrow(mergedFamWilcoxDegsAdrian %>% filter(cluster == "Plasma" & purity_correlation < -0.5)))
candPlasma <- mergedFamWilcoxDegsAdrian %>% 
  filter(cluster == "Plasma" & purity_correlation < -0.5) %>% 
  select(Symbol, medianEh, purity_correlation, avg_log2FC)

plasmaSet <- candPlasma$Symbol

# 9 genes are negatively correlated with tumor purity and a Tcells DEG
print(nrow(mergedFamWilcoxDegsAdrian %>% filter(cluster == "Tcells" & purity_correlation < -0.5)))
candTcells <- mergedFamWilcoxDegsAdrian %>% 
  filter(cluster == "Tcells" & purity_correlation < -0.5) %>% 
  select(Symbol, medianEh, purity_correlation, avg_log2FC)

tcellSet <- candTcells$Symbol

geneSets <- list("ctvt"=ctvtSet,
             "bcell"=bcellSet,
             "myeloid"=candMyeloid$Symbol,
             "caf.endo"=cafendoSet,
             "plasma"=plasmaSet,
             "tcell"=tcellSet)


genes <- sort(unique(unlist(geneSets)))


sets <- t(+sapply(geneSets, "%in%", x = genes)) %>%
  t() %>%
  as.data.frame()

rownames(sets) <- genes

sets <- sets %>%
  rownames_to_column("Gene")


#         Gene1 Gene2 Gene3
#pathwayX     0     0     1
#pathwayY     0     1     1
#pathwayz     1     1     1

data.frame(mat)  ## data.frame output

