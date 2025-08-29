library(Seurat)
library(zellkonverter)
library(SingleCellExperiment)
library(dplyr)
# test <- readH5AD("~/Downloads/atlas_dataset/all_GSE132509.h5ad")
# print(test)
# print(colnames(colData(test)))
# table(test$Cancer.type)
# table(test$Celltype)
# class: SingleCellExperiment 
# dim: 36601 25637 
# metadata(0):
#   assays(1): X
# rownames(36601): MIR1302-2HG FAM138A ... AC007325.4 AC007325.2
# rowData names(0):
#   colnames(25637): all_GSE132509_HHD_1_AAACCTGAGCCAGAAC-1 all_GSE132509_HHD_1_AAACCTGAGTGTGGCA-1 ...
# all_GSE132509_ETV6_RUNX1_1_TTTGTCAAGGACGAAA-1 all_GSE132509_ETV6_RUNX1_1_TTTGTCAGTTCGGCAC-1
# colData names(8): Dataset Organ_origin ... cnv_status Celltype
# reducedDimNames(0):
#   mainExpName: NULL
# altExpNames(0):

# test2 <- readH5AD("~/Downloads/atlas_dataset/br_tnbc_pan_GSE148673.h5ad")
# print(test2)
# print(colnames(colData(test2)))
# table(test2$Cancer.type)
# table(test2$Celltype)
# 
# coldata <- colData(test2) %>% as.data.frame()
# 
# table(test2$cnv_status)
# 
# sc <- list.files("~/Downloads/atlas_dataset", pattern = ".h5ad$", full.names = T)
# scNames <- list.files("~/Downloads/atlas_dataset", pattern = ".h5ad$", full.names = F)

# for(file in sc){
  temp <- readH5AD("~/Downloads/atlas_dataset/skin_GSE144236.h5ad")
  coldata <- colData(temp) %>% as.data.frame()
  assay(temp, "counts") <- as.matrix(assay(temp, "X"))
  # Files have already been processed! --> log counts
  assay(temp, "logcounts") <- as.matrix(assay(temp, "X"))
  temp <- as.Seurat(temp)
  temp@meta.data <- coldata %>% as.data.frame()
  # > table(temp$Celltype)
  # B cell Dendritic cell    Endothelial     Epithelial     Fibroblast     Macrophage           Mast        NK cell    Plasma cell         T cell 
  # 141           4339            444          29802            900           9396             23            223            102           1480 
  # pDC 
  # 266 
# }

gc()

# > length(intersect(rownames(obj), rownames(temp)))
# [1] 13096

#### Orthologs ####
library(biomaRt)

# Human genes
human <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
human = useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", mirror = "useast")

# Map human gene symbols -> Ensembl IDs
gene_map_human <- getBM(
  attributes = c("external_gene_name", "ensembl_gene_id"),
  filters = "external_gene_name",
  values = rownames(temp),
  mart = human
)

# Dog genes
dog = useEnsembl(biomart = "ensembl", dataset = "clfamiliaris_gene_ensembl", mirror = "useast")
dog <- useMart("ensembl", dataset = "clfamiliaris_gene_ensembl")

gene_map_dog <- getBM(
  attributes = c("external_gene_name", "ensembl_gene_id"),
  filters = "external_gene_name",
  values = rownames(obj),
  mart = dog
)

rm(dog, human, ensembl, hub)
gc()

mapping <- read.delim("~/Downloads/mart_export.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

obj <- readRDS(paste0("~/gitClones/01_Deconvolution/03_Deconvolution_CTVT/031_Data/SCE_concat_CTVT_includingAllGenes_logCorrecDone_vstCorrecDone_ClustDone2_ManAnot.rds"))
rowdata <- rowData(obj) %>%
  as.data.frame()
rm(obj)
gc()

# Of the 15,684 unique ID's in rowdata, none are in mapping table?
# Between the Symbols, 12657 are in the mapping table. 

rowdataFiltered <- rowdata %>%
  filter(Symbol %in% mapping$Dog.gene.name)

# [1] 12447
length(intersect(gene_map_human$external_gene_name, rowdataFiltered$Symbol))

# [1] 13104
length(intersect(rowdata$Symbol, rownames(temp)))

# obj <- obj[rowdataFiltered$Symbol]

# What genes are in our sn data that aren't found in the human object? Includes many ENSCAF ID's and DLA etc. 
# > length(unique(differenceHuman$unionHumanDog))
# [1] 2575


differenceHuman <- rowdata %>%
  filter(!Symbol %in% rownames(temp)) %>%
  left_join(mapping %>% dplyr::rename("Symbol"="Dog.gene.name")) %>%
  left_join(gene_map_human %>% rename("Gene.stable.ID"="ensembl_gene_id")) %>%
  # If there is a dog Symbol but no exteral_gene_name, take the dog Symbol. 
  mutate(unionHumanDog = ifelse(is.na(Symbol) == F & is.na(external_gene_name) == T,
                                Symbol,
                                external_gene_name)
         ) %>%
  drop_na(unionHumanDog)


length(intersect(differenceHuman$unionHumanDog, rownames(temp)))

# > length(intersect(differenceHuman$external_gene_name, rownames(temp)))
# [1] 86 overlap. 
# > length(intersect(differenceHuman$unionHumanDog, rownames(temp)))
# [1] 86 overlap. 

write.xlsx(differenceHuman, 
           "~/gitClones/01_Deconvolution/03_Deconvolution_CTVT/033_Results/03_Deconvolution_CTVT_orthologs/03_Deconvolution_CTVT_orthologs_differenceHuman.xlsx")



#### Final mapping ####
# [1] 13104
# length(intersect(rowdata$Symbol, rownames(temp)))
# + 86 genes that have mappings between Dog and Human genome (based on Gene ID)

intersecting <- rowdata %>%
  filter(Symbol %in% rownames(temp)) %>%
  dplyr::select(Symbol) %>%
  rename("Dog"="Symbol") %>%
  mutate(Human = Dog,
         Type = "Intersecting")

nonintersecting <- differenceHuman %>%
  dplyr::filter(unionHumanDog %in% rownames(temp)) %>%
  rename("Dog"="Symbol",
         "Human"="external_gene_name") %>%
  mutate(Type = "Mapped") %>%
  dplyr::select(Dog, Human, Type) %>%
  distinct()

mapping <- rbind(intersecting,
                 nonintersecting)

length(unique(mapping$Dog))
length(unique(mapping$Human))
# > length(unique(mapping$Dog))
# [1] 13184
# > length(unique(mapping$Human))
# [1] 13189 --> 5 genes multiple Dog matches???
# > nrow(mapping)
# [1] 13228 --> multiple many-many mappings?


mappingEdit <- mapping %>%
  # count how many humans per dog
  left_join(mapping %>%
              count(Dog, name = "HumanCount"), by = "Dog") %>%
  # count how many dogs per human
  left_join(mapping %>%
              count(Human, name = "DogCount"), by = "Human") %>%
  mutate(Mapping = case_when(
    HumanCount == 1 & DogCount == 1 ~ "unique",
    HumanCount > 1 & DogCount == 1 ~ "humanManyDogOne",   # one dog → many humans
    HumanCount == 1 & DogCount > 1 ~ "dogManyHumanOne",   # many dogs → one human
    HumanCount > 1 & DogCount > 1 ~ "manyToMany"
  ))

# table(mappingEdit$Mapping)
# > table(mappingEdit$Mapping)
# 
# dogManyHumanOne humanManyDogOne      manyToMany          unique 
# 2               7              53           13166 


write.xlsx(mappingEdit, 
           "~/gitClones/01_Deconvolution/03_Deconvolution_CTVT/033_Results/03_Deconvolution_CTVT_orthologs/03_Deconvolution_CTVT_orthologs_mappingEdit.xlsx")

# Take only unique mappings for now!!

mappingUnique <- mappingEdit %>%
  filter(Mapping == "unique")

#### Integration ####
# Using subset counts, by gene Id's, integrate objects. 

# length(unique(mappingUnique$Dog))
# [1] 13166
# length(intersect(rownames(obj), mappingUnique$Dog))
# [1] 13166

objFiltered <- obj[mappingUnique$Dog,]
rm(obj)
head(rownames(objFiltered))

# length(unique(mappingUnique$Human))
# [1] 13166
# length(intersect(rownames(temp), mappingUnique$Human))
# [1] 13166

humanFiltered <- temp[mappingUnique$Human, ]
rm(temp)
head(rownames(humanFiltered))


objFilteredHuman <- objFiltered

mappingUnique <- mappingUnique[]


