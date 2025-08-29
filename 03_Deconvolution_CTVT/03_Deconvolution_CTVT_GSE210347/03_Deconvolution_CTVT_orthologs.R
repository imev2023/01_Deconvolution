library(Seurat)
library(zellkonverter)
library(SingleCellExperiment)
library(dplyr)

temp <- readRDS(paste0("~/gitClones/01_Deconvolution/03_Deconvolution_CTVT/031_Data/GSE210347/GSE210347_counts.Rds"))

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

obj <- readRDS(paste0("~/gitClones/01_Deconvolution/03_Deconvolution_CTVT/031_Data/SCE_concat_CTVT_includingAllGenes_logCorrecDone_vstCorrecDone_ClustDone2_ManAnot.rds"))

gene_map_dog <- getBM(
  attributes = c("external_gene_name", "ensembl_gene_id"),
  filters = "external_gene_name",
  values = rownames(obj),
  mart = dog
)

rm(dog, human, ensembl, hub)
gc()

##############################################
# ## File generated from ensembl BioMart webpage
# - https://www.ensembl.org/info/data/biomart/index.html. navigate to the Data Mining Tool. 
# - Choose "Ensembl Genes 114" and "Human genes (GRCh38.p14)"
# - On left-hand side click on "Attributes" and select "Homologues". Open "Gene" and also select "Gene name" (in addition to the four default columns)
# - Open "ORTHOLOGUES A-E" and select "Dog Orthologues" columns "Dog gene stable ID", "Dog gene name", and "Dog protein or transcript stable id".
# - Export as XML file. 
# ** Of note, the dog genes come from ROS Cfam 1.0. This will need to be updated to match CANFAM soon **
##############################################

mapping <- read.delim("~/Downloads/mart_export.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

rowdata <- rowData(obj) %>%
  as.data.frame()

rm(obj)
gc()

# Between the Symbols, 12657 are in the mapping table. 
rowdataFiltered <- rowdata %>%
  filter(Symbol %in% mapping$Dog.gene.name)

# [1] 12412
length(intersect(gene_map_human$external_gene_name, rowdataFiltered$Symbol))

# [1] 13064
length(intersect(rowdata$Symbol, rownames(temp)))

##############################################
### Initial mapping ####


# What genes are in our snRNAseq data that aren't found in the human object? Includes many ENSCAF ID's and DLA etc. 
# > length(unique(differenceHuman$unionHumanDog))
# [1] 2575

differenceHuman <- rowdata %>%
  filter(!Symbol %in% rownames(temp)) %>%
  left_join(mapping %>% dplyr::rename("Symbol"="Dog.gene.name")) %>%
  left_join(gene_map_human %>% rename("ensembl_gene_id"="Gene.stable.ID")) %>%
  # If there is a dog Symbol but no exteral_gene_name, take the dog Symbol. 
  mutate(unionHumanDog = ifelse(is.na(Symbol) == F & is.na(external_gene_name) == T,
                                Symbol,
                                external_gene_name)
         ) %>%
  drop_na(unionHumanDog)

length(intersect(differenceHuman$unionHumanDog, rownames(temp)))

# > length(intersect(differenceHuman$external_gene_name, rownames(temp)))
# [1] 81 overlap. 
# > length(intersect(differenceHuman$unionHumanDog, rownames(temp)))
# [1] 81 overlap. 

write.xlsx(differenceHuman, 
           "~/gitClones/01_Deconvolution/03_Deconvolution_CTVT/033_Results/03_Deconvolution_CTVT_GSE210347/03_Deconvolution_CTVT_orthologs_differenceHuman.xlsx")

##############################################


##############################################
#### Final mapping ####

# [1] 13104
# length(intersect(rowdata$Symbol, rownames(temp)))
# + 81 genes that have mappings between Dog and Human genome (based on Gene ID)


### Intersecting = genes that have matching Symbols in dogs and humans. 
intersecting <- rowdata %>%
  filter(Symbol %in% rownames(temp)) %>%
  dplyr::select(Symbol) %>%
  rename("Symbol"="Dog") %>%
  mutate(Human = Dog,
         Type = "Intersecting")

### Non-Intersecting = genes that have mapped gene id's but not matching symbols. (ex. HLA- vs DLA-)
nonintersecting <- differenceHuman %>%
  dplyr::filter(unionHumanDog %in% rownames(temp)) %>%
  rename("Symbol"="Dog",
         "external_gene_name"="Human") %>%
  mutate(Type = "Mapped") %>%
  dplyr::select(Dog, Human, Type) %>%
  distinct()

mapping <- rbind(intersecting,
                 nonintersecting)

length(unique(mapping$Dog))
length(unique(mapping$Human))

mappingEdit <- mapping %>%
  # count how many humans per dog
  left_join(mapping %>%
              dplyr::count(Dog, name = "HumanCount"), by = "Dog") %>%
  # count how many dogs per human
  left_join(mapping %>%
              dplyr::count(Human, name = "DogCount"), by = "Human") %>%
  mutate(Mapping = case_when(
    HumanCount == 1 & DogCount == 1 ~ "unique",
    HumanCount > 1 & DogCount == 1 ~ "humanManyDogOne",   # one dog --> many human genes
    HumanCount == 1 & DogCount > 1 ~ "dogManyHumanOne",   # many dog genes --> one human
    HumanCount > 1 & DogCount > 1 ~ "manyToMany"
  ))

# table(mappingEdit$Mapping)
# > table(mappingEdit$Mapping)
# 
# dogManyHumanOne humanManyDogOne      manyToMany          unique 
# 2               7              53           13166 


write.xlsx(mappingEdit, 
           "~/gitClones/01_Deconvolution/03_Deconvolution_CTVT/033_Results/03_Deconvolution_CTVT_GSE210347/03_Deconvolution_CTVT_orthologs_mappingEdit.xlsx")
##############################################

# Taking only unique mappings for integration (for now)!!

mappingUnique <- mappingEdit %>%
  filter(Mapping == "unique")

#### Integration ####
# Using subset counts, by gene Id's, integrate objects. 

# length(unique(mappingUnique$Dog))
# [1] 13166
# length(intersect(rownames(obj), mappingUnique$Dog))
# [1] 13166



