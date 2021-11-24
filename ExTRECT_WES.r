setwd("/home/jupyter/notebooks/Melanoma_Combi_D_TK/")
dirC <- function(fv){invisible(sapply(fv,function(s) if(!dir.exists(s)){dir.create(s)}))}
fileCopy <- function(bpath){
    oPath <- gsub(".*\\/Exome.*\\/v[0-9]\\/(.*)$","\\1",bpath)
    system(glue("gsutil cp {bpath} data/",
                "{oPath}"),intern = T)
    oFile <- paste0("data/",oPath)
    return(oFile)
}  
#system("rm data/*.gstmp")                                     
#system("ls -lr data/| grep txt |wc -l")    
#ls -lr | grep txt |wc -l

setwd("/home/jupyter/notebooks/Melanoma_Combi_D_TK/")
options(quietly=T)
require(dplyr)
require(glue)
require(TcellExTRECT)
import::from(magrittr,"%<>%","set_colnames","set_rownames","set_names")
import::from(tibble,"column_to_rownames","rownames_to_column")
import::from(tidyr,"spread","unite","gather","separate")
import::from(stringr,"str_split")

SampleTable <- read.table("data/sampleTable125.tsv", sep = '\t', header = T)
SampleTable%<>%
    select(Entity = entity.sample_id,
           BamPath = clean_bam_file_capture,
           BaiPath = BAM_IDX,
           ID = squid_sample_id_capture)%>%
    filter(BamPath!="" & BaiPath!="")

SampleTable

CovWrapper <- function(s){
    bampath <- s[[2]]
    baipath <- s[[3]]
    bamFile <- fileCopy(bampath)
    baiFile <- fileCopy(baipath)
    data(TCRA_exons_hg19)
    cov.file <- getCovFromBam(bamPath = bamFile,
                          outPath = 'data/',
                          vdj.seg = tcra_seg_hg19)
    file.remove(c(bamFile,baiFile)) 
}

for (i in 1:nrow(SampleTable)){
  s <-   SampleTable[i,]%>%as.list()
  txtOut <- paste0("data/",gsub("\\s","_",s[[4]]),"_TCRA.txt")
  if(!file.exists(txtOut)){
      message("Processing Sample",s[[4]],"...")
      CovWrapper(s)
  } 
}

data(TCRA_exons_hg19)
estD <- lapply(list.files("data/",".txt",full.names = T),function(s){
    cov_df <- loadCov(s)
    TCRA.out <- runTcellExTRECT(cov_df, TCRA_exons_hg19, tcra_seg_hg19, 'hg19',sample_name=s)
    return(TCRA.out)
})

estD[!is.na(estD)]%>%bind_rows()%>%openxlsx::write.xlsx(.,paste0("data/EstTCF_",Sys.Date(),".xlsx"),overwrite=T)

estD[!is.na(estD)]%>%bind_rows()%>%dim()

bucket <- Sys.getenv('WORKSPACE_BUCKET')
system(glue('gsutil cp -r data/EstTCF_2021-11-04.xlsx {bucket}/'), intern = T)

options(repr.plot.width=12, repr.plot.height=6)
plotTcellExTRECT(cov_df, TCRA_exons_hg19,
                 tcra_seg_hg19,'hg19', sample_name = 'TEST')

fullSampleBucket <- "gs://fc-ada443cd-d676-4953-8be6-71e755301ad1/20192015Cohort.xlsx"
system(glue('gsutil cp -r {fullSampleBucket} data/'), intern = T)

subsetPat <- openxlsx::read.xlsx("data/20192015Cohort.xlsx", sep = '\t')

head(subsetPat)

CovWrapper2 <- function(s){
    bampath <- s[[1]]
    baipath <- s[[2]]
    bamFile <- fileCopy(bampath)
    baiFile <- fileCopy(baipath)
    data("tcra_seg_hg19")
    cov.file <- getCovFromBam(bamPath = bamFile,
                          outPath = 'data/',
                          vdj.seg = tcra_seg_hg19)
    file.remove(c(bamFile,baiFile)) 
}

for (i in 1:nrow(subsetPat)){
  s <-   subsetPat[i,]%>%as.list()
  txtOut <- paste0("data/",gsub("\\s","_",s[[3]]),"_TCRA.txt")
  if(!file.exists(txtOut)){
      message("Processing Sample",s[[4]],"...")
      CovWrapper2(s)
  } 
}

allOutput <- list.files("data",pattern = "TCRA.txt",full.names = T)
CombinedOut <- pbmcapply::pbmclapply(allOutput,function(t) read.delim(t,header = F)%>%
    set_colnames(c("chr",'pos','reads'))%>%
    mutate(sample = gsub("data\\/(.*)\\_TCRA.txt","\\1",t)),mc.cores=15)

readr::write_delim(CombinedOut,"data/combinedTCRA_Coverage.txt")

bucket <- Sys.getenv('WORKSPACE_BUCKET')
system(glue('gsutil cp -r data/combinedTCRA_Coverage.txt {bucket}/'), intern = T)

system(glue('gsutil cp -r {bucket}/Agilent_WES_hg19.rds data'), intern = T)

allOutput <- list.files("data",pattern = "TCRA.txt",full.names = T)
AgilengExon <- readRDS("data/Agilent_WES_hg19.rds")
data("tcra_seg_hg19")

CombinedTCF_thresD <- list()
N <- length(allOutput)
for (i in 1:N){
    message("Estimating for ",i," ...")
    CombinedTCF_thresD[[i]] <- read.delim(allOutput[i],header = F)[,2:3]%>%
        set_colnames(c('pos','reads'))%>%
        runTcellExTRECT(., AgilengExon, tcra_seg_hg19, 'hg19',
                       sample_name = gsub("data\\/(.*)\\_TCRA.txt","\\1",allOutput[i])) 
}

CombinedTCF_thresD%<>%.[!is.na(.)]%>% bind_rows()
saveRDS(CombinedTCF_thresD,"data/CombinedTCF_thresD.rds")

bucket <- Sys.getenv('WORKSPACE_BUCKET')
system(glue('gsutil cp -r data/CombinedTCF_thresD.rds {bucket}/'), intern = T)

bucket <- Sys.getenv('WORKSPACE_BUCKET')
system(glue('gsutil cp -r {bucket}/adj14_CNfetch.rds  data/'), intern = T)

adj14_CN <- readRDS("data/adj14_CNfetch.rds")

absoluteCN <- lapply(adj14_CN$cn,function(s){
  system(glue('gsutil cp -r {s} data/'),intern=T)
    f <- gsub(".*/","",s)
    tcraCN <- data.table::fread(glue("data/{f}"))%>%
       .[,c(1,2,3,c(10,11,14,15))]%>%
        filter(Chromosome==14 &Start.bp <22090057)
  return(tcraCN)})%>%
    bind_rows()

saveRDS(absoluteCN,"data/absoluteCN.rds")

system(glue('gsutil cp -r data/absoluteCN.rds {bucket}/'), intern = T)