#===== Library  loading and set-up =====
pacman::p_load(glue,TcellExTRECT,mgsub,
               ggplot2,RColorBrewer,ggpubr,dplyr)
import::from(magrittr,"%<>%","set_colnames","set_rownames","set_names")
import::from(tibble,"column_to_rownames","rownames_to_column")
import::from(tidyr,"spread","unite","gather","separate")
setwd("~/ZoneIn/2021-11-24_BioPhy_Draft/Progress/")
adjTCF <- function(TCRA.out, trustPurity = TRUE){
  # solve visible binding issue
  TCRA.tcell.fraction <- TCRA.tcell.fraction.lwr <- TCRA.tcell.fraction.upr <- NULL
  rawRatio <- rawRatio.lwr <- rawRatio.upr <- maxPossible <- highTcellFlag <- NULL
  TCRA.tcell.fraction.adj <- TCRA.tcell.fraction.adj.lwr <-TCRA.tcell.fraction.adj.upr <- NULL
  
  TCRA.out <- TCRA.out %>%
    dplyr::mutate(rawRatio = 1-TCRA.tcell.fraction) %>%
    dplyr::mutate(rawRatio.lwr = 1-TCRA.tcell.fraction.lwr) %>%
    dplyr::mutate(rawRatio.upr = 1-TCRA.tcell.fraction.upr) %>%
    dplyr::mutate(TCRA.tcell.fraction.adj= 1 - ((1-purity+(purity*TCRAcn)/2)*rawRatio) - purity + ((purity*TCRAcn)/2)) %>%
    dplyr::mutate(TCRA.tcell.fraction.adj.lwr= 1 - ((1-purity+(purity*TCRAcn)/2)*rawRatio.lwr) - purity + ((purity*TCRAcn)/2)) %>%
    dplyr::mutate(TCRA.tcell.fraction.adj.upr= 1 - ((1-purity+(purity*TCRAcn)/2)*rawRatio.upr) - purity + ((purity*TCRAcn)/2)) %>%
    dplyr::mutate(maxPossible = 1- purity) %>%
    dplyr::mutate(highTcellFlag = TCRA.tcell.fraction.adj.lwr > maxPossible)
  
  if(trustPurity){
    TCRA.out <- TCRA.out %>%
      dplyr::mutate(TCRA.tcell.fraction.adj = ifelse(highTcellFlag, TCRA.tcell.fraction,TCRA.tcell.fraction.adj)) %>%
      dplyr::mutate(TCRA.tcell.fraction.adj.lwr = ifelse(highTcellFlag, TCRA.tcell.fraction.lwr,TCRA.tcell.fraction.adj.lwr)) %>%
      dplyr::mutate(TCRA.tcell.fraction.adj.upr = ifelse(highTcellFlag, TCRA.tcell.fraction.upr, TCRA.tcell.fraction.adj.upr))
  }
  return(TCRA.out)
  
}
#===== Customized WES capture kit =====
Agilent_WES<- "Data/ExTRECT/whole_exome_agilent_1.1_refseq_plus_3_boosters.Homo_sapiens_assembly19.targets.bed"
TCRA_exons_hg19_custom <- createExonDFBed(Agilent_WES, 'hg19')
saveRDS(TCRA_exons_hg19_custom,"Agilent_WES_hg19.rds")
#===== Matrix Input for Timer2 =====
tpmMatrix <- data.table::fread("~/Downloads/TPM_matrix_rows.txt")%>%t()
colnames(tpmMatrix) <- tpmMatrix[1,]
tpmMatrix%<>%.[-1,]
sampleInput <- data.table::fread("exampleForLUAD.csv")
Matrix2Timer <- tpmMatrix[intersect(rownames(tpmMatrix),sampleInput$V1),]%>%
  as.data.frame()%>%
  mutate_all(as.numeric)
RNAseqcohorts <- openxlsx::read.xlsx("20192015Cohort.xlsx")
FullSampleTable <- read.table("FullSample_Melanoma.tsv", sep = '\t', header = T)
FullSampleTable%<>%
  select(Entity = entity.sample_id,
         BamPath = clean_bam_file_capture,
         BaiPath = BAM_IDX,
         ID = squid_sample_id_capture)%>%
  filter(BamPath!="" & BaiPath!="")%>%
  filter(grepl("MEL-Patient|MEL-IPI",Entity))%>%
  filter(ID %in% RNAseqcohorts$ID)
rnaseq_Matrix <- pbmcapply::pbmclapply(FullSampleTable$Entity,function(s) grep(s,colnames(Matrix2Timer),value=T),mc.cores = 7)%>%
  .[lengths(.)==1]%>%unlist()
Matrix2Timer%<>%.[,rnaseq_Matrix]
colnames(Matrix2Timer)%<>%gsub("\\_1","",.)
write.csv(Matrix2Timer,file = "TimerExprRNAseq.csv") # File uploaded to Timer2 server: http://timer.cistrome.org/
#===== Raw TCF =====
rawTCF <- readRDS("Data/ExTRECT/CombinedTCF_thresD.rds")%>%
  mutate(sample = case_when(grepl("BRAF|[Pp]re",sample)~mgsub(sample,c("Pat ","Pat_"),c("Pat","Pat")),TRUE~sample))%>%
  mutate(sample = case_when(grepl("BRAF|[Pp]re",sample)~gsub("_",".",sample),TRUE~sample))%>%
  mutate(sample = gsub(".BRAF.Pre",".pre",sample))
#===== Purity and TCRA cn Collection =====
adj14 <- read.delim("~/Downloads/pair.tsv")%>%
  select(Entity = entity.pair_id,
         cn=absolute_seg_file,
         purity = purity)
saveRDS(adj14,"Data/ExTRECT/adj14_CNfetch.rds")
adjCN_14_abs <- readRDS("Data/ExTRECT/absoluteCN.rds")
adj15 <- openxlsx::read.xlsx("Data/SampleInfo/2015_IO_pairs_census.xlsx")%>%
  select(TumorEntity = case_sample,
         purity = facets_purity)
adj19 <- openxlsx::read.xlsx("Data/SampleInfo/ClinicalChar_2019.xlsx")[-1,c(1,45)]%>%
  set_colnames(c("PatientID","purity"))%>%
  filter(grepl("Patient",PatientID))%>%
  mutate(purity =round(as.numeric(purity),2))
adjCN_1519_facets <- pbmcapply::pbmclapply(list.files(".","cncf.txt",recursive = T),
       function(cnf)
        read.delim(cnf)%>%
         filter(chrom==14 & start <22090057 & end >22090057)%>%
         select(chrom,tcn.em,lcn.em,start,end)%>%
         mutate(sample = gsub(".*\\/(.*)\\.facets_cncf.txt","\\1",cnf)),mc.cores = 7)%>%
  bind_rows()
adjPurityCN_19 <- adjCN_1519_facets%>%
  filter(grepl("Patient",sample))%>%
  mutate(PatientID = gsub("MEL-(.*)\\-T[PM].*","\\1",sample))%>%
  select(PatientID,TCRAcn = tcn.em)%>%
  inner_join(adj19,by="PatientID")%>%
  mutate(PubYr = "Y2019")%>%
  distinct()
adjPurityCN_15 <- adjCN_1519_facets%>%
  filter(grepl("IPI",sample))%>%
  mutate(sampleV1 = gsub("TP-N[BT]","Tumor",sample))%>%
  mutate(sampleV2 = gsub("\\-SM\\-[1-z0-9]{5}$","",sampleV1))%>%
  select(TumorEntity = sampleV2,
         TCRAcn = tcn.em)%>%
  inner_join(adj15,by="TumorEntity")%>%
  mutate(PubYr = "Y2015")%>%
  distinct()
adjPurityCN_14 <- adjCN_14_abs%>%
  filter(rescaled_total_cn>0)%>%
  select(Entity=sample,TCRAcn = rescaled_total_cn)%>%
  inner_join(adj14%>%select(-cn),by="Entity")%>%
  mutate(PubYr = "Y2014")%>%
  distinct()%>%
  group_by(Entity,purity,PubYr)%>%
  summarise(TCRAcn = mean(TCRAcn))%>%
  as.data.frame()%>%
  mutate(Entity =mgsub(Entity,
                         c("MEL-Pat_40-TM-NB","MEL-Pat_40-TR-NB","MEL-Pat_41-TM-NB","MEL-Pat_41-TR-NB"),
                         c("MEL-Pat_40-Tumor-SM-38ZZ6","MEL-Pat_40-Tumor-SM-39113",
                           "MEL-Pat_41-Tumor-SM-38ZZ7","MEL-Pat_41-Tumor-SM-39114")))%>%
  mutate(sampleV1 = gsub("T[MR]-N[BT]","Tumor",Entity))%>%
  mutate(TumorEntity = case_when(grepl("-Pat_4[01]-",sampleV1)~sampleV1,
                                 TRUE~gsub("\\-SM\\-[1-z0-9]{5}$","",sampleV1)))%>%
  select(TumorEntity,TCRAcn,purity,PubYr)%>%
  mutate(TumorEntity = gsub("mel","MEL",TumorEntity))
#===== TCF Adjustment =====
FullSampleTable <- read.table("Data/SampleInfo/FullSample_Melanoma.tsv", sep = '\t', header = T)%>%
  filter(clean_bam_file_capture!="" & BAM_IDX!="")%>%
  select(Entity = entity.sample_id,
         ID = squid_sample_id_capture)%>%
  mutate(Entity = gsub("mel","MEL",Entity))%>%
  filter(grepl("MEL-Patient|MEL-IPI|MEL-Pat\\_",Entity,ignore.case = T))%>%
  mutate(ID = mgsub(ID,c("_re$"," pre"," BRAF","Pat ","Pat_"),c("",".pre",".BRAF","Pat","Pat"),ignore.case=T))
TCFsubSample <- FullSampleTable%>%
  mutate(PubYr = case_when(grepl("MEL-Patient",Entity)~"Y2019",
                           grepl("MEL-IPI",Entity)~"Y2015",
                           grepl("MEL-Pat_",Entity)~"Y2014"))%>%
  mutate(ID = case_when(PubYr=="Y2014"~(gsub("_",".",ID)%>%gsub(".BRAF.Pre",".pre",.)),TRUE~ID))%>%
  inner_join(rawTCF,by=c("ID"="sample"))%>%
  mutate(PatientID = case_when(PubYr=="Y2014"~gsub("\\..*","",ID),
                               PubYr=="Y2015"~gsub(".*\\_","",ID),
                               PubYr=="Y2019"~gsub("\\_.*","",ID)))
Combined_ExTRECT <- TCFsubSample%>%
  inner_join(rbind(adjPurityCN_14,adjPurityCN_15),by=c("Entity"="TumorEntity","PubYr"))%>%
  replace(is.na(.),0.1)%>%
  adjTCF()%>%
  select(PubYr,TumorEntity = Entity,ID,PatientID,
         RAW_TCF = TCRA.tcell.fraction,ADJ_TCF=TCRA.tcell.fraction.adj,
         TCRAcn,Purity = purity,maxPossible)%>%
  rbind(TCFsubSample%>%
          filter(grepl("Tumor",Entity) & PubYr=="Y2019")%>%
          inner_join(adjPurityCN_19,by=c("PatientID","PubYr"))%>%
          adjTCF()%>%
          select(PubYr,TumorEntity = Entity,ID,PatientID,
                 RAW_TCF = TCRA.tcell.fraction,ADJ_TCF=TCRA.tcell.fraction.adj,
                 TCRAcn,Purity = purity,maxPossible))%>%
  mutate(RAW_TCF = ifelse(RAW_TCF<0.01,0.01,round(RAW_TCF,2)))%>%
  mutate(ADJ_TCF = ifelse(ADJ_TCF<0.01,0.01,round(ADJ_TCF,2)))
#===== Clinical Characteristics (3 Pubs) =====
Clin14 <- openxlsx::read.xlsx("Data/SampleInfo/ClinicalChar_2014.xlsx")[-c(1),1:8]%>%
  set_colnames(c("PatientID","Gender","StartAge","Medication",
                 "EarlyResistance","TherapyDuration",
                 "Response","NormalTissue"))
Clin15 <- openxlsx::read.xlsx("Data/SampleInfo/ClinicalChar_2015.xlsx")
Clin19 <- openxlsx::read.xlsx("Data/SampleInfo/ClinicalChar_2019.xlsx",startRow = 2)
colnames(Clin19)[1] <- "PatientID"
Clin19%<>%.[c(1,3,9,12,13,29,35)]%>%
  set_colnames(c("PatientID","nonsynonymous","gender","BR","PFS","progression","preMAPK"))
#===== Q: ExTRECT versus Timer2 =====
timer2Est <- read.csv("Data/Timer2/Timer2Ested.csv")%>%t()
colnames(timer2Est) <- timer2Est[1,]
combInfil <- timer2Est%>%
  .[-1,]%>%as.data.frame()%>%
  rownames_to_column("PatientID")%>%
  filter(!grepl("Normal",PatientID))%>%
  mutate(PatientID = (gsub("MEL\\.(.*)\\.Tumor.*","\\1",PatientID)%>%
                        gsub("IPI\\_","",.)))%>%
  inner_join(Combined_ExTRECT%>%filter(PubYr!="Y2014")%>%select(PatientID,ADJ_TCF),by="PatientID")
cor2metrics <- apply(combInfil[,-1],2,as.numeric)%>%cor(.,method = "pearson")%>%.[,120]
apply((combInfil%>%filter(grepl("Patient",PatientID)))[,-1],2,as.numeric)%>%cor(.,method = "pearson")%>%.[,120]%>%View()
apply((combInfil%>%filter(grepl("Pat[0-9]",PatientID)))[,-1],2,as.numeric)%>%cor(.,method = "pearson")%>%.[,120]%>%View()
#===== Figure 1: Correlation Bar Plot =====
tcellTimerNames <- grep("T cell CD8\\+\\_",colnames(timer2Est),value=T)
require(forcats)
cor2metrics[tcellTimerNames]%>%
  as.data.frame()%>%
  set_colnames("Pearson")%>%
  rownames_to_column("cellEst")%>%
  mutate(Pearson = as.numeric(Pearson))%>%
  mutate(cellEst = gsub(".*\\_","",cellEst))%>%
  mutate(cellEst =fct_reorder(cellEst,-Pearson))%>%
  ggplot(., aes(x=cellEst, y=Pearson, fill=cellEst)) +
  viridis::scale_fill_viridis(discrete = T,alpha = 0.95,option='E')+
  geom_bar(stat="identity")+
  theme(plot.title = element_text(color='black', hjust = 0.5),
        plot.background = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(color = "black", size = rel(1), fill = NA),
        panel.grid = element_blank(),
        axis.text.x = element_text(color='black',family = "Helvetica",size = rel(0.7),angle = 30,
                                   vjust = 1,hjust = 1),
        axis.title.y =  element_text(color='black',family = "Helvetica",size = rel(1)),
        axis.title.x = element_blank(),
        legend.key = element_rect(fill = NA),
        legend.position ="none")+
  labs(y="Pearson Correlation with TCF")
ggsave(filename = "draftPlot/Figure1_Corplot.png",width = 5,height = 5)
correlations <- cor(apply(combInfil[,-1]%>%.[,c(tcellTimerNames,"ADJ_TCF")],2,as.numeric))
corrplot(correlations, method="circle",cl.lim = c(-1,1))
#===== Figure 2: Pre-Post Treatment for 2014 Cohort =====
PrePost <- Combined_ExTRECT%>%
  filter(PubYr=="Y2014")%>%
  select(PatientID,ID,ADJ_TCF)%>%
  mutate(condition = gsub(".*\\.","",ID))%>%
  mutate(condition = factor(condition,levels = c("BRAF","pre")))
CliniTCF <- Clin14%>%
  mutate(PatientID = gsub(" ","",PatientID))%>%
  mutate(Response = as.character(Response))%>%
  inner_join(PrePost,by="PatientID")%>%
  mutate(Response = gsub(" ","",Response))

pbox <- ggplot(CliniTCF%>%filter(Response %in% c('SD',"PD")), aes(x=condition, y= ADJ_TCF,fill=condition)) +
  geom_boxplot()+
  theme(plot.title = element_text(color='black', hjust = 0.5),
        panel.background = element_blank(),
        panel.border = element_rect(color = "black", size = 1, fill = NA),
        panel.grid = element_blank(),
        axis.text = element_text(color='black'),
        legend.key = element_rect(fill = NA),
        legend.position ="none")+
  labs(x=NULL,y="Adjusted TCF",title = "Non-Progressed")+
  stat_compare_means()
pscatter <- ggplot(CliniTCF%>%filter(Response %in% c('SD',"PD")), aes(x=condition, y= ADJ_TCF,group=PatientID)) +
  geom_point(aes(color=condition)) +
  geom_line()+
  theme(plot.title = element_text(color='black', hjust = 0.5),
        panel.background = element_blank(),
        panel.border = element_rect(color = "black", size = 1, fill = NA),
        panel.grid = element_blank(),
        axis.text = element_text(color='black'),
        legend.key = element_rect(fill = NA),
        legend.position ="none")+
  labs(x=NULL,y="Adjusted TCF")
p1 <- ggarrange(pbox,pscatter,nrow = 1)

pbox <- ggplot(CliniTCF%>%filter(Response %in% c("PR")), aes(x=condition, y= ADJ_TCF,fill=condition)) +
  geom_boxplot()+
  theme(plot.title = element_text(color='black', hjust = 0.5),
        panel.background = element_blank(),
        panel.border = element_rect(color = "black", size = 1, fill = NA),
        panel.grid = element_blank(),
        axis.text = element_text(color='black'),
        legend.key = element_rect(fill = NA),
        legend.position ="none")+
  labs(x=NULL,y="Adjusted TCF",title = "Progressed")+
  stat_compare_means()
pscatter <- ggplot(CliniTCF%>%filter(Response %in% c("PR")), aes(x=condition, y= ADJ_TCF,group=PatientID)) +
  geom_point(aes(color=condition)) +
  geom_line()+
  theme(plot.title = element_text(color='black', hjust = 0.5),
        panel.background = element_blank(),
        panel.border = element_rect(color = "black", size = 1, fill = NA),
        panel.grid = element_blank(),
        axis.text = element_text(color='black'),
        legend.key = element_rect(fill = NA),
        legend.position ="none")+
  labs(x=NULL,y="Adjusted TCF")
p2 <- ggarrange(pbox,pscatter,nrow = 1)

ggsave(plot = ggarrange(p1,p2,
          ncol = 1,
          labels = c("Non-Progressed","Progressed")),
       filename = "draftPlot/Figure2_PrePost.png")


p3 <- ggplot(CliniTCF, aes(x=ADJ_TCF))+
  geom_histogram(position="identity", alpha=0.5,fill="#E69F00",color = "#E69F00",binwidth = 0.01)+
  #annotate("text", x = 75, y = 320,label = "Median Pick Order: 7")+
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(color='black',family = "Helvetica"),
        axis.title = element_text(color='black',family = "Helvetica"),
        panel.border = element_rect(color = "black", size = 1, fill = NA))+
  labs(x="Adjusted TCF",y="Frequency")
pPrePost <- ggarrange(p1,p2,
                      ncol = 1)
ggsave(plot = ggarrange(p3,pPrePost,nrow=1,widths = c(1,1.5)),
         filename = "draftPlot/Figure2_PrePost.png",width = 8,height = 6)

Clin14%>%
  mutate(PatientID = gsub(" ","",PatientID))%>%
  mutate(Response = as.character(Response))%>%
  inner_join(PrePost,by="PatientID")%>%
  dplyr::select(-ID)%>%
  spread(condition,ADJ_TCF)%>%
  mutate(diffTCF = BRAF-pre)%>%
  filter(PatientID %in% subPat)%>%
  filter(grepl("PR",Response))%>%
  mutate(f = ifelse(diffTCF>0,1,0))%>%
  group_by(f)%>%
  summarise(n=n())
#===== Q: With and without BRAF for 2015 and 2019 =====
C15D <- Clin15%>%
  rename(PatientID = patient)%>%
  inner_join(Combined_ExTRECT%>%
               filter(PubYr=="Y2015")%>%
               dplyr::select(ADJ_TCF,PatientID),by="PatientID")%>%
  mutate(preBRAF = case_when(pre_BRAF==1~"Yes",TRUE~"No"))

ggplot(C15D, aes(x=preBRAF, y= ADJ_TCF,fill=preBRAF)) +
  geom_boxplot()+
  theme(plot.title = element_text(color='black', hjust = 0.5),
        plot.background = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(color = "black", size = 1, fill = NA),
        panel.grid = element_blank(),
        axis.text = element_text(color='black'),
        legend.key = element_rect(fill = NA),
        legend.position ="none")+
  labs(x=NULL,y="Adjusted TCF")+
  stat_compare_means()
C15D%<>%
  mutate(group = case_when(group!="nonresponse"~"response",
                           TRUE~group)%>%factor(.,levels=c("nonresponse","response")))%>%
  mutate(preBRAF = factor(preBRAF,levels=c("Yes","No")))

ggplot(C15D, aes(x=ADJ_TCF, color=preBRAF,fill=preBRAF)) +
  geom_histogram(aes(y=..density..,color=preBRAF,fill=preBRAF), alpha=0.5, 
                 position="identity")+
  geom_density(alpha=.2) +
  theme(plot.title = element_text(color='black', hjust = 0.5),
        plot.background = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(color = "black", size = 1, fill = NA),
        panel.grid = element_blank(),
        axis.text = element_text(color='black'),
        legend.key = element_rect(fill = NA),
        legend.position ="bottom")

glm(ADJ_TCF~as.factor(preBRAF)+nonsynonymous+age_start+gender,data=C15D)%>%
  summary()
glm(as.factor(group)~ADJ_TCF+neos50+age_start+gender+as.factor(preBRAF),data=C15D,family =binomial)%>%
  summary()

C19D <- Clin19%>%
  inner_join(Combined_ExTRECT%>%
               filter(PubYr=="Y2019")%>%
               dplyr::select(ADJ_TCF,PatientID),by="PatientID")%>%
  mutate(preMAPK = case_when(preMAPK==1~"Yes",TRUE~"No"))
#===== Figure 3: Survival Analaysis =====
# https://rviews.rstudio.com/2017/09/25/survival-analysis-with-r/
require(survival)
cox15 <- coxph(Surv(overall_survivial, progression) ~ pre_BRAF + age_start + gender+ADJ_TCF+nonsynonymous, data = C15D)
summary(cox15)
library("survminer")
C15SuvD <- Clin15%>%
  mutate(preBRAF = case_when(pre_BRAF==1~"Yes",TRUE~"No"))%>%
  mutate(preBRAF = factor(preBRAF,levels=c("Yes","No")))
p1_15 <- ggsurvplot(survfit(Surv(progression_free, progression) ~ pre_BRAF,data=C15SuvD), data = C15SuvD)
p2_15 <- ggsurvplot(survfit(Surv(overall_survival, dead) ~ pre_BRAF,data=C15SuvD), data = C15SuvD)
comb15 <- arrange_ggsurvplots(list(p1_15,p2_15),ncol = 1, nrow = 2,print = F)
ggsave("draftPlot/Figure3a_Survival15.pdf", comb15)
C19SuvD <- Clin19%>%
  mutate(pretrt = case_when(priorMAPKTx==1~"Yes",TRUE~"No"))%>%
  mutate(pretrt = factor(pretrt,levels=c("Yes","No")))
p1_19 <- ggsurvplot(survfit(Surv(PFS, progressed) ~ pretrt,data=C19SuvD), data = C19SuvD)
p2_19 <- ggsurvplot(survfit(Surv(OS, dead) ~ pretrt,data=C19SuvD), data = C19SuvD)
comb19 <- arrange_ggsurvplots(list(p1_19,p2_19),ncol = 1, nrow = 2,print = F)
ggsave("draftPlot/Figure3b_Survival19.pdf", comb19)

ggplot(Combined_ExTRECT%>%
         filter(PubYr!="Y2014"), aes(x=ADJ_TCF))+
  geom_histogram(position="identity", alpha=0.5,fill="#E69F00",color = "#E69F00",binwidth = 0.01)+
  facet_wrap(~PubYr,nrow = 1,scales = "free")+
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(color='black',family = "Helvetica"),
        axis.title = element_text(color='black',family = "Helvetica"),
        panel.border = element_rect(color = "black", size = 1, fill = NA))+
  labs(x="Adjusted TCF",y="Frequency")
 ggsave("draftPlot/Figure3c_Distribution.pdf",width = 7,height = 4)
#===== Q: T cell CDK4+, CDK8+ TIMER (MAPK pre-treated or not, 2015 & 2019) =====
tcellNames <- grep("T cell",colnames(timer2Est),value=T)
cor(timer2Est[,tcellNames]%>%as.data.frame()%>%mutate_all(as.numeric))%>%View()
Timer2Result <- timer2Est[,c("T cell CD4+_TIMER", "T cell CD8+_TIMER")]%>%
  as.data.frame()%>%
  set_colnames(c("TimerCDK4","TimerCDK8"))%>%
  rownames_to_column("Entity")%>%
  filter(!grepl("Normal",Entity))%>%
  mutate(PatientID = (gsub("MEL\\.(.*)\\.Tumor.*","\\1",Entity)%>%
                        gsub("IPI\\_","",.)))
RNAseqTimer<- read.csv("Data/Timer2/TimerExprRNAseq.csv")%>%
  as.data.frame()%>%t()
colnames(RNAseqTimer) <- RNAseqTimer[1,]
RNAseqTimer%<>%.[-1,]%>%as.data.frame()%>%mutate_all(as.numeric)
RNAseqTimer%<>%
  as.data.frame()%>%
  .[,c("PDCD1","CD274")]%>%
  rownames_to_column("Entity")%>%
  inner_join(Timer2Result,by="Entity")
RNAseqC15D <- C15D%>%
  dplyr::select(PatientID,age_start,gender,group,preBRAF,nonsynonymous)%>%
  rename(MutLoad = nonsynonymous)%>%
  inner_join(RNAseqTimer%>%filter(grepl("IPI",Entity)),by="PatientID")%>%
  mutate(TimerCDK4 = as.numeric(TimerCDK4)+1e-5,
         TimerCDK8 = as.numeric(TimerCDK8)+1e-5)%>%
  mutate(group = case_when(group!="nonresponse"~"response",
                           TRUE~group)%>%
           factor(.,levels=c("nonresponse","response")))
lm(TimerCDK4~as.factor(preBRAF)+MutLoad +age_start+gender+PDCD1+CD274,
    data=RNAseqC15D)%>%summary()
lm(TimerCDK8~as.factor(preBRAF)+MutLoad +age_start+gender+PDCD1+CD274,
   data=RNAseqC15D)%>%summary()
glm(as.factor(group)~TimerCDK4+MutLoad+age_start+gender+PDCD1+CD274+as.factor(preBRAF),
    data=RNAseqC15D,family=binomial)%>%summary() 
glm(as.factor(group)~TimerCDK8+MutLoad+age_start+gender+PDCD1+CD274+as.factor(preBRAF),
    data=RNAseqC15D,family=binomial)%>%summary()  

Clin19_2nd <- openxlsx::read.xlsx("Data/SampleInfo/2019_clinical_data.xlsx")
Clin19_2nd%<>%
  mutate(PatientID = gsub("MEL-","",patient_id))%>%
  dplyr::select(PatientID,age,sex)%>%
  inner_join(Clin19%>%dplyr::select(-gender),by="PatientID")
RNAseqC19D <- Clin19_2nd%>%
  inner_join(RNAseqTimer%>%filter(grepl("Patient",Entity)),by="PatientID")
RNAseqC19D%<>%
  mutate(group = case_when(BR %in% c("SD","PD")~"Nonresponse",
                           TRUE~"Response"))%>%
  mutate(TimerCDK4 = as.numeric(TimerCDK4)+1e-5,
         TimerCDK8 = as.numeric(TimerCDK8)+1e-5)
glm(as.factor(group)~TimerCDK4+nonsynonymous+age+as.factor(sex)+PDCD1+CD274+as.factor(preMAPK),
    data=RNAseqC19D,family=binomial)%>%summary() 
glm(as.factor(group)~TimerCDK8+nonsynonymous+age+as.factor(sex)+PDCD1+CD274+as.factor(preMAPK),
    data=RNAseqC19D,family=binomial)%>%summary() 
lm(TimerCDK4~as.factor(preMAPK)+nonsynonymous +age+as.factor(sex)+PDCD1+CD274,
   data=RNAseqC19D)%>%summary()
lm(TimerCDK8~as.factor(preMAPK)+nonsynonymous +age+as.factor(sex)+PDCD1+CD274,
   data=RNAseqC19D)%>%summary()



#===== Figure 4: 2015 Combined 2019, Cell Type & GeneMarkers =====
tcellNames <- grep("XCELL",colnames(timer2Est),value=T, ignore.case = T)
CellTypes <- c("T cell CD8+_XCELL",
  "T cell CD4+ Th1_XCELL","T cell CD4+ Th2_XCELL",
  "Macrophage M1_XCELL","Macrophage M2_XCELL",
  "T cell regulatory (Tregs)_XCELL","NK cell_XCELL",
  "Myeloid dendritic cell_XCELL","B cell_XCELL")
GeneMarkers <-c("HHLA1","HHLA2","HHLA3","PDCD1","CD274","ADRB1","ADRB2","ADRB3")
RNAseqTimer<- read.csv("Data/Timer2/TimerExprRNAseq.csv")%>%
  as.data.frame()%>%t()
colnames(RNAseqTimer) <- RNAseqTimer[1,]
RNAseqTimer%<>%.[-1,]%>%as.data.frame()%>%mutate_all(as.numeric)
Timer2Result <- timer2Est[,CellTypes]%>%
  as.data.frame()%>%
  rownames_to_column("Entity")%>%
  filter(!grepl("Normal",Entity))%>%
  mutate(PatientID = (gsub("MEL\\.(.*)\\.Tumor.*","\\1",Entity)%>%
                        gsub("IPI\\_","",.)))
RNAseqTimer%<>%
  as.data.frame()%>%
  .[,GeneMarkers]%>%
  rownames_to_column("Entity")%>%
  inner_join(Timer2Result,by="Entity")
RNAseqC1519D <- C15D%>%
  dplyr::select(PatientID,preTRT = preBRAF)%>%
  rbind(C19D%>%
           dplyr::select(PatientID,preTRT = preMAPK))%>%
  inner_join(RNAseqTimer,by="PatientID")
RNAseqC1519D%>%
  dplyr::select(-Entity)%>%
  gather(features,value,-PatientID,-preTRT)%>%
  mutate(features = gsub("\\_XCELL","",features))%>%
  mutate(value = as.numeric(value))%>%
  ggplot(., aes(x=preTRT, y=value, fill=preTRT)) +
  geom_boxplot()+
  facet_wrap(~features,nrow = 3,scales="free")+
  theme(plot.title = element_text(color='black', hjust = 0.5),
        panel.background = element_blank(),
        panel.border = element_rect(color = "black", size = 1, fill = NA),
        panel.grid = element_blank(),
        axis.text = element_text(color='black'),
        legend.key = element_rect(fill = NA))
ggsave("draftPlot/Figure4_GeneMarker_ImmuneCellEst.png",width = 11.5,height = 7)
for(i in colnames(RNAseqC1519D)[-c(1:3)]){
  p <- RNAseqC1519D%>%
    dplyr::select(preTRT,i)%>%
    rename("marker" = i)%>%
    mutate(marker = as.numeric(marker))%>%
    ggplot(., aes(x=preTRT, y=marker, fill=preTRT)) +
    geom_boxplot()+
    theme(plot.title = element_text(color='black', hjust = 0.5),
          plot.background = element_blank(),
          panel.background = element_blank(),
          panel.border = element_rect(color = "black", size = 1, fill = NA),
          panel.grid = element_blank(),
          axis.text = element_text(color='black'),
          legend.key = element_rect(fill = NA))+
    labs(title = i)+
    stat_compare_means()
  ggsave(plot = p,filename = glue("~/Desktop/{i}.png"))
}
#===== Q: 
#===== Q: Pre-Post treated and clinical characteristics =====
PPH15 <- openxlsx::read.xlsx("Data/SampleInfo/2015_IO_pairs_census.xlsx")%>%
  dplyr::select(facets_purity,facets_ploidy,participant)%>%
  mutate(PatientID = gsub(".*\\_","",participant))
C15D%>%
  dplyr::select(PatientID,preTRT = preBRAF)%>%
  inner_join(PPH15,by="PatientID")%>%
  ggplot(., aes(x=preTRT, y=facets_purity, fill=preTRT)) +
  geom_boxplot()+
  theme(plot.title = element_text(color='black', hjust = 0.5),
        plot.background = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(color = "black", size = 1, fill = NA),
        panel.grid = element_blank(),
        axis.text = element_text(color='black'),
        legend.key = element_rect(fill = NA))

C15D%>%
  dplyr::select(PatientID,preTRT = preBRAF)%>%
  inner_join(PPH15,by="PatientID")%>%
  ggplot(., aes(x=preTRT, y=facets_ploidy, fill=preTRT)) +
  geom_boxplot()+
  theme(plot.title = element_text(color='black', hjust = 0.5),
        plot.background = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(color = "black", size = 1, fill = NA),
        panel.grid = element_blank(),
        axis.text = element_text(color='black'),
        legend.key = element_rect(fill = NA))


#===== BackUp Code ======
lm(TimerCDK4~as.factor(preBRAF)+MutLoad +age_start+gender+PDCD1+CD274,
   data=RNAseqC15D)%>%summary()
lm(TimerCDK8~as.factor(preBRAF)+MutLoad +age_start+gender+PDCD1+CD274,
   data=RNAseqC15D)%>%summary()
glm(as.factor(group)~TimerCDK4+MutLoad+age_start+gender+PDCD1+CD274+as.factor(preBRAF),
    data=RNAseqC15D,family=binomial)%>%summary() 
glm(as.factor(group)~TimerCDK8+MutLoad+age_start+gender+PDCD1+CD274+as.factor(preBRAF),
    data=RNAseqC15D,family=binomial)%>%summary()  

Clin19_2nd <- openxlsx::read.xlsx("Data/SampleInfo/2019_clinical_data.xlsx")
Clin19_2nd%<>%
  mutate(PatientID = gsub("MEL-","",patient_id))%>%
  dplyr::select(PatientID,age,sex)%>%
  inner_join(Clin19%>%dplyr::select(-gender),by="PatientID")
RNAseqC19D <- Clin19_2nd%>%
  inner_join(RNAseqTimer%>%filter(grepl("Patient",Entity)),by="PatientID")
RNAseqC19D%<>%
  mutate(group = case_when(BR %in% c("SD","PD")~"Nonresponse",
                           TRUE~"Response"))%>%
  mutate(TimerCDK4 = as.numeric(TimerCDK4)+1e-5,
         TimerCDK8 = as.numeric(TimerCDK8)+1e-5)
glm(as.factor(group)~TimerCDK4+nonsynonymous+age+as.factor(sex)+PDCD1+CD274+as.factor(preMAPK),
    data=RNAseqC19D,family=binomial)%>%summary() 
glm(as.factor(group)~TimerCDK8+nonsynonymous+age+as.factor(sex)+PDCD1+CD274+as.factor(preMAPK),
    data=RNAseqC19D,family=binomial)%>%summary() 
lm(TimerCDK4~as.factor(preMAPK)+nonsynonymous +age+as.factor(sex)+PDCD1+CD274,
   data=RNAseqC19D)%>%summary()
lm(TimerCDK8~as.factor(preMAPK)+nonsynonymous +age+as.factor(sex)+PDCD1+CD274,
   data=RNAseqC19D)%>%summary()


