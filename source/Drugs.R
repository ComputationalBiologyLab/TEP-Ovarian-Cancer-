#########################  Drugs for Genes with AUC >0.6 in both datasets  #####################

#Get working directory 

getwd()

#Set working directory to the location of Data/drugs through setwd()

library(readr)
library(tidygraph)
library(tidyr)
library(tidyselect)
library(tidytree)
library(tidyverse)
library(dplyr)
library(dbplyr)
library(dtplyr)
library(magrittr)
library(readxl)


AUC_defined_DEGs=read_excel("Supplementary table.xlsx", sheet=4, skip=1)

AUC_Stage1=AUC_defined_DEGs%>%
  filter(Stage=="Stage1")

AUC_Stage2=AUC_defined_DEGs%>%
  filter(Stage=="Stage2")

AUC_Stage3=AUC_defined_DEGs%>%
  filter(Stage=="Stage3")

AUC_Stage4=AUC_defined_DEGs%>%
  filter(Stage=="Stage4")

######## Stage 1#######
#UP DEGs

AUC1_UP=AUC_Stage1%>%
  filter(Status=="Up")

AUC1_UP
gene_expression_MDM4_ovary_All_information <- read_csv("S1 drugs/gene_expression_MDM4_ovary_All_information.csv")
gene_expression_KIF13B_ovary_All_information <- read_csv("S1 drugs/gene_expression_KIF13B_ovary_All_information.csv")
gene_expression_DNAJC3_ovary_All_information <- read_csv("S1 drugs/gene_expression_DNAJC3_ovary_All_information.csv")
gene_expression_SLC38A1_ovary_All_information <- read_csv("S1 drugs/gene_expression_SLC38A1_ovary_All_information.csv")
gene_expression_NEDD9_ovary_All_information <- read_csv("S1 drugs/gene_expression_NEDD9_ovary_All_information.csv")
gene_expression_FZR1_ovary_All_information <- read_csv("S1 drugs/gene_expression_FZR1_ovary_All_information.csv")
gene_expression_CLEC4E_ovary_All_information <- read_csv("S1 drugs/gene_expression_CLEC4E_ovary_All_information.csv")
gene_expression_BTAF1_ovary_All_information <- read_csv("S1 drugs/gene_expression_BTAF1_ovary_All_information.csv")
gene_expression_FNDC3B_ovary_All_information <- read_csv("S1 drugs/gene_expression_FNDC3B_ovary_All_information.csv")
gene_expression_CPD_ovary_All_information <- read_csv("S1 drugs/gene_expression_CPD_ovary_All_information.csv")
gene_expression_MGAT4A_ovary_All_information <- read_csv("S1 drugs/gene_expression_MGAT4A_ovary_All_information.csv")
gene_expression_S100A12_ovary_All_information <- read_csv("S1 drugs/gene_expression_S100A12_ovary_All_information.csv")
gene_expression_IL1R2_ovary_All_information <- read_csv("S1 drugs/gene_expression_IL1R2_ovary_All_information.csv")


###########
AUC1_DN=AUC_Stage1%>%
  filter(Status=="Down")

#DOWN DEGs
gene_expression_GTPBP2_ovary_All_information <- read_csv("S1 drugs/gene_expression_GTPBP2_ovary_All_information.csv")
gene_expression_PARP1_ovary_All_information <- read_csv("S1 drugs/gene_expression_PARP1_ovary_All_information.csv")
gene_expression_CASP3_ovary_All_information <- read_csv("S1 drugs/gene_expression_CASP3_ovary_All_information.csv")
gene_expression_MCM3_ovary_All_information <- read_csv("S1 drugs/gene_expression_MCM3_ovary_All_information.csv")
gene_expression_PREP_ovary_All_information <- read_csv("S1 drugs/gene_expression_PREP_ovary_All_information.csv")
gene_expression_MRPL24_ovary_All_information <- read_csv("S1 drugs/gene_expression_MRPL24_ovary_All_information.csv")
gene_expression_LCK_ovary_All_information <- read_csv("S1 drugs/gene_expression_LCK_ovary_All_information.csv")
gene_expression_SLC27A4_ovary_All_information <- read_csv("S1 drugs/gene_expression_SLC27A4_ovary_All_information.csv")
gene_expression_IARS2_ovary_All_information <- read_csv("S1 drugs/gene_expression_IARS2_ovary_All_information.csv")
gene_expression_MSH2_ovary_All_information <- read_csv("S1 drugs/gene_expression_MSH2_ovary_All_information.csv")
gene_expression_GLT8D1_ovary_All_information <- read_csv("S1 drugs/gene_expression_GLT8D1_ovary_All_information.csv")
gene_expression_DIS3L_ovary_All_information <- read_csv("S1 drugs/gene_expression_DIS3L_ovary_All_information.csv")
gene_expression_VWA5A_ovary_All_information <- read_csv("S1 drugs/gene_expression_VWA5A_ovary_All_information.csv")
gene_expression_CUL1_ovary_All_information <- read_csv("S1 drugs/gene_expression_CUL1_ovary_All_information.csv")
gene_expression_SEC23IP_ovary_All_information <- read_csv("S1 drugs/gene_expression_SEC23IP_ovary_All_information.csv")
gene_expression_FARSA_ovary_All_information <- read_csv("S1 drugs/gene_expression_FARSA_ovary_All_information.csv")
gene_expression_PIK3AP1_ovary_All_information <- read_csv("S1 drugs/gene_expression_PIK3AP1_ovary_All_information.csv")

S1_drugs=rbind(gene_expression_MDM4_ovary_All_information, 
               gene_expression_KIF13B_ovary_All_information, 
               gene_expression_DNAJC3_ovary_All_information, 
               gene_expression_SLC38A1_ovary_All_information, 
               gene_expression_NEDD9_ovary_All_information, 
               gene_expression_FZR1_ovary_All_information, 
               gene_expression_CLEC4E_ovary_All_information, 
               gene_expression_BTAF1_ovary_All_information, 
               gene_expression_FNDC3B_ovary_All_information, 
               gene_expression_CPD_ovary_All_information, 
               gene_expression_MGAT4A_ovary_All_information, 
               gene_expression_S100A12_ovary_All_information, 
               gene_expression_IL1R2_ovary_All_information,
               gene_expression_GTPBP2_ovary_All_information, 
               gene_expression_PARP1_ovary_All_information, 
               gene_expression_CASP3_ovary_All_information, 
               gene_expression_MCM3_ovary_All_information, 
               gene_expression_PREP_ovary_All_information, 
               gene_expression_MRPL24_ovary_All_information, 
               gene_expression_LCK_ovary_All_information, 
               gene_expression_SLC27A4_ovary_All_information, 
               gene_expression_IARS2_ovary_All_information, 
               gene_expression_MSH2_ovary_All_information, 
               gene_expression_GLT8D1_ovary_All_information, 
               gene_expression_DIS3L_ovary_All_information, 
               gene_expression_VWA5A_ovary_All_information, 
               gene_expression_CUL1_ovary_All_information, 
               gene_expression_SEC23IP_ovary_All_information, 
               gene_expression_FARSA_ovary_All_information, 
               gene_expression_PIK3AP1_ovary_All_information)


S1_drugs <- S1_drugs %>%
  mutate(Status = ifelse(row_number() <= 4458, "Up", "Down"))

S1_drugs_up=S1_drugs%>%
  filter(Status=="Up")

S1_drugs_up=S1_drugs_up%>%
  filter(correlation < -0.6)

S1_drugs_up$F="favorable for up"

S1_drugs_down=S1_drugs%>%
  filter(Status=="Down")

S1_drugs_down=S1_drugs_down%>%
  filter(correlation>0.6)

S1_drugs_down$F="favorable for down"

S1_drugs=rbind(S1_drugs_up, S1_drugs_down)

S1_drugs=S1_drugs%>%
  filter(pvalue<0.05)

S1_drugs=S1_drugs[,-c(3,4,6,7)]

length(unique(S1_drugs$drug_name))
length(unique(S1_drugs$gene))


all(S1_drugs$gene %in% AUC_Stage1$hgnc_symbol)

######Stage2#####
#UP DEGs

AUC2_UP=AUC_Stage2%>%
  filter(Status=="Up")

gene_expression_CALR_ovary_All_information <- read_csv("S2 drugs/gene_expression_CALR_ovary_All_information.csv")
gene_expression_CTSW_ovary_All_information <- read_csv("S2 drugs/gene_expression_CTSW_ovary_All_information.csv")
gene_expression_KLHDC8B_ovary_All_information <- read_csv("S2 drugs/gene_expression_KLHDC8B_ovary_All_information.csv")
gene_expression_ARL2_ovary_All_information <- read_csv("S2 drugs/gene_expression_ARL2_ovary_All_information.csv")
gene_expression_NME4_ovary_All_information <- read_csv("S2 drugs/gene_expression_NME4_ovary_All_information.csv")
gene_expression_LY6E_ovary_All_information <- read_csv("S2 drugs/gene_expression_LY6E_ovary_All_information.csv")
gene_expression_QSOX1_ovary_All_information <- read_csv("S2 drugs/gene_expression_QSOX1_ovary_All_information.csv")
gene_expression_UPP1_ovary_All_information <- read_csv("S2 drugs/gene_expression_UPP1_ovary_All_information.csv")
gene_expression_ARID5A_ovary_All_information <- read_csv("S2 drugs/gene_expression_ARID5A_ovary_All_information.csv")
gene_expression_TPM2_ovary_All_information <- read_csv("S2 drugs/gene_expression_TPM2_ovary_All_information.csv")
gene_expression_COL6A3_ovary_All_information <- read_csv("S2 drugs/gene_expression_COL6A3_ovary_All_information.csv")
gene_expression_LPCAT1_ovary_All_information <- read_csv("S2 drugs/gene_expression_LPCAT1_ovary_All_information.csv")
gene_expression_SMURF2_ovary_All_information <- read_csv("S2 drugs/gene_expression_SMURF2_ovary_All_information.csv")
gene_expression_DGKQ_ovary_All_information <- read_csv("S2 drugs/gene_expression_DGKQ_ovary_All_information.csv")
gene_expression_NFKB1_ovary_All_information <- read_csv("S2 drugs/gene_expression_NFKB1_ovary_All_information.csv")

#DOWN DEGs
AUC2_DN=AUC_Stage2%>%
  filter(Status=="Down")

AUC2_DN
gene_expression_CDYL_ovary_All_information <- read_csv("S2 drugs/gene_expression_CDYL_ovary_All_information.csv")
gene_expression_EIF4ENIF1_ovary_All_information <- read_csv("S2 drugs/gene_expression_EIF4ENIF1_ovary_All_information.csv")
gene_expression_METAP1_ovary_All_information <- read_csv("S2 drugs/gene_expression_METAP1_ovary_All_information.csv")
gene_expression_WDR44_ovary_All_information <- read_csv("S2 drugs/gene_expression_WDR44_ovary_All_information.csv")
gene_expression_IGSF6_ovary_All_information <- read_csv("S2 drugs/gene_expression_IGSF6_ovary_All_information.csv")
gene_expression_MOB3C_ovary_All_information <- read_csv("S2 drugs/gene_expression_MOB3C_ovary_All_information.csv")
gene_expression_G3BP1_ovary_All_information <- read_csv("S2 drugs/gene_expression_G3BP1_ovary_All_information.csv")
gene_expression_TBC1D14_ovary_All_information <- read_csv("S2 drugs/gene_expression_TBC1D14_ovary_All_information.csv")
gene_expression_CLEC2D_ovary_All_information <- read_csv("S2 drugs/gene_expression_CLEC2D_ovary_All_information.csv")
gene_expression_NSF_ovary_All_information <- read_csv("S2 drugs/gene_expression_NSF_ovary_All_information.csv")
gene_expression_STAG1_ovary_All_information <- read_csv("S2 drugs/gene_expression_STAG1_ovary_All_information.csv")
gene_expression_NUDT5_ovary_All_information <- read_csv("S2 drugs/gene_expression_NUDT5_ovary_All_information.csv")
gene_expression_AMICA1_ovary_All_information <- read_csv("S2 drugs/gene_expression_AMICA1_ovary_All_information.csv")
gene_expression_SLC25A20_ovary_All_information <- read_csv("S2 drugs/gene_expression_SLC25A20_ovary_All_information.csv")
gene_expression_ZNF542_ovary_All_information <- read_csv("S2 drugs/gene_expression_ZNF542P_ovary_All_information.csv")
gene_expression_MAT2A_ovary_All_information <- read_csv("S2 drugs/gene_expression_MAT2A_ovary_All_information.csv")
gene_expression_CASP1_ovary_All_information <- read_csv("S2 drugs/gene_expression_CASP1_ovary_All_information.csv")
gene_expression_OAS1_ovary_All_information <- read_csv("S2 drugs/gene_expression_OAS1_ovary_All_information.csv")
gene_expression_SF3B3_ovary_All_information <- read_csv("S2 drugs/gene_expression_SF3B3_ovary_All_information.csv")
gene_expression_WARS_ovary_All_information <- read_csv("S2 drugs/gene_expression_WARS1_ovary_All_information.csv")

S2_drugs=rbind(gene_expression_CALR_ovary_All_information, 
               gene_expression_CTSW_ovary_All_information, 
               gene_expression_KLHDC8B_ovary_All_information, 
               gene_expression_ARL2_ovary_All_information, 
               gene_expression_NME4_ovary_All_information, 
               gene_expression_LY6E_ovary_All_information, 
               gene_expression_QSOX1_ovary_All_information, 
               gene_expression_UPP1_ovary_All_information, 
               gene_expression_ARID5A_ovary_All_information, 
               gene_expression_TPM2_ovary_All_information, 
               gene_expression_COL6A3_ovary_All_information, 
               gene_expression_LPCAT1_ovary_All_information, 
               gene_expression_SMURF2_ovary_All_information, 
               gene_expression_DGKQ_ovary_All_information, 
               gene_expression_NFKB1_ovary_All_information, 
               gene_expression_CDYL_ovary_All_information, 
               gene_expression_EIF4ENIF1_ovary_All_information, 
               gene_expression_METAP1_ovary_All_information, 
               gene_expression_WDR44_ovary_All_information, 
               gene_expression_IGSF6_ovary_All_information, 
               gene_expression_MOB3C_ovary_All_information, 
               gene_expression_G3BP1_ovary_All_information, 
               gene_expression_TBC1D14_ovary_All_information, 
               gene_expression_CLEC2D_ovary_All_information, 
               gene_expression_NSF_ovary_All_information, 
               gene_expression_STAG1_ovary_All_information, 
               gene_expression_NUDT5_ovary_All_information, 
               gene_expression_AMICA1_ovary_All_information, 
               gene_expression_SLC25A20_ovary_All_information, 
               gene_expression_ZNF542_ovary_All_information, 
               gene_expression_MAT2A_ovary_All_information, 
               gene_expression_CASP1_ovary_All_information, 
               gene_expression_OAS1_ovary_All_information, 
               gene_expression_SF3B3_ovary_All_information, 
               gene_expression_WARS_ovary_All_information
)

#Aliases
S2_drugs$gene <- str_replace(S2_drugs$gene, "JAML", "AMICA1")
S2_drugs$gene <- str_replace(S2_drugs$gene, "ZNF542P", "ZNF542")
S2_drugs$gene <- str_replace(S2_drugs$gene, "WARS1", "WARS")


S2_drugs <- S2_drugs %>%
  mutate(Status = ifelse(row_number() <= 5074, "Up", "Down"))

S2_drugs_up=S2_drugs%>%
  filter(Status=="Up")

S2_drugs_up=S2_drugs_up%>%
  filter(correlation < -0.6)

S2_drugs_up$F="favorable for up"

S2_drugs_down=S2_drugs%>%
  filter(Status=="Down")

S2_drugs_down=S2_drugs_down%>%
  filter(correlation>0.6)

S2_drugs_down$F="favorable for down"

S2_drugs=rbind(S2_drugs_up, S2_drugs_down)

S2_drugs=S2_drugs%>%
  filter(pvalue<0.05)

S2_drugs=S2_drugs[,-c(3,4,6,7)]

length(unique(S2_drugs$drug_name))
length(unique(S2_drugs$gene))

all(S2_drugs$gene %in% AUC_Stage2$hgnc_symbol)


#####Stage3####
#UP DEGs
gene_expression_RANBP2_ovary_All_information <- read_csv("S3 drugs/gene_expression_RANBP2_ovary_All_information.csv")
gene_expression_IRF1_ovary_All_information <- read_csv("S3 drugs/gene_expression_IRF1_ovary_All_information.csv")
gene_expression_GZMB_ovary_All_information <- read_csv("S3 drugs/gene_expression_GZMB_ovary_All_information.csv")
gene_expression_PIM3_ovary_All_information <- read_csv("S3 drugs/gene_expression_PIM3_ovary_All_information.csv")
gene_expression_SYTL3_ovary_All_information <- read_csv("S3 drugs/gene_expression_SYTL3_ovary_All_information.csv")
gene_expression_PELI1_ovary_All_information <- read_csv("S3 drugs/gene_expression_PELI1_ovary_All_information.csv")
gene_expression_PRSS50_ovary_All_information <- read_csv("S3 drugs/gene_expression_PRSS50_ovary_All_information.csv")
#Down DEGs
gene_expression_ATP13A4_ovary_All_information <- read_csv("S3 drugs/gene_expression_ATP13A4_ovary_All_information.csv")

S3_drugs=rbind(gene_expression_RANBP2_ovary_All_information, 
               gene_expression_IRF1_ovary_All_information, 
               gene_expression_GZMB_ovary_All_information, 
               gene_expression_PIM3_ovary_All_information, 
               gene_expression_SYTL3_ovary_All_information, 
               gene_expression_PELI1_ovary_All_information, 
               gene_expression_PRSS50_ovary_All_information, 
               gene_expression_ATP13A4_ovary_All_information)


S3_drugs <- S3_drugs %>%
  mutate(Status = ifelse(row_number() <= 2267, "Up", "Down"))

S3_drugs_up=S3_drugs%>%
  filter(Status=="Up")

S3_drugs_up=S3_drugs_up%>%
  filter(correlation < -0.6)

S3_drugs_up$F="favorable for up"

S3_drugs_down=S3_drugs%>%
  filter(Status=="Down")

S3_drugs_down=S3_drugs_down%>%
  filter(correlation>0.6)

S3_drugs_down$F="favorable for down"

S3_drugs=rbind(S3_drugs_up, S3_drugs_down)

S3_drugs=S3_drugs%>%
  filter(pvalue<0.05)

S3_drugs=S3_drugs[,-c(3,4,6,7)]

length(unique(S3_drugs$drug_name))
length(unique(S3_drugs$gene))

######################### Stage4 Drugs for AUC >0.6 genes  #####################
#UP DEGs
AUC4_UP=AUC_Stage4%>%
  filter(Status=="Up")

AUC4_UP
gene_expression_CA1_ovary_All_information <- read_csv("S4 drugs/gene_expression_CA1_ovary_All_information.csv")
gene_expression_TEX9_ovary_All_information <- read_csv("S4 drugs/gene_expression_TEX9_ovary_All_information.csv")
gene_expression_HBB_ovary_All_information <- read_csv("S4 drugs/gene_expression_HBB_ovary_All_information.csv")
gene_expression_HPSE_ovary_All_information <- read_csv("S4 drugs/gene_expression_HPSE_ovary_All_information.csv")
gene_expression_FCGR2A_ovary_All_information <- read_csv("S4 drugs/gene_expression_FCGR2A_ovary_All_information.csv")
gene_expression_VPS8_ovary_All_information <- read_csv("S4 drugs/gene_expression_VPS8_ovary_All_information.csv")
gene_expression_SLC25A39_ovary_All_information <- read_csv("S4 drugs/gene_expression_SLC25A39_ovary_All_information.csv")
gene_expression_SLC25A37_ovary_All_information <- read_csv("S4 drugs/gene_expression_SLC25A37_ovary_All_information.csv")
gene_expression_HBD_ovary_All_information <- read_csv("S4 drugs/gene_expression_HBD_ovary_All_information.csv")
gene_expression_AHR_ovary_All_information <- read_csv("S4 drugs/gene_expression_AHR_ovary_All_information.csv")
gene_expression_CD96_ovary_All_information <- read_csv("S4 drugs/gene_expression_CD96_ovary_All_information.csv")
gene_expression_ENTPD5_ovary_All_information <- read_csv("S4 drugs/gene_expression_ENTPD5_ovary_All_information.csv")
gene_expression_S100A10_ovary_All_information <- read_csv("S4 drugs/gene_expression_S100A10_ovary_All_information.csv")
gene_expression_LXN_ovary_All_information <- read_csv("S4 drugs/gene_expression_LXN_ovary_All_information.csv")
gene_expression_MXI1_ovary_All_information <- read_csv("S4 drugs/gene_expression_MXI1_ovary_All_information.csv")
gene_expression_TMSB4XP8_ovary_All_information <- read_csv("S4 drugs/gene_expression_TMSB4XP8_ovary_All_information.csv")
gene_expression_FECH_ovary_All_information <- read_csv("S4 drugs/gene_expression_FECH_ovary_All_information.csv")
gene_expression_HBM_ovary_All_information <- read_csv("S4 drugs/gene_expression_HBM_ovary_All_information.csv")
gene_expression_NUDT4_ovary_All_information <- read_csv("S4 drugs/gene_expression_NUDT4_ovary_All_information.csv")
gene_expression_IFI6_ovary_All_information <- read_csv("S4 drugs/gene_expression_IFI6_ovary_All_information.csv")
gene_expression_MTRNR2L12_ovary_All_information <- read_csv("S4 drugs/gene_expression_MTRNR2L12_ovary_All_information.csv")
gene_expression_MTRNR2L2_ovary_All_information <- read_csv("S4 drugs/gene_expression_MTRNR2L2_ovary_All_information.csv")
gene_expression_SOS1_ovary_All_information <- read_csv("S4 drugs/gene_expression_SOS1_ovary_All_information.csv")
gene_expression_MTRNR2L4_ovary_All_information <- read_csv("S4 drugs/gene_expression_MTRNR2L4_ovary_All_information.csv")
gene_expression_IFI27_ovary_All_information <- read_csv("S4 drugs/gene_expression_IFI27_ovary_All_information.csv")
gene_expression_MTRNR2L9_ovary_All_information <- read_csv("S4 drugs/gene_expression_MTRNR2L9_ovary_All_information.csv")
gene_expression_MTRNR2L8_ovary_All_information <- read_csv("S4 drugs/gene_expression_MTRNR2L8_ovary_All_information.csv")
gene_expression_HERC5_ovary_All_information <- read_csv("S4 drugs/gene_expression_HERC5_ovary_All_information.csv")
gene_expression_HBG2_ovary_All_information <- read_csv("S4 drugs/gene_expression_HBG2_ovary_All_information.csv")
gene_expression_NISCH_ovary_All_information <- read_csv("S4 drugs/gene_expression_NISCH_ovary_All_information.csv")
gene_expression_TFRC_ovary_All_information <- read_csv("S4 drugs/gene_expression_TFRC_ovary_All_information.csv")
gene_expression_MCM2_ovary_All_information <- read_csv("S4 drugs/gene_expression_MCM2_ovary_All_information.csv")
gene_expression_RALGAPB_ovary_All_information <- read_csv("S4 drugs/gene_expression_RALGAPB_ovary_All_information.csv")
gene_expression_QPCT_ovary_All_information <- read_csv("S4 drugs/gene_expression_QPCT_ovary_All_information.csv")
gene_expression_GNS_ovary_All_information <- read_csv("S4 drugs/gene_expression_GNS_ovary_All_information.csv")
gene_expression_LTF_ovary_All_information <- read_csv("S4 drugs/gene_expression_LTF_ovary_All_information.csv")
gene_expression_DEFA3_ovary_All_information <- read_csv("S4 drugs/gene_expression_DEFA3_ovary_All_information.csv")
gene_expression_PDE4B_ovary_All_information <- read_csv("S4 drugs/gene_expression_PDE4B_ovary_All_information.csv")
gene_expression_DYSF_ovary_All_information <- read_csv("S4 drugs/gene_expression_DYSF_ovary_All_information.csv")
gene_expression_RPL23AP82_ovary_All_information <- read_csv("S4 drugs/gene_expression_RPL23AP82_ovary_All_information.csv")
gene_expression_PPP2R1B_ovary_All_information <- read_csv("S4 drugs/gene_expression_PPP2R1B_ovary_All_information.csv")
gene_expression_PADI2_ovary_All_information <- read_csv("S4 drugs/gene_expression_PADI2_ovary_All_information.csv")
gene_expression_RABL2A_ovary_All_information <- read_csv("S4 drugs/gene_expression_RABL2A_ovary_All_information.csv")

#DOWN DEGs
AUC4_DN=AUC_Stage4%>%
  filter(Status=="Down")

AUC4_DN

gene_expression_HSD17B3_ovary_All_information <- read_csv("S4 drugs/gene_expression_HSD17B3_ovary_All_information.csv")
gene_expression_ZNF385D_ovary_All_information <- read_csv("S4 drugs/gene_expression_ZNF385D_ovary_All_information.csv")
gene_expression_MPL_ovary_All_information <- read_csv("S4 drugs/gene_expression_MPL_ovary_All_information.csv")
gene_expression_PTK7_ovary_All_information <- read_csv("S4 drugs/gene_expression_PTK7_ovary_All_information.csv")
gene_expression_ANKH_ovary_All_information <- read_csv("S4 drugs/gene_expression_ANKH_ovary_All_information.csv")
gene_expression_CADM2_ovary_All_information <- read_csv("S4 drugs/gene_expression_CADM2_ovary_All_information.csv")
gene_expression_NOMO1_ovary_All_information <- read_csv("S4 drugs/gene_expression_NOMO1_ovary_All_information.csv")
gene_expression_SLC25A43_ovary_All_information <- read_csv("S4 drugs/gene_expression_SLC25A43_ovary_All_information.csv")
gene_expression_FCRL6_ovary_All_information <- read_csv("S4 drugs/gene_expression_FCRL6_ovary_All_information.csv")
gene_expression_TCF7_ovary_All_information <- read_csv("S4 drugs/gene_expression_TCF7_ovary_All_information.csv")
gene_expression_DDX11L5_ovary_All_information <- read_csv("S4 drugs/gene_expression_DDX11L5_ovary_All_information.csv")
gene_expression_EXOC3L2_ovary_All_information <- read_csv("S4 drugs/gene_expression_EXOC3L2_ovary_All_information.csv")
gene_expression_DDX11L10_ovary_All_information <- read_csv("S4 drugs/gene_expression_DDX11L10_ovary_All_information.csv")
gene_expression_EML2_ovary_All_information <- read_csv("S4 drugs/gene_expression_EML2_ovary_All_information.csv")
gene_expression_FCMR_ovary_All_information <- read_csv("S4 drugs/gene_expression_FCMR_ovary_All_information.csv")
gene_expression_C1orf87_ovary_All_information <- read_csv("S4 drugs/gene_expression_C1orf87_ovary_All_information.csv")
gene_expression_PGM2_ovary_All_information <- read_csv("S4 drugs/gene_expression_PGM2_ovary_All_information.csv")
gene_expression_LINS1_ovary_All_information <- read_csv("S4 drugs/gene_expression_LINS1_ovary_All_information.csv")
gene_expression_PDE3B_ovary_All_information <- read_csv("S4 drugs/gene_expression_PDE3B_ovary_All_information.csv")


S4_drugs=rbind(
  gene_expression_CA1_ovary_All_information, 
  gene_expression_TEX9_ovary_All_information, 
  gene_expression_HBB_ovary_All_information, 
  gene_expression_HPSE_ovary_All_information, 
  gene_expression_FCGR2A_ovary_All_information, 
  gene_expression_VPS8_ovary_All_information, 
  gene_expression_SLC25A39_ovary_All_information, 
  gene_expression_SLC25A37_ovary_All_information, 
  gene_expression_HBD_ovary_All_information, 
  gene_expression_AHR_ovary_All_information, 
  gene_expression_CD96_ovary_All_information, 
  gene_expression_ENTPD5_ovary_All_information, 
  gene_expression_S100A10_ovary_All_information, 
  gene_expression_LXN_ovary_All_information, 
  gene_expression_MXI1_ovary_All_information, 
  gene_expression_TMSB4XP8_ovary_All_information, 
  gene_expression_FECH_ovary_All_information, 
  gene_expression_HBM_ovary_All_information, 
  gene_expression_NUDT4_ovary_All_information, 
  gene_expression_IFI6_ovary_All_information, 
  gene_expression_MTRNR2L12_ovary_All_information, 
  gene_expression_MTRNR2L2_ovary_All_information, 
  gene_expression_SOS1_ovary_All_information, 
  gene_expression_MTRNR2L4_ovary_All_information, 
  gene_expression_IFI27_ovary_All_information, 
  gene_expression_MTRNR2L9_ovary_All_information, 
  gene_expression_MTRNR2L8_ovary_All_information, 
  gene_expression_HERC5_ovary_All_information, 
  gene_expression_HBG2_ovary_All_information, 
  gene_expression_NISCH_ovary_All_information, 
  gene_expression_TFRC_ovary_All_information, 
  gene_expression_MCM2_ovary_All_information, 
  gene_expression_RALGAPB_ovary_All_information, 
  gene_expression_QPCT_ovary_All_information, 
  gene_expression_GNS_ovary_All_information, 
  gene_expression_LTF_ovary_All_information, 
  gene_expression_DEFA3_ovary_All_information, 
  gene_expression_PDE4B_ovary_All_information, 
  gene_expression_DYSF_ovary_All_information, 
  gene_expression_RPL23AP82_ovary_All_information, 
  gene_expression_PPP2R1B_ovary_All_information, 
  gene_expression_PADI2_ovary_All_information, 
  gene_expression_RABL2A_ovary_All_information, 
  gene_expression_HSD17B3_ovary_All_information, 
  gene_expression_ZNF385D_ovary_All_information, 
  gene_expression_MPL_ovary_All_information, 
  gene_expression_PTK7_ovary_All_information, 
  gene_expression_ANKH_ovary_All_information, 
  gene_expression_CADM2_ovary_All_information, 
  gene_expression_NOMO1_ovary_All_information, 
  gene_expression_SLC25A43_ovary_All_information, 
  gene_expression_FCRL6_ovary_All_information, 
  gene_expression_TCF7_ovary_All_information, 
  gene_expression_DDX11L5_ovary_All_information, 
  gene_expression_EXOC3L2_ovary_All_information, 
  gene_expression_DDX11L10_ovary_All_information, 
  gene_expression_EML2_ovary_All_information, 
  gene_expression_FCMR_ovary_All_information, 
  gene_expression_C1orf87_ovary_All_information, 
  gene_expression_PGM2_ovary_All_information, 
  gene_expression_LINS1_ovary_All_information, 
  gene_expression_PDE3B_ovary_All_information
)

#One gene deprecated ENSEMBL ID missing ENSG00000237805 and does not have drug info on CREMMIST
#https://www.ensembl.org/Homo_sapiens/Gene/Idhistory?g=ENSG00000237805, date 19/8/2024
#One gene that does not have drug info on CREAMMIST ENSG00000267243 lncRNA
length(unique(S4_drugs$gene))

S4_drugs <- S4_drugs %>%
  mutate(Status = ifelse(row_number() <= 13900, "Up", "Down"))

S4_drugs_up=S4_drugs%>%
  filter(Status=="Up")

S4_drugs_up=S4_drugs_up%>%
  filter(correlation < -0.6)

S4_drugs_up$F="favorable for up"

S4_drugs_down=S4_drugs%>%
  filter(Status=="Down")

S4_drugs_down=S4_drugs_down%>%
  filter(correlation>0.6)

S4_drugs_down$F="favorable for down"

S4_drugs=rbind(S4_drugs_up, S4_drugs_down)

S4_drugs=S4_drugs%>%
  filter(pvalue<0.05)

S4_drugs=S4_drugs[,-c(3,4,6,7)]

S4_drugs$gene <- str_replace(S4_drugs$gene, "FCMR", "FAIM3")
S4_drugs$gene <- str_replace(S4_drugs$gene, "LINS1", "LINS")
length(unique(S4_drugs$gene))

length(unique(S1_drugs$drug_name))
length(unique(S2_drugs$drug_name))
length(unique(S3_drugs$drug_name))
length(unique(S4_drugs$drug_name))


S1_drugs$Stage="Stage1"
S2_drugs$Stage="Stage2"
S3_drugs$Stage="Stage3"
S4_drugs$Stage="Stage4"


drugs=rbind(S1_drugs, S2_drugs, S3_drugs, S4_drugs)

unique_degs=rbind(AUC_Stage1, AUC_Stage2, AUC_Stage3, AUC_Stage4)

length(unique(drugs$gene))
length(unique(drugs$drug_name))

write.csv(S1_drugs, "drugs1.csv", row.names = FALSE)
write.csv(S2_drugs, "drugs2.csv", row.names = FALSE)
write.csv(S3_drugs, "drugs3.csv", row.names = FALSE)
write.csv(S4_drugs, "drugs4.csv", row.names = FALSE)






