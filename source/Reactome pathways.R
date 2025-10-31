#Get working directory 

getwd()

#Set working directory to the location of Data through setwd() then navigate to reactomegsa folder through uncommenting the next line and replacing X with the location of Data folder
#setwd("X/Reactomegsa/Output data for reactomegsa")

# Loading libraries ----


library(dplyr)
library(ggplotify)
library(ggVennDiagram)
library(ggvenn)
library(magrittr)
library(openxlsx)
library(patchwork)
library(readr)
library(readxl)
library(ReactomeGSA)
library(tibble)
library(tidyr)
library(VennDetail)
library(writexl)


#### Reading the csv file obtained from Reactome API using the analysis token #### 

Pathways_Genes1 <- read_csv("Pathways_Genes1.csv")
Pathways_Genes2 <- read_csv("Pathways_Genes2.csv")
Pathways_Genes3 <- read_csv("Pathways_Genes3.csv")
Pathways_Genes4 <- read_csv("Pathways_Genes4.csv")


### Reading ReactomeGSA results for every Stage ####

Reactome_S1=read_excel("S1.xlsx")
colnames(Reactome_S1)=c("Name","Pathway Status", "FDR", "Pval", "NGenes", "avFC", "sig")
Reactome_S1=Reactome_S1%>%
  filter(FDR<0.05)

Reactome_S2=read_excel("S2.xlsx")
colnames(Reactome_S2)=c("Name","Pathway Status", "FDR", "Pval", "NGenes", "avFC", "sig")
Reactome_S2=Reactome_S2%>%
  filter(FDR<0.05)

Reactome_S3=read_excel("S3.xlsx")
colnames(Reactome_S3)=c("Name","Pathway Status", "FDR", "Pval", "NGenes", "avFC", "sig")
Reactome_S3=Reactome_S3%>%
  filter(FDR<0.05)

Reactome_S4=read_excel("S4.xlsx")
colnames(Reactome_S4)=c("Name","Pathway Status", "FDR", "Pval", "NGenes", "avFC", "sig")
Reactome_S4=Reactome_S4%>%
  filter(FDR<0.05)

Reactome_S1$Stage="Stage1"
Reactome_S2$Stage="Stage2"
Reactome_S3$Stage="Stage3"
Reactome_S4$Stage="Stage4"

Reactome_full=rbind(Reactome_S1, Reactome_S2, Reactome_S3, Reactome_S4)

length(unique(Reactome_full$Name)) ###Need to add this as supplementary

write.xlsx(Reactome_full, "Reactome full.xlsx")

#### Align pathways with their IDs and associated Genes ####

# Stage 1
Reactome_S1 <- Reactome_S1 %>%
  left_join(Pathways_Genes1 %>% select(`Pathway name`, `Submitted entities found`),
            by = c("Name" = "Pathway name"))

# Stage 2
Reactome_S2 <- Reactome_S2 %>%
  left_join(Pathways_Genes2 %>% select(`Pathway name`, `Submitted entities found`),
            by = c("Name" = "Pathway name")) 

# Stage 3
Reactome_S3 <- Reactome_S3 %>%
  left_join(Pathways_Genes3 %>% select(`Pathway name`, `Submitted entities found`),
            by = c("Name" = "Pathway name")) 

# Stage 4
Reactome_S4 <- Reactome_S4 %>%
  left_join(Pathways_Genes4 %>% select(`Pathway name`, `Submitted entities found`),
            by = c("Name" = "Pathway name")) 

#### Get AUC-Defined DEGs ####


AUC_defined_DEGs <- read_excel("C:/Users/shaim/Desktop/Data/Differential_exp_results.xlsx", sheet = 4)
AUC1=AUC_defined_DEGs%>%
  filter(Stage=="Stage1")

AUC2=AUC_defined_DEGs%>%
  filter(Stage=="Stage2")

AUC3=AUC_defined_DEGs%>%
  filter(Stage=="Stage3")

AUC4=AUC_defined_DEGs%>%
  filter(Stage=="Stage4")

#### Create temporary vectors with ensembl IDs of the AUC-Defined DEGs to use them in the filteration function ####

TMP1=AUC1$ensembl_gene_id
TMP2=AUC2$ensembl_gene_id
TMP3=AUC3$ensembl_gene_id
TMP4=AUC4$ensembl_gene_id

#### Venn diagram #####
R <-list('Stage 1'=Reactome_S1$Name,'Stage 2'=Reactome_S2$Name, 
         'Stage 3'=Reactome_S3$Name,'Stage 4'=Reactome_S4$Name) 

P=ggvenn(
  R,
  show_percentage = FALSE, 
  fill_color = c("yellow2", "dodgerblue", "green3", "red3"),  
  stroke_color = "black",
  stroke_size = 0.2,  
  set_name_color = "black",
  set_name_size = 0,
  text_color = "black",
  text_size = 4
) +
  ggtitle("a)") +  
  theme_void() +
  annotate("text", x = -2, y = 0.8, label = "Stage 1", size = 3) +
  annotate("text", x = -1,   y = 1.2, label = "Stage 2", size = 3) +
  annotate("text", x =  0.9, y = 1.2, label = "Stage 3", size = 3) +
  annotate("text", x =  1.65, y = 0.8, label = "Stage 4", size = 3) +
  theme(
    plot.title = element_text(hjust = 0.01, vjust = 0.7, size = 10),
    plot.margin = margin(1, 1, 1, 1)
  ) +
  coord_cartesian(clip = "off")


ggsave("SSPs.png", P,
       width = 1800, height = 900, units = "px", dpi = 300, bg = "white")

# Find unique pathways in Stage 1
unique_path_stage1 <- setdiff(R$`Stage 1`, union(union(R$`Stage 2`, R$`Stage 3`), R$`Stage 4`))
unique_path_stage1 <- Reactome_S1[Reactome_S1$Name %in% unique_path_stage1, ]
length(unique(unique_path_stage1$Name))

# Find unique pathways in Stage 2
unique_path_stage2 <- setdiff(R$`Stage 2`, union(union(R$`Stage 1`, R$`Stage 3`), R$`Stage 4`))
unique_path_stage2 <- Reactome_S2[Reactome_S2$Name%in% unique_path_stage2, ]
length(unique(unique_path_stage2$Name))

# Find unique pathways in Stage 3
unique_path_stage3 <- setdiff(R$`Stage 3`, union(union(R$`Stage 1`, R$`Stage 2`), R$`Stage 4`))
unique_path_stage3 <- Reactome_S3[Reactome_S3$Name%in% unique_path_stage3, ]
length(unique(unique_path_stage3$Name))

# Find unique pathways in Stage 4
unique_path_stage4 <- setdiff(R$`Stage 4`, union(union(R$`Stage 1`, R$`Stage 2`), R$`Stage 3`))
unique_path_stage4 <- Reactome_S4[Reactome_S4$Name %in% unique_path_stage4, ]
length(unique(unique_path_stage4$Name))

Unique_S1=unique_path_stage1
Unique_S2=unique_path_stage2
Unique_S3=unique_path_stage3
Unique_S4=unique_path_stage4

table(Unique_S1$`Pathway Status`)
table(Unique_S2$`Pathway Status`)
table(Unique_S3$`Pathway Status`)
table(Unique_S4$`Pathway Status`)


#### filter genes and align them with their associated pathways #####
#Stage1 
unique_path_stage1=unique_path_stage1%>%
  separate_rows(`Submitted entities found`, sep = ";")

unique_path_stage1=unique_path_stage1%>%
  filter(`Submitted entities found` %in% TMP1)

#Stage2

unique_path_stage2=unique_path_stage2%>%
  separate_rows(`Submitted entities found`, sep = ";")

unique_path_stage2=unique_path_stage2%>%
  filter(`Submitted entities found` %in% TMP2)

#Stage3
unique_path_stage3=unique_path_stage3%>%
  separate_rows(`Submitted entities found`, sep = ";")

unique_path_stage3=unique_path_stage3%>%
  filter(`Submitted entities found` %in% TMP3)

#Stage4

unique_path_stage4=unique_path_stage4%>%
  separate_rows(`Submitted entities found`, sep = ";")

unique_path_stage4=unique_path_stage4%>%
  filter(`Submitted entities found` %in% TMP4)


significant_pathways=rbind(Reactome_S1, Reactome_S2, Reactome_S3, Reactome_S4)
stage_specific_pathways=rbind(Unique_S1, Unique_S2, Unique_S3, Unique_S4)

pathways=createWorkbook()

addWorksheet(pathways, sheetName = "significant_pathways")
writeData(pathways, sheet = 1, x = significant_pathways)

addWorksheet(pathways, sheetName = "stage_specific_pathways")
writeData(pathways, sheet = 2, x = stage_specific_pathways)

##### Attach corresponding gene symbols and their status ####
#Stage 1
unique_path_stage1 <- unique_path_stage1 %>%
  left_join(AUC1 %>% select(ensembl_gene_id, hgnc_symbol),
            by = c("Submitted entities found" = "ensembl_gene_id"))

unique_path_stage1 <- unique_path_stage1 %>%
  left_join(AUC1 %>% select(ensembl_gene_id, Status),
            by = c("Submitted entities found" = "ensembl_gene_id"))

colnames(unique_path_stage1)[11]="Gene Status"


#Stage 2
unique_path_stage2 <- unique_path_stage2 %>%
  left_join(AUC2 %>% select(ensembl_gene_id, hgnc_symbol),
            by = c("Submitted entities found" = "ensembl_gene_id"))

unique_path_stage2 <- unique_path_stage2 %>%
  left_join(AUC2 %>% select(ensembl_gene_id, Status),
            by = c("Submitted entities found" = "ensembl_gene_id"))

colnames(unique_path_stage2)[11]="Gene Status"

#Stage 3
unique_path_stage3 <- unique_path_stage3 %>%
  left_join(AUC3 %>% select(ensembl_gene_id, hgnc_symbol),
            by = c("Submitted entities found" = "ensembl_gene_id"))

unique_path_stage3 <- unique_path_stage3 %>%
  left_join(AUC3 %>% select(ensembl_gene_id, Status),
            by = c("Submitted entities found" = "ensembl_gene_id"))

colnames(unique_path_stage3)[11]="Gene Status"


#Stage4
unique_path_stage4 <- unique_path_stage4 %>%
  left_join(AUC4 %>% select(ensembl_gene_id, hgnc_symbol),
            by = c("Submitted entities found" = "ensembl_gene_id"))

unique_path_stage4 <- unique_path_stage4 %>%
  left_join(AUC4 %>% select(ensembl_gene_id, Status),
            by = c("Submitted entities found" = "ensembl_gene_id"))

colnames(unique_path_stage4)[11]="Gene Status"


#### Venn diagram #####
V <-list('Stage 1'=unique_path_stage1$Name,'Stage 2'=unique_path_stage2$Name, 
         'Stage 3'=unique_path_stage3$Name,'Stage 4'=unique_path_stage4$Name) 


U=ggvenn(
  V,
  show_percentage = FALSE, 
  fill_color = c("darkorange", "mediumorchid", "deepskyblue", "limegreen"), 
  stroke_color = "black",
  stroke_size = 0.2,  
  set_name_color = "black",
  set_name_size = 0,
  text_color = "black",
  text_size = 4
) +
  ggtitle("b)") + 
  theme_void() +
  annotate("text", x = -1.7, y = 0.8, label = "Stage 1", size = 3) +
  annotate("text", x = -1,   y = 1.2, label = "Stage 2", size = 3) +
  annotate("text", x =  0.9, y = 1.2, label = "Stage 3", size = 3) +
  annotate("text", x =  1.9, y = 0.8, label = "Stage 4", size = 3) +
  theme(
    plot.title = element_text(hjust = 0.01, vjust = 0.7, size = 10),
    plot.margin = margin(1, 1, 1, 1)
  ) +
  coord_cartesian(clip = "off")


P_plot <- as.ggplot(P)
U_plot <- as.ggplot(U)


combined_plot <- P_plot + U_plot

ggsave("Pathways combined.png", combined_plot,
       width = 1800, height = 900, units = "px", dpi = 300, bg = "white")



stage_specific_path_DEGs=rbind(unique_path_stage1, unique_path_stage2, unique_path_stage3, unique_path_stage4)

addWorksheet(pathways, sheetName = "specific_path associated DEGs")
writeData(pathways, sheet = 3, x = stage_specific_path_DEGs)

saveWorkbook(pathways, file = "pathways.xlsx", overwrite = TRUE)

sessionInfo()

#End of code


