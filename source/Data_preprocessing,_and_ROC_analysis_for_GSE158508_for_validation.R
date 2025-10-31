#Get working directory 

getwd()

#Set working directory to the location of Data through setwd()

# Loading libraries ----

library(magrittr)
library(dplyr)
library(edgeR)
library(ggvenn)
library(ggVennDiagram)
library(EDASeq)
library(dplyr)
library(Biobase)
library(scater)
library(ggplot2)
library(tidyr)
library(RUVSeq)
library(SingleCellExperiment)
library(scater)
library(ggplot2)
library(foreach)
library(readr)
library(limma)
library(pROC)
library(ROCR)
library(metaSeq)
library(readxl)
library(openxlsx)

# Load GSE158508 data ----
load("dgeGenesEnsembl75.rdata")

ImPlatelet<-read.delim("ImPlatelet_samples.tsv",sep="\t")

ImPlatelet_OC=ImPlatelet%>%
  filter(OriginalGroup=="OC")

ImPlatelet_OC <- ImPlatelet_OC[-((nrow(ImPlatelet_OC) - 4):nrow(ImPlatelet_OC)), ]

#remove samples with missing stage and males

ImPlatelet_HC=ImPlatelet%>%
  filter(OriginalGroup=="HC")

ImPlatelet_HC=ImPlatelet%>%
  filter(Gender=="2")

ImPlatelet_samples=rbind(ImPlatelet_HC, ImPlatelet_OC)

# Exclude Y-chromosome linked genes ----

ImPlatelet_counts<-read.delim("ImPlatelet_counts_raw.tsv",sep="\t")

Y_chr=filter(genes, genes$chromosome_name=='Y')
Y_rows=rownames(Y_chr)

Y_genes <- ImPlatelet_counts[rownames(ImPlatelet_counts) %in% Y_rows, ]
Y_genes =rownames(Y_genes)

# Arrange samples and stages

ImPlatelet_samples$EarlyStage[is.na(ImPlatelet_samples$EarlyStage)] <- "Control"

ImPlatelet_samples$Stage[is.na(ImPlatelet_samples$Stage)] <- "Control"

ImPlatelet_samples$Stage[ImPlatelet_samples$Stage %in% c("IC", "IA")] <- "Stage1"

ImPlatelet_samples$Stage[ImPlatelet_samples$Stage %in% c("IIA", "IIB", "IIC")] <- "Stage2"

ImPlatelet_samples$Stage[ImPlatelet_samples$Stage %in% c("IIIA", "IIIB", "IIIC")] <- "Stage3"

ImPlatelet_samples$Stage[ImPlatelet_samples$Stage %in% c("IVB")] <- "Stage4"

ImPlatelet_samples=ImPlatelet_samples[, -c(8,9,10)]

selected_samples=ImPlatelet_samples$Id

ImPlatelet_counts <- ImPlatelet_counts[, selected_samples]

table(ImPlatelet_samples$Stage)

###Create DGE object ####

dge <- DGEList(counts = ImPlatelet_counts,
               group = ImPlatelet_samples$Stage,
               genes = genes[which(rownames(genes) %in% rownames(ImPlatelet_counts)),]
)           

ImPlatelet_samples=ImPlatelet_samples[,-c(2,3,4,9,6,7)]

dge$samples <- cbind(dge$samples, ImPlatelet_samples)

dge$samples$Stage<- factor(dge$samples$Stage, levels = c(levels(dge$samples$Stage),"Control","Stage1","Stage2","Stage3","Stage4"))

dge$samples$Stage[which(dge$samples$Stage %in% c("Control"))] <- "Control"

dge$samples$Stage[which(dge$samples$Stage %in% c("Stage1"))] <- "Stage1"

dge$samples$Stage[which(dge$samples$Stage %in% c("Stage2"))] <- "Stage2"

dge$samples$Stage[which(dge$samples$Stage %in% c("Stage3"))] <- "Stage3"

dge$samples$Stage[which(dge$samples$Stage %in% c("Stage4"))] <- "Stage4"

dim(dge)

summary(dge$samples$lib.size)
summary(dge$samples$group)
summary(dge$samples$Stage)

# Filter lowly expressed genes ----

keep.exprs <- filterByExpr(dge, group=dge$samples$Stage)
dge <- dge[keep.exprs,, keep.lib.sizes=FALSE]
dim(dge)
dge=calcNormFactors(dge)

#RLE plot before RUVg correction ----

sce1 <- SingleCellExperiment(assays = list(counts =dge$counts))

plotRLE(sce1, exprs_values = "counts", exprs_logged=FALSE, style = "minimal")+
  ggtitle("Before RUVg correction")

# Make RUVg correction ----

ord_count=dge$counts

filtered_cts <- ord_count

factor=as.factor(dge$samples$Stage)

set1 <- newSeqExpressionSet(as.matrix(ord_count),
                            phenoData = data.frame(factor, 
                                                   row.names=colnames(filtered_cts)))

design = model.matrix(~0+factor)
design=as.matrix(design)
design <- design[, c(2,3,4,5,1)]
y <- DGEList(counts=counts(set1), group=factor)
y <- calcNormFactors(y, method="TMM")
y <- estimateGLMCommonDisp(y, design)
y <- estimateGLMTagwiseDisp(y, design)
fit <- glmFit(y, design)
lrt <- glmLRT(fit, coef=ncol(fit$design))
top <- topTags(lrt, n=nrow(set1))$table
empirical <- rownames(set1)[which(!(rownames(set1) %in% rownames(top)[1:5000]))]

set2 <- RUVg(set1, empirical, k=3, isLog = FALSE)

W=pData(set2)

ruv_counts=as.data.frame(set2@assayData[["normalizedCounts"]])

rownames(W)=NULL

W=W[,-1]

W=as.matrix(W)

sce2 <- SingleCellExperiment(assays = list(counts = ruv_counts))

colData(sce2)$Stage <- dge$samples$Stage

plotRLE(sce2, exprs_values = "counts", exprs_logged = FALSE, style = "minimal") +
  ggtitle("After RUVg correction")


design <- model.matrix(~0+dge$samples$Stage)

colnames(design)[1] ="Control"
colnames(design)[2] ="Stage1"
colnames(design)[3] ="Stage2"
colnames(design)[4] ="Stage3"
colnames(design)[5] ="Stage4"

design <- design[, c(2,3,4,5,1)]

IMP_lcpm=cpm.DGEList(dge, log = TRUE)

dge$counts=IMP_lcpm

IMP_lcpm=removeBatchEffect(dge, design=design, covariates = as.matrix(W))

IMP_lcpm=as.data.frame(IMP_lcpm)

cpm=IMP_lcpm

AUC_original_dataset=read_excel("Differential_exp_results.xlsx", sheet=3)

UniqueS1 <- AUC_original_dataset%>%
  filter(Stage=="Stage1")

rownames(UniqueS1)=UniqueS1$ensembl_gene_id

UniqueS2 <- AUC_original_dataset%>%
  filter(Stage=="Stage2")

rownames(UniqueS2)=UniqueS2$ensembl_gene_id

UniqueS3 <- AUC_original_dataset%>%
  filter(Stage=="Stage3")

rownames(UniqueS3)=UniqueS3$ensembl_gene_id

UniqueS4 <- AUC_original_dataset%>%
  filter(Stage=="Stage4")

rownames(UniqueS4)=UniqueS4$ensembl_gene_id

####Stage 1 ----

Samples_Control=ImPlatelet_samples %>%
  filter(Stage == "Control")

S1_samples <- ImPlatelet_samples %>%
  filter(Stage == "Stage1")

S1_samples <- rbind(S1_samples, Samples_Control)

S1_samples <- S1_samples$Id
label1 <- c(rep(1, 4), rep(0, 123 - 4))

Stage1_genes <- UniqueS1$ensembl_gene_id
Stage1_counts <- cpm[Stage1_genes, ]
Stage1 <- Stage1_counts[, S1_samples]
Stage1 <- rbind(Stage1, label1)

rownames(Stage1)[nrow(Stage1)] <- "label"

Stage1 <- as.matrix(Stage1)

genes1 <- rownames(Stage1)[-nrow(Stage1)]

num_genes1 <- length(genes1)

roc_objects_stage1 <- vector("list", length = num_genes1)
auc_values_stage1 <- numeric(0)
genes_above_threshold_stage1 <- character(0)
for (i in 1:num_genes1) {
  gene_expression <- as.numeric(Stage1[i, -ncol(Stage1)])
  labels <- as.numeric(Stage1[nrow(Stage1), -ncol(Stage1)])
  model <- glm(labels ~ gene_expression, family = binomial(link = "logit"))
  predictions <- predict(model, type = "response")
  prediction_object <- prediction(predictions, labels)
  roc_obj <- roc(labels, predictions)
  roc_objects_stage1[[i]] <- roc_obj
  auc_values_stage1 <- c(auc_values_stage1, auc(roc_obj))
  genes_above_threshold_stage1 <- c(genes_above_threshold_stage1, genes1[i])
}

auc_stage1 <- data.frame(Gene = genes_above_threshold_stage1, AUC = auc_values_stage1)

for (i in seq_along(roc_objects_stage1)) {
  gene_name <- genes1[i]
  plot(roc_objects_stage1[[i]], main = paste("ROC Curve -", gene_name))
  text(0.5, 0.4, paste("AUC =", round(auc(roc_objects_stage1[[i]]), 2)), adj = c(0.5, 0.5), cex = 1)
}

auc_stage1 <- auc_stage1 %>% arrange(desc(AUC))

auc_stage1 <- auc_stage1 %>%
  left_join(UniqueS1 %>% dplyr::select(ensembl_gene_id, Status), 
            by = c("Gene" = "ensembl_gene_id"))

auc_stage1 <- auc_stage1 %>%
  left_join(UniqueS1 %>% dplyr::select(ensembl_gene_id, Gene), 
            by = c("Gene" = "ensembl_gene_id"))

colnames(auc_stage1)[4]="hgnc_symbol"


####Stage 2 ----

S2_samples <- ImPlatelet_samples %>%
  filter(Stage == "Stage2")

S2_samples <- rbind(S2_samples, Samples_Control)

S2_samples <- S2_samples$Id
label2 <- c(rep(1, 5), rep(0, 124 - 5))

Stage2_genes <- UniqueS2$ensembl_gene_id
Stage2_counts <- cpm[Stage2_genes, ]
Stage2 <- Stage2_counts[, S2_samples]
Stage2 <- rbind(Stage2, label2)

rownames(Stage2)[nrow(Stage2)] <- "label"

Stage2 <- as.matrix(Stage2)

genes2 <- rownames(Stage2)[-nrow(Stage2)]

num_genes2 <- length(genes2)


roc_objects_stage2 <- vector("list", length = num_genes2)
auc_values_stage2 <- numeric(0)
genes_above_threshold_stage2 <- character(0)
for (i in 1:num_genes2) {
  gene_expression <- as.numeric(Stage2[i, -ncol(Stage2)])
  labels <- as.numeric(Stage2[nrow(Stage2), -ncol(Stage2)])
  model <- glm(labels ~ gene_expression, family = binomial(link = "logit"))
  predictions <- predict(model, type = "response")
  prediction_object <- prediction(predictions, labels)
  roc_obj <- roc(labels, predictions)
  roc_objects_stage2[[i]] <- roc_obj
  auc_values_stage2 <- c(auc_values_stage2, auc(roc_obj))
  genes_above_threshold_stage2 <- c(genes_above_threshold_stage2, genes2[i])
}

auc_stage2 <- data.frame(Gene = genes_above_threshold_stage2, AUC = auc_values_stage2)

for (i in seq_along(roc_objects_stage2)) {
  gene_name <- genes2[i]
  plot(roc_objects_stage2[[i]], main = paste("ROC Curve -", gene_name))
  text(0.5, 0.4, paste("AUC =", round(auc(roc_objects_stage2[[i]]), 2)), adj = c(0.5, 0.5), cex = 1)
}

auc_stage2 <- auc_stage2 %>% arrange(desc(AUC))


auc_stage2 <- auc_stage2 %>%
  left_join(UniqueS2 %>% dplyr::select(ensembl_gene_id, Status), 
            by = c("Gene" = "ensembl_gene_id"))

auc_stage2 <- auc_stage2 %>%
  left_join(UniqueS2 %>% dplyr::select(ensembl_gene_id, Gene), 
            by = c("Gene" = "ensembl_gene_id"))

colnames(auc_stage2)[4]="hgnc_symbol"

####Stage 3 ----

S3_samples <- ImPlatelet_samples %>%
  filter(Stage == "Stage3")

S3_samples <- rbind(S3_samples, Samples_Control)

S3_samples <- S3_samples$Id
label3 <- c(rep(1, 14), rep(0, 133 - 14))


Stage3_genes <- UniqueS3$ensembl_gene_id
Stage3_counts <- cpm[Stage3_genes, ]
Stage3 <- Stage3_counts[, S3_samples]
Stage3 <- rbind(Stage3, label3)

rownames(Stage3)[nrow(Stage3)] <- "label"

Stage3 <- as.matrix(Stage3)

genes3 <- rownames(Stage3)[-nrow(Stage3)]

num_genes3 <- length(genes3)

roc_objects_stage3 <- vector("list", length = num_genes3)
auc_values_stage3 <- numeric(0)
genes_above_threshold_stage3 <- character(0)
for (i in 1:num_genes3) {
  gene_expression <- as.numeric(Stage3[i, -ncol(Stage3)])
  labels <- as.numeric(Stage3[nrow(Stage3), -ncol(Stage3)])
  model <- glm(labels ~ gene_expression, family = binomial(link = "logit"))
  predictions <- predict(model, type = "response")
  prediction_object <- prediction(predictions, labels)
  roc_obj <- roc(labels, predictions)
  roc_objects_stage3[[i]] <- roc_obj
  auc_values_stage3 <- c(auc_values_stage3, auc(roc_obj))
  genes_above_threshold_stage3 <- c(genes_above_threshold_stage3, genes3[i])
}

auc_stage3 <- data.frame(Gene = genes_above_threshold_stage3, AUC = auc_values_stage3)

for (i in seq_along(roc_objects_stage3)) {
  gene_name <- genes3[i]
  plot(roc_objects_stage3[[i]], main = paste("ROC Curve -", gene_name))
  text(0.5, 0.4, paste("AUC =", round(auc(roc_objects_stage3[[i]]), 2)), adj = c(0.5, 0.5), cex = 1)
}

auc_stage3 <- auc_stage3 %>% arrange(desc(AUC))

auc_stage3 <- auc_stage3 %>%
  left_join(UniqueS3 %>% dplyr::select(ensembl_gene_id, Status), 
            by = c("Gene" = "ensembl_gene_id"))

auc_stage3 <- auc_stage3 %>%
  left_join(UniqueS3 %>% dplyr::select(ensembl_gene_id, Gene), 
            by = c("Gene" = "ensembl_gene_id"))

colnames(auc_stage3)[4]="hgnc_symbol"

####Stage 4 ----

S4_samples <- ImPlatelet_samples %>%
  filter(Stage == "Stage4")

S4_samples <- rbind(S4_samples, Samples_Control)

S4_samples <- S4_samples$Id
label4 <- c(rep(1, 2), rep(0, 121 - 2))


Stage4_genes <- UniqueS4$ensembl_gene_id
Stage4_counts <- cpm[Stage4_genes, ]
Stage4 <- Stage4_counts[, S4_samples]
Stage4 <- rbind(Stage4, label4)

rownames(Stage4)[nrow(Stage4)] <- "label"

Stage4 <- as.matrix(Stage4)

genes4 <- rownames(Stage4)[-nrow(Stage4)]

num_genes4 <- length(genes4)

roc_objects_stage4 <- vector("list", length = num_genes4)
auc_values_stage4 <- numeric(0)
genes_above_threshold_stage4 <- character(0)
for (i in 1:num_genes4) {
  gene_expression <- as.numeric(Stage4[i, -ncol(Stage4)])
  labels <- as.numeric(Stage4[nrow(Stage4), -ncol(Stage4)])
  model <- glm(labels ~ gene_expression, family = binomial(link = "logit"))
  predictions <- predict(model, type = "response")
  prediction_object <- prediction(predictions, labels)
  roc_obj <- roc(labels, predictions)
  roc_objects_stage4[[i]] <- roc_obj
  auc_values_stage4 <- c(auc_values_stage4, auc(roc_obj))
  genes_above_threshold_stage4 <- c(genes_above_threshold_stage4, genes4[i])
}

auc_stage4 <- data.frame(Gene = genes_above_threshold_stage4, AUC = auc_values_stage4)

for (i in seq_along(roc_objects_stage4)) {
  gene_name <- genes4[i]
  plot(roc_objects_stage4[[i]], main = paste("ROC Curve -", gene_name))
  text(0.5, 0.4, paste("AUC =", round(auc(roc_objects_stage4[[i]]), 2)), adj = c(0.5, 0.5), cex = 1)
}

auc_stage4 <- auc_stage4 %>% arrange(desc(AUC))


auc_stage4 <- auc_stage4 %>%
  left_join(UniqueS4 %>% dplyr::select(ensembl_gene_id, Status), 
            by = c("Gene" = "ensembl_gene_id"))

auc_stage4 <- auc_stage4 %>%
  left_join(UniqueS4 %>% dplyr::select(ensembl_gene_id, Gene), 
            by = c("Gene" = "ensembl_gene_id"))

colnames(auc_stage4)[4]="hgnc_symbol"

####AUC comparison ----

auc_stage1_original=UniqueS1

auc_stage2_original=UniqueS2

auc_stage3_original=UniqueS3

auc_stage4_original=UniqueS4

colnames(auc_stage1_original)[2]="AUC original dataset"
colnames(auc_stage2_original)[2]="AUC original dataset"
colnames(auc_stage3_original)[2]="AUC original dataset"
colnames(auc_stage4_original)[2]="AUC original dataset"

colnames(auc_stage4)[2]="AUC validation dataset"
colnames(auc_stage3)[2]="AUC validation dataset"
colnames(auc_stage2)[2]="AUC validation dataset"
colnames(auc_stage1)[2]="AUC validation dataset"


AUC1=auc_stage1_original$ensembl_gene_id
AUC2=auc_stage2_original$ensembl_gene_id
AUC3=auc_stage3_original$ensembl_gene_id
AUC4=auc_stage4_original$ensembl_gene_id

auc_stage1 <- auc_stage1[match(AUC1, auc_stage1$Gene), ]
auc_stage2 <- auc_stage2[match(AUC2, auc_stage2$Gene), ]
auc_stage3 <- auc_stage3[match(AUC3, auc_stage3$Gene), ]
auc_stage4 <- auc_stage4[match(AUC4, auc_stage4$Gene), ]

####Stage 1----

Stage1_AUC=cbind(auc_stage1_original, auc_stage1)
Stage1_AUC=Stage1_AUC[,-c(1,8,10)]

Stage1_AUC <- Stage1_AUC %>%
  mutate(`AUC in both > 0.6` = ifelse(
    `AUC original dataset` > 0.6 & `AUC validation dataset` > 0.6,
    TRUE,
    FALSE
  ))

Stage1_AUC <- Stage1_AUC%>%
  filter(`AUC in both > 0.6`=="TRUE")

table(Stage1_AUC$Status)

####Stage 2----

Stage2_AUC <- cbind(auc_stage2_original, auc_stage2)
Stage2_AUC=Stage2_AUC[,-c(1,8,10)]


Stage2_AUC <- Stage2_AUC %>%
  mutate(`AUC in both > 0.6` = ifelse(
    `AUC original dataset` > 0.6 & `AUC validation dataset` > 0.6,
    TRUE,
    FALSE
  ))


table(Stage2_AUC$`AUC in both > 0.6`)

Stage2_AUC <- Stage2_AUC%>%
  filter(`AUC in both > 0.6`=="TRUE")

table(Stage2_AUC$Status)

####Stage 3----

Stage3_AUC <- cbind(auc_stage3_original, auc_stage3)
Stage3_AUC=Stage3_AUC[,-c(1,8,10)]


Stage3_AUC <- Stage3_AUC %>%
  mutate(`AUC in both > 0.6` = ifelse(
    `AUC original dataset` > 0.6 & `AUC validation dataset` > 0.6,
    TRUE,
    FALSE
  ))

table(Stage3_AUC$`AUC in both > 0.6`)

Stage3_AUC <- Stage3_AUC%>%
  filter(`AUC in both > 0.6`=="TRUE")

table(Stage3_AUC$Status)

####Stage 4----

Stage4_AUC <- cbind(auc_stage4_original, auc_stage4)
Stage4_AUC=Stage4_AUC[,-c(1,8,10)]

Stage4_AUC <- Stage4_AUC %>%
  mutate(`AUC in both > 0.6` = ifelse(
    `AUC original dataset` > 0.6 & `AUC validation dataset` > 0.6,
    TRUE,
    FALSE
  ))


table(Stage4_AUC$`AUC in both > 0.6`)

Stage4_AUC <- Stage4_AUC%>%
  filter(`AUC in both > 0.6`=="TRUE")

table(Stage4_AUC$Status)

AUC_validation_dataset=rbind(Stage1_AUC, Stage2_AUC, Stage3_AUC, Stage4_AUC)

AUC_validation_dataset=AUC_validation_dataset[c(3,8,4,5,6,2,1,7,9)]


Differential_exp_result=loadWorkbook("Differential_exp_results.xlsx")

addWorksheet(Differential_exp_result, "AUC_validation_dataset")
writeData(Differential_exp_result, "AUC_validation_dataset", AUC_validation_dataset)

saveWorkbook(Differential_exp_result, 
             file = "Differential_exp_results.xlsx", overwrite = TRUE)

#### VennDiagram ####

R <-list('Stage 1'=Stage1_AUC$ensembl_gene_id,'Stage 2'=Stage2_AUC$ensembl_gene_id, 
         'Stage 3'=Stage3_AUC$ensembl_gene_id,'Stage 4'=Stage4_AUC$ensembl_gene_id) 


ggvenn(
  R,
  show_percentage = FALSE,
  fill_color = c("#E69F00", "#56B4E9", "#009E73", "red"),
  stroke_color = "black",
  stroke_size = 0.5,  
  set_name_color = "black",
  set_name_size = 5,
  text_color = "black",
  text_size = 4
) +
  theme_void() + 
  theme(
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5, size = 12, face = "bold", vjust = -1)
  ) +
  ggtitle("") +  
  annotate("text", x = 0, y = -1.9, label = "Stage-specific DEGs with ROC-AUC>0.6", 
           size = 4.5, fontface = "bold")

sessionInfo()

#End of code

