#Get working directory 

getwd()

#Set working directory to the location of Data through setwd()

# Loading libraries ----
library(Biobase)
library(cowplot)
library(dplyr)
library(EDASeq)
library(edgeR)
library(foreach)
library(Glimma)
library(ggplot2)
library(ggvenn)
library(grid)
library(limma)
library(magrittr)
library(metaSeq)
library(openxlsx)
library(patchwork)
library(pROC)
library(plotly)
library(plyr)
library(readr)
library(RUVSeq)
library(scater)
library(SingleCellExperiment)
library(tidyselect)
library(tidyverse)
library(tidyr)


# Loading data ----

url <- "https://raw.githubusercontent.com/MyronBest/InTVeld_Pancancer_TSOO/main/TableS2.csv"
data <- read_csv(url)
sampleInfo=as.data.frame(data)
sampleInfo <- sampleInfo[, -1]



Samples_Control=filter(sampleInfo, sampleInfo$Group=="Asymptomatic controls")
Samples_Control=filter(Samples_Control, Samples_Control$Sex=="F")

Samples_Control$Stage <- replace(Samples_Control$Stage, Samples_Control$Stage == "n.a.", "Control")

Samples_Cancer=filter(sampleInfo, sampleInfo$Group=="Ovarian cancer")
Samples_Cancer <- replace(Samples_Cancer, Samples_Cancer == "I", "Stage1")
Samples_Cancer <- replace(Samples_Cancer, Samples_Cancer == "II", "Stage2")
Samples_Cancer <- replace(Samples_Cancer, Samples_Cancer == "III", "Stage3")
Samples_Cancer <- replace(Samples_Cancer, Samples_Cancer == "IV", "Stage4")
Samples_Cancer <- replace(Samples_Cancer, Samples_Cancer == "n.a.", "Unknown Stage")

# Exclude samples with unknown stages ----
unknown_stage=Samples_Cancer%>%
  filter(Stage=="Unknown Stage")
unknown_stage=unknown_stage$Sample.ID

Samples_Cancer <- Samples_Cancer[!(Samples_Cancer$Stage == "Unknown Stage"), ]

load("GSE183635_TEP_Count_Matrix.rdata")

CountMatrix=as.data.frame(TEP_Count_Matrix)
colnames(CountMatrix) <- sub("-", "#", colnames(CountMatrix))
colnames(CountMatrix) <- sub(".*#", "", colnames(CountMatrix))

CountMatrix <- CountMatrix[, !(names(CountMatrix) %in% unknown_stage)]

sampleInfo=rbind(Samples_Control, Samples_Cancer)

table(sampleInfo$Group)

table(sampleInfo$Stage)


row_names=sampleInfo$Sample.ID



Counts=CountMatrix %>%
  select(any_of(row_names))

load("dgeGenesEnsembl75.rdata")

# Exclude Y-chr genes ----

Y_chr=filter(genes, genes$chromosome_name=='Y')
Y_rows=rownames(Y_chr)

Y_genes <- Counts[rownames(Counts) %in% Y_rows, ]
Y_genes =rownames(Y_genes)

Counts=Counts[!(row.names(Counts) %in% Y_genes),]

# Order counts by institute ----

Inst1=sampleInfo%>%
  filter(Sample.supplying.institution=="Institute 1")%>%
  pull(Sample.ID)

Inst2 <- sampleInfo %>%
  filter(Sample.supplying.institution == "Institute 2") %>%
  pull(Sample.ID)

Inst3 <- sampleInfo %>%
  filter(Sample.supplying.institution == "Institute 3") %>%
  pull(Sample.ID)

Inst4 <- sampleInfo %>%
  filter(Sample.supplying.institution == "Institute 4") %>%
  pull(Sample.ID)

Inst5 <- sampleInfo %>%
  filter(Sample.supplying.institution == "Institute 5") %>%
  pull(Sample.ID)

Inst10 <- sampleInfo %>%
  filter(Sample.supplying.institution == "Institute 10") %>%
  pull(Sample.ID)

Inst12 <- sampleInfo %>%
  filter(Sample.supplying.institution == "Institute 12") %>%
  pull(Sample.ID)

Inst13 <- sampleInfo %>%
  filter(Sample.supplying.institution == "Institute 13") %>%
  pull(Sample.ID)

combined_ids <- c(Inst1, Inst2, Inst3, Inst4, Inst5, Inst10, Inst12, Inst13)


sampleInfo <- sampleInfo %>%
  arrange(match(Sample.ID, combined_ids))


ord_count=Counts %>%
  select(any_of(combined_ids))

#### Extract Control and Stage 1 samples to be used later for ReactomeGSA####

Control_S1 <- sampleInfo %>%
  filter(Stage %in% c("Control", "Stage1"))

write.csv(Control_S1, "Control_S1.csv")

Control_S1=Control_S1$Sample.ID

Control_S1_counts=ord_count[,Control_S1]

Control_S1_counts=as.data.frame(Control_S1_counts)

write.csv(Control_S1_counts, "Control_S1_counts.csv")

#### Extract Control and Stage 2 samples to be used later for ReactomeGSA####

Control_S2 <- sampleInfo %>%
  filter(Stage %in% c("Control", "Stage2"))

write.csv(Control_S2, "Control_S2.csv")

Control_S2=Control_S2$Sample.ID

Control_S2_counts=ord_count[,Control_S2]

Control_S2_counts=as.data.frame(Control_S2_counts)

write.csv(Control_S2_counts, "Control_S2_counts.csv")

#### Extract Control and Stage 3 samples to be used later for ReactomeGSA####

Control_S3 <- sampleInfo %>%
  filter(Stage %in% c("Control", "Stage3"))

write.csv(Control_S3, "Control_S3.csv")

Control_S3=Control_S3$Sample.ID

Control_S3_counts=ord_count[,Control_S3]

Control_S3_counts=as.data.frame(Control_S3_counts)

write.csv(Control_S3_counts, "Control_S3_counts.csv")

#### Extract Control and Stage 4 samples to be used later for ReactomeGSA####

Control_S4 <- sampleInfo %>%
  filter(Stage %in% c("Control", "Stage4"))

write.csv(Control_S4, "Control_S4.csv")

Control_S4=Control_S4$Sample.ID

Control_S4_counts=ord_count[,Control_S4]

Control_S4_counts=as.data.frame(Control_S4_counts)

write.csv(Control_S4_counts, "Control_S4_counts.csv")


# Create a DGE object ----

dge <- DGEList(counts = ord_count,
               group = sampleInfo$Stage,
               genes = genes[which(rownames(genes) %in% rownames(ord_count)),]
)           

#dge$samples$lib.size <- sampleInfo$lib.size

sampleInfo=sampleInfo[,-c(2,7,8,9,10,11)]

dge$samples <- cbind(dge$samples, sampleInfo)

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

sce1 <- SingleCellExperiment(assays = list(counts =ord_count))

colData(sce1)$Institute=sampleInfo$Sample.supplying.institution
colData(sce1)

before=plotRLE(sce1, exprs_values = "counts", exprs_logged = FALSE, color_by = colData(sce1), style = "minimal") +
  ggtitle("A)                       
                                   Before batch correction - GSE183635")+
  theme(
    plot.title = element_text(hjust = 0.001, vjust = -0.5, size = 5), 
    plot.margin = margin(t = 100, r = 100, b = 100, l = 100), 
    axis.title.x = element_blank(), 
    axis.title.y = element_text(size=5)
  )

# Make RUVg correction ----

ord_count=dge$counts

filtered_cts <- ord_count

factor=as.factor(sampleInfo$Stage)
batch=as.factor(sampleInfo$Sample.supplying.institution)

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
empirical <- rownames(set1)[which(!(rownames(set1) %in% rownames(top)[1:1000]))]


set2 <- RUVg(set1, empirical, k=3, isLog = FALSE)

W=pData(set2)

ruv_counts=as.data.frame(set2@assayData[["normalizedCounts"]])

rownames(W)=NULL

W=W[,-1]

W=as.matrix(W)


ruv_counts=ruv_counts %>%
  select(any_of(combined_ids))

#https://bioconductor.org/packages/release/bioc/vignettes/RUVSeq/inst/doc/RUVSeq.html
#https://support.bioconductor.org/p/67791/#9155690


sce2 <- SingleCellExperiment(assays = list(counts = ruv_counts))

colData(sce2)$Institute <- sampleInfo$Sample.supplying.institution
colData(sce2)$Stage <- sampleInfo$Stage


before <- plotRLE(sce1, exprs_values = "counts", exprs_logged = FALSE, 
                  color_by = "Institute", style = "minimal") +
  ggtitle("A) Before batch correction - GSE183635") +
  theme(
    plot.title = element_text(hjust = 0, vjust = 0, size = 8),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 6),
    plot.margin = margin(5, 5, 5, 5)
  )

after <- plotRLE(sce2, exprs_values = "counts", exprs_logged = FALSE, 
                 color_by = "Institute", style = "minimal") +
  ggtitle("B) After batch correction - GSE183635") +
  theme(
    plot.title = element_text(hjust = 0, vjust = 0, size = 8),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 6),
    plot.margin = margin(5, 5, 5, 5)
  )

combined_plot <- (before / after) +
  plot_layout(guides = "collect", heights = c(1, 1)) &
  theme(
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 6),
    legend.key.size = unit(0.4, "cm"),
    legend.spacing.y = unit(0.2, "cm"),
    legend.direction = "vertical",
    legend.justification = "left",
    legend.position = "right", ,
    plot.margin = margin(1, 1, 1, 1)
  )

ggsave("combined_plot.png", combined_plot,
       width = 1800, height = 900, dpi = 300, units = "px", bg = "white")


sce2 <- runPCA(sce2, ncomponents = 3, exprs_values = "counts")


plotReducedDim(sce2, dimred = "PCA", colour_by = "Stage") +
  ggtitle("PCA with RUVg correction")


pca_data <- reducedDim(sce2, "PCA")
pca_data <- as.data.frame(pca_data)
pca_data$Stage <- colData(sce2)$Stage
pca_data$Institute <- colData(sce2)$Institute


fig <- plot_ly(pca_data, x = ~PC1, y = ~PC2, z = ~PC3,
               color = ~Stage, colors = c("Control" = "blue", "Stage1" = "red", "Stage2" = "green", "Stage3" = "purple", "Stage4" = "orange"),
               marker = list(size = 5),
               text = ~paste("Sample:", rownames(pca_data), "<br>Institution:", Institute)) %>%
  add_markers() %>%
  layout(scene = list(xaxis = list(title = "PC1"),
                      yaxis = list(title = "PC2"),
                      zaxis = list(title = "PC3")),
         title = "3D PCA of Samples with RUVg Correction")

fig


design <- model.matrix(~0+dge$samples$Stage + W)

colnames(design)[1] ="Control"
colnames(design)[2] ="Stage1"
colnames(design)[3] ="Stage2"
colnames(design)[4] ="Stage3"
colnames(design)[5] ="Stage4"
colnames(design)[6] ="W_1"
colnames(design)[7] ="W_2"
colnames(design)[8] ="W_3"


design <- design[, c(2,3,4,5,1,6,7,8)]


v=voom(dge, design = design, plot = TRUE) #transform counts to logcpm
vfit <- lmFit(v, design) #fit linear model for each gene
contrast.matrix <- makeContrasts(Stage1vsControl=Stage1-Control,
                                 Stage2vsControl=Stage2-Control,
                                 Stage3vsControl=Stage3-Control,
                                 Stage4vsControl=Stage4-Control, 
                                 levels = colnames(design))

vfit <- contrasts.fit(vfit, contrasts=contrast.matrix)
efit <- eBayes(vfit) #calculate the linear model statistics
plotSA(efit,main = "Final model: Mean-variance trend")
# Extract limma results ----
summary(decideTests(efit, method = "separate", adjust.method = "BH", lfc=1)) 
#"separate" applies multiple testing adjustments separately to 
#each column of a matrix of p-values or t-statistics, 
#effectively treating each contrast independently. This is because the same 
#t-statistic cutoff cannot correspond to statistical significance 
#for different contrasts (stages).
finalDEGs=(decideTests(efit, method = "separate", adjust.method = "BH", lfc=1))

#https://support.bioconductor.org/p/9156456/

vennDiagram(finalDEGs, 
            include=c("up", "down"),
            counts.col=c("red", "blue"),
            circle.col = c("red", "blue", "green3"))+title(main = "DEGs Ovarian Cancer")


#Standard RNAseq analysis approaches like DESeq2 and limma typically rely on 
#the Benjamini–Hochberg (BH) procedure to control the FDR 
#A typical RNAseq analysis involves applying a differential expression 
#analysis method to an RNAseq count matrix, where rows are genes, columns are 
#replicate samples for each condition and the values in the matrix are read counts, 
#generating a gene-p-value matrix (Fig. 1B, setting 1). 
#Often, an FDR correction (e.g. BH) is then applied to this 
#gene-p-value matrix to generate multiplicity-adjusted p-values, 
#after which a research decision may be made 
#https://academic.oup.com/bioinformatics/article/39/1/btac718/6795009?login=false



#A p-value and FDR-value <0.05 was considered as statistically significant.
#https://www.cell.com/cancer-cell/fulltext/S1535-6108(22)00370-1#sectitle0030


#false discovery rate (FDR) <0.05 and an absolute log2 fold change >1 were set 
#as the threshold for a significant differential expression.
#https://www.frontiersin.org/journals/oncology/articles/10.3389/fonc.2022.844520/full


#The genes with adjusted p value (padj)<0.05 
#and fold change |log2FC|≥1 were screened as significant DEGs.
#https://link.springer.com/article/10.1007/s12033-022-00611-z


#Log2 fold change of 1.0 and a cutoff P value of 0.05 were used
#https://www.tandfonline.com/doi/full/10.1080/09537104.2023.2212071


result=summary(decideTests(efit, method = "separate", adjust.method = "BH", lfc=1))
result=as.data.frame(result)

colnames(result)[1] ="Status"
colnames(result)[2] ="Stage"
colnames(result)[3] ="No. of DEGs"
result=result[-c(2, 5, 8, 11),]

# Stage 1 DEGs ----

summary(decideTests(efit, method = "separate", adjust.method = "BH", lfc = 1))
finalDEGs=decideTests(efit, method = "separate", adjust.method = "BH", lfc = 1)

vennDiagram(finalDEGs,
             include=c("up"), counts.col= "black", circle.col = c("coral", "blue", "deeppink3", "steelblue3"))


S1=topTable(efit, coef = "Stage1vsControl", number = 5000, genelist = efit$genes,
            adjust.method = "BH", sort.by = "AveExpr", p.value = 0.05, lfc=1)

S1['Stage']="Stage 1"
S1['Status']=""
S1$Status <- ifelse(S1$logFC > 1, "Up", "Down")
S1_UP=S1%>%
  filter(Status=="Up")

S1_down=S1%>%
  filter(Status=="Down")


summary(decideTests(efit, method = "separate", adjust.method = "BH", p.value= 0.05,lfc = 1))

# Stage 2 DEGs ----

S2=topTable(efit, coef = "Stage2vsControl", number = 5000, genelist = efit$genes,
            adjust.method = "BH", sort.by = "AveExpr",
            p.value = 0.05, lfc = 1)

S2['Stage']="Stage 2"
S2['Status']=""
S2$Status <- ifelse(S2$logFC > 1, "Up", "Down")
S2_UP=S2%>%
  filter(Status=="Up")

S2_down=S2%>%
  filter(Status=="Down")

# Stage 3 DEGs ----

S3=topTable(efit, coef = "Stage3vsControl", number = 5000, genelist = efit$genes,
            adjust.method = "BH", sort.by = "AveExpr",
            p.value = 0.05, lfc = 1)

S3['Stage']="Stage 3"
S3['Status']=""
S3$Status <- ifelse(S3$logFC > 1, "Up", "Down")
S3_UP=S3%>%
  filter(Status=="Up")

S3_down=S3%>%
  filter(Status=="Down")


# Stage 4 DEGs ----

S4=topTable(efit, coef = "Stage4vsControl", number = 5000, genelist = efit$genes,
            adjust.method = "BH", sort.by = "AveExpr",
            p.value = 0.05, lfc = 1)
S4['Stage']="Stage 4"
S4['Status']=""
S4$Status <- ifelse(S4$logFC > 1, "Up", "Down")
S4_UP=S4%>%
  filter(Status=="Up")

S4_down=S4%>%
  filter(Status=="Down")



Non_stage_specific_DEGs=rbind(S1, S2, S3, S4)

vlist <-list('Stage 1'=S1$ensembl_gene_id,'Stage 2'=S2$ensembl_gene_id, 
             'Stage 3'=S3$ensembl_gene_id,'Stage 4'=S4$ensembl_gene_id) 


p <- ggvenn(
  vlist,
  show_percentage = FALSE, 
  fill_color = c("coral", "blue", "deeppink3", "steelblue3"),
  stroke_color = "black",
  stroke_size = 0.1,  
  set_name_size = 0,  
  text_color = "black",
  text_size = 3
) +
  ggtitle("a)") +  # Add title here
  theme_void() +
  annotate("text", x = -2.1, y = 0.6, label = "Stage 1", size = 2) +
  annotate("text", x = -1,   y = 1.1, label = "Stage 2", size = 2) +
  annotate("text", x =  0.9, y = 1.1, label = "Stage 3", size = 2) +
  annotate("text", x =  1.9, y = 0.7, label = "Stage 4", size = 2) +
  theme(
    plot.title = element_text(hjust = 0, vjust = 1, size = 10),
    plot.margin = margin(1, 1, 1, 1)
  ) +
  coord_cartesian(clip = "off")

ggsave("Figure 2A.png", plot = p,
       width = 1200, height = 600, units = "px", dpi = 300,
       bg = "white")



# Extract Unique DEGs ----

dgs=data.frame(finalDEGs)

# Extract Stage 1 Unique DEGs ----

UniqueS1=dgs[dgs$Stage1vsControl !=0 & dgs$Stage2vsControl ==0 & dgs$Stage3vsControl ==0 & dgs$Stage4vsControl ==0, ]
UniqueS1=rownames(UniqueS1)
UniqueS1=as.data.frame(UniqueS1)
rownames(UniqueS1)=UniqueS1$UniqueS1
colnames(UniqueS1)[1] = "ensembl_gene_id"

UniqueS1=match_df(S1, UniqueS1, on="ensembl_gene_id")

# Extract Stage 2 Unique DEGs ----

UniqueS2=dgs[dgs$Stage1vsControl ==0 & dgs$Stage2vsControl !=0 & dgs$Stage3vsControl ==0 & dgs$Stage4vsControl ==0, ]
UniqueS2=rownames(UniqueS2)
UniqueS2=as.data.frame(UniqueS2)
rownames(UniqueS2)=UniqueS2$UniqueS2
colnames(UniqueS2)[1] = "ensembl_gene_id"

UniqueS2=match_df(S2, UniqueS2, on="ensembl_gene_id")

# Extract Stage 3 Unique DEGs ----

UniqueS3=dgs[dgs$Stage1vsControl ==0 & dgs$Stage2vsControl ==0 & dgs$Stage3vsControl !=0 & dgs$Stage4vsControl ==0, ]
UniqueS3=rownames(UniqueS3)
UniqueS3=as.data.frame(UniqueS3)
rownames(UniqueS3)=UniqueS3$UniqueS3
colnames(UniqueS3)[1] = "ensembl_gene_id"

UniqueS3=match_df(S3, UniqueS3, on="ensembl_gene_id")

# Extract Stage 4 Unique DEGs ----

UniqueS4 <- dgs[dgs$Stage1vsControl == 0 & dgs$Stage2vsControl == 0 & dgs$Stage3vsControl == 0 & dgs$Stage4vsControl != 0, ]
UniqueS4=rownames(UniqueS4)
UniqueS4=as.data.frame(UniqueS4)
colnames(UniqueS4)[1] = "ensembl_gene_id"

UniqueS4=match_df(S4, UniqueS4, on="ensembl_gene_id")

#Save stage-specific DEGs

stage_specific_DEGs=rbind(UniqueS1, UniqueS2, UniqueS3, UniqueS4)

#####################    Visualization   ###################
#####################   MA PLOT   #######################
####################   Stage 1 Visualization ###################


S1_coef="Stage1vsControl"
S1_coef <- match(S1_coef, colnames(efit$coefficients))

S2_coef="Stage2vsControl"
S2_coef <- match(S2_coef, colnames(efit$coefficients))

S3_coef="Stage3vsControl"
S3_coef <- match(S3_coef, colnames(efit$coefficients))

S4_coef="Stage4vsControl"
S4_coef <- match(S4_coef, colnames(efit$coefficients))

# ######################### VOLCANO PLOT #####################
# ######################### Stage 1 VOLCANO PLOT #########################



volcano_df1 <- data.frame(
  Gene = rownames(efit),
  logFC = efit$coefficients[, 1],   
  pval = efit$p.value[, 1]          
) %>%
  mutate(
    negLogP = -log10(pval),
    Significance = case_when(
      logFC > 1 & pval < 0.05 ~ "Upregulated",
      logFC < -1 & pval < 0.05 ~ "Downregulated",
      TRUE ~ "Not significant"
    )
  )

volcano_plot1 <- ggplot(volcano_df1, aes(x = logFC, y = negLogP)) +
  geom_point(aes(color = Significance), alpha = 0.6, size = 0.5) +
  scale_color_manual(values = c("Upregulated" = "red", 
                                "Downregulated" = "blue", 
                                "Not significant" = "black")) +
  labs(
    title = "a)    ",
    x = "Log2 Fold Change",
    y = "-log10(P-value)",
    color = " "
  ) +
  theme_minimal(base_size = 7) +
  theme(
    plot.title = element_text(hjust = 0.001, size = 7),
    legend.position = c(0.9, 1),   
    legend.background = element_rect(fill = "white", 
                                     color = "white", size = 0.1),
    legend.key.size = unit(0.5, "lines"),
    legend.text = element_text(size = 4)
  )


# ######################### Stage 2 VOLCANO PLOT #########################

volcano_df2 <- data.frame(
  Gene = rownames(efit),
  logFC = efit$coefficients[, 2],   # Stage 2 comparison
  pval = efit$p.value[, 2]
) %>%
  mutate(
    negLogP = -log10(pval),
    Significance = case_when(
      logFC > 1 & pval < 0.05 ~ "Upregulated",
      logFC < -1 & pval < 0.05 ~ "Downregulated",
      TRUE ~ "Not significant"
    )
  )

volcano_plot2 <- ggplot(volcano_df2, aes(x = logFC, y = negLogP)) +
  geom_point(aes(color = Significance), alpha = 0.6, size = 0.5) +
  scale_color_manual(values = c("Upregulated" = "red", 
                                "Downregulated" = "blue", 
                                "Not significant" = "black")) +
  labs(
    title = "b) ",
    x = "Log2 Fold Change",
    y = "-log10(P-value)",
    color = " "
  ) +
  theme_minimal(base_size = 7) +
  theme(
    plot.title = element_text(hjust = 0.001, size = 7),
    legend.position = c(0.9, 1),
    legend.background = element_rect(fill = "white", 
                                     color = "white", size = 0.1),
    legend.key.size = unit(0.5, "lines"),
    legend.text = element_text(size = 4)
  )

ggsave("Volcano2.png", plot = volcano_plot2, 
       width = 4, 
       height = 2, 
       dpi = 300, 
       bg = "white")

print(volcano_plot2)

# ######################### Stage 3 VOLCANO PLOT #########################

volcano_df3 <- data.frame(
  Gene = rownames(efit),
  logFC = efit$coefficients[, 3],   # Stage 3 comparison
  pval = efit$p.value[, 3]
) %>%
  mutate(
    negLogP = -log10(pval),
    Significance = case_when(
      logFC > 1 & pval < 0.05 ~ "Upregulated",
      logFC < -1 & pval < 0.05 ~ "Downregulated",
      TRUE ~ "Not significant"
    )
  )

volcano_plot3 <- ggplot(volcano_df3, aes(x = logFC, y = negLogP)) +
  geom_point(aes(color = Significance), alpha = 0.6, size = 0.5) +
  scale_color_manual(values = c("Upregulated" = "red", 
                                "Downregulated" = "blue", 
                                "Not significant" = "black")) +
  labs(
    title = "c) ",
    x = "Log2 Fold Change",
    y = "-log10(P-value)",
    color = " "
  ) +
  theme_minimal(base_size = 7) +
  theme(
    plot.title = element_text(hjust = 0.001, size = 7),
    legend.position = c(0.9, 1),
    legend.background = element_rect(fill = "white", 
                                     color = "white", size = 0.1),
    legend.key.size = unit(0.5, "lines"),
    legend.text = element_text(size = 4)
  )

ggsave("Volcano3.png", plot = volcano_plot3, 
       width = 4, 
       height = 2, 
       dpi = 300, 
       bg = "white")

print(volcano_plot3)

# ######################### Stage 4 VOLCANO PLOT #########################

volcano_df4 <- data.frame(
  Gene = rownames(efit),
  logFC = efit$coefficients[, 4],   # Stage 4 comparison
  pval = efit$p.value[, 4]
) %>%
  mutate(
    negLogP = -log10(pval),
    Significance = case_when(
      logFC > 1 & pval < 0.05 ~ "Upregulated",
      logFC < -1 & pval < 0.05 ~ "Downregulated",
      TRUE ~ "Not significant"
    )
  )

volcano_plot4 <- ggplot(volcano_df4, aes(x = logFC, y = negLogP)) +
  geom_point(aes(color = Significance), alpha = 0.6, size = 0.5) +
  scale_color_manual(values = c("Upregulated" = "red", 
                                "Downregulated" = "blue", 
                                "Not significant" = "black")) +
  labs(
    title = "d) ",
    x = "Log2 Fold Change",
    y = "-log10(P-value)",
    color = " "
  ) +
  theme_minimal(base_size = 7) +
  theme(
    plot.title = element_text(hjust = 0.001, size = 7),
    legend.position = c(0.9, 1),
    legend.background = element_rect(fill = "white", 
                                     color = "white", size = 0.1),
    legend.key.size = unit(0.5, "lines"),
    legend.text = element_text(size = 4)
  )

combined_volcano <- (volcano_plot1 + volcano_plot2) /
  (volcano_plot3 + volcano_plot4) +
  plot_layout(guides = "collect") &
  theme(legend.position = "right")  


ggsave("combined_volcano.png", combined_volcano,
       width = 1800, height = 1200, dpi = 300, units = "px", bg = "white")



# ROC curve for unique DEGs ----

lcpm=cpm.DGEList(dge, log=TRUE)
dge$counts=lcpm

design <- model.matrix(~0+dge$samples$Stage)

colnames(design)[1] ="Control"
colnames(design)[2] ="Stage1"
colnames(design)[3] ="Stage2"
colnames(design)[4] ="Stage3"
colnames(design)[5] ="Stage4"

design <- design[, c(2,3,4,5,1)]

lcpm=removeBatchEffect(dge, design=design, covariates = as.matrix(W))

cpm=lcpm


# ROC curve for Stage 1 unique DEGs ----

S1_samples <- Samples_Cancer %>%
  filter(Stage == "Stage1")

S1_samples <- rbind(S1_samples, Samples_Control)

S1_samples <- S1_samples$Sample.ID
label1 <- c(rep(1, 36), rep(0, 253 - 36))

Stage1 <- UniqueS1 

Stage1_genes <- UniqueS1$ensembl_gene_id
Stage1_counts <- cpm[rownames(cpm) %in% Stage1_genes, ]
Stage1 <- Stage1_counts[, S1_samples]
Stage1 <- rbind(Stage1, label1)

rownames(Stage1)[nrow(Stage1)] <- "label"

Stage1 <- as.matrix(Stage1)

genes1 <- rownames(Stage1)[-nrow(Stage1)] 

UniqueS1 <- UniqueS1[match(genes1, UniqueS1$ensembl_gene_id), ]

genes1=as.character(UniqueS1$hgnc_symbol)

num_genes1 <- length(genes1)

roc_objects_stage1 <- vector("list", length = num_genes1)
auc_values_stage1 <- numeric(num_genes1)
genes_above_threshold_stage1 <- character(num_genes1)

for (i in 1:num_genes1) {
  gene_expression <- as.numeric(Stage1[i, -ncol(Stage1)])
  labels <- as.numeric(Stage1[nrow(Stage1), -ncol(Stage1)])
  model <- glm(labels ~ gene_expression, family = binomial(link = "logit"))
  predictions <- predict(model, type = "response")
  roc_obj <- roc(labels, predictions)
  roc_objects_stage1[[i]] <- roc_obj
  auc_values_stage1[i] <- auc(roc_obj)
  genes_above_threshold_stage1[i] <- genes1[i]
}


auc_stage1 <- data.frame(Gene = genes_above_threshold_stage1, AUC = auc_values_stage1)

for (i in seq_along(roc_objects_stage1)) {
  if (!is.null(roc_objects_stage1[[i]])) {
    gene_name <- genes1[i]
    plot(roc_objects_stage1[[i]], main = paste("ROC Curve -", gene_name))
    text(0.5, 0.4, paste("AUC =", round(auc(roc_objects_stage1[[i]]), 2)), adj = c(0.5, 0.5), cex = 1)
  }
}

auc_stage1 <- auc_stage1 %>% arrange(desc(AUC))

auc_stage1 <- auc_stage1 %>%
  left_join(UniqueS1 %>% dplyr::select(hgnc_symbol, Status), 
            by = c("Gene" = "hgnc_symbol"))

auc_stage1 <- auc_stage1 %>%
  left_join(UniqueS1 %>% dplyr::select(hgnc_symbol, ensembl_gene_id), 
            by = c("Gene" = "hgnc_symbol"))

auc_stage1 <- auc_stage1 %>%
  left_join(UniqueS1 %>% dplyr::select(hgnc_symbol, logFC), 
            by = c("Gene" = "hgnc_symbol"))

auc_stage1 <- auc_stage1 %>%
  left_join(UniqueS1 %>% dplyr::select(hgnc_symbol, adj.P.Val), 
            by = c("Gene" = "hgnc_symbol"))


# ROC curve for Stage 2 unique DEGs ----

S2_samples <- Samples_Cancer %>%
  filter(Stage == "Stage2")

S2_samples <- rbind(S2_samples, Samples_Control)

S2_samples <- S2_samples$Sample.ID
label2 <- c(rep(1, 18), rep(0, 235 - 18))

Stage2 <- UniqueS2 

Stage2_genes <- Stage2$ensembl_gene_id
Stage2_counts <- cpm[rownames(cpm) %in% Stage2_genes, ]
Stage2 <- Stage2_counts[, S2_samples]
Stage2 <- rbind(Stage2, label2)

rownames(Stage2)[nrow(Stage2)] <- "label"

Stage2 <- as.matrix(Stage2)

genes2 <- rownames(Stage2)[-nrow(Stage2)] 

UniqueS2 <- UniqueS2[match(genes2, UniqueS2$ensembl_gene_id), ]

genes2 = as.character(UniqueS2$hgnc_symbol)

num_genes2 <- length(genes2)

roc_objects_stage2 <- vector("list", length = num_genes2)
auc_values_stage2 <- numeric(num_genes2)
genes_above_threshold_stage2 <- character(num_genes2)

for (i in 1:num_genes2) {
  gene_expression <- as.numeric(Stage2[i, -ncol(Stage2)])
  labels <- as.numeric(Stage2[nrow(Stage2), -ncol(Stage2)])
  model <- glm(labels ~ gene_expression, family = binomial(link = "logit"))
  predictions <- predict(model, type = "response")
  roc_obj <- roc(labels, predictions)
  roc_objects_stage2[[i]] <- roc_obj
  auc_values_stage2[i] <- auc(roc_obj)
  genes_above_threshold_stage2[i] <- genes2[i]
}

auc_stage2 <- data.frame(Gene = genes_above_threshold_stage2, AUC = auc_values_stage2)

for (i in seq_along(roc_objects_stage2)) {
  if (!is.null(roc_objects_stage2[[i]])) {
    gene_name <- genes2[i]
    plot(roc_objects_stage2[[i]], main = paste("ROC Curve -", gene_name))
    text(0.5, 0.4, paste("AUC =", round(auc(roc_objects_stage2[[i]]), 2)), adj = c(0.5, 0.5), cex = 1)
  }
}

auc_stage2 <- auc_stage2 %>% arrange(desc(AUC))

auc_stage2 <- auc_stage2 %>%
  left_join(UniqueS2 %>% dplyr::select(hgnc_symbol, Status), 
            by = c("Gene" = "hgnc_symbol"))

auc_stage2 <- auc_stage2 %>%
  left_join(UniqueS2 %>% dplyr::select(hgnc_symbol, ensembl_gene_id), 
            by = c("Gene" = "hgnc_symbol"))

auc_stage2 <- auc_stage2 %>%
  left_join(UniqueS2 %>% dplyr::select(hgnc_symbol, logFC), 
            by = c("Gene" = "hgnc_symbol"))

auc_stage2 <- auc_stage2 %>%
  left_join(UniqueS2 %>% dplyr::select(hgnc_symbol, adj.P.Val), 
            by = c("Gene" = "hgnc_symbol"))


# ROC curve for Stage 3 unique DEGs ----

S3_samples=Samples_Cancer%>%
  filter(Stage=="Stage3")

S3_samples=rbind(S3_samples, Samples_Control)

S3_samples=S3_samples$Sample.ID

label3=c(rep(1, 49), rep(0, 266 - 49))

Stage3 <- UniqueS3 

Stage3_genes <- Stage3$ensembl_gene_id
Stage3_counts <- cpm[rownames(cpm) %in% Stage3_genes, ]
Stage3 <- Stage3_counts[, S3_samples]
Stage3 <- rbind(Stage3, label3)

rownames(Stage3)[nrow(Stage3)] <- "label"

Stage3 <- as.matrix(Stage3)

genes3 <- rownames(Stage3)[-nrow(Stage3)] 

UniqueS3 <- UniqueS3[match(genes3, UniqueS3$ensembl_gene_id), ]

genes3 = as.character(UniqueS3$hgnc_symbol)

num_genes3 <- length(genes3)

roc_objects_stage3 <- vector("list", length = num_genes3)
auc_values_stage3 <- numeric(num_genes3)
genes_above_threshold_stage3 <- character(num_genes3)

for (i in 1:num_genes3) {
  gene_expression <- as.numeric(Stage3[i, -ncol(Stage3)])
  labels <- as.numeric(Stage3[nrow(Stage3), -ncol(Stage3)])
  model <- glm(labels ~ gene_expression, family = binomial(link = "logit"))
  predictions <- predict(model, type = "response")
  roc_obj <- roc(labels, predictions)
  roc_objects_stage3[[i]] <- roc_obj
  auc_values_stage3[i] <- auc(roc_obj)
  genes_above_threshold_stage3[i] <- genes3[i]
}

auc_stage3 <- data.frame(Gene = genes_above_threshold_stage3, AUC = auc_values_stage3)

for (i in seq_along(roc_objects_stage3)) {
  if (!is.null(roc_objects_stage3[[i]])) {
    gene_name <- genes3[i]
    plot(roc_objects_stage3[[i]], main = paste("ROC Curve -", gene_name))
    text(0.5, 0.4, paste("AUC =", round(auc(roc_objects_stage3[[i]]), 2)), adj = c(0.5, 0.5), cex = 1)
  }
}

auc_stage3 <- auc_stage3 %>% arrange(desc(AUC))

auc_stage3 <- auc_stage3 %>%
  left_join(UniqueS3 %>% dplyr::select(hgnc_symbol, Status), 
            by = c("Gene" = "hgnc_symbol"))

auc_stage3 <- auc_stage3 %>%
  left_join(UniqueS3 %>% dplyr::select(hgnc_symbol, ensembl_gene_id), 
            by = c("Gene" = "hgnc_symbol"))

auc_stage3 <- auc_stage3 %>%
  left_join(UniqueS3 %>% dplyr::select(hgnc_symbol, logFC), 
            by = c("Gene" = "hgnc_symbol"))

auc_stage3 <- auc_stage3 %>%
  left_join(UniqueS3 %>% dplyr::select(hgnc_symbol, adj.P.Val), 
            by = c("Gene" = "hgnc_symbol"))


# ROC curve for Stage 4 unique DEGs ----

S4_samples=Samples_Cancer%>%
  filter(Stage=="Stage4")

S4_samples=rbind(S4_samples, Samples_Control)

S4_samples=S4_samples$Sample.ID

label4=c(rep(1, 35), rep(0, 252 - 35))

Stage4 <- UniqueS4 

Stage4_genes <- UniqueS4$ensembl_gene_id
Stage4_counts <- cpm[rownames(cpm) %in% Stage4_genes, ]
Stage4 <- Stage4_counts[, S4_samples]
Stage4 <- rbind(Stage4, label4)

rownames(Stage4)[nrow(Stage4)] <- "label"

Stage4 <- as.matrix(Stage4)

genes4 <- rownames(Stage4)[-nrow(Stage4)] 



num_genes4 <- length(genes4)

roc_objects_stage4 <- vector("list", length = num_genes4)
auc_values_stage4 <- numeric(num_genes4)
genes_above_threshold_stage4 <- character(num_genes4)

for (i in 1:num_genes4) {
  gene_expression <- as.numeric(Stage4[i, -ncol(Stage4)])
  labels <- as.numeric(Stage4[nrow(Stage4), -ncol(Stage4)])
  model <- glm(labels ~ gene_expression, family = binomial(link = "logit"))
  predictions <- predict(model, type = "response")
  roc_obj <- roc(labels, predictions)
  roc_objects_stage4[[i]] <- roc_obj
  auc_values_stage4[i] <- auc(roc_obj)
  genes_above_threshold_stage4[i] <- genes4[i]
}

auc_stage4 <- data.frame(Gene = genes_above_threshold_stage4, AUC = auc_values_stage4)

for (i in seq_along(roc_objects_stage4)) {
  if (!is.null(roc_objects_stage4[[i]])) {
    gene_name <- genes4[i]
    plot(roc_objects_stage4[[i]], main = paste("ROC Curve -", gene_name))
    text(0.5, 0.4, paste("AUC =", round(auc(roc_objects_stage4[[i]]), 2)), adj = c(0.5, 0.5), cex = 1)
  }
}

auc_stage4 <- auc_stage4 %>% arrange(desc(AUC))

auc_stage4 <- auc_stage4 %>%
  left_join(UniqueS4 %>% dplyr::select(ensembl_gene_id, Status), 
            by = c("Gene" = "ensembl_gene_id"))

auc_stage4 <- auc_stage4 %>%
  left_join(UniqueS4 %>% dplyr::select(ensembl_gene_id, hgnc_symbol), 
            by = c("Gene" = "ensembl_gene_id"))

auc_stage4 <- auc_stage4 %>%
  left_join(UniqueS4 %>% dplyr::select(ensembl_gene_id, logFC), 
            by = c("Gene" = "ensembl_gene_id"))

auc_stage4 <- auc_stage4 %>%
  left_join(UniqueS4 %>% dplyr::select(ensembl_gene_id, adj.P.Val), 
            by = c("Gene" = "ensembl_gene_id"))

colnames(auc_stage4)[1]="ensembl_gene_id"
colnames(auc_stage4)[4]="Gene"

auc_stage4=auc_stage4[c(4,2,3,1,5,6)]


auc_stage1['Stage']="Stage1"
auc_stage2['Stage']="Stage2"
auc_stage3['Stage']="Stage3"
auc_stage4['Stage']="Stage4"

DEGs=rbind(S1,S2,S3,S4)
AUC_original_dataset=rbind(auc_stage1, auc_stage2, auc_stage3, auc_stage4)


results=createWorkbook()
addWorksheet(results, "Differential_exp_result")
addWorksheet(results, "stage_specific_DEGs")
addWorksheet(results, "AUC_original_dataset")

writeData(results, "Differential_exp_result", DEGs)
writeData(results, "stage_specific_DEGs", stage_specific_DEGs)
writeData(results, "AUC_original_dataset", AUC_original_dataset)

saveWorkbook(results, "Differential_exp_results.xlsx", overwrite = TRUE)

sessionInfo()

#End of code
#Next step is validation with dataset GSE158508




