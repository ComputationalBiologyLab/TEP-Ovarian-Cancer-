#Get working directory 

getwd()

#Set working directory to the location of Data through setwd()

# Loading libraries ----

library(BSgenome)
library(BSgenome.Hsapiens.UCSC.hg19)
library(caret)
library(DBI)            
library(dbplyr)
library(dplyr)
library(DT)
library(dtplyr)
library(edgeR)
library(factoextra)
library(FactoMineR)
library(genefilter)
library(ggplot2)         
library(glmnet)
library(gplots)
library(gProfileR)
library(limma)
library(magrittr)
library(openxlsx)
library(patchwork)
library(readxl)
library(RColorBrewer)
library(segmented)
library(SummarizedExperiment)
library(survival)
library(survminer)
library(tidygraph)
library(tidyr)
library(tidyselect)
library(tidytree)
library(tidyverse)
library(TCGAbiolinks)
library(TCGAWorkflowData)
library(Trendy)

# Loading clinical data ----

clinical=read.delim("ov_tcga_clinical_data.tsv", sep="\t")

tcga_clinical=GDCquery_clinic("TCGA-OV", type = "clinical")

table(clinical$Neoplasm.American.Joint.Committee.on.Cancer.Clinical.Group.Stage)
sum(is.na(clinical$Overall.Survival..Months.))
sum(is.na(clinical$Neoplasm.American.Joint.Committee.on.Cancer.Clinical.Group.Stage.))
sum(is.na(clinical$Overall.Survival.Status))

clinical <- clinical %>%
  filter(
    !is.na(Overall.Survival..Months.),
    !is.na(Overall.Survival.Status)
  )

clinical <- clinical %>%
  filter(grepl("01$", Sample.ID))

clinical <- clinical %>%
  filter(!is.na(Neoplasm.American.Joint.Committee.on.Cancer.Clinical.Group.Stage))


samples=clinical$Patient.ID


colnames(tcga_clinical)
colnames(clinical)


clinical <- clinical %>%
  left_join(tcga_clinical %>% dplyr::select(submitter_id, days_to_last_follow_up),
            by = c("Patient.ID" = "submitter_id"))

clinical <- clinical %>%
  left_join(tcga_clinical %>% dplyr::select(submitter_id, figo_stage),
            by = c("Patient.ID" = "submitter_id"))

data = readRDS(file = "tcga_data.RDS")

batch=get_IDs(data)
gene_info=rowRanges(data)
gene_info=gene_info@elementMetadata@listData

clinical=clinical%>%
     filter(Patient.ID %in% batch$patient)

colnames(clinical)

###### RNA Seq preprocessing ######

exp_ov_preprocessed <- TCGAanalyze_Preprocessing(
  object = data,
  cor.cut = 0.6,    
  datatype = "unstranded",
  filename = "OV_IlluminaHiSeq_RNASeqV2.png"
)

exp_ov_normalized <- TCGAanalyze_Normalization(
  tabDF = exp_ov_preprocessed,
  geneInfo = TCGAbiolinks::geneInfoHT,
  method = "gcContent"
)

exp_ov_filtered <- TCGAanalyze_Filtering(
  tabDF = exp_ov_normalized,
  method = "quantile",
  qnt.cut =  0.25
) 

rm("exp_ov_normalized")
rm("exp_ov_preprocessed")

counts=as.data.frame(exp_ov_filtered)

clinical$deceased = clinical$Overall.Survival.Status == "1:DECEASED"

#https://doi.org/10.3389/fonc.2020.00162
#https://pmc.ncbi.nlm.nih.gov/articles/PMC11782890/
#https://doi.org/10.1016/j.mehy.2013.01.027

clinical$figo_stage = gsub("Stage IC", "1", clinical$figo_stage)
clinical$figo_stage <- gsub("Stage IIA|Stage IIB|Stage IIC", "2", clinical$figo_stage)
clinical$figo_stage <- gsub("Stage IIIA|Stage IIIB|Stage IIIC", "3", clinical$figo_stage)
clinical$figo_stage <- gsub("Stage IV", "4", clinical$figo_stage)

clinical$figo_stage = gsub("1|2", "early", clinical$figo_stage)
clinical$figo_stage = gsub("3|4", "late", clinical$figo_stage)

table(clinical$figo_stage)

clinical$figo_stage <- as.factor(clinical$figo_stage)

clinical <- clinical %>%
      left_join(batch %>% dplyr::select(patient, barcode),
                by = c("Patient.ID" = "patient"))


table(clinical$Patient.ID==batch$patient)
 
counts <- counts[, clinical$barcode[clinical$barcode %in% colnames(counts)]]

table(colnames(counts)==clinical$barcode)

colnames(clinical)[47]="overall_survival"

counts <- counts[rowSums(counts[, -1]) >= 50, ]

#### Process all samples then subset early and late ####

dge = DGEList(
  counts=counts,
  samples=clinical)

dge = calcNormFactors(dge,method="TMM")

dim(dge)

####Get candidate biomarker genes (CBGs) ----

CBGs=read_excel("C:/Users/shaim/Desktop/Data/Supplementary table.xlsx", sheet = 4, skip = 1)

norm_counts=cpm.DGEList(dge)

norm_counts=as.data.frame(norm_counts)

clinical <- clinical %>%
  relocate(days_to_last_follow_up, .after = overall_survival)

colnames(clinical)[46:50]

##### EXTRACT EARLY COUNTS ####

early_clinical=clinical%>%
  filter(figo_stage=="early")

early_exp <- norm_counts[, early_clinical$barcode[early_clinical$barcode %in% colnames(norm_counts)]]

##### Survival function ####

perform_survival_analysis <- function(gene_list, norm_counts, survival_data, sig_level = 0.05) {

  results <- list()
  
  for (gene in gene_list) {
    if (!(gene %in% rownames(norm_counts))) {
      warning(paste("Gene", gene, "not found in norm_counts. Skipping."))
      next
    }
    
    gene_expression <- norm_counts[gene, , drop = FALSE]
    gene_expression <- t(gene_expression)
    gene_expression <- as.data.frame(gene_expression)
    colnames(gene_expression)[1] <- gene
    
    if (any(is.na(gene_expression))) {
      warning(paste("Gene expression data for", gene, "contains missing values. Skipping."))
      next
    }
    
    gene_data <- cbind(gene_expression, survival_data)
    
    if (nrow(gene_data) == 0 ||
        any(is.na(gene_data$overall_survival)) ||
        any(is.na(gene_data$deceased)) ||
        any(is.na(gene_data[[gene]]))) {
      warning(paste("No valid data or missing values detected for gene", gene, ". Skipping."))
      next
    }
    
    median_value <- median(gene_data[[gene]], na.rm = TRUE)
    gene_data$gene_class <- ifelse(gene_data[[gene]] >= median_value, "High", "Low")
    
    fit <- tryCatch({
      survfit(Surv(overall_survival, deceased) ~ gene_class, data = gene_data)
    }, error = function(e) {
      warning(paste("Error in survfit for gene", gene, ":", e$message))
      return(NULL)
    })
    
    if (is.null(fit)) next
    
    pval <- tryCatch({
      surv_pvalue(fit, data = gene_data)$pval
    }, error = function(e) {
      warning(paste("Error in surv_pvalue for gene", gene, ":", e$message))
      return(NA)
    })
    
    
    km_plot <- ggsurvplot(
      fit,
      data = gene_data,
      pval = paste0("p = ", signif(pval, 3)),
      risk.table = FALSE,   
      ggtheme = theme_minimal(),
      palette = c("red", "blue"),
      title = paste("Survival Curve for", gene)
    )
    
    print(km_plot)
    
    results[[gene]] <- list(
      fit = fit,
      pval = pval,
      km_plot = km_plot
    )
  }
  
  return(results)
}


CBG1=CBGs%>%filter(Stage=="Stage1")
CBG2=CBGs%>%filter(Stage=="Stage2")
CBG3=CBGs%>%filter(Stage=="Stage3")
CBG4=CBGs%>%filter(Stage=="Stage4")


result1 <- perform_survival_analysis(CBG1$ensembl_gene_id, early_exp, early_clinical)

result2 <- perform_survival_analysis(CBG2$ensembl_gene_id, early_exp, early_clinical)

####Extract late counts#####

late_clinical=clinical%>%
  filter(figo_stage=="late")

late_exp <- norm_counts[, late_clinical$barcode[late_clinical$barcode %in% colnames(norm_counts)]]

result3 <- perform_survival_analysis(CBG3$ensembl_gene_id, late_exp, late_clinical)

result4 <- perform_survival_analysis(CBG4$ensembl_gene_id, late_exp, late_clinical)


################Customized plot for ENSG00000155629 (PIK3AP1) #############################

PIK3AP1 <- early_exp[c("ENSG00000155629"),]

PIK3AP1 = t(PIK3AP1)

PIK3AP1 = as.data.frame(PIK3AP1)
colnames(PIK3AP1)[1]="PIK3AP1"
PIK3AP1$PIK3AP1=as.numeric(PIK3AP1$PIK3AP1)
PIK3AP1_median_value = median(PIK3AP1$PIK3AP1)
print(PIK3AP1_median_value)

PIK3AP1 = cbind(PIK3AP1, early_clinical)

PIK3AP1$group = ifelse(PIK3AP1$PIK3AP1 >= PIK3AP1_median_value, "High Expression", "Low Expression")

PIK3AP1_fit = survfit(Surv(overall_survival, deceased) ~ group, data=PIK3AP1)

pval = surv_pvalue(PIK3AP1_fit, data=PIK3AP1)$pval
print(pval)

plot1 <- ggsurvplot(
  PIK3AP1_fit, 
  data = PIK3AP1, 
  pval = TRUE, 
  risk.table = FALSE, 
  title = "                                                         
a)                                 PIK3AP1",
  palette = c("red", "blue"),
  xlab = "Time (months)",
  xscale = 1,
  size = 0.2,
  censor.size = 1,
  ggtheme = theme_minimal(base_size = 3) + theme(
    aspect.ratio = 0.3,
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.3),  
    axis.ticks = element_line(size = 0.3),
    axis.text = element_text(size = 4),
    axis.title = element_text(size = 3),
    plot.title = element_text(size = 4, hjust = 0.01),
    plot.margin = margin(1, 1, 1, 1),
    legend.title = element_text(size = 3),
    legend.text = element_text(size = 3),
    legend.key = element_blank()
  ),
  pval.size = 1
)

ggsave("PIK3AP1_survival_plot_small.png", plot = plot1$plot,
       width = 1200, height = 600, units = "px", dpi = 300, bg = "white")

################Customized plot for ENSG00000137752 (CASP1) ######################

CASP1 <- early_exp[c("ENSG00000137752"),]

CASP1 = t(CASP1)

CASP1 = as.data.frame(CASP1)
colnames(CASP1)[1]="CASP1"
CASP1$CASP1=as.numeric(CASP1$CASP1)
CASP1_median_value = median(CASP1$CASP1)
print(CASP1_median_value)

CASP1 = cbind(CASP1, early_clinical)

CASP1$group = ifelse(CASP1$CASP1 >= CASP1_median_value, "High Expression", "Low Expression")

CASP1_fit = survfit(Surv(overall_survival, deceased) ~ group, data=CASP1)

pval = surv_pvalue(CASP1_fit, data=CASP1)$pval
print(pval)

plot2 <- ggsurvplot(
  CASP1_fit, 
  data = CASP1, 
  pval = TRUE, 
  risk.table = FALSE, 
  title = "                                                          
b)                                 CASP1",
  palette = c("red", "blue"),
  xlab = "Time (months)",
  xscale = 1,
  size = 0.2,
  censor.size = 1,
  ggtheme = theme_minimal(base_size = 3) + theme(
    aspect.ratio = 0.3,
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black", size = 0.3),  
    axis.ticks = element_line(size = 0.3),
    axis.text = element_text(size = 4),
    axis.title = element_text(size = 3),
    plot.title = element_text(size = 4, hjust = 0.01),
    plot.margin = margin(1, 1, 1, 1),
    legend.title = element_text(size = 3),
    legend.text = element_text(size = 3),
    legend.key = element_blank()
  ),
  pval.size = 1
)


ggsave("CASP1_survival_plot_small.png", plot = plot2$plot,
       width = 1200, height = 600, units = "px", dpi = 300, bg = "white")

################Customized plot for ENSG00000116260 (QSOX1) ######################

QSOX1 <- early_exp[c("ENSG00000116260"),]

QSOX1 = t(QSOX1)

QSOX1 = as.data.frame(QSOX1)
colnames(QSOX1)[1]="QSOX1"
QSOX1$QSOX1=as.numeric(QSOX1$QSOX1)
QSOX1_median_value = median(QSOX1$QSOX1)
print(QSOX1_median_value)

QSOX1 = cbind(QSOX1, early_clinical)

QSOX1$group = ifelse(QSOX1$QSOX1 >= QSOX1_median_value, "High Expression", "Low Expression")

QSOX1_fit = survfit(Surv(overall_survival, deceased) ~ group, data=QSOX1)

pval = surv_pvalue(QSOX1_fit, data=QSOX1)$pval
print(pval)

plot3 <- ggsurvplot(
  QSOX1_fit, 
  data = QSOX1, 
  pval = TRUE, 
  risk.table = FALSE, 
  title = "                                                          
c)                                 QSOX1",
  palette = c("red", "blue"),
  xlab = "Time (months)",
  xscale = 1,
  size = 0.2,
  censor.size = 1,
  ggtheme = theme_minimal(base_size = 3) + theme(
    aspect.ratio = 0.3,
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black", size = 0.3),  
    axis.ticks = element_line(size = 0.3),
    axis.text = element_text(size = 4),
    axis.title = element_text(size = 3),
    plot.title = element_text(size = 4, hjust = 0.01),
    plot.margin = margin(1, 1, 1, 1),
    legend.title = element_text(size = 3),
    legend.text = element_text(size = 3),
    legend.key = element_blank()
  ),
  pval.size = 1
)


ggsave("QSOX1_survival_plot_small.png", plot = plot3$plot,
       width = 1200, height = 600, units = "px", dpi = 300, bg = "white")

################Customized plot for ENSG00000172543 (CTSW) ######################

CTSW <- early_exp[c("ENSG00000172543"),]

CTSW = t(CTSW)

CTSW = as.data.frame(CTSW)
colnames(CTSW)[1]="CTSW"
CTSW$CTSW=as.numeric(CTSW$CTSW)
CTSW_median_value = median(CTSW$CTSW)
print(CTSW_median_value)

CTSW = cbind(CTSW, early_clinical)

CTSW$group = ifelse(CTSW$CTSW >= CTSW_median_value, "High Expression", "Low Expression")

CTSW_fit = survfit(Surv(overall_survival, deceased) ~ group, data=CTSW)

pval = surv_pvalue(CTSW_fit, data=CTSW)$pval
print(pval)

plot4 <- ggsurvplot(
  CTSW_fit, 
  data = CTSW, 
  pval = TRUE, 
  risk.table = FALSE, 
  title = "                                                          
d)                                 CTSW",
  palette = c("red", "blue"),
  xlab = "Time (months)",
  xscale = 1,
  size = 0.2,
  censor.size = 1,
  ggtheme = theme_minimal(base_size = 3) + theme(
    aspect.ratio = 0.3,
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black", size = 0.3),  
    axis.ticks = element_line(size = 0.3),
    axis.text = element_text(size = 4),
    axis.title = element_text(size = 3),
    plot.title = element_text(size = 4, hjust = 0.01),
    plot.margin = margin(1, 1, 1, 1),
    legend.title = element_text(size = 3),
    legend.text = element_text(size = 3),
    legend.key = element_blank()
  ),
  pval.size = 1
)


ggsave("CTSW_survival_plot_small.png", plot = plot4$plot,
       width = 1200, height = 600, units = "px", dpi = 300, bg = "white")

################Customized plot for ENSG00000100453 (GZMB) ######################

GZMB <- late_exp[c("ENSG00000100453"),]

GZMB = t(GZMB)

GZMB = as.data.frame(GZMB)
colnames(GZMB)[1]="GZMB"
GZMB$GZMB=as.numeric(GZMB$GZMB)
GZMB_median_value = median(GZMB$GZMB)
print(GZMB_median_value)

GZMB = cbind(GZMB, late_clinical)

GZMB$group = ifelse(GZMB$GZMB >= GZMB_median_value, "High Expression", "Low Expression")

GZMB_fit = survfit(Surv(overall_survival, deceased) ~ group, data=GZMB)

pval = surv_pvalue(GZMB_fit, data=GZMB)$pval
print(pval)

plot5 <- ggsurvplot(
  GZMB_fit, 
  data = GZMB, 
  pval = TRUE, 
  risk.table = FALSE, 
  title = "                                                          
e)                                 GZMB",
  palette = c("red", "blue"),
  xlab = "Time (months)",
  xscale = 1,
  size = 0.2,
  censor.size = 1,
  ggtheme = theme_minimal(base_size = 3) + theme(
    aspect.ratio = 0.3,
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black", size = 0.3),  
    axis.ticks = element_line(size = 0.3),
    axis.text = element_text(size = 4),
    axis.title = element_text(size = 3),
    plot.title = element_text(size = 4, hjust = 0.01),
    plot.margin = margin(1, 1, 1, 1),
    legend.title = element_text(size = 3),
    legend.text = element_text(size = 3),
    legend.key = element_blank()
  ),
  pval.size = 1
)


ggsave("GZMB_survival_plot_small.png", plot = plot5$plot,
       width = 1200, height = 600, units = "px", dpi = 300, bg = "white")

################Customized plot for ENSG00000133742 (CA1) ######################

CA1 <- late_exp[c("ENSG00000133742"), ]  

CA1 <- t(CA1)
CA1 <- as.data.frame(CA1)
colnames(CA1)[1] <- "CA1"
CA1$CA1 <- as.numeric(CA1$CA1)

CA1_median_value <- median(CA1$CA1, na.rm = TRUE)
print(CA1_median_value)

CA1 <- cbind(CA1, late_clinical)

CA1$group = ifelse(CA1$CA1 >= CA1_median_value, "High Expression", "Low Expression")

CA1_fit <- survfit(Surv(overall_survival, deceased) ~ group, data = CA1)

pval <- surv_pvalue(CA1_fit, data = CA1)$pval
print(pval)

plot6 <- ggsurvplot(
  CA1_fit, 
  data = CA1, 
  pval = TRUE, 
  risk.table = FALSE, 
  title = "                                                          
f)                                 CA1",
  palette = c("red", "blue"),
  xlab = "Time (months)",
  xscale = 1,
  size = 0.2,
  censor.size = 1,
  ggtheme = theme_minimal(base_size = 3) + theme(
    aspect.ratio = 0.3,
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black", size = 0.3),  
    axis.ticks = element_line(size = 0.3),
    axis.text = element_text(size = 4),
    axis.title = element_text(size = 3),
    plot.title = element_text(size = 4, hjust = 0.01),
    plot.margin = margin(1, 1, 1, 1),
    legend.title = element_text(size = 3),
    legend.text = element_text(size = 3),
    legend.key = element_blank()
  ),
  pval.size = 1
)

ggsave("CA1_survival_plot_small.png", plot = plot6$plot,
       width = 1200, height = 600, units = "px", dpi = 300, bg = "white")


plots <- list(
  plot1$plot + theme(legend.position = "none"),
  plot2$plot + theme(legend.position = "none"),
  plot3$plot + theme(legend.position = "none"), 
  plot4$plot + theme(legend.position = "none"),
  plot5$plot + theme(legend.position = "none"),
  plot6$plot + theme(legend.position = "right",
                     legend.direction = "vertical",
                     legend.box = "vertical")
)

spacer <- plot_spacer() + plot_layout(widths = unit(0.5, "cm"))
combined_plot <- (
  plots[[1]] | spacer | plots[[2]]
) /
  (
    plots[[3]] | spacer | plots[[4]]
  ) /
  (
    plots[[5]] | spacer | plots[[6]]
  ) +
  plot_layout(guides = "collect") &
  theme(plot.margin = margin(3, 0, 0, 0))


# Save the final plot
ggsave("combined_survival 2.png", plot = combined_plot,
       width = 1200, height = 600, units = "px", dpi = 300, bg = "white")



sessionInfo()

#End of code 

