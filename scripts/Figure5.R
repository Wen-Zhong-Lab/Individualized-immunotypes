
pack_R <- c("dplyr","ggplot2","ggrepel","umap", "edgeR",
            "RColorBrewer", "pheatmap", "tidyverse",
            "igraph", "biomaRt", 
            "org.Hs.eg.db", "ggpubr", 
            "clusterProfiler", "vcfR",
            "igraph", "readxl")

for (i in 1:length(pack_R)) {
  library(pack_R[i], character.only = TRUE)
}

set.seed(1)

# Load ####
source('C:/Users/albze08/Desktop/postDoc/functions/my_fisher_test.R')

vcf <- read.vcfR("original.genotype.vcf")
clinical <- read.table("original.clinical.txt", sep="\t", header = T)
metadata <- read.table("original.metadata.txt", sep="\t", header = T)

rna.log <- read.table("original.rna.txt", sep="\t", header = T)
cytof <- read.table("original.cytof.txt", sep="\t", header = T)
protein <- read.table("original.protein.txt", sep="\t", header = T)

# Group populations into categories ####

# Natural Killer (NK) Cells
nk_cells <- c(
  "CD56neg_CD57neg_HLADRneg", "CD56neg_CD57neg_HLADRpos", "CD56neg_CD57pos_HLADRneg",
  "CD56neg_CD57pos_HLADRpos", "CD56_bright", "CD56dim_CD57neg_HLADRneg",
  "CD56dim_CD57neg_HLADRpos", "CD56dim_CD57pos_HLADRneg", "CD56dim_CD57pos_HLADRpos"
)

# Monocytes
monocytes_classical <- c("CD16neg_Classical_monocytes")
monocytes_intermediate <- c("Intermediate_Monocytes")
monocytes_nonclassical <- c("Non_classical_monocytes")

# Dendritic Cells (DCs)
mDC <- "mDC"
pDC <- "pDC"

# Basophils
basophils <- c("basophils")

# Innate Lymphoid Cells (ILCs)
ilcs <- c("Non_NK_ILC")

# Naive B Cells
naive_b_cells <- c("Naive_B_cells")

# B Cells
b_cells <- c(
  "Unswitched_memory_B_cells",
  "CD27neg_switched_memory_B_cells_CD22neg_CD24neg", "CD27neg_switched_memory_B_cells_CD22neg_CD24pos",
  "CD27neg_switched_memory_B_cells_CD22pos_CD24neg", "CD27neg_switched_memory_B_cells_CD22pos_CD24pos",
  "plasmablasts", "CD27pos_switched_memory_B_cells_CD22neg_CD24neg",
  "CD27pos_switched_memory_B_cells_CD22neg_CD24pos", "CD27pos_switched_memory_B_cells_CD22pos_CD24neg",
  "CD27pos_switched_memory_B_cells_CD22pos_CD24pos"
)

# Naive T Cells
naive_cd4 <- "CD4pos_naive"
naive_cd8 <- "CD8pos_naive"

# Central Memory (T_CM)
central_memory_cd4 <- c(
  "CD4pos_central_memory_CD57pos_CD39pos", "CD4pos_central_memory_CD57pos_CD39neg",
  "CD4pos_central_memory_CD57neg_CD39neg", "CD4pos_central_memory_CD57neg_CD39pos")

central_memory_cd8 <- c(
  "CD8pos_Central_memory_CD57neg_CD39neg", "CD8pos_Central_memory_CD57neg_CD39pos",
  "CD8pos_Central_memory_CD57pos_CD39neg", "CD8pos_Central_memory_CD57pos_CD39pos"
)

# Effector Memory (T_EM)
effector_memory_cd4 <- c(
  "CD4pos_Effector_memory_CD57pos_CD39neg", "CD4pos_Effector_memory_CD57neg_CD39pos",
  "CD4pos_Effector_memory_CD57neg_CD39neg", "CD4pos_Effector_memory_CD57pos_CD39pos")

effector_memory_cd8 <- c(
  "CD8pos_Effector_memory_CD57pos_CD39pos", "CD8pos_Effector_memory_CD57neg_CD39neg",
  "CD8pos_Effector_memory_CD57neg_CD39pos", "CD8pos_Effector_memory_CD57pos_CD39neg"
)

# TEMRA
temra_cd4 <- c(
  "CD4pos_TEMRA_CD57pos_CD39neg", "CD4pos_TEMRA_CD57neg_CD39pos",
  "CD4pos_TEMRA_CD57neg_CD39neg", "CD4pos_TEMRA_CD57pos_CD39pos"
)

temra_cd8 <- c(
  "CD8pos_TEMRA_CD57pos_CD39pos", "CD8pos_TEMRA_CD57pos_CD39neg",
  "CD8pos_TEMRA_CD57neg_CD39pos", "CD8pos_TEMRA_CD57neg_CD39neg"
)

# Combine into a list
pop.group <- list(
  NK_Cells = nk_cells,
  Monocytes_classical = monocytes_classical,
  Monocytes_intermediate = monocytes_intermediate,
  Monocytes_nonclassical = monocytes_nonclassical,
  Plasmacytoid_DC = pDC,
  Myeloid_DC = mDC,
  Basophils = basophils,
  ILCs = ilcs,
  Naive_B_Cells = naive_b_cells,
  B_Cells = b_cells,
  Naive_CD4 = naive_cd4,
  Naive_CD8 = naive_cd8,
  Central_Memory_CD4 = central_memory_cd4,
  Central_Memory_CD8 = central_memory_cd8,
  Effector_Memory_CD4 = effector_memory_cd4,
  Effector_Memory_CD8 = effector_memory_cd8,
  TEMRA_CD4 = temra_cd4,
  TEMRA_CD8 = temra_cd8
)
pop.ord <- names(pop.group)


cytof.nonnegative <- cytof
cytof.nonnegative[cytof.nonnegative<0] <- 2.2e-16
for (n in 1:nrow(cytof.nonnegative)){
  cytof.nonnegative[n,] <- cytof.nonnegative[n,] / sum(cytof.nonnegative[n,])
}

cytof.group <- matrix(0, nrow(cytof), length(pop.group)) %>% as.data.frame()
rownames(cytof.group) <- rownames(cytof)
colnames(cytof.group) <- names(pop.group)
for (n in 1:ncol(cytof.group)){
  if (length(pop.group[[n]])==1){
    cytof.group[,n] <- cytof.nonnegative[,pop.group[[n]]]
  } else {
    cytof.group[,n] <- rowSums(cytof.nonnegative[,pop.group[[n]]])
  }
}

macro.anno.my <- setNames(rep("", ncol(cytof)), colnames(cytof))
for (pop in names(pop.group)){
  macro.anno.my[names(macro.anno.my) %in% pop.group[[pop]]] <- pop
}


# Extract GWAS results of specific immune families ####

pval.rna <- 5e-8/ncol(rna.log)
pval.cytof <- 5e-8/ncol(cytof)

# gather cQTLs
gather_assoc <- function(folder_QTL, pval_loose=5e-8){
  
  QTL_files <- grep("linear", list.files(folder_QTL), value=T)
  df_summary <- data.frame()
  for (n in 1:length(QTL_files)){
    f <- read.table(paste0(folder_QTL, QTL_files[n]), header=F) 
    colnames(f) <- c("CHR", "BP", "SNP", paste0("V", 1:8), "beta", "V9", "V10", "P", "V11")
    f <- f[f$P<pval_loose,]
    if (nrow(f)>0){
      f$group <- gsub("out\\.", "", QTL_files[n])
      f$group <- gsub("\\.glm\\.linear", "", f$group)
      
      df_summary <- rbind(df_summary, f)
    }
  }
  df_summary$id <- paste0(df_summary$CHR, ":", df_summary$BP)
  df_summary$BP <- as.numeric(df_summary$BP)
  
  return(df_summary)
}

gather_conditional_assoc <- function(folder_QTL, pval_loose=5e-8){
  var <- grep("linear", list.files(paste0(folder_QTL, "conditional")), value=T)
  var <-  gsub("\\..*", "", var)
  QTL_files <- grep(paste(var, collapse="|"), list.files(folder_QTL), value=T)
  QTL_files <- QTL_files[grep("linear", QTL_files)]
  df_summary <- data.frame()
  for (n in 1:length(QTL_files)){
    f <- read.table(paste0(folder_QTL, QTL_files[n]), header=F) 
    colnames(f) <- c("CHR", "BP", "SNP", paste0("V", 1:8), "beta", "V9", "V10", "P", "V11")
    f <- f[f$P<pval_loose,]
    if (nrow(f)>0){
      f$group <- gsub("out\\.", "", QTL_files[n])
      f$group <- gsub("\\.glm\\.linear", "", f$group)
      
      df_summary <- rbind(df_summary, f)
    }
  }
  df_summary$family <- macro.anno.my[match(df_summary$group, names(macro.anno.my))]
  GWAS.cytof <- df_summary
  GWAS.cytof$id <- paste0(GWAS.cytof$CHR, ":", GWAS.cytof$BP)
  GWAS.cytof$BP <- as.numeric(GWAS.cytof$BP)
  
  return(GWAS.cytof)
}

#LD
gather_assoc_LD <- function(folder_QTL, pval.ths){
  QTL_files <- grep("clumps", list.files(paste0(folder_QTL, "LD")), value=T)
  df_summary.QTL.clump <- data.frame()
  for (n in 1:length(QTL_files)){
    var <- gsub("LD\\.out\\.", "", QTL_files[n])
    var <- gsub("\\.glm\\.linear\\.clumps", "", var)
    
    f <- read.table(paste0(folder_QTL, "LD/", QTL_files[n]), header=F)
    colnames(f) <- c("CHR", "BP", "SNP", "P", paste0("V", 1:6))
    f <- f[f$P<pval.ths,]
    if (nrow(f)>0){
      f$group <- var
      df_summary.QTL.clump <- rbind(df_summary.QTL.clump, f[which.min(f$P),c("SNP", "P", "group")])
    }
    if (nrow(f)>1){
      conditional_file <- grep(paste0(var, "\\."), list.files(paste0(folder_QTL, "conditional")), value=T)
      conditional_file <- grep("linear", conditional_file, value=T)
      if (length(conditional_file)==1){
        f_conditional <- read.table(paste0(folder_QTL, "conditional/", conditional_file), header=F)
        colnames(f_conditional) <- c("CHR", "BP", "SNP", paste0("V", 1:11), "P", "ERRCODE")
        f_conditional$group <- var
        
        idx <- which(f_conditional$P<0.1) #this is the conditional p-value, i.e. when doing linear model that includes the dominant variant
        f_conditional <- f_conditional[idx, ]
        f_conditional$P <- f$P[match(f_conditional$SNP, f$SNP)] #store the actual p-value from GWAS, not the conditional
        df_summary.QTL.clump <- rbind(df_summary.QTL.clump, f_conditional[,c("SNP", "P", "group")])
      } else {
        print(paste0("[WARNING] Number of files found for ", var, " : ", length(conditional_file)))
      }
    }
  }
  GWAS.cytof.LD <- df_summary.QTL.clump
  
  str <- strsplit(GWAS.cytof.LD$SNP, "\\:") %>% unlist()
  GWAS.cytof.LD$CHR <- str[seq(1,length(str),4)] %>% as.numeric()
  GWAS.cytof.LD$BP <- str[seq(2,length(str),4)] %>% as.numeric()
  GWAS.cytof.LD$id <- paste0(GWAS.cytof.LD$CHR, ":", GWAS.cytof.LD$BP)
  
  GWAS.cytof.LD <- GWAS.cytof.LD %>% arrange(CHR,BP,P)
  return(GWAS.cytof.LD)
}

closest_snp_fun <- function(df1, df2){
  df.out <- data.frame()
  for (n in 1:nrow(df1)){
    idx <- which(df2$CHR == df1$CHR[n])
    d <- abs(df1$BP[n] - df2$BP[idx])
    
    dist.min <- min(d) #to account for multiple matches with the same distance
    idx.min.vec <- idx[which(d==dist.min)]
    
    for (idx.min in idx.min.vec){
      df.out <- rbind(df.out, 
                      data.frame(df1[n,],
                                 idx=idx.min,
                                 distance=dist.min,
                                 CHR_other=df2$CHR[idx.min],
                                 BP_other=df2$BP[idx.min]) )
    }
  }
  return(df.out %>% arrange(distance))
}

GWAS.lmm <- gather_assoc("cQTL-mean-53/", pval_loose = 1e-6)
GWAS.lmm$family <- macro.anno.my[match(GWAS.lmm$group, names(macro.anno.my))]

GWAS.lmm.LD <- gather_assoc_LD("cQTL-mean-53/", pval.cytof)
GWAS.lmm.LD$family <- macro.anno.my[match(GWAS.lmm.LD$group, names(macro.anno.my))]
GWAS.lmm.LD$beta <- GWAS.lmm$beta[match(GWAS.lmm.LD$SNP, GWAS.lmm$SNP)]

SNP_LD <- GWAS.lmm$SNP %>% unique()

GWAS.lmm.eqtl <- gather_assoc("eQTL-mean-8K/", pval_loose = 5e-8)
GWAS.lmm.LD.eqtl <- gather_assoc_LD("eQTL-mean-8K/", pval.rna)
GWAS.lmm.LD.eqtl$beta <- GWAS.lmm.eqtl$beta[match(GWAS.lmm.LD.eqtl$SNP, GWAS.lmm.eqtl$SNP)]

GWAS.lmm.eqtl.filt <- GWAS.lmm.eqtl %>% filter(P<5e-8)
GWAS.shared <- closest_snp_fun(GWAS.lmm.LD %>% filter(beta>0), GWAS.lmm.eqtl.filt) 
GWAS.shared$gene <- GWAS.lmm.eqtl.filt$group[GWAS.shared$idx]
GWAS.shared$gene[GWAS.shared$distance==0] %>% unique()

# add gene column
mart <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl", version=115)
genes <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol",
                 "chromosome_name", "start_position", "end_position", "strand"),
  mart = mart
)

gene.anno <- genes %>% filter(chromosome_name %in% 1:22) %>%
  filter(hgnc_symbol!="") %>%
  mutate(CHR=chromosome_name, BP=(start_position+end_position)/2)

snp_gene_anno <- closest_snp_fun(GWAS.lmm.LD, gene.anno)
snp_gene_anno$gene <- gene.anno$hgnc_symbol[snp_gene_anno$idx]

GWAS.lmm.LD$gene <- snp_gene_anno$gene[match(GWAS.lmm.LD$SNP, snp_gene_anno$SNP)]


GWAS.lmm$gene <- snp_gene_anno$gene[match(GWAS.lmm$SNP, snp_gene_anno$SNP)]
GWAS.lmm$gene[GWAS.lmm$P > pval.cytof] <- NA

GWAS.lmm.eqtl$gene <- ifelse(GWAS.lmm.eqtl$P<pval.rna, GWAS.lmm.eqtl$group, NA)


# Supplementary Table 4 ####
df <- GWAS.lmm.LD[, c("SNP", "beta", "P", "group", "family", "gene")]
colnames(df) <- c("SNP", "Coefficient", "P-value", "Population", "Major population", "Closest gene")

df2 <- GWAS.lmm.LD.eqtl[, c("SNP", "beta", "P", "group")]
colnames(df2) <- c("SNP", "Coefficient", "P-value", "Gene")

wb <- createWorkbook()
addWorksheet(wb, "Immune frequencies GWAS")
addWorksheet(wb, "RNA-seq GWAS")
writeData(wb, "Immune frequencies GWAS", df)
writeData(wb, "RNA-seq GWAS", df2)

saveWorkbook(wb, "Supplementary Table 4.xlsx", overwrite = TRUE)


#Manhattan plot function ####
manhattan_plot <- function(df, pval.ths=5e-8, MAX_SNP=5e4){
  
  #Prepare GWAS CyTOF
  df_plot <- df[order(df$P, decreasing=F), c("CHR", "BP", "P", "group", "family", "SNP", "gene")]
  df_plot <- df_plot[!df_plot$CHR %in% c("PAR1", "PAR2"),]
  df_plot$CHR <- as.numeric(df_plot$CHR)
  
  #merge
  chr.df <- data.frame(CHR=1:22)
  chr.df$BP_max <- 0
  chr.df$BP_mean <- 0
  df_plot$BP_cum <- 0
  max.previous.chr <- 0
  
  stopifnot(is.numeric(df_plot$CHR))
  for (n in 1:22){ #chromosomes
    idx <- which(df_plot$CHR==n)
    
    #new BP in plot df
    df_plot$BP_cum[idx] <- df_plot$BP[idx] + max.previous.chr 
    
    #chromosome charac
    chr.df$BP_max[n] <- max(df_plot$BP_cum[idx])
    chr.df$BP_mean[n] <- mean(df_plot$BP_cum[idx])
    
    #update baseline
    max.previous.chr <- max(df_plot$BP_cum[idx])
  }
  df_plot$color <- ifelse(df_plot$CHR %% 2 ==0, "even", "odd")

  p <- ggplot(df_plot, aes(x=BP_cum, y=-log10(P), fill=color)) +
    geom_point(shape=21) + 
    geom_hline(yintercept=-log10(pval.ths), linetype=2) +
    theme_bw() +
    geom_text(aes(label=gene), cex=2, color="black") +
    scale_x_continuous(label = chr.df$CHR, breaks = chr.df$BP_mean) +
    scale_fill_manual(values=c(even="gray77", odd="gray23")) + 
    scale_size_continuous(range = c(0.5,3)) + 
    labs(x = NULL, y = "-log<sub>10</sub>(p)") +
    theme(panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.x = element_text(angle=0, size=8, vjust=0.5),
          legend.position="none",
          text=element_text(size=12)) +
    ylab("-Log10 P-value")
  #p  
  
  return(list(df=df_plot, p=p, chr.df=chr.df))
}

# Combined frequency-expression Manhattan ####
GWAS.both <- rbind(GWAS.lmm %>% mutate(omic="CyTOF"),
                   GWAS.lmm.eqtl %>% mutate(family="", omic="RNA")) %>% arrange(P)
idx.cytof <- which(GWAS.both$omic=="CyTOF")
idx.rna <- which(GWAS.both$omic=="RNA")

manhattan.both <- manhattan_plot(GWAS.both, MAX_SNP=Inf)
manhattan.cytof <- manhattan.both$df[idx.cytof,] %>% arrange(P)

manhattan.rna <- manhattan.both$df[idx.rna,] %>% arrange(P)
manhattan.rna <- manhattan.rna[1:3e4,]

chr.df <- manhattan.both$chr.df

# CyTOF
manhattan.cytof$gene <- ifelse(manhattan.cytof$P<pval.cytof, manhattan.cytof$gene, NA)

pop.palette <- c(
  "Monocytes_classical" = "#C24100", "Monocytes_intermediate"  = "#E66101", "Monocytes_nonclassical"  = "#FDB863",  
  "Basophils" = "#FFD92F", "Myeloid_DC" = "#E6AB02", "Plasmacytoid_DC" = "#B8A200",   
  "B_Cells" = "dodgerblue4", "Naive_B_Cells" = "dodgerblue2",  
  "Naive_CD4" = "#33A02C", "Central_Memory_CD4" = "#66C2A5","Effector_Memory_CD4" = "#1B9E77", "TEMRA_CD4" = "#41AE76",  
  "Naive_CD8" = "#FB9A99", "Central_Memory_CD8" = "#E31A1C", "Effector_Memory_CD8" = "#B2182B", "TEMRA_CD8" = "#67001F",  
  "NK_Cells" = "#CC78BC",  
  "ILCs" = "#B3B3B3"   
)

manhattan.cytof$color <- ifelse(manhattan.cytof$P<pval.cytof, manhattan.cytof$family, manhattan.cytof$color)

p1 <- ggplot(manhattan.cytof, aes(x=BP_cum, y=-log10(P), fill=color)) +
  geom_point(shape=21) + 
  geom_hline(yintercept=-log10(pval.cytof), linetype=2) +
  theme_bw() +
  geom_text(aes(label=gene), cex=2, color="black") +
  scale_x_continuous(label = chr.df$CHR, breaks = chr.df$BP_mean) +
  scale_fill_manual(values=c(pop.palette, even="white", odd="gray23")) + 
  scale_size_continuous(range = c(0.5,3)) + 
  labs(x = NULL, y = "-log<sub>10</sub>(p)") +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(angle=0, size=8, vjust=0.5),
        legend.position="none",
        text=element_text(size=12)) +
  ylab("-Log10 P-value")

#RNA-seq
df_min_p <- manhattan.rna %>%
  filter(P<1e-30) %>%
  group_by(group) %>%
  slice_min(P, n = 1, with_ties = FALSE) %>%
  ungroup()

cell.gene.best <- cell.gene %>% group_by(gene) %>% slice_min(pval, n=1) %>% ungroup()
manhattan.rna$pop.network <- cell.gene.best$pop[match(manhattan.rna$group, cell.gene.best$gene)]
table(manhattan.rna$pop.network) %>% sort()

manhattan.rna$color[!is.na(manhattan.rna$pop.network)] <- manhattan.rna$pop.network[!is.na(manhattan.rna$pop.network)]


p2 <- ggplot(manhattan.rna, aes(x=BP_cum, y=-log10(P), fill=color)) +
  geom_point(data=manhattan.rna %>% filter(!color %in% names(pop.palette)), shape=21) + 
  geom_point(data=manhattan.rna %>% filter(color %in% names(pop.palette)), shape=21) + 
  geom_hline(yintercept=-log10(pval.rna), linetype=2) +
  theme_bw() +
  geom_text(data=df_min_p, aes(label=group), cex=2, color="black") +
  scale_x_continuous(label = chr.df$CHR, breaks = chr.df$BP_mean) +
  scale_fill_manual(values=c(pop.palette, even="white", odd="gray23")) + 
  scale_size_continuous(range = c(0.5,3)) + 
  labs(x = NULL, y = "-log<sub>10</sub>(p)") +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(angle=0, size=8, vjust=0.5),
        legend.position="none",
        text=element_text(size=12)) +
  ylab("-Log10 P-value")

#Plot
pdf("manhattan_CyTOF.pdf", width=7, height=4)
p1
dev.off()

pdf("manhattan_RNA.pdf", width=7, height=4)
p2
dev.off()


# Load Genotype ####
snp <- vcf@fix %>% as.data.frame()
geno <- vcf@gt %>% as.data.frame()
rm(vcf)
geno <- geno[,2:ncol(geno)]
colnames(geno) <- anno$Subject[match(colnames(geno), paste0("P", anno$Barcode))]

stopifnot(all(SNP_LD %in% snp$ID))

idx <- match(SNP_LD, snp$ID) 
snp <- snp[idx,]
geno <- geno[idx,]


# Transform geno in a matrix ####
geno.X <- matrix(NA, nrow(geno), ncol(geno))
rownames(geno.X) <- rownames(geno)
colnames(geno.X) <- colnames(geno)
for (n in 1:nrow(geno.X)){
  vec <- geno[n,]
  vec[vec=="0/0"] <- 0
  vec[vec=="0/1"] <- 1
  vec[vec=="1/1"] <- 2
  vec <- as.numeric(vec)
  geno.X[n,] <- vec
}


# Init ####
var.group <- data.frame(var=c("subject_id",
                              "Age_at_Visit", 
                              "Gender", 
                              "Height", "Weight", "Waist", "Hip", "BMI", "Bioimp_fat", "Bioimp_muscle", "Bioimp_bone", 
                              "SBP", "DBP", 
                              "Cap_Gluc", "Gluc", "HbA1c", 
                              "Chol", "LDL", "HDL", "TG", "ApoB", "ApoA1", "ApoB.apoA1",   
                              "CRP", 
                              "ALAT", "GGT", 
                              "ProBNP", "TNT", 
                              "Urate", 
                              "CystC", "Crea",
                              "Smoking", 
                              "MaritalStatus_change", "Housing_change", "Employment_change",
                              "Stress",
                              "PhysicalActivity", "SedentaryTime_hours",
                              "PerceivedHealth",
                              "Common_cold_influenzae",
                              "NSAID_painmed", "Bp_med", "Lipid_med", "Antibiotics_med"),
                        group=c("subject_id",
                                "Age", 
                                "Sex", 
                                "Body composition", "Body composition", "Body composition", "Body composition", "Body composition", "Body composition", "Body composition", "Body composition", 
                                "Blood pressure", "Blood pressure", 
                                "Glucose homeostasis", "Glucose homeostasis", "Glucose homeostasis", 
                                "Lipid profile", "Lipid profile", "Lipid profile", "Lipid profile", "Lipid profile", "Lipid profile", "Lipid profile",   
                                "CRP", 
                                "Organ biomarker", "Organ biomarker", 
                                "Organ biomarker", "Organ biomarker", 
                                "Organ biomarker", 
                                "Organ biomarker", "Organ biomarker",
                                "Lifestyle", 
                                "Lifestyle", "Lifestyle", "Lifestyle",
                                "Lifestyle",
                                "Lifestyle", "Lifestyle",
                                "Lifestyle",
                                "Lifestyle",
                                "Lifestyle", "Lifestyle", "Lifestyle", "Lifestyle"))



var.group.SNP <- data.frame(var=unique(SNP_LD),
                            group="Genetics")
var.group.SNP$var <- paste0("SNP:", var.group.SNP$var)
var.group.SNP$var <- gsub(":", "\\.", var.group.SNP$var)

var.group <- rbind(var.group, var.group.SNP)


get_mcfadden <- function(df){
  null_model <- glm(y ~ 1, data = df, family = binomial)
  ll_null <- logLik(null_model)
  
  # Fit full model
  full_model <- glm(formula(df), data = df, family = binomial)
  ll_full <- logLik(full_model)
  
  # McFadden R^2
  R2_mcfadden <- 1 - ll_full/ll_null
  
  return(R2_mcfadden)
}

get_r2 <- function(df){
  lm.out <- lm(formula(df), df) %>% summary()
  r2 <- lm.out$r.squared
  return(r2)
}

color_palette <- c("darkseagreen",	"#A30059",	"#9e5a28ff",		
                   "tomato",
                   "#4a6fe3",	"lightblue",	"steelblue" ,
                   "bisque",	"#1CE6FF",
                   "goldenrod")
names(color_palette) <- c("Genetics", "Sex", "Age", 
                          "Blood pressure", 
                          "Body composition", "Lipid profile","Glucose homeostasis",
                          "Organ biomarker", "CRP",
                          "Lifestyle")


clinical <- metadata %>%
  subset(select=-c(SCAPISorIGT_id, subject, Study, Visitdate, id, visit, Birthdate, subject_id,
                   Hb, WBC, Plt, RBC, Hct, MCV, MCH, MCHC, Neut, Lymph, Mono, Eos, Baso,
                   Housing_change, Housing_current, Address_change, MaritalStatus_change, MaritalStatus_current, ShareHousehold_change, 
                   ShareHousehold_current, Employment_change, Employment_current, Tobacco_change, Tobacco_current, PerceivedHealth,       
                   Stress, PhysicalActivity, SedentaryTime_hours, SedentaryTime_minutes, SedentaryTime_unknown, TravelAbroad,          
                   TravelAbroad_country, Animals, Animals_type,
                   Calcitriol, Calcidiol, Health_Other, Common_cold_influenzae, NSAID_painmed, Bp_med, Lipid_med,             
                   Antibiotics_med, Diab_med, Med_details, Smoking, Smoking_hours)) %>%
  mutate(Gender=ifelse(Gender=="m",1,0))

# Plot contribution of immune composition and proteomics to the clustering ####
get_auc_crossfold <- function(df, nfold=5){
  
  idx <- sample(rep(1:nfold, length.out = nrow(df)))
  mean_auc <- rep(NA, nfold)
  for (k in 1:nfold){
    idx.valid <- which(idx==k)
    idx.train <- setdiff(1:nrow(df), idx.valid)
    
    df.train <- df[idx.train,]
    df.valid <- df[idx.valid,]
    
    model.train <- glm(formula(df.train), data = df.train, family = binomial)
    pred.valid <- predict(model.train, df.valid[,-1], type = "response")
    pr <- prediction(pred.valid, df.valid$y)
    auc <- performance(pr, measure = "auc")
    mean_auc[k] <- auc@y.values[[1]]
  }
  auc <- mean(mean_auc)
  
  return(auc)
}

  
  df.plot <- shap_values %>% melt()
  
  p <- ggplot(df.plot, aes(x=value, y=Var2)) + 
    geom_quasirandom(aes(color=value)) +
    scale_color_gradient(low="royalblue", high="tomato") +
    geom_vline(xintercept=0, linetype=2) +
    theme_classic() + theme(legend.position="bottom") +
    ylab("") + xlab("SHAP value")
  p
  return(list(p=p, df=df.plot))
}


df.out <- data.frame()
common.samples <- intersect(rownames(cytof.group), rownames(protein)) %>% intersect(primary_cluster$sample)
cluster_id <- unique(primary_cluster$primary_cluster)
for (id in cluster_id){

  X.full <- cbind(
    clinical[common.samples, c("Gender", "Age_at_Visit", "BMI")],
    cytof.group[common.samples,], 
    protein[common.samples,] %>% makeX(na.impute = T) %>% 
      condense_PCA(Npca = 20) %>% as.data.frame() |> rename_with(~ paste0(.x, "protein")) %>% scale()
  ) %>% as.data.frame()
  
  
  df.y <- ifelse(primary_cluster[common.samples, "primary_cluster"] == id, 1, 0) %>% as.data.frame() |> rename_with(~"y")
  df <- cbind(df.y, X.full)
  
  #prefiltering
  lasso_model <- cv.glmnet(X.full %>% as.matrix(), df$y, alpha = 1, family = "binomial")
  out.lasso <- coef(lasso_model, s = "lambda.min") 
  selected_var <- rownames(out.lasso)[as.numeric(out.lasso) != 0] %>% setdiff("(Intercept)" ) 
  
  selected_var <- c("y", selected_var %>% union(c("Gender", "Age_at_Visit", "BMI")))
  stopifnot(all(selected_var %in% colnames(df)))
  df <- df[, selected_var]
  
  # do McFadden
  Niter <- 100
  for (iter in 1:Niter){
    
    new_order <- sample(setdiff(colnames(df), "y"))
    df.shuffle <- df[, c("y", new_order)]
    
    df.n <- df.shuffle[,1,drop=F]
    r2_prev <- 0
    for (n in 2:ncol(df.shuffle)){
      df.n <- cbind(df.n, df.shuffle[,n,drop=F])
      r2 <- get_mcfadden(df.n)
      if (r2-r2_prev<0){
        #stopifnot(1<0)
        print(paste0(id, ": ", colnames(df.shuffle)[n]))
        df.n <- df.n[,1:(ncol(df.n)-1)]
      } else {
        df.out <- rbind(df.out,
                        data.frame(cluster=id, r2_cum=r2, r2=r2-r2_prev, var=colnames(df.shuffle)[n]))
        r2_prev <- r2
      }
    }
  }
}
df.out.summ <- df.out %>% group_by(var, cluster) %>% summarise(r2=mean(r2))

df.out.summ$group <- var.group$group[match(df.out.summ$var, var.group$var)]
df.out.summ$group <- NA
df.out.summ$group[df.out.summ$var == "Gender"] <- "Sex"
df.out.summ$group[df.out.summ$var == "Age_at_Visit"] <- "Age"
df.out.summ$group[df.out.summ$var == "BMI"] <- "BMI"
df.out.summ$group[df.out.summ$var %in% colnames(cytof.group)] <- "CyTOF"
df.out.summ$group[grep("protein", df.out.summ$var)] <- "Protein"

df.plot <- df.out.summ %>% group_by(group, cluster) %>% summarise(r2=sum(r2)) 
ord.group <- df.plot %>% group_by(group) %>% summarise(mean_r2=mean(r2)) %>% arrange(mean_r2) %>% pull(group)

df.plot$group <- factor(df.plot$group, levels=ord.group)
p <- ggplot(df.plot, aes(x=cluster, y=r2, fill=group)) +
  geom_bar(stat="identity", position="stack", color="black") +
  scale_fill_manual(values=c(CyTOF="royalblue", Protein="tomato", Sex="#A30059", Age="#9e5a28ff", BMI="green")) + 
  theme_classic() +
  theme(legend.position="none") +
  theme(plot.title = element_text(hjust = 0.5)) +
  ylab("Cumulative McFadden R2") + xlab("") + ggtitle("Clusters")
p

pdf("clusters_contributions_omics.pdf", height=5, width=3)
print(p)
dev.off()


# Contribution of genetic and non-genetic to clusters ####
common.samples <- intersect(rownames(primary_cluster), rownames(clinical)) 
X <- clinical[common.samples,] %>% scale()
snp_n <- GWAS.lmm.LD$SNP %>% unique()

get_geno <- function(snp_n, X){
  idx <- match(snp_n, snp$ID)
  geno_n <- geno[idx,] %>% t() %>% as.data.frame()
  colnames(geno_n) <- paste0("SNP:", snp_n)
  colnames(geno_n) <- gsub(":", "\\.", colnames(geno_n) ) #formula transforms them into . anyway
  
  subject <- gsub("\\:.*", "", rownames(X))
  X.geno <- geno_n[subject,,drop=F] 
  X.geno[X.geno=="0/0"] <- 0
  X.geno[X.geno=="0/1"] <- 1
  X.geno[X.geno=="1/1"] <- 2
  for (n in 1:ncol(X.geno)){
    X.geno[,n] <- as.numeric(X.geno[,n])
  }
  return(X.geno)
}
X <- cbind(X, get_geno(snp_n, X))


Niter <- 100
df.out <- data.frame()
cluster_id <- unique(primary_cluster$primary_cluster)
for (id in cluster_id){
  
  df <- data.frame(y=ifelse(primary_cluster[common.samples, "primary_cluster"]==id, 1, 0), X)
  df$y <- as.factor(df$y)
  
  #prefiltering
  lasso_model <- cv.glmnet(X %>% makeX(na.impute = T), df$y, alpha = 1, family = "binomial")
  out.lasso <- coef(lasso_model, s = "lambda.min") 
  selected_var <- rownames(out.lasso)[as.numeric(out.lasso)!=0] %>% setdiff("(Intercept)" )
  
  stopifnot(all(selected_var %in% colnames(df)))
  df <- df[, c("y", selected_var)]
  for (iter in 1:Niter){
    new_order <- sample(setdiff(colnames(df), c("y", "Gender", "Age_at_Visit")))
    df.shuffle <- df[, c("y", "Gender", "Age_at_Visit", new_order)]
    
    r2_prev <- 0
    for (n in 4:ncol(df.shuffle)){
      df.n <- df.shuffle[,1:n]
      r2 <- get_mcfadden(df.n)
      if (r2-r2_prev<0){
        stopifnot(1<0)
      }
      
      df.out <- rbind(df.out, 
                      data.frame(cluster=id, r2=r2-r2_prev, var=colnames(df.shuffle)[n]))
      r2_prev <- r2
    }
  }
}
df.out.summ <- df.out %>% group_by(var, cluster) %>% summarise(r2=mean(r2))

df.out.summ$group <- var.group$group[match(df.out.summ$var, var.group$var)]

df.plot <- df.out.summ %>% group_by(group, cluster) %>% summarise(r2=sum(r2)) 
ord.group <- df.plot %>% group_by(group) %>% summarise(mean_r2=mean(r2)) %>% arrange(mean_r2) %>% pull(group)

df.plot$group <- factor(df.plot$group, levels=ord.group)
p <- ggplot(df.plot, aes(x=cluster, y=r2, fill=group)) +
  geom_bar(stat="identity", position="stack", color="black") +
  scale_fill_manual(values=color_palette) + 
  theme_classic() +
  theme(legend.position="none") +
  theme(plot.title = element_text(hjust = 0.5)) +
  ylab("Cumulative McFadden R2") + xlab("") + ggtitle("Clusters")
p

pdf("clusters_contributions_2.pdf", height=5, width=3)
print(p)
dev.off()


# Contribution of genetic and non-genetic to major pops ####
GWAS.cytof.group.LD <- GWAS.lmm.LD
for (g in names(pop.group)){
  GWAS.cytof.group.LD$group[GWAS.lmm.LD$group %in% pop.group[[g]]] <- g
}

common.samples <- intersect(rownames(cytof.group), rownames(clinical))
X <- clinical[common.samples,] %>% scale() %>% as.data.frame()

Niter <- 100
major.pop <- colnames(cytof.group)
df.out <- data.frame()
for (pop in major.pop){
  
  snp_n <- GWAS.cytof.group.LD$SNP[GWAS.cytof.group.LD$group == pop]
  if (length(snp_n)>0){
    geno_pop <- get_geno(snp_n, X)
    X_pop <- cbind(X, geno_pop)
  } else {
    X_pop <- X
  }
  df <- data.frame(y=cytof.group[common.samples, pop], X_pop)
  
  for (iter in 1:Niter){
    
    new_order <- sample(setdiff(colnames(df), c("y", "Gender", "Age_at_Visit")))
    df.shuffle <- df[, c("y", "Gender", "Age_at_Visit", new_order)]
    
    r2_prev <- 0
    for (n in 4:ncol(df.shuffle)){
      df.n <- df.shuffle[,1:n]
      r2 <- get_r2(df.n)
      if (r2-r2_prev<0){
        #stopifnot(1<0)
        print(paste0(pop, ": ", colnames(df.shuffle)[n]))
      }
      
      df.out <- rbind(df.out, 
                      data.frame(pop=pop, r2=r2-r2_prev, var=colnames(df.shuffle)[n]))
      r2_prev <- r2
    }
  }
}
df.out.summ <- df.out %>% group_by(var, pop) %>% summarise(r2=mean(r2))
df.out.summ$group <- var.group$group[match(df.out.summ$var, var.group$var)]

df.plot <- df.out.summ %>% group_by(group, pop) %>% summarise(r2=sum(r2)) 

ord.group <- df.plot %>% group_by(group) %>% summarise(mean_r2=mean(r2)) %>% arrange(mean_r2) %>% pull(group)
ord.pop <- df.plot %>% group_by(pop) %>% summarise(tot_r2=sum(r2)) %>% arrange(tot_r2) %>% pull(pop) %>% rev()

df.plot$pop <- factor(df.plot$pop, levels=ord.pop)
df.plot$group <- factor(df.plot$group, levels=ord.group)

p <- ggplot(df.plot, aes(x=pop, y=r2, fill=group)) +
  geom_bar(stat="identity", position="stack", color="black") +
  scale_fill_manual(values=color_palette) + 
  theme_classic() +
  theme(legend.position="none") +
  theme(plot.title = element_text(hjust = 0.5)) +
  ylab("Cumulative R2") + xlab("") + ggtitle("Major populations")
p

pdf("contribution_major_pop_2.pdf", width=4, height=3)
print(p)
dev.off()
    
# Plot examples of genotype--> frequency ####
plot_SNP_frequency <- function(snp_id, pop, major.pop=NULL){
  idx.snp <- which(snp$ID==snp_id)
  if (is.null(major.pop)){
    df.plot <- data.frame(sample=rownames(cytof.nonnegative), value=dfcytof.nonnegative[,pop], omic="CyTOF")
    pop.name <- pop
  } else {
    df.plot <- data.frame(sample=rownames(cytof.nonnegative), value=cytof.group[,major.pop], omic="CyTOF")
    pop.name <- major.pop
  }
  df.plot$rank <- 100*(df.plot$value %>% rank())/nrow(df.plot)
  df.plot$ind <- gsub("\\:.*", "", df.plot$sample)
  df.plot$visit <- gsub(".*\\:", "", df.plot$sample)
  df.plot$geno <- geno[idx.snp, df.plot$ind] %>% as.character()
  df.plot <- na.omit(df.plot)
  df.plot$geno <- ifelse(df.plot$geno =="0/0", 0, 1) %>% as.character()
  
  df.plot.summ <- df.plot %>% group_by(geno, visit) %>% 
    summarise(min_value=min(value), max_value=max(value), median_value=median(value), low_05=quantile(value, 0.025), high_05=quantile(value, 0.975)) %>% 
    as.data.frame()
  
  p4 <- ggplot(df.plot.summ, aes(x = visit, fill = geno, group = geno)) +
    geom_ribbon(aes(ymin = low_05, ymax = high_05), alpha = 0.3)  +
    geom_line(aes(group=geno, y=median_value, color=geno)) +
    geom_point(size=4, aes(y=median_value), shape=21) +
    theme_classic() +
    theme(legend.position = "bottom") +
    ylab("Percentile in the cohort") + xlab(snp_id) +
    ggtitle(pop.name)
  p4
  
  return(list(df=df.plot, p4=p4))
}

GWAS.lmm.LD <- GWAS.lmm.LD %>% arrange(desc(beta))

GWAS.lmm.LD.B <- GWAS.lmm.LD %>% filter(family=="B_Cells") %>% arrange(beta)
pdf(paste0("pdfs/2610-Bcells-major", ".pdf"), height=3, width=4)
for (n in 1:nrow(GWAS.lmm.LD.B)){
  
  snp_n <- GWAS.lmm.LD.B$SNP[n]
  pop_n <- GWAS.lmm.LD.B$group[n]
  
  out <- plot_SNP_frequency(snp_n, pop_n, major.pop=GWAS.lmm.LD.B$family[n])
  p4 <- out$p4
  print(p4)
}
dev.off()


# B cell freq based on PGS ####
major.pop <- "B_Cells"
df.snp <- GWAS.lmm.LD %>% filter(family == major.pop)
df.snp$beta_weighted <- NA

pop.spec <- unique(df.snp$group)
for (pop in pop.spec){
  idx <- which(df.snp$group == pop)
  relative.freq <- mean(cytof.nonnegative[,pop])/mean(cytof.group[,major.pop])
  
  df.snp$beta_weighted[idx] <- df.snp$beta[idx]*relative.freq
}

idx.geno <- match(df.snp$SNP, snp$ID)
PRS <- apply(geno.X[idx.geno,], 2, function(x){sum(x*df.snp$beta_weighted, na.rm=T)}) # scale by relative freq of that specific pop?

pheno <- data.frame(y=cytof.group[,major.pop], sample=rownames(cytof.group))
pheno$ind <- gsub("\\:.*", "", pheno$sample)
pheno$PRS <- PRS[pheno$ind]
pheno$PRS_rank <- rank(pheno$PRS)
pheno$PRS_norm <- scale(pheno$PRS) %>% as.numeric()

ind.uniq <- unique(pheno$ind)
color.UMAP <- grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)] %>% sample(length(ind.uniq))
color.UMAP <- setNames(color.UMAP, ind.uniq)

p <- ggplot(pheno, aes(x=PRS, y=y)) +
  geom_point(aes(fill=ind), shape=21, size=2) +
  geom_smooth(method="lm", color="red", size=1) +
  geom_text(label=paste0("Spearman's r= ", cor(pheno$PRS, pheno$y, method="spearman") %>% round(digits=2)), x=5, y=0.03) +
  scale_fill_manual(values=color.UMAP) +
  theme_classic() + theme(legend.position = "none") +
  xlab("Polygenic score") + ylab("Memory B cells frequency")
p

pdf("pdfs/PGS-memory B-mean.pdf", width=3, height=3)
print(p)
dev.off()

cor(pheno$y, pheno$PRS_norm, method="spearman")

#PGS based on cluster
pheno$cluster <- primary_cluster$primary_cluster[match(pheno$sample, primary_cluster$sample)]

comparisons <- list(c(1,2), c(1,3), c(2,3))
p <- ggplot(pheno %>% na.omit(), aes(x=cluster, y=PRS_rank)) +
  geom_boxplot() +
  geom_quasirandom(shape=21, aes(fill=cluster), size=2) +
  scale_fill_manual(values=c("1"="#E69F00", "2"="#56B4E9", "3"="#009E73")) + 
  stat_compare_means(method="wilcox.test", comparisons=comparisons, size=3) +
  theme_classic()
p
pdf("PGS-cluster.pdf", width=4, height=3)
print(p)
dev.off()


# PRS and gene expression ####
Bcell.g <- cell.gene %>% filter(pop %in% c("Naive_B_Cells", "B_Cells")) %>% pull(gene) %>% unique()
df.DEG <- data.frame()
for (g in Bcell.g){
  
  df <- pheno
  df$gene <- rna.log[match(df$sample, rownames(rna.log)),g] 
  
  lmFit <- lm(gene~PRS_norm, df) %>% summary() %>% coefficients() %>% as.data.frame() %>% mutate(gene=g)
  
  df.DEG <- rbind(df.DEG,
                  lmFit["PRS",])
}
df.DEG$adj.pval <- p.adjust(df.DEG$`Pr(>|t|)`, method="BH")
df.DEG$fill <- ifelse(df.DEG$adj.pval<0.05, "notDEG", "DEG")

df.DEG <- df.DEG %>% arrange(adj.pval)

df.DEG$label <- NA
df.DEG$label[1:50] <- df.DEG$gene
p <- ggplot(df.DEG, aes(x=Estimate, y=-log10(adj.pval), fill=fill)) +
  geom_point(shape=21, size=2) +
  geom_text(aes(label=label), size=2) +
  geom_vline(xintercept=0, linetype=2) +
  geom_hline(yintercept=0, linetype=2) +
  geom_hline(yintercept=-log10(0.05), linetype=2) +
  scale_fill_manual(values=c(notDEG="gray77", DEG="gray35")) +
  theme_classic() + theme(legend.position = "none")
p

pdf("pdfs/PGS-memory B-genes-mean.pdf", width=3, height=3)
print(p)
dev.off()
    

