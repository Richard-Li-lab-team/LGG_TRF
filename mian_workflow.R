#
##
###   LGG(Normal, pLGG, rLGG, sGBM) DEGs

setwd('/export3/zhuyw/Atask/20211119_CBTeng-sGBM/runtime')
save.image("working.RData") #保存工作空间

par(mfrow = c(1,1))
par(pin=c(3,3))

## 基本配置
# 设置颜色
bioCol <- c('skyblue3', "deeppink3","deepskyblue3",'orange3', 
            'tomato3', 'turquoise3', 'seagreen3')
#  19种颜色
bioCol =c("#999999","#FF0099", "#E69F00", "#56B4E9", "#009E73", "#F0E442", 
          "#0072B2", "#D55E00", "#CC79A7","#990000","#9900cc","#66FF66",
          "#663300","#0000FF","#CC0033","#FF0000","#000099","#660066","#333333")
barplot(1:20,col=bioCol)


###
##########      数据的读取与保存      ##########
#
##   基因列表
DNA5mC <- read.table('/export2/zhuyw/database/gencode/DNA_5mC_gene.list', 
                     header = T, sep = '\t')

TNF_Gene <- read.table('/export2/zhuyw/database/gencode/TNF_family_gene.list', 
                       header = F,  sep = ' ')

#
##   DataSet01-GTEX_brain (正常脑组织)
GTEX_brain_exp <- read.table('/export2/zhuyw/database/GTEX/Brain_sample_exp.tsv', 
                             header = T, sep = '\t')
tmp_exp <- GTEX_brain_exp
tmp_exp <- avereps(tmp_exp, ID = tmp_exp[,1])
gene_symbol <- gsub(tmp_exp[,1], pattern = '-', replacement = '_')
tmp_exp <- apply(tmp_exp[,c(-1)], 2, as.numeric)
rownames(tmp_exp) <- gene_symbol
GTEX_brain_exp <- as.data.frame(tmp_exp)


#
##   DataSet02-TCGA-GBMLGG (TCGA胶质瘤)
TCGA_GBMLGG_exp <- read.table('/export2/zhuyw/database/TCGA/GBMLGG/TCGA_GBMLGG_HiSeqV2.gz', 
                              header = T, sep = '\t')
tmp_exp <- TCGA_GBMLGG_exp
tmp_exp <- avereps(tmp_exp, ID = tmp_exp[,1])
gene_symbol <- gsub(tmp_exp[,1], pattern = '-', replacement = '_')
tmp_exp <- apply(tmp_exp[,c(-1)], 2, as.numeric)
rownames(tmp_exp) <- gene_symbol
#colnames(tmp_exp) <- substr(colnames(tmp_exp), 1, 12)
TCGA_GBMLGG_exp <- tmp_exp


TCGA_GBMLGG_trait <- read.table('/export2/zhuyw/database/TCGA/GBMLGG/TCGA_GBMLGG_traits.tsv', 
                             header = T, row.names = 1, sep = '\t')
tmp_traits <- TCGA_GBMLGG_trait
colnames(tmp_traits) <- c("Histology", "Grade", "Age", "Gender",
                         "OS.time", "OS", "IDH.status", "X1p.19q.codel", "IDH.codel.subtype")
tmp_traits <- as.data.frame(tmp_traits)
tmp_traits$OS.time <- 30*abs(tmp_traits$OS.time)
tmp_traits$Age.l <- ifelse(tmp_traits$Age>60, 'Age>60', 'Age<=60')
rownames(tmp_traits) <- gsub(rownames(tmp_traits), pattern = '-', replacement = '.')
TCGA_GBMLGG_trait <- tmp_traits


#
##    DataSet03-CGGA325 & CGGA693 (CGGA胶质瘤)
#   CGGA325
CGGA325_exp <- read.table('/export2/zhuyw/database/TCGA/GBM/CGGA.mRNAseq_325.RSEM-genes.20200506.txt', 
                          header = T, row.names = 1, sep = '\t')
tmp_exp <- CGGA325_exp
rownames(tmp_exp) <- gsub(rownames(tmp_exp), pattern = '-', replacement = '_')
CGGA325_exp <- tmp_exp

CGGA325_trait <- read.table('/export2/zhuyw/database/TCGA/GBM/CGGA.mRNAseq_325_clinical.20200506.txt', 
                            header = T, row.names = 1, sep = '\t')
tmp_traits <- CGGA325_trait
colnames(tmp_traits) <- c('PRS.type', 'Histology', 'Grade', 'Gender', 'Age', 'OS.time', 'OS',
                         'RadioT', 'ChemoT', 'IDH.mutation', 'X1p19q.status', 'MGMTp.status')
tmp_traits$Age.l <- ifelse(tmp_traits$Age>60, 'Age>60', 'Age<=60')
CGGA325_trait <- tmp_traits


#   CGGA693
CGGA693_exp <- read.table('/export2/zhuyw/database/TCGA/GBM/CGGA.mRNAseq_693.RSEM-genes.20200506.txt', 
                          header = T, row.names = 1, sep = '\t')
tmp_exp <- CGGA693_exp
rownames(tmp_exp) <- gsub(rownames(tmp_exp), pattern = '-', replacement = '_')
CGGA693_exp <- tmp_exp

CGGA693_trait <- read.table('/export2/zhuyw/database/TCGA/GBM/CGGA.mRNAseq_693_clinical.20200506.txt', 
                            header = T, row.names = 1, sep = '\t')
tmp_traits <- CGGA693_trait
colnames(tmp_traits) <- c('PRS.type', 'Histology', 'Grade', 'Gender', 'Age', 'OS.time', 'OS',
                         'RadioT', 'ChemoT', 'IDH.mutation', 'X1p19q.status', 'MGMTp.status')
tmp_traits$Age.l <- ifelse(tmp_traits$Age>50, 'Age>60', 'Age<=60')
CGGA693_trait <- tmp_traits


##
###   DataSet04-GEO数据集





##
###     数据保存
write.table(data.frame(Gene=rownames(DEGs.VS1), DEGs.VS1),
            "./output/tableS1.1_CGGA325_DEGs_VS1.tsv",
            row.names = F, sep = "\t", quote = F)
write.table(data.frame(Gene=rownames(DEGs.VS2), DEGs.VS2),
            "./output/tableS1.2_CGGA325_DEGs_VS2.tsv",
            row.names = F, sep = "\t", quote = F)
write.table(data.frame(Gene=rownames(DEGs.VS3), DEGs.VS3),
            "./output/tableS1.3_CGGA325_DEGs_VS3.tsv",
            row.names = F, sep = "\t", quote = F)
write.table(data.frame(Gene=rownames(DEGs.VS4), DEGs.VS4),
            "./output/tableS1.4_CGGA325_DEGs_VS4.tsv",
            row.names = F, sep = "\t", quote = F)
write.table(data.frame(Gene=rownames(DEGs.VS5), DEGs.VS5[, 2:ncol(DEGs.VS5)]),
            "./output/tableS1.6_CGGA325_DEGs_VS5.tsv",
            row.names = F, sep = "\t", quote = F)
write.table(data.frame(Gene=rownames(DEGs.Union), DEGs.Union[, 2:ncol(DEGs.Union)]),
            "./output/tableS1.5_CGGA325_DEGs_Union.tsv",
            row.names = F, sep = "\t", quote = F)

write.table(ego,
            "./output/tableS1.7_CGGA325_TRF_DEG.DOWN_GO.tsv",
            row.names = F, sep = "\t", quote = F)
write.table(ego,
            "./output/tableS1.8_CGGA325_TRF_DEG.UP_GO.tsv",
            row.names = F, sep = "\t", quote = F)

write.table(uni_cox_df,
            "./output/tableS2.1_CGGA325_TRF_UniCox.tsv",
            row.names = F, sep = "\t", quote = F)
write.table(uni_cox_df.filter2,
            "./output/tableS2.1_CGGA693_TRF_UniCox.tsv",
            row.names = F, sep = "\t", quote = F)
write.table(uni_cox_df.filter3,
            "./output/tableS2.1_TCGA_LGG_TRF_UniCox.tsv",
            row.names = F, sep = "\t", quote = F)

write.table(data.frame(SID=rownames(x), x),
            "./output/tableS3.1_TCGA_LGG_RS_model.tsv",
            row.names = F, sep = "\t", quote = F)
write.table(data.frame(SID=rownames(x), x),
            "./output/tableS3.2_CGGA325_LGG_RS_model.tsv",
            row.names = F, sep = "\t", quote = F)
write.table(data.frame(SID=rownames(x), x),
            "./output/tableS3.3_CGGA693_LGG_RS_model.tsv",
            row.names = F, sep = "\t", quote = F)

write.table(data.frame(top.table),
            "./output/tableS4.1_RS_Low_vs_High_DEGs.tsv",
            row.names = F, sep = "\t", quote = F)
write.table(data.frame(KEGG@result),
            "./output/tableS4.2_RS_Low_vs_High_DEGs_ssGSEA.tsv",
            row.names = F, sep = "\t", quote = F)

write.table(data.frame(top.table),
            "./output/tableS5.1_RS_Low_vs_High_ssGSEA_DEPs.tsv",
            row.names = F, sep = "\t", quote = F)
write.table(data.frame(Pathway=rownames(GSEA_matrix), GSEA_matrix),
            "./output/tableS5.2_RS_Low_vs_High_ssGSEA_hallmark.ES.tsv",
            row.names = F, sep = "\t", quote = F)
write.table(data.frame(top.table),
            "./output/tableS5.3_RS_Low_vs_High_ssGSEA_hallmark.DEPs.tsv",
            row.names = F, sep = "\t", quote = F)

write.table(data.frame(uni_cox_df),
            "./output/tableS6.1_Nomogram_UniCox.tsv",
            row.names = F, sep = "\t", quote = F)
write.table(data.frame(SID=rownames(tmpdata), tmpdata),
            "./output/tableS6.2_Nomogram_model.tsv",
            row.names = F, sep = "\t", quote = F)



##   临时读取
DEGs.Union <- read.table('./output/tableS1.5_CGGA325_DEGs_Union.tsv',
                       header = T, row.names = 1, sep = '\t')
rownames(DEGs.VS5) <- DEGs.VS5$Gene



###
######      1-LGG中复发因子的鉴定      ######
#
##   1.1-差异表达分析
#  DEGs_VS1: CGGA325_rLGG__VS__GTEX
#  DEGs_VS2: CGGA325_rLGG__VS__CGGA325_pLGG
#  DEGs_VS3: CGGA325_sGBM__VS__GTEX
#  DEGs_VS4: CGGA325_sGBM__VS__CGGA325_pLGG
#  DEGs_VS5: sGBM__VS__rLGG

gene_list <- intersect(rownames(CGGA325_exp), rownames(GTEX_brain_exp))
tmp_exp <- data.frame(log2(CGGA325_exp[gene_list, ]+1),
                      GTEX_brain_exp[gene_list, ])
tmp_traits <- data.frame(Sample.type=c(CGGA325_trait$Histology,
                                      rep('GTEX', ncol(GTEX_brain_exp))))
rownames(tmp_traits) <- colnames(tmp_exp)
tmp_traits$Tumor.type <- ifelse(tmp_traits$Sample.type %in% c('A', 'AA', 'AO', 'O'), 'pLGG',
                               ifelse(tmp_traits$Sample.type %in% c('rA', 'rO', 'rAA'), 'rLGG',
                                      ifelse(tmp_traits$Sample.type %in% c('sGBM'), 'sGBM',
                                             ifelse(tmp_traits$Sample.type %in% c('rGBM'), 'rGBM',
                                                    ifelse(tmp_traits$Sample.type %in% c('GBM'), 'pGBM',
                                                           ifelse(tmp_traits$Sample.type %in% c('GTEX'), 'GTEX', 'UN'))))))

##  DEGs_VS1
sample_list <- rownames(tmp_traits[tmp_traits$Tumor.type %in% c('GTEX', 'rLGG'), ])

##  DEGs_VS2
sample_list <- rownames(tmp_traits[tmp_traits$Tumor.type %in% c('pLGG', 'rLGG'), ])

##  DEGs_VS3
sample_list <- rownames(tmp_traits[tmp_traits$Tumor.type %in% c('GTEX', 'sGBM'), ])

##  DEGs_VS4
sample_list <- rownames(tmp_traits[tmp_traits$Tumor.type %in% c('pLGG', 'sGBM'), ])


counts <- floor(2^tmp_exp[, sample_list]-1)
group <- tmp_traits[sample_list, 'Tumor.type']

top.table$Gene <- rownames(top.table)
top.table <- top.table[,c("Gene", names(top.table)[1:6])]
top.table$change = as.factor(ifelse(abs(top.table$logFC) > 0.6, ifelse(top.table$logFC > 0.6, 'UP','DOWN'),'NOT'))
top.table$symbol <- gsub(top.table$Gene, pattern = '_', replacement = '-')

top.table.filter <- top.table[top.table$adj.P.Val < 0.05 &
                                abs(top.table$logFC) > 1, ]
table(top.table.filter$change)

DEGs.VS1 <- top.table.filter
DEGs.VS2 <- top.table.filter
DEGs.VS3 <- top.table.filter
DEGs.VS4 <- top.table.filter

#   取交集
Core_TEF <- Reduce(intersect,  list(v1 = paste(DEGs.VS1$Gene, DEGs.VS1$change, sep = '|'),
                        v2 = paste(DEGs.VS2$Gene, DEGs.VS2$change, sep = '|'),
                        v3 = paste(DEGs.VS3$Gene, DEGs.VS3$change, sep = '|'),
                        V4 = paste(DEGs.VS4$Gene, DEGs.VS4$change, sep = '|')))
Core_TEF <- data.frame(Genes=gsub('(.*)\\|(.*)', '\\1', Core_TEF),
                       change=gsub('(.*)\\|(.*)', '\\2', Core_TEF))


#   取并集
gene_list <- Reduce(union,  list(v1 = DEGs.VS1$Gene,
                                 v2 = DEGs.VS2$Gene,
                                 v3 = DEGs.VS3$Gene,
                                 V4 = DEGs.VS4$Gene))

DEGs.Union <- data.frame(Genes=gene_list)
rownames(DEGs.Union) <- DEGs.Union$Genes
DEGs.Union$VS1 <- DEGs.VS1[rownames(DEGs.Union), 'change']
DEGs.Union$VS2 <- DEGs.VS2[rownames(DEGs.Union), 'change']
DEGs.Union$VS3 <- DEGs.VS3[rownames(DEGs.Union), 'change']
DEGs.Union$VS4 <- DEGs.VS4[rownames(DEGs.Union), 'change']

tmp <- DEGs.Union[,2:5]
for(i in 1:ncol(tmp))
{
  tmp[ ,i] <- ifelse(is.na(tmp[,i]), 0,
                     ifelse(tmp[,i] %in% "UP", 1, -1))
}
colnames(tmp) <- paste(colnames(tmp), 'l', sep = '.')
tmp$change.sum <- apply(tmp, 1, sum)

DEGs.Union <- data.frame(DEGs.Union, tmp)


#
##   GO  富集分析
gene_list <- rownames(DEGs.Union[DEGs.Union$change.sum %in% c(3,4), ])
gene_list <- gsub(gene_list, pattern = "_", replacement = "-")

#
##   1.2-单因素分析（UniCox）
tmp_traits <- CGGA325_trait
tmp_exp <- CGGA325_exp

tmp_traits <- tmp_traits[tmp_traits$Histology %in% c("A", "AA", "AO", "O", "AOA", "OA") &
                           !is.na(tmp_traits$OS) &
                           !is.na(tmp_traits$OS.time) &
                           tmp_traits$OS.time>0, ]

tmp_traits <- CGGA693_trait
tmp_exp <- CGGA693_exp

tmp_traits <- TCGA_GBMLGG_traits
tmp_exp <- TCGA_GBMLGG_exp

tmp_traits <- tmp_traits[tmp_traits$Histology %in% c("astrocytoma", "oligoastrocytoma", "oligodendroglioma") & 
                           substr(rownames(tmp_traits),14,15)<10 &
                           !is.na(tmp_traits$OS) &
                           !is.na(tmp_traits$OS.time) &
                           tmp_traits$OS.time>0, ]

sample_list <- intersect(rownames(tmp_traits),
                         colnames(tmp_exp))
gene_list <- rownames(DEGs.Union[DEGs.Union$change.sum %in% c(-4,-3,3,4), ])
gene_list <- intersect(gene_list, rownames(tmp_exp))

tmp <- data.frame(tmp_traits[sample_list, c('OS', 'OS.time')],
                  t(tmp_exp[gene_list, sample_list]))

uni_cox_df <- uni_cox_in_bulk(gene_list = gene_list, 
                              survival_info_df = as.data.frame(tmp))
rownames(uni_cox_df) <- uni_cox_df$uni_cox_sig_gene


##   UniCox的结果也需要利用两个外部数据集进行验证
#  CGGA325(pLGG)
uni_cox_df.filter1 <- uni_cox_df[uni_cox_df$Likelihood_pvalue<0.05 &
                                  abs(uni_cox_df$bate)>0.1, ]
#  CGGA693(pLGG)
uni_cox_df.filter2 <- uni_cox_df[uni_cox_df$Likelihood_pvalue<0.05 &
                                   abs(uni_cox_df$bate)>0.1, ]
#  TCGA_LGG(pLGG)
uni_cox_df.filter3 <- uni_cox_df[uni_cox_df$Likelihood_pvalue<0.05 &
                                   abs(uni_cox_df$bate)>0.1, ]

gene_list <- intersect(uni_cox_df.filter1$uni_cox_sig_gene,
                       uni_cox_df.filter2$uni_cox_sig_gene)
gene_list <- intersect(gene_list, uni_cox_df.filter3$uni_cox_sig_gene)


#   以 TCGA_LGG 为训练集的情况
#   after LASSO
gene_list <-  c("IGF2BP3", "HOXA1", "KIF18A", "CALN1", "RAB42", "KIF4A", "FAM133A", "HGF", "MN1")
#  after stepwise
gene_list <- c("HOXA1", "KIF18A", "FAM133A", "HGF", "MN1")
#   TCGA_LGG  AUC: 0.756, 0.877, 0.88
#   CGGA_325  AUC: 0.799, 0.839, 0.844
#   CGGA_693  AUC: 0.77, 0.826, 0.858
#Call:
#  coxph(formula = mutil_formulas, data = as.data.frame(tmp))

#n= 465, number of events= 113 

#coef exp(coef) se(coef)      z Pr(>|z|)    
#HOXA1    0.17207   1.18776  0.06842  2.515 0.011902 *  
#  KIF18A   0.29861   1.34799  0.07197  4.149 3.33e-05 ***
#  FAM133A -0.13995   0.86940  0.06405 -2.185 0.028894 *  
#  HGF      0.16022   1.17377  0.06224  2.574 0.010047 *  
#  MN1     -0.27572   0.75902  0.08267 -3.335 0.000852 ***
#  ---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#exp(coef) exp(-coef) lower .95 upper .95
#HOXA1      1.1878     0.8419    1.0387    1.3582
#KIF18A     1.3480     0.7418    1.1707    1.5522
#FAM133A    0.8694     1.1502    0.7668    0.9857
#HGF        1.1738     0.8520    1.0390    1.3261
#MN1        0.7590     1.3175    0.6455    0.8925

#Concordance= 0.82  (se = 0.024 )
#Likelihood ratio test= 124.8  on 5 df,   p=<2e-16
#Wald test            = 130.7  on 5 df,   p=<2e-16
#Score (logrank) test = 169.2  on 5 df,   p=<2e-16



#   以 CGGA_325 为训练集的情况
#   after LASSO
gene_list <-  c("OIP5", "KSR2", "HMGCLL1", "GNAL", "WSCD2", "FAM133A", "KCNJ11", "AMZ1", "FAM72D", "CACNG2")
#  after stepwise
gene_list <- c("KSR2", "HMGCLL1", "AMZ1", "FAM72D")
#   TCGA_LGG  AUC: 0.696032628151244, 0.744486394110743, 0.772689545815242
#   CGGA_325  AUC: 0.844330797144937, 0.903790513791139, 0.908965473747312
#   CGGA_693  AUC: 0.726464063858832, 0.766214449485494, 0.768411895220038



##
######      2-RS.model与多组学数据的关系      ######
tmp_traits <- TCGA_GBMLGG_traits
tmp_exp <- TCGA_GBMLGG_exp

tmp_traits <- tmp_traits[tmp_traits$Histology %in% c("astrocytoma", "oligoastrocytoma", "oligodendroglioma") & 
                           substr(rownames(tmp_traits),14,15)<10 &
                           !is.na(tmp_traits$OS) &
                           !is.na(tmp_traits$OS.time) &
                           tmp_traits$OS.time>0, ]

sample_list <- intersect(rownames(tmp_traits),
                         colnames(tmp_exp))

#
##   2.1  PCA分析，验证分组的一致性
tmp <- t(tmp_exp[gene_list, sample_list])
res.pca <- PCA(tmp, graph = FALSE)
fviz_pca_ind(res.pca,
             geom.ind = "point", # show points only (nbut not "text")
             col.ind = as.factor(tmp_traits[sample_list, "TRF.RS.l"]), # color by groups
             palette = bioCol[4:2], #c("#00AFBB", "#E7B800", "#FC4E07", ""),
             addEllipses = TRUE, # Concentration ellipses
             legend.title = "RiskScore"
)


#
###   2.2   ICI(免疫细胞浸润)特征
#   cibersort
ICI_cibersort <- read.table("/export2/zhuyw/database/TCGA/TCGA_ICI/CIBERSORT/TCGA.GBMLGG.cibersort.relative.tsv", 
                            header = T, row.names = 1, sep = "\t")
ICI_cibersort$SID <- substr(rownames(ICI_cibersort),1,15)

tmpdata <- ICI_cibersort[ICI_cibersort$SID %in% sample_list, 2:23]
tmpdata$SID <- substr(rownames(tmpdata),1,15)
tmpdata <- avereps(tmpdata[, 1:22], ID = tmpdata[ ,23])

#   timer
ICI_timer <- read.table("/export2/zhuyw/database/TCGA/TCGA_ICI/TIMER/TCGA.all.timer.scores.tsv", 
                            header = T, sep = "\t")
ICI_timer$SID <- substr(ICI_timer$X,1,15)

tmpdata <- ICI_timer[ICI_timer$SID %in% sample_list, ]
rownames(tmpdata) <- tmpdata$SID
tmpdata <- avereps(tmpdata[, 2:7], ID = tmpdata[ ,8])

#   xcell
ICI_xcell <- read.table("/export2/zhuyw/database/TCGA/TCGA_ICI/xCell/TCGA.all.ICI.xCell.tsv", 
                            header = T, row.names = 1, sep = "\t")
tmpdata <- t(ICI_xcell[, colnames(ICI_xcell) %in% sample_list])


tmpdata <- as.data.frame(tmpdata)
tmpdata$Type <- tmp_traits[rownames(tmpdata), "TRF.RS.l"]
tmpdata <- melt(tmpdata, id.vars=c("Type"))
colnames(tmpdata)=c("Type", "Immune", "Fraction")
tmpdata$Type <- as.factor(tmpdata$Type)
tmpdata$Fraction <- as.numeric(tmpdata$Fraction)
tmpdata$Fraction <- scale(tmpdata[,3], center = F)
#data$Fraction[data$Fraction>12]=12
ggboxplot(tmpdata, x="Immune", y="Fraction", fill = "Type", palette = bioCol[4:3], #fill="ICIcluster",
          ylab="Scale of Fraction",
          xlab="",
          legend.title="RiskScore") + 
  #palette=bioCol)
  rotate_x_text(50) + #scale_fill_brewer(palette = bioCol[4:3]) +#"Set3") +
  stat_compare_means(aes(group=Type),
                     symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), 
                                      symbols = c("***", "**", "*", "ns")),label = "p.signif")


##
######      3-ssGSEA分析      ######
#
##   3.1-KEGG(186 GeneSets) ssGSEA
#   DEGs 分析 RS_Low_vs_High_DEGs
counts <- floor(2^tmp_exp[, sample_list]-1)
group <- factor(tmp_traits[sample_list, "TRF.RS.l"])

genelist <- top.table$logFC[top.table$change %in% c("UP", "DOWN")]
names(genelist) <- top.table$Gene[top.table$change %in% c("UP", "DOWN")]
genelist <- sort(genelist, decreasing = T)

gmtFile <- read.gmt("/export2/zhuyw/database/GSEA/c2.cp.kegg.v7.2.symbols.gmt")

KEGG <- GSEA(genelist, TERM2GENE = gmtFile, pvalueCutoff = 0.05)  #  GSEA分析
gseaplot2(KEGG, 1:6,
          pvalue_table = T)
gseaplot2(KEGG, 7:13, 
          pvalue_table = T)


#
##   3.2-TCGA ssGSEA (1387 GeneSets)  RS_Low_vs_High
GSEA_matrix <- read.table("/export2/zhuyw/database/TCGA/LGG/TCGA.LGG_ssGSEA_1387GeneSets_ZScore.tsv",
                          header = T, row.names = 1, sep = "\t")
tmp <- GSEA_matrix[, intersect(colnames(GSEA_matrix), sample_list)]

group <- factor(tmp_traits[colnames(tmp), "TRF.RS.l"])

mm <- model.matrix(~0 + group)
fit <- lmFit(tmp[, colnames(tmp)], mm)
head(coef(fit))
contr <- makeContrasts(groupLow - groupHigh, levels = colnames(coef(fit)))
tmp <- contrasts.fit(fit, contr)
tmp <- eBayes(tmp)

top.table <- topTable(tmp, sort.by = "P", number=Inf)
head(top.table, 10)

top.table$Pathway <- rownames(top.table)
top.table <- top.table[,c("Pathway", names(top.table)[1:6])]
top.table$change = as.factor(ifelse(abs(top.table$logFC) > 0.15, ifelse(top.table$logFC > 0.15, 'UP','DOWN'),'NOT'))
top.table$change[top.table$P.Value>0.05] <- "NOT"
table(top.table$change)

gene_list <- top.table[top.table$change %in% c("UP", "DOWN"), "Pathway"]

#   热图可视化
tmpdata <- tmp_traits[intersect(colnames(GSEA_matrix), sample_list),
                      c("OS", "Histology", "Grade", "Gender", "IDH.status", "X1p.19q.codel", "IDH.codel.subtype", "Age.l", "TRF.RS.l")]
tmpdata <- tmpdata[order(tmpdata$TRF.RS.l), ]

x <- list()
for(i in 1:length(colnames(tmpdata))){
  x[[i]] <- as.factor(tmpdata[, i])
}
anno_col <- as.data.frame(x)
colnames(anno_col) <- colnames(tmpdata)
rownames(anno_col) <- rownames(tmpdata)

x <- list()
for(i in 1:length(colnames(tmpdata))){
  n = length(unique(sort(tmpdata[, i])))+2
  x[[i]] <- cols[3:n]
  names(x[[i]]) <- unique(sort(tmpdata[, i]))
}
names(x) <- colnames(tmpdata)
anno_colors <- x

pheatmap(scale(GSEA_matrix[gene_list, rownames(tmpdata)], center = F),
         cluster_rows = T,
         cluster_cols = FALSE,
         color=colorRampPalette(colors = c("purple", "white", "green"))(10),
         #display_numbers=data_mark, 
         fontsize_number=3.6,
         #gaps_row = c(2,13),
         #legend_breaks = c(-4,4),
         annotation_col = anno_col,
         annotation_colors = anno_colors,
         show_colnames = F)


#
##   3.3-HallMark (50 GeneSets) ssGAEA
#  读取基因集
gmtFile <- read.gmt("/export2/zhuyw/database/GSEA/h.all.v7.5.symbols.gmt")  #读gmt文件
#geneset <- GSEABase::getGmt(gmtFile) 
list <- split(as.matrix(gmtFile)[,2], gmtFile[,1])

#  计算富集评分es
tmp <- floor(2**tmp_exp[, sample_list]-1)
rownames(tmp) <- gsub(rownames(tmp), pattern = '_', replacement = '-')
GSEA_matrix <- GSVA::gsva(as.matrix(tmp), list,
                 min.sz=10, max.sz=500, verbose=TRUE)




###
######      4-TMB.mut and TMB.met characteristic between high and low RS group      ######
#
##   4.1-TMB.mut
TMB.mut <- read.table("/export2/zhuyw/database/TCGA/LGG/TCGA.LGG.varscan.DR-10.0.somatic.maf.TMB",
                      header = T, row.names = 1, sep = "\t")
TMB.mut$SID <- gsub(substr(rownames(TMB.mut),1,15), pattern = "-", replacement = ".")
rownames(TMB.mut) <- TMB.mut$SID

tmp <- tmp_traits[intersect(sample_list, rownames(TMB.mut)), ]
tmp$TMB.mut <- TMB.mut[rownames(tmp), "TMB.met"]

res.cut <- surv_cutpoint(tmp, time = "OS.time", event = "OS",
                         variables = "TMB.mut")
summary(res.cut)
plot(res.cut)

tmp$TMB.mut.l <- ifelse(tmp$TMB.mut>0.36, "High", "Low")

fit_km <- survfit(Surv(OS.time, OS) ~ TMB.mut.l , data = tmp)
ggsurvplot(fit_km, pval=TRUE,
           conf.int =TRUE,
           xlab ="Time of days",
           ggtheme =theme_light(),
           risk.table = T, 
           palette = bioCol[c(4:3)],
           legend.labs=c('High', 'Low'),
           legend.title='TMB.mut')


#   线性回归
x <- tmp[, c("TRF.RS", "TRF.RS.l", "TMB.mut", "TMB.mut.l")]
x <- na.omit(x)
x$TRF.RS.l <- as.factor(x$TRF.RS.l)
x$TMB.mut <- log(x$TMB.mut)
x$TMB.mut[x$TMB.mut>4] <- 4

ggplot(x, aes(TMB.mut, TRF.RS)) + 
  xlab("Tumor Burden Mutation")+ylab("Risk Score")+
  geom_point(aes(colour=TRF.RS.l))+
  scale_color_manual(values=bioCol[4:3])+ 
  geom_smooth(method="lm",formula = y ~ x) + theme_bw()+
  stat_cor(method = 'spearman', aes(x =TMB.mut, y =TRF.RS))


#   比例boxplot
ggplot(x, aes(TMB.mut.l, fill=factor(TRF.RS.l))) +
  geom_bar(position='fill') +
  scale_fill_manual(values = bioCol[4:3]) #c("deeppink3", "skyblue3")) #bioCol[5:6])

table(tmp[tmp$TMB.mut.l=='High', 'TRF.RS.l'])/length(tmp[tmp$TMB.mut.l=='High', 'TRF.RS.l'])
#   High       Low 
#   0.48       0.52 
table(tmp[tmp$TMB.mut.l=='Low', 'TRF.RS.l'])/length(tmp[tmp$TMB.mut.l=='Low', 'TRF.RS.l'])
#   High       Low 
#   0.050847 0.949153



#   小提琴图
group=levels(factor(x$TMB.mut.l))
comp=combn(group, 2)
my_comparisons=list()
for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}

ggviolin(x, x="TMB.mut.l", y="TRF.RS", fill = "TMB.mut.l",
         ylab=c("Risk Score"),
         xlab="",
         legend.title="TMB level",
         palette=bioCol[4:3], #bioCol[1:3], 
         add="boxplot", add.params = list(fill="white")) + 
  stat_compare_means(comparisons = my_comparisons)  ##, method = "t.test")

#
###   4.2-TMB.met
TMB.mut <- read.table("/export2/zhuyw/database/TCGA/LGG/TCGA.LGG.HumanMethylation450.TMBmet",
                      header = T, row.names = 1, sep = "\t")
#  count  TMB.met_High
#High        Low 
#0.06865672 0.93134328
#   count  TMB.met_Low
#High       Low 
#0.5230769 0.4769231 



##
###   4.3-瀑布图
library('maftools')
maf=read.maf(maf="./RS_high_new.maf")
maf=read.maf(maf="./RS_low_new.maf")#, clinicalData="ann.txt")

oncoplot(maf=maf, 
         #clinicalFeatures=colnames(data), 
         top = 30,
         #genes=gene_list[1:30, 2], 
         #annotationColor=ann_colors, 
         #keepGeneOrder=T
)



###
######      5-Nomogram with clinical traits      ######
gene_list <- c("Histology", "Grade", "Age", "Gender", "IDH.status", "X1p.19q.codel",
               "IDH.codel.subtype", "Age.l", "TRF.RS", "TRF.RS.l", "TMB.mut", "TMB.mut.l", "TMB.met", "TMB.met.l")
tmp <- tmp_traits[, c("OS", "OS.time", gene_list)]
tmp <- na.omit(tmp)

#  UniCox
gene_list <- c("Histology", "Grade", "Gender", "IDH.status", "X1p.19q.codel",
               "IDH.codel.subtype", "Age.l", "TRF.RS.l", "TMB.mut.l", "TMB.met.l")

#  MutilCox
gene_list <- c("Histology", "Grade", "IDH.status", "X1p.19q.codel",
               "Age.l", "TRF.RS.l")

#  Nomogram
gene_list <- c("IDH.status", "Age.l", "TRF.RS.l")





###
######      6-Depmap 5-signature genes      ######
tmp <- read.table("./data/Depmap_5gene.tmp",
                  header = T, row.names = 1, sep = "\t")

#   箱线图
ggboxplot(tmp, x="Primary.Disease", y="CRISPR..DepMap.22Q1.Public.Score..Chronos.", #fill = "Type", #color="ICIcluster", fill="ICIcluster",
          ylab="Gene Effect (chronos)", xlab="", title = "MN1") +
  geom_hline(yintercept=0, linetype=3.5) +
  coord_flip()

#   散点图
x <- tmp[tmp$Primary.Disease %in% c("Brain Cancer"), c(1,6)]
x <- na.omit(x)
colnames(x) <- c("Geff", "Gexp")

ggplot(x, aes(Geff, Gexp)) + xlab("Gene Effect (chronos)")+ylab("Gene Expression Scale")+
  geom_point() + ggtitle("MN1") + 
  #scale_color_manual(values=bioCol[4:3])+ 
  geom_smooth(method="lm",formula = y ~ x) + theme_bw()+
  stat_cor(method = 'spearman', aes(x =Geff, y =Gexp))





###
######      7-model 5-signature genes prediction in immunothreapy      ######
IMvigor210_exp <- read.table('/export2/zhuyw/database/IMvigor210/IMvigor210_counts_symbol_uni.tsv', 
                             header = T, row.names = 1,sep = '\t')
IMvigor210_traits <- read.table('/export2/zhuyw/database/IMvigor210/IMvigor210_traits.tsv',
                                header = T, row.names = 1,sep = '\t')

#  after stepwise
gene_list <- c("HOXA1", "KIF18A", "FAM133A", "HGF", "MN1")
gene_list <- intersect(gene_list, rownames(IMvigor210_exp))
gene_list <- c("KIF18A")

tmp <- IMvigor210_traits[!is.na(IMvigor210_traits$OS) &
                         IMvigor210_traits$OS.time>0, ]
#tmp$OS.time <- floor(tmp$OS.time*7)

tmp <- data.frame(tmp, log2(t(IMvigor210_exp[gene_list, rownames(tmp)])+1))

#  比例
#  High_CR/PR : 0.1658,  High_SD/PD: 0.8342
#  Low_CR/PR : 0.3333,  Low_SD/PD: 0.6667


#   violin

tmpdata <- na.omit(tmp[,c("RS", "Response")])
colnames(tmpdata) <- c('RS', 'Group')

group=levels(factor(tmpdata$Group))
#group=levels(factor(data$anatomic.level))
comp=combn(group, 2)
my_comparisons=list()
for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}

ggviolin(tmpdata, x="Group", y="RS", fill = "Group",
         ylab=paste0("Risk Score"),
         xlab="",
         legend.title="DMS",
         palette=bioCol[6:2], 
         add="boxplot", add.params = list(fill="white"))+ 
  stat_compare_means(comparisons = my_comparisons,
                     symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), 
                                      symbols = c("***", "**", "*", "ns")),label = "p.signif")





###
######      调试      ######
#
##   1.1-按照北哥的要求进行调试
DEGs.VS5$VS1.change <- DEGs.VS1[rownames(DEGs.VS5), 'change']
DEGs.VS5$VS3.change <- DEGs.VS3[rownames(DEGs.VS5), 'change']

# UniCox
gene_list <- DEGs.VS5[!is.na(DEGs.VS5$VS1.change) & 
                        !is.na(DEGs.VS5$VS3.change), 'Gene']

tmp <- data.frame(tmp_traits[sample_list, c('OS', 'OS.time')],
                  t(tmp_exp[gene_list, sample_list]))

uni_cox_df <- uni_cox_in_bulk(gene_list = gene_list, 
                              survival_info_df = as.data.frame(tmp))
rownames(uni_cox_df) <- uni_cox_df$uni_cox_sig_gene
uni_cox_df.filter <- uni_cox_df[uni_cox_df$Likelihood_pvalue<0.05 &
                                  abs(uni_cox_df$bate)>0.3, ]

DEGs.VS5$UniCox.HR <- uni_cox_df.filter[rownames(DEGs.VS5), 'Hazard_Ratio']
DEGs.VS5$UniCox.pvalue <- uni_cox_df.filter[rownames(DEGs.VS5), 'Likelihood_pvalue']
























