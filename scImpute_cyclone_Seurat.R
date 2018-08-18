#use scImpute to impute the dropouts before downstream analysis.Assumed that all R packages needed are installed successfully.
library(scImpute)
library(Seurat)
setwd("E:/scRNA-Seq/HP1_outs/scImpute_k5")
HP1.data<-Read10X(data.dir="E:/scRNA-Seq/HP1_outs/outs/filtered_gene_bc_matrices/hg19")
write.table(as.matrix(HP1.data),"E:/scRNA-Seq/HP1_outs/HP1_data.txt",sep="\t")
scimpute("E:/scRNA-Seq/HP1_outs/HP1_data.txt",infile = "txt", outfile = "txt", out_dir="E:/scRNA-Seq/HP1_outs/scImpute_k5/",labeled = FALSE,drop_thre =0.5, Kcluster =5,labels = NULL,ncores = 1)

#cyclone to evaluate the stage of cells 
library(scran)
library(org.Hs.eg.db)
hs.pairs <- readRDS(system.file("exdata", "human_cycle_markers.rds", package="scran"))
HSMM_expr_matrix <- read.table("E:/scRNA-Seq/HP1_outs/scImpute_k5/scImpute/scimpute_count.txt")
HSMM_expr_matrix_ok=HSMM_expr_matrix
dim(HSMM_expr_matrix_ok)
genes=read.table("E:/scRNA-Seq/genes.tsv",header=T)
HSMM_expr_matrix_ok$genename=rownames(HSMM_expr_matrix_ok)
HSMM_expr_matrix_ok=merge(genes,HSMM_expr_matrix_ok,by="genename")
dim(HSMM_expr_matrix_ok)
rownames(HSMM_expr_matrix_ok)=HSMM_expr_matrix_ok$ensembl
HSMM_expr_matrix_ok=subset.data.frame(HSMM_expr_matrix_ok,select = -c(ensembl,genename))
dim(HSMM_expr_matrix_ok)
assigned<- cyclone(as.matrix(HSMM_expr_matrix_ok), hs.pairs, gene.names=rownames(HSMM_expr_matrix_ok))
write.table(assigned,"assigned.txt",sep="\t")
assigned=read.delim("assigned.txt",sep="\t",header=T,row.names=1)
cells=colnames(HSMM_expr_matrix)
rownames(assigned)=cells
write.table(assigned,"cell_sample_sheet.txt",sep="\t",row.names=F)

col <- character(ncol(HSMM_expr_matrix_ok))
is.G1=grep("G1",assigned$phases)
is.G2M=grep("G2M",assigned$phases)
is.S=grep("S",assigned$phases)
col[is.G1] <- "red"
col[is.G2M] <- "blue"
col[is.S] <- "darkgreen"
plot(assigned$score$G1, assigned$score$G2M, col=col, pch=16)

#Data quality control with Seurat package
library(Seurat)
setwd("E:/scRNA-Seq/HP1_outs/scImpute_k5")
HSMM_expr_matrix <- read.table("E:/scRNA-Seq/HP1_outs/scImpute_k5/scImpute/scimpute_count.txt")
cell_sample_sheet <- read.table("E:/scRNA-Seq/HP1_outs/scImpute_k5/cell_sample_sheet.txt",header=T,row.names = 1,sep="\t")
HP1<-CreateSeuratObject(raw.data=HSMM_expr_matrix ,min.cells=10,min.genes=200,project="HP1",is.expr=0,meta.data=cell_sample_sheet)
mito.genes<-grep(pattern="^MT-",x=rownames(x=HP1@data),value=TRUE)
percent.mito<-Matrix::colSums(HP1@raw.data[mito.genes,])/Matrix::colSums(HP1@raw.data)
HP1<-AddMetaData(HP1,metadata=percent.mito,col.name="percent.mito")
Total_mRNAs <- Matrix::colSums(HP1@raw.data)
HP1<-AddMetaData(HP1,metadata=Total_mRNAs,col.name="Total_mRNAs")
upper_bound <- 10^(mean(log10(HP1@meta.data$Total_mRNAs)) +2*sd(log10(HP1@meta.data$Total_mRNAs)))
lower_bound <- 10^(mean(log10(HP1@meta.data$Total_mRNAs)) -2*sd(log10(HP1@meta.data$Total_mRNAs)))
HP1<-FilterCells(HP1,subset.names=c("percent.mito","Total_mRNAs"),low.thresholds=c(-Inf,lower_bound),high.thresholds=c(0.05,upper_bound))
HP1<-NormalizeData(HP1,normalization.method="LogNormalize",scale.factor=10000,display.progress=T)
HP1<-AddMetaData(HP1,metadata=HP1@data["HBB",],col.name="HBB_exprs")
HP1<-AddMetaData(HP1,metadata=HP1@data["HBG2",],col.name="HBG2_exprs")
HP1<-AddMetaData(HP1,metadata=HP1@data["HBG1",],col.name="HBG1_exprs")
HP1<-AddMetaData(HP1,metadata=HP1@data["HBA2",],col.name="HBA2_exprs")
HP1<-AddMetaData(HP1,metadata=HP1@data["HBA1",],col.name="HBA1_exprs")
HP1<-AddMetaData(HP1,metadata=HP1@data["GYPA",],col.name="GYPA_exprs")
HP1<-ScaleData(HP1,vars.to.regress="phases")

write.table(HP1@meta.data,"HP1.meta.data.txt",sep="\t")
save(HP1,file="HP1.Rdata")

plot(HP1@meta.data$HBG2_exprs,HP1@meta.data$HBB_exprs,pch=20)
plot(HP1@meta.data$HBG1_exprs,HP1@meta.data$HBB_exprs,pch=20)
plot(HP1@meta.data$HBA2_exprs,HP1@meta.data$HBB_exprs,pch=20)
plot(HP1@meta.data$HBA1_exprs,HP1@meta.data$HBB_exprs,pch=20)
plot(HP1@meta.data$GYPA_exprs,HP1@meta.data$HBB_exprs,pch=20)

#N:globin-Negtive cells.P:globin-Positive cells
HP1_N_cells=rownames(HP1@meta.data[which(HP1@meta.data$HBB_exprs <= 5 & HP1@meta.data$HBG2_exprs <= 6 & HP1@meta.data$HBG1_exprs <= 6 & HP1@meta.data$HBA2_exprs <= 5  & HP1@meta.data$HBA1_exprs <= 5 & HP1@meta.data$GYPA_exprs ==0 ),])
HP1_P_cells=setdiff(HP1@cell.names,HP1_N_cells) 
HP1P=SubsetData(HP1,cells.use=HP1_P_cells,subset.raw=T)
HP1P<-NormalizeData(HP1P,normalization.method="LogNormalize",scale.factor=10000,display.progress=T)
HP1P<-ScaleData(HP1P,vars.to.regress="phases")
HP1P<-FindVariableGenes(HP1P, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)

write.table(HP1P@meta.data,"HP1P.meta.data.txt",sep="\t")
save(HP1P,file="HP1P.Rdata")

###prepare normalized expression for Cellrouter input
write.table(as.matrix(HP1P@data),"HP1P.normalized_expression.txt",sep="\t")
write.table(colnames(HP1P@data),"HP1P.cell_names.csv",sep="\t")
write.table(rownames(HP1P@data),"HP1P.gene_names.csv",sep="\t")