#run this script in linux system and oracle java
source('/mnt/md1200/2/scRNA_Seq/Cleandata/cellrouter-master/CellRouter_Class.R')
libdir <- '/mnt/md1200/2/scRNA_Seq/Cleandata/cellrouter-master/CellRouter/'
library(dplyr)
library(plotrix)
matrix=read.table("/mnt/md1200/2/scRNA_Seq/Cleandata/HP1/cellrangerRkit/HP1P_outs/outs/analysis/tsne/2_components/projection.csv",sep=",",header=T,row.names=1)
colnames(matrix) <- c('tSNE1','tSNE2')
rownames(matrix)=gsub("-1","",rownames(matrix))
ndata <- read.table('../HP1P.normalized_expression.txt',sep="\t",header=T,row.names=1)
genes <-as.vector(rownames(ndata))
map <- data.frame(id=rownames(ndata),symbol=genes,stringsAsFactors = FALSE)
ndata <- averageIds(ndata,map,'symbol')

#Remove genes with zero variance across all cells
var <- apply(ndata,1,var)
var <- var[which(var > 0)]
ndata <- ndata[names(var),]

### selecting genes to use as regulated along developmental trajectories.
#pca <- prcomp(t(ndata),scale=TRUE,center=TRUE)
#loadings <- pca$rotation
#num_pc <- 5
#quantile <- 0.975
#genes2use <- unique(as.vector(unlist(apply(loadings[,1:num_pc],2,function(x){names(x[which(abs(x) >= quantile(x,quantile))])}))))
genes2use=rownames(ndata)
ggrn <- buildGRN('Hs',ndata,genes2use,2,'results/GRN.R') #original 5

### Subpopulation identification and gene signatures with CellRouter
cellrouter <- CellRouter(expdata=ndata,annotations=colnames(ndata))
cellrouter@rdimension <- matrix
cellrouter <- findsubpopulations(cellrouter,90,'jaccard','results/kNN_network.gml')

df=read.table("HP1P.meta.data.txt",header=T,row.names=1,sep="\t")
df$sample_id=rownames(df)
df=merge(cellrouter@sampTab,df,by="sample_id",all=T)
write.table(df,"HP1P.meta.data.withSP.txt",sep="\t")

lengths(cellrouter@graph$subpopulation)
cellrouter <- diffexpr(cellrouter,column='population',pvalue = 0.05)
markers <- findmarkers(cellrouter)
write.table(markers,"results/HP1P.markers.txt",sep="\t")
plotReducedDimension(cellrouter,5,5,filename='results/HP1P.tSNE.pdf')
table(cellrouter@sampTab$population)
write.table(cellrouter@sampTab,"results/HP1P.cellrouter_sampTab.txt",sep="\t")

######## Trajectory Detection using CellRouter ###
cellrouter <- createKNN(cellrouter,90,'jaccard','results/paths/kNN_network_trajectory.gml') #10 before this 90
filename <- "results/paths/cell_edge_weighted_network.txt"
write.table(cellrouter@graph$edges,file=filename,sep='\t',row.names=FALSE,col.names = FALSE,quote=FALSE) #input network

##select starting subpopulation,all other subpopulations are targets
sources <- c('SP_10') #from SP_10 to SP_8
targets <- setdiff(as.vector(cellrouter@sampTab$population),sources)
methods <- c("euclidean","maximum","manhattan","canberra","binary",'graph') #graph for distances in KNN
cellrouter <- findpaths(cellrouter,libdir,paste(getwd(),'results/paths',sep='/'),method="graph")
ranks <- c('path_cost','path_flow','rank','length')
cellrouter <- processtrajectories(cellrouter,genes2use,path.rank=ranks[3],num.cells = 3,neighs = 1)
names <- unique(names(cellrouter@pathsinfo$distr))
clusters.show <- names
cellrouter <- correlationpseudotime(cellrouter,type='spearman')
cellrouter <- topgenes(cellrouter,0.85,0.15)
cellrouter <- smoothdynamics(cellrouter,names)
cellrouter <- clusterGenesPseudotime(cellrouter,10)
save(cellrouter,file='results/CellRouter_StemID_Processed.R')

##plot begins####
###positive and negative controls
p <- c('SP_10.SP_8') 
cellrouter@signatures$SP_1$subpopulation="SP_1"
cellrouter@signatures$SP_2$subpopulation="SP_2"
cellrouter@signatures$SP_3$subpopulation="SP_3"
cellrouter@signatures$SP_4$subpopulation="SP_4"
cellrouter@signatures$SP_5$subpopulation="SP_5"
cellrouter@signatures$SP_6$subpopulation="SP_6"
cellrouter@signatures$SP_7$subpopulation="SP_7"
cellrouter@signatures$SP_8$subpopulation="SP_8"
cellrouter@signatures$SP_9$subpopulation="SP_9"
cellrouter@signatures$SP_10$subpopulation="SP_10"
cellrouter@signatures$SP_11$subpopulation="SP_11"
cellrouter@signatures$SP_12$subpopulation="SP_12"
data=rbind(cellrouter@signatures$SP_1,cellrouter@signatures$SP_2,cellrouter@signatures$SP_3,cellrouter@signatures$SP_4,cellrouter@signatures$SP_5,cellrouter@signatures$SP_6,cellrouter@signatures$SP_7,cellrouter@signatures$SP_8,cellrouter@signatures$SP_9,cellrouter@signatures$SP_10,cellrouter@signatures$SP_11,cellrouter@signatures$SP_12)
data$gene_names=rownames(data)
write.table(data,"results/HP1P.dif_genes.txt",sep="\t")
write.table(as.matrix(unlist(cellrouter@top.correlations$up)),"results/SP_10.SP_8_up_top_correlations.genes.txt",sep="\t")
write.table(as.matrix(unlist(cellrouter@top.correlations$down)),"results/SP_10.SP_8_down_top_correlations.genes.txt",sep="\t")

genlist=c("ALAS2","BNIP3L","SEC62","CA1","HISTIH4C","PTMA","TMCC2","ARL4A","BPGM","HIST1H4C","HMGB2","GPX1","KRT1","EIF1AY","MALAT1","MGST3","FTL","HBD")
plotPathHeatmap2(cellrouter,p,genelist,TRUE,2,2,10,10,paste('results/',"heatmap_along_trajectory__",sep=''))

## GRN score for selected transitions
tfs <- find_tfs(species = 'Hs')
save(tfs,file="results/tfs.R")
x <- grnscores(cellrouter,tfs,p,direction='both',dir.targets='up',q.up=0.9,q.down=0.1,columns=1,width=15,height=10,flip=T,filename=paste('results/',p,sep=''))#barplot,considering this for more endpoints
plottr(cellrouter,p,x[[p]]$scores,cluster=TRUE,1,15,10,paste('results/',p,'_up_diff_dynamics.pdf',sep='')) #two panel heatmap

pdf("results/grndynamics.cor.SP_10.SP_8.pdf")
grndynamics(cellrouter, tfs,p, 100)
dev.off()
#scores <- x[[p]]$scores
#plottrajectories(cellrouter, p, names(scores), rescale = TRUE, columns=1, width=20, height=10, filename='results/GRNscores_dynamics_curve.SP_10.SP_8.pdf') #trend lines

#transitions=names(cellrouter@pathsinfo$path)
#grntransition(cellrouter, tfs, transitions, dir.targets='up',q.up=0.95, q.down=0.05, columns=2, 50, 50, "results/grntransition")

cluster1_genes=names(cellrouter@clusters$SP_10.SP_8$clustering[which(cellrouter@clusters$SP_10.SP_8$clustering==1)])
cluster2_genes=names(cellrouter@clusters$SP_10.SP_8$clustering[which(cellrouter@clusters$SP_10.SP_8$clustering==2)])
cluster3_genes=names(cellrouter@clusters$SP_10.SP_8$clustering[which(cellrouter@clusters$SP_10.SP_8$clustering==3)])
cluster4_genes=names(cellrouter@clusters$SP_10.SP_8$clustering[which(cellrouter@clusters$SP_10.SP_8$clustering==4)])
cluster5_genes=names(cellrouter@clusters$SP_10.SP_8$clustering[which(cellrouter@clusters$SP_10.SP_8$clustering==5)])
cluster6_genes=names(cellrouter@clusters$SP_10.SP_8$clustering[which(cellrouter@clusters$SP_10.SP_8$clustering==6)])
cluster7_genes=names(cellrouter@clusters$SP_10.SP_8$clustering[which(cellrouter@clusters$SP_10.SP_8$clustering==7)])
cluster8_genes=names(cellrouter@clusters$SP_10.SP_8$clustering[which(cellrouter@clusters$SP_10.SP_8$clustering==8)])
cluster9_genes=names(cellrouter@clusters$SP_10.SP_8$clustering[which(cellrouter@clusters$SP_10.SP_8$clustering==9)])
cluster10_genes=names(cellrouter@clusters$SP_10.SP_8$clustering[which(cellrouter@clusters$SP_10.SP_8$clustering==10)])

genelist <- c("GYPA","EPOR","SNCA")
plotDRExpression(cellrouter,genelist,TRUE,2,10,10,paste('results/',p,"_some_genes_DRE.pdf",sep='')) #set less ploting genes,show marker exressing region 
plotheatmap(cellrouter, names(cellrouter@dynamics), 0.7, crows=5, ccols=1, width=10, height=10, filename="results/cor.path_heatmap.SP_10.SP_8.pdf") 
dotplot2(cellrouter,cellrouter@sampTab,markers_select,"population",1,logtransform=TRUE,15,15,filename="results/dotplot2.markers_select.SP_10.SP_8..pdf")
plotclusters(cellrouter, p,2,10,10, filename="results/plotclusters.SP_10.SP_8")
plotpaths(cellrouter, p, genelist, columns=7, width=23, height=10, file_prefix="results/plotpaths_1") #curve plus scatter means expression
plottrajectories(cellrouter,p,genelist,rescale = TRUE,columns=1,width=10,height=6,filename=paste('results/',p,'.genelist_1_dynamics_curve.SP_10.SP_8.pdf',sep=''))

#### Pathway enrichment analysis on selected trajectories
clustergenes=c(cluster1_genes,cluster2_genes,cluster3_genes,cluster4_genes,cluster5_genes,cluster6_genes,cluster7_genes,cluster8_genes,cluster9_genes,cluster10_genes)
write.table(clustergenes,"results/clustergenes.SP_10.SP_8.txt",sep="\t")
length(clustergenes)
ids <- read.table("results/ids.SP_10.SP_8.txt",header=T,sep="\t",stringsAsFactors=F) #cook right ids object of 10 cluster genes with geneIDannotation function
dim(ids)
colnames(ids)=c('entrezgene','external_gene_name','description','cytogenetic_location')
paths <- names(cellrouter@pathsinfo$path)
cellrouter <- pathwayenrichment(cellrouter,paths,cc=NULL,'human','org.Hs.eg.db',ids)
enr_UP <- pathwaycluster(cellrouter,cellrouter@pathwayenrichment$UP$GOBP,30,TRUE,30,20,'results/UP_GOBP.SP_10.SP_8.pdf')
enr_DOWN <- pathwaycluster(cellrouter,cellrouter@pathwayenrichment$DOWN$GOBP,30,TRUE,30,20,'results/DOWN_GOBP.SP_10.SP_8.pdf')
##regulatornetwork for several TFs
library('ggnetwork')
library('GGally')
library('geomnet')
library('network')
library('sna')
#regulators=c("FOXO3","RNF10","SNCA","ZNF737","PER1","IRF3")
regulators=c("TERF2IP","NFIX","RNF10","GTF2B","SNCA","PBX1")
regulatornetwork(x, regulators, 4, 2, 10, 15, 'results/regulator_networks.TERF2IP.pdf')