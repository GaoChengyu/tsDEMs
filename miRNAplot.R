#miRNA plot
###PLOT pheatmap tpm
library(pheatmap)
#将相对丰度log2标准化
#scale_test <- apply(plot_mat, 2, function(x) log2(x+1))
plot_mat<-read.csv('d:/data/miRNA/vm-data/quantifier4/IL_VM/IL差异tpm.CSV',header = T,row.names = 1)
scale_test <- apply(plot_mat, 2, function(x) log2(x+1))
#添加列注释
matadata2<-data.frame(Sample=c('IB1','IB2','IB3','IL1','IL2','IL3','VM1','VM2','VM3'),
                      Group=c('IB:infected apple branch tissue','IB:infected apple branch tissue','IB:infected apple branch tissue','IL:infected apple leaf tissue','IL:infected apple leaf tissue','IL:infected apple leaf tissue','MVm:in vitro mycelium','MVm:in vitro mycelium','MVm:in vitro mycelium'))


anol_col<- data.frame(row.names = colnames(plot_mat), 
                      SampleGroup = matadata2$Group)

pheatmap(mat =scale_test,cluster_cols = F,annotation_col = anol_col,angle_col = "45")


######
#Correlation analysis and plot
#loading data
library(ggplot2)
matrix_data<-read.table('d:/data/miRNA/vm-data/quantifier4/VM_all_miRNA_expression.tsv',header = T,row.names = 1)
colnames(matrix_data)<-c('IB1','IB2','IB3','IL1','IL2','IL3','MVm1','MVm2','MVm3')
matrix.cor<-cor(matrix_data,method ='pearson' )
library(Hmisc)#add p-value
res <- rcorr(as.matrix(matrix.cor))
res$r<-round(res$r, 2)
CorMatrix <- function(cor,p) {
  ut <- upper.tri(cor,diag = T)
  data.frame(row = rownames(cor)[row(cor)[ut]],
             column = rownames(cor)[col(cor)[ut]], 
             cor =(cor)[ut], 
             p = p[ut] )
}            
df<-CorMatrix (res$r, res$P)
df[is.na(df)]<-0
colnames(df)<-c('var_x','var_y','value','signi')

library(dplyr)
df%>%
  mutate(label=case_when(
    signi<0.001 ~ "***",
    signi>0.001&signi<0.01 ~ "**",
    signi>0.01&signi<0.05 ~ "*",
    TRUE ~ ""
  )
  ) -> df1
df1%>%
  filter(var_x!=var_y) -> df2

ggplot(data=df1,aes(x=var_x,y=var_y))+
  geom_point(aes(size=value,color=value))+
  scale_color_gradient(low = "#80fcfe",high = "#ff80fc",
                       breaks=seq(-1,1,0.2))+
  scale_size_continuous(range = c(5,10))+
  guides(size='none')+
  theme_bw()+
  geom_text(data=df2,aes(x=var_y,y=var_x,
                         label=paste0(value,label)))+
  theme(legend.key.height = unit(4,'cm'),
        legend.justification = c(0,0),
        legend.title = element_blank())+
  xlab('')+
  ylab('')


###volcanic plot for DESeq2
library(ggplot2)
library(ggrepel)
library(tidyverse)
dataset<-read.csv(file = 'd:/data/miRNA/vm-data/quantifier4/IL_VM/All_results.csv',header = T, row.names=1, stringsAsFactors=FALSE, as.is = TRUE)
dataset<-na.omit(dataset)
#output <- 'E:/rnaseq/rnadiff_diffparts/live_dead/'
#output <- 'C:/Users/86183/Documents/Tencent Files/565811128/FileRecv/'
#set FDR and logFC cut_off
cut_off_fdr = 0.05
cut_off_logFC= 1


dataset$change = ifelse(dataset$padj < cut_off_fdr & abs(dataset$log2FoldChange) >= cut_off_logFC, 
                        ifelse(dataset$log2FoldChange> cut_off_logFC ,'Up','Down'),
                        'noSig')

head(dataset)

p <- ggplot(
  # draw plot
  dataset, aes(x = log2FoldChange, y = -log10(padj), colour=change)) +
  geom_point(alpha=1, size=2) +
  scale_color_manual(values=c("red", "#00B2FF","orange"))+
  # draw line
  geom_vline(xintercept=c(-1,1),lty=4,col="#990000",lwd=0.8) +
  geom_hline(yintercept = -log10(cut_off_fdr),lty=4,col="#990000",lwd=0.8) +
  # change labs
  labs(x="log2(fold change)",
       y="-log10 (FDR)")+
  theme_bw()+
  # set theme
  theme(plot.title = element_text(hjust = 0.5), 
        legend.position="right", 
        legend.title = element_blank())
  #ggtitle("Different expression milRNAs between IL and MVm")
p

#p+geom_text_repel(data =gene,aes(log2FoldChange,-log10(padj),label=rownames(gene)),max.overlaps=100)

#gene<-dataset[which(dataset$change=='Up'|dataset$change=='Down'),]


#####GO & KEGG
library(clusterProfiler)
library(tidyr)
gene1<-read.delim('D:/data/miRNA/vm-data/degradationData/DEG/DEGIB特异下调vmTARGET.txt',header = F)
gene1 <- gene1$V1[1:nrow(gene1)]#genelist必须为vector

term2gene <-read.delim('d:/db/genome_db/Valsa_mali/annotation/GO/term2gene.txt',header = T)
term2name<-read.delim('d:/db/genome_db/Valsa_mali/annotation/GO/term2name.txt',header = T)


df <- enricher(gene = gene1, TERM2GENE = term2gene, TERM2NAME = term2name, pvalueCutoff = 1, qvalueCutoff = 1)
input<-df@result
input<-separate(input,GeneRatio,c('Study.term','Study.total'),sep = '\\/')
input<-separate(input,BgRatio,c('Pop.term','Pop.total'),sep = '\\/')
go2ont <- go2ont(term2gene$GO)
mergedf<-merge(input,go2ont,by.x = 'ID',by.y = 'go_id',all.x = T)
write.table(mergedf,'D:/data/miRNA/vm-data/degradationData/DEG/DEGIL特异上调mdTARGETGOEA.txt',quote = F,row.names = F,sep = '\t')


library(clusterProfiler)
library(tidyr)
gene1<-read.delim('D:/data/miRNA/md-data/degradationData/DEM/tsIBDEM/tsIBDEM_target.txt',header = F)
gene1 <- gene1$V1[1:nrow(gene1)]#genelist必须为vector

ko2gene <-read.delim('d:/db/genome_db/Malus_domestica/annotation/KEGG/ko2gene.txt',header = T,colClasses ="character")
ko2name<-read.delim('d:/db/genome_db/Malus_domestica/annotation/KEGG/ko2name.txt',header = T,colClasses ="character")


df <- enricher(gene = gene1, TERM2GENE = ko2gene, TERM2NAME = ko2name, pvalueCutoff = 0.05, qvalueCutoff = 0.2)
input<-df@result
input<-separate(input,GeneRatio,c('Study.term','Study.total'),sep = '\\/')
input<-separate(input,BgRatio,c('Pop.term','Pop.total'),sep = '\\/')


write.table(input,'D:/data/miRNA/vm-data/degradationData/vmtargetgeneU/DEGvmTarget/DEGILvmTARGETas10KEGGEA.txt',quote = F,row.names = F,sep = '\t')


######milRNA靶基因GO气泡图
#4.1 bubble plot
library(squash)
library(graphics)
pdf("go_up_pcu.pdf",width = 6, height = 6)

layout(matrix(c(1,1,1,4,2,3),2,3,byrow = TRUE),c(1.5,1,1),c(10,1))#布局
#m<-matrix(1:4,2,2);m #建立矩阵m,2列2行
#layout(m,widths=c(1,3),heights=c(3,1)) #将当时装置按照m进行划分，宽度之比为1:3，高度之比为3:1
go<-read.delim("D:/data/miRNA/vm-data/degradationData/DEG/DEGIL特异上调mdTARGETGOEA.txt", header = T, sep = "\t", check.names = T)
#go <- go[go$p.adjust<0.05,]
go <- go[go$pvalue<0.05,]
go <- go[order(go$Ontology),]
gofold<-(go$Study.term/go$Study.total)/(go$Pop.term/go$Pop.total)
gof<--log10(go$pvalue)

par(mar=c(0,1,0,1))  #mar：以数值向量表示边界大小，顺序为"下、左、上、右"，单位为英分，默认值c(5, 4, 4, 2)+0.1

plot(gof,seq(1,length(go$Description)),   #plot(x,y)
     xlab="-log10(p)",
     ylab="",
     axes=F,
     xlim=c(0,46),  #设置图宽刻度分化
     col="white",
     frame.plot = T,
     adj=0.15
)

axis(1,at=seq(0,max(gof)+1,4))

ngo <- length(go$Ontology)
nmf <- sum(go$Ontology == "MF")
ncc <- sum(go$Ontology == "CC")
nbp <- sum(go$Ontology == "BP")

#abline(h=seq(0.5,50.5,1),lty=3,col="grey75",lwd=0.6)
#abline(v=seq(-4,max(gof)+1,1),lty=3,col="grey75",lwd=0.6)
#rect(xleft = 1, ybottom = 1, xright = 5, ytop = 5)
rect(max(gof)+1.2,0.8,47.63,0.8+nbp,col="#FFFFCC", border = "white")
rect(max(gof)+1.2,0.8+nbp,47.63,0.8+nbp+ncc,col="#CCCCFF", border = "white")
if (nmf==0){
  rect(max(gof)+1.2,0.8+nbp+ncc,47.63,0.8+nbp+ncc+nmf,col="#FFCCCC", border = "white")
}else{
  rect(max(gof)+1.2,0.8+nbp+ncc,47.63,0.2+nbp+ncc+nmf,col="#FFCCCC", border = "white")
}




points(gof,seq(1,length(go$Description)),
       col=jet(ngo)[rank(go$Study.term)],
       pch=16,
       cex=0.2+log(gofold,2)
)

text(max(gof)+2,seq(1,length(go$Description)),paste(go$ID,go$Description),cex=1.2,adj = 0)

abline(v=45.5,col="white", lwd=2)

abline(v=max(gof)+1)

if (nbp!=0){
  text(46.5,nbp/2+0.5,"BP",cex=1.2) 
}
if (ncc!=0){
  text(46.5,nbp+ncc/2+0.5,"CC",cex=1.2)  
}
if (nmf!=0){
  text(46.5,nbp+ncc+nmf/2+0.5,"MF",cex=1.2)
}




# legend 1
par(mar=c(1.7,1,1.1,2),xpd = T)
barplot(rep(1,11),border = "NA", space = 0,
        ylab="",
        xlab="",
        xlim=c(0.40,10.54),
        axes = F, col=jet(11) )
box()
text(1,1.6,"Less", cex = 0.9)
text(10,1.6,"More", cex = 0.9)
text(6,-0.8,"Number of enriched genes", adj=0.55, cex = 0.9)

# legend 2
par(mar=c(1.7,1,1.1,2),xpd = T)
plot(c(1:5),rep(1,5),
     ylab="",
     xlab="",
     xlim=c(1,5.2),
     axes = F,pch=16, cex=0.75*(0.2+fivenum(log(gofold,2))),col="green")
box()
text(1.15,1.9,"Low", cex = 0.9)
text(5,1.9,"High", cex = 0.9)
text(3,-0.1,"Enrichment factor (fold)", adj=0.5, cex = 0.9)

# legend main
plot(0,0, axes = F,col="NA")
text(-0.3,0.5,"-log10(p-value)", adj=0.5,pos = 1)

dev.off()


#########KEGG
library(ggplot2)


go<-read.delim("d:/data/miRNA/vm-data/degradationData/vmtargetgeneU/DEGvmTarget/DEGILvmTARGETas10KEGGEA.txt", header = T, sep = "\t", check.names = T)
go <- go[go$pvalue <0.05,]
go<-go[order(go$pvalue,decreasing = T),]


#go <- go[go$p<0.05,]
gofold<-(as.numeric(go$Study.term)/as.numeric(go$Study.total))/(as.numeric(go$Pop.term)/as.numeric(go$Pop.total))
gof<--log10(go$pvalue)

#Category <- go$name.space
Number <- go$Study.term

y<-factor(go$Description,levels = go$Description)


ggplot(go,aes(gofold,y))+
  geom_point(aes(size=gof,color=Number))+
  scale_color_gradient(low = "#6a82fb", high = "#fc5c7d")+ 
  labs(color='Counts',alpha="Number of enriched genes",size="−log10(pvalue)",x="Enrichment factor (fold)",y="",title="")+
  theme_bw()+theme(axis.line = element_line(colour = "black"), axis.text = element_text(color = "black",size = 14),legend.text = element_text(size = 14),legend.title=element_text(size=14),axis.title.x = element_text(size = 14))+
  scale_size_continuous(range=c(4,8))

#自定义气泡大小，以防太小看不见





##########################################
######################################
######################miRNA mRNA 共表达图
library(ggraph)
library(igraph)
library(tidyverse)
library(RColorBrewer)
######################
#prepare file 
datasetIB24<-read.csv('d:/data/冯老师组学/CleanData.outdir/24b_c/VM_IB.deseq2.outdir/All_results.csv',header = T)
cut_off_fdr = 0.05
cut_off_logFC= 1
datasetIB24$changeIB24 = ifelse(datasetIB24$padj < cut_off_fdr & abs(datasetIB24$log2FoldChange) >= cut_off_logFC, 
                                ifelse(datasetIB24$log2FoldChange> cut_off_logFC ,'Up','Down'),
                                'Stable')
datasetIB24<-datasetIB24[,c(1,8)]






datasetIL24<-read.csv('d:/data/冯老师组学/CleanData.outdir/24l_c/VM_IL.deseq2.outdir/All_results.csv',header = T)
cut_off_fdr = 0.05
cut_off_logFC= 1
datasetIL24$changeIL24 = ifelse(datasetIL24$padj < cut_off_fdr & abs(datasetIL24$log2FoldChange) >= cut_off_logFC, 
                                ifelse(datasetIL24$log2FoldChange> cut_off_logFC ,'Up','Down'),
                                'Stable')
datasetIL24<-datasetIL24[,c(1,8)]


dataset<-merge(datasetIB24,datasetIL24,by="X",all = T)

dataset[is.na(dataset)]<-"Stable"
write.csv(dataset,"d:/data/冯老师组学/CleanData.outdir/featureCounts/Vm_DEGs_message.csv",row.names = F)


peplist<-read.table("d:/data/miRNA/vm-data/degradationData/DEG/DEGIB特异下调vmTARGET.txt",header = F)
pep<-merge(peplist,dataset,by.x = "V1",by.y = "X", all.x = T)
pep<-pep[,c(1,2)]
colnames(pep)<-c("GeneID","change")
write.csv(pep,"d:/data/miRNA/vm-data/degradationData/DEG/DEGIB特异下调vmTARGETinRNAseq.txt",quote = F,row.names = F)
############################
#draw plot

#创建一个数据框，给出你个人的层次结构
d1<-data.frame(from='origin',to=paste('group',seq(1,2),sep=""))
milR<-read.table('d:/data/miRNA/vm-data/degradationData/DEG/PPI/tsIBDOWN_target/tsIBdownvmMIR.txt',header = F)
vmTar<-read.table('d:/data/miRNA/vm-data/degradationData/DEG/PPI/tsIBDOWN_target/DEGIB特异下调vmTARGET.txt',header = F)
#mdTar<-read.table('d:/data/miRNA/vm-data/degradationData/mdtargetgeneU/DEGmdTarget/DEGIBmdTARGET.txt',header = F)
d2<-data.frame(from=rep("group2",each=33),to=vmTar$V1)
d3<-data.frame(from=rep('group1',each=5),to=milR$V1)

#d4<-data.frame(from=rep('group3',each=479),to=mdTar$V1)
edges<-rbind(d1,d2,d3)


####load miRNA和靶基因对应关系表以及靶基因间PPI表
connect<-read.table('d:/data/miRNA/vm-data/degradationData/DEG/PPI/tsIBDOWN_target/tsIBdownPPI.txt',header = F)
colnames(connect)<-c('from','to')
connect$value<-runif(nrow(connect))

#创建顶点数据框，层次结构中的每个对象一行
vertices<-data.frame(
  name=unique(c(as.character(edges$from),as.character(edges$to)))
)
#用每个名称的组添加一个列
vertices$group<-edges$from[match(vertices$name,edges$to)]

#添加关于我们将要添加的标签的信息，角度，水平调整和潜在翻转计算标签的角度
vertices$id<-NA
myleaves<-which(is.na(match(vertices$name,edges$from)))
nleaves<-length(myleaves)
vertices$id[myleaves]<-seq(1:nleaves)
vertices$angle<-155-360*vertices$id/nleaves
#计算标签的对齐方式：左对齐还是右对齐
#如果在图的左侧，标签当前的角度为< -90
vertices$hjust<-ifelse(vertices$angle< -90,1,0)
#翻转角度BY使它们可读
vertices$angle<-ifelse(vertices$angle< -90,vertices$angle+180,vertices$angle)

#load data
#小RNA的靶基因在转录水平上的差异表达情况
#添加靶基因的基因表达情况 up down stable

vertices<-merge(vertices,pep,by.x = "name",by.y = "GeneID", all.x = T)

#创建graph对象
mygraph<-graph_from_data_frame(edges,vertices = vertices)
#连接对象必须引用叶节点的id
from<- match(connect$from,vertices$name)
to<- match(connect$to,vertices$name)

ggraph(mygraph,layout = 'dendrogram',circular=T)+
  geom_conn_bundle(data=get_con(from=from,to=to),alpha=0.2,width=0.9,aes(colour="skyblue"))+
  geom_node_text(aes(x=x*1.4,y=y*1.4,filter=leaf,label=name,angle=angle,hjust=hjust,colour=change),size=2,alpha=1)+
  geom_node_point(aes(filter=leaf,x=x*1.05,y=y*1.05,colour=change,alpha=0.2))+
  scale_color_manual(values = c('#7f00ff','#4ac29a','#FFA500','black','yellow'))+
  scale_size_continuous(range=c(0.1,10))+
  theme_void()+
  theme(
    legend.position = "none",
    plot.margin = unit(c(0,0,0,0),"cm"),
  )+
  expand_limits(x=c(-1.3,1.3),y=c(-1.3,1.3))