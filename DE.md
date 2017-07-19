# Differential expression analysis with DESeq2

## Upgrade R to the very latest (3.4.x)

```
sudo echo "deb http://cran.rstudio.com/bin/linux/ubuntu xenial/" | sudo tee -a /etc/apt/sources.list
gpg --keyserver keyserver.ubuntu.com --recv-key E084DAB9
gpg -a --export E084DAB9 | sudo apt-key add -

sudo apt-get update && sudo apt-get install -y r-base r-base-dev gdebi-core
```

## Make sure you're running RStudio

For this, we will again be working exclusively in RStudio!  Try to connect to a
running RStudio Web server instance -- you can get the Web address by
running this command:

```
echo My RStudio Web server is running at: http://$(hostname):8787/
```

## Install RStudio Web server

If you cannot connect, you'll need to install it:

```
wget https://download2.rstudio.org/rstudio-server-1.0.143-amd64.deb
sudo gdebi -n rstudio-server-1.0.143-amd64.deb
```

And, finally, change the password to something you can remember:
```
sudo passwd username
```

## Install DESeq2 prereqs

```
sudo apt-get install -y libxml2 libxml2-dev libcurl4-gnutls-dev libssl-dev
```

and then install DESeq2:

```
curl -O -L https://github.com/ngs-docs/angus/raw/2017/_static/install-deseq2.R

sudo Rscript --no-save install-deseq2.R
```

## Move salmon output quant files to their own directory

```
export PROJECT=/mnt/work
cd $PROJECT/quant
mkdir salmon_out

```

## Move the gene names to your home directory (to easily access it)
```
cp /mnt/work/annotation/trinity.nema.fasta.dammit/nema_gene_name_id.csv ~/
```

## Grab a special script plotPCAWithSampleNames.R

```
cd
wget https://raw.githubusercontent.com/ngs-docs/2017-dibsi-rnaseq/master/plotPCAWithSampleNames.R
```

## RStudio!

From this point on, we will be typing these commands into R Studio.

Load libraries
```
library(DESeq2)
library("lattice")
library(tximport)
library(readr)
library(gplots)
library(RColorBrewer)
source('~/plotPCAWithSampleNames.R')
```

Tell RStudio where your files are and ask whether they exist:
```
setwd("/mnt/work/quant/salmon_out/")
dir<-"/mnt/work/quant/"
files_list = list.files()
files <- file.path(dir, "salmon_out",files_list, "quant.sf")
names(files) <- c("0Hour_1","0Hour_2","0Hour_3","0Hour_4","0Hour_5","6Hour_1","6Hour_2","6Hour_3","6Hour_4","6Hour_5")
files
print(file.exists(files))
```

Grab the gene name and transcript ID file we generated yesterday: 
```
gene_names <- read.csv("~/nema_gene_name_id.csv")
cols<-c("row","transcript_id","gene_id")
colnames(gene_names)<-cols
tx2gene<-gene_names[,2:3]
head(tx2gene)
txi.salmon <- tximport(files, type = "salmon", tx2gene = tx2gene)
head(txi.salmon$counts)
dim(txi.salmon$counts)
```

Assign experimental variables:

```
condition = factor(c("0Hour","0Hour","0Hour","0Hour","0Hour","6Hour","6Hour","6Hour","6Hour","6Hour"))
ExpDesign <- data.frame(row.names=colnames(txi.salmon$counts), condition = condition)
ExpDesign
```

Run DESeq2:

```
dds <- DESeqDataSetFromTximport(txi.salmon, ExpDesign, ~condition)
levels(dds$condition) <- c("6Hour","0Hour")
```

Get counts:
```
counts_table = counts( dds, normalized=TRUE )
```

Filter out low expression transcripts:

```
filtered_norm_counts<-counts_table[!rowSums(counts_table==0)>=1, ]
filtered_norm_counts<-as.data.frame(filtered_norm_counts)
GeneID<-rownames(filtered_norm_counts)
filtered_norm_counts<-cbind(filtered_norm_counts,GeneID)
dim(filtered_norm_counts)
head(filtered_norm_counts)
```

Estimate dispersion:

```
plotDispEsts(dds)
```

PCA:
```
log_dds<-rlog(dds)
plotPCAWithSampleNames(log_dds, intgroup="condition", ntop=40000)
```

Get DE results:

```
res<-results(dds,contrast=c("condition","6Hour","0Hour"))
head(res)
res_ordered<-res[order(res$padj),]
GeneID<-rownames(res_ordered)
res_ordered<-as.data.frame(res_ordered)
res_genes<-cbind(res_ordered,GeneID)
dim(res_genes)
head(res_genes)
dim(res_genes)
res_genes_merged <- merge(res_genes,filtered_norm_counts,by=unique("GeneID"))
dim(res_genes_merged)
head(res_genes_merged)
res_ordered<-res_genes_merged_names_unique[order(res_genes_merged_names_unique$padj),]
write.csv(res_ordered, file="nema_DESeq_all.csv" )
```

Set a threshold cutoff of padj<0.05 and ± log2FC 1:

```
resSig = res_ordered[res_ordered$padj < 0.05, ]
resSig = resSig[resSig$log2FoldChange > 1 | resSig$log2FoldChange < -1,]
write.csv(resSig,file="nema_DESeq_padj0.05_log2FC1.csv")
```


MA plot with gene names:

```
resSig = res_ordered[res_ordered$padj < 0.05, ]
resSig = resSig[resSig$log2FoldChange > 1 | resSig$log2FoldChange < -1,]
plot(log2(res_ordered$baseMean), res_ordered$log2FoldChange, col=ifelse(res_ordered$padj < 0.05, "red","gray67"),main="nema (padj<0.05, log2FC = ±1)",xlim=c(1,20),pch=20,cex=1,ylim=c(-12,12))
abline(h=c(-1,1), col="blue")
genes<-resSig$GeneID
mygenes <- resSig[,]
baseMean_mygenes <- mygenes[,"baseMean"]
log2FoldChange_mygenes <- mygenes[,"log2FoldChange"]
text(log2(baseMean_mygenes),log2FoldChange_mygenes,labels=genes,pos=2,cex=0.60)
```

Heatmap

```
up_down<-resSig
dim(up_down)
up_down_FC<-subset(up_down,up_down$log2FoldChange>1 | up_down$log2FoldChange< -1)
dim(up_down_FC)
d<-up_down_FC
d<-na.omit(d)
dim(d)
head(d)
colnames(d)
d<-up_down_FC[,c(8:17)]
d<-as.matrix(d)
d<-as.data.frame(d)
d<-as.matrix(d)
rownames(d) <- up_down_FC[,1]
head(d)

hr <- hclust(as.dist(1-cor(t(d), method="pearson")), method="complete")
mycl <- cutree(hr, h=max(hr$height/1.5))
clusterCols <- rainbow(length(unique(mycl)))
myClusterSideBar <- clusterCols[mycl]
myheatcol <- greenred(75)
heatmap.2(d, main="nema (padj<0.05, log2FC = ±1)", 
          Rowv=as.dendrogram(hr),
          cexRow=0.75,cexCol=0.8,srtCol= 90,
          adjCol = c(NA,0),offsetCol=2.5, 
          Colv=NA, dendrogram="row", 
          scale="row", col=myheatcol, 
          density.info="none", 
          trace="none", RowSideColors= myClusterSideBar)
```


