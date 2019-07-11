setwd("/mnt/bigdata/Genetic/Projects/shg047/methylation/Pancancer")
newimp<-head(RLT,n=200)
TSG<-read.table("https://raw.githubusercontent.com/Shicheng-Guo/AnnotationDatabase/master/TSGene2.0.txt",head=T,sep="\t")
TSGTARGET<-newimp[newimp$Symbol %in% TSG[,2],]
write.table(TSGTARGET,file=paste("../../RF_LUAD_LUSC.BUR.PAN.Top6000VIP.TSG.txt",sep=""),sep="\t",quote=F,col.names = NA,row.names = T)

newimp<-head(RLT,n=170)
GWAS<-read.table("https://raw.githubusercontent.com/Shicheng-Guo/HowtoBook/master/TCGA/lung/lung_GWAS.Genelist.txt")
DMRGWAS<-na.omit(data.frame(newimp,GWAS[match(newimp$Symbol,GWAS[,1]),]))
unique(DMRGWAS$Symbol)
write.table(DMRGWAS,file=paste("../../RF_LUAD_LUSC.BUR.PAN.Top6000VIP.WGBS.txt",sep=""),sep="\t",quote=F,col.names = NA,row.names = T)

newimp<-head(RLT,n=3000)
DEG<-read.table("https://raw.githubusercontent.com/Shicheng-Guo/HowtoBook/master/TCGA/Pancancer_FPKM_UQ/TCGA-Pancancer-RNAseq-FPKM-UQ.Meta.diff.Symbol.txt",sep="\t",head=T)
DMRDEG<-data.frame(newimp,DEG[match(newimp$Symbol,DEG$hgnc_symbol),])
DMRDEGTarget<-subset(DMRDEG,beta<0 & pval<10^-5)
unique(DMRDEGTarget$hgnc_symbol)
ADOS<-read.table("~/hpc/methylation/TCGA_LUAD_FPKM-UQ.DGE_OS_HR_PanDiff.All.txt",head=T,sep="\t")
SCOS<-read.table("~/hpc/methylation/TCGA_LUSC_FPKM-UQ.DGE_OS_HR_PanDiff.All.txt",head=T,sep="\t")
ADOSS<-subset(ADOS,Estimate<0 & beta<0 & Pr...t..<10^-5 & Pr...z..<0.01)
SCOSS<-subset(SCOS,Estimate<0 & beta<0 & Pr...t..<10^-5 & Pr...z..<0.01)
xx<-newimp[which(newimp$Symbol %in% ADOSS$Symbol),]
yy<-newimp[which(newimp$Symbol %in% SCOSS$Symbol),]
MethPANExpOS<-unique(rbind(xx,yy))
unique(MethPANExpOS$Symbol)
write.table(MethPANExpOS,file=paste("../../RF_LUAD_LUSC.BUR.PAN.Top6000VIP.OS.txt",sep=""),sep="\t",quote=F,col.names = NA,row.names = T)
