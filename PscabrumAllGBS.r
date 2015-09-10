# TASSEL GBS HMN post-processing script by Justin Borevitz, with additional modifications by Chuah Aaron. (c) Australian National University, 2012-2015
## for scabrum
# Genotyping By Sequencing

# HapMapFilter.R
# to run from your computer first
# download R from http://r-project.org

# Adjusted for Nic Dussex - UOtago
library(compiler) 
R_COMPILE_PKGS=TRUE
R_ENABLE_JIT=3
enableJIT(3)
library(ape)
library(RColorBrewer)
library(NCBI2R)

plotPCoA <- function(pcoa,x=1,y=2,dataset,labels,pcoa.color,pcoa.glyph,legend,legend.color,dist.method) {
  plot(pcoa[,y]~pcoa[,x], xlab=paste0("PC",x," var ",per.var[x],"%"),
       ylab=paste0("PC",y," var ",per.var[y],"%"),
       main=paste0(dataset," ",dist.method," Principal Coordinate Analysis (PC",x," vs PC",y,")"),col=pcoa.color,bg=adjustcolor(pcoa.color,0.2),pch=pcoa.glyph,cex=4)
  #	text(pcoa[,y]~pcoa[,x], labels=labels, cex=1.5, col=pcoa.color)
  legend("topright", legend, col=legend.color, lty=1)
}

colLab <- function(n) { #helper function to color dendrogram labels
  if (is.leaf(n)) {
    a <- attributes(n)
    labCol<-color.pal[1]
    for(i in 1:length(categories)) {
      if(!is.na(pmatch(categories[i],a$label))) {
        labCol<-color.pal[i+1]
      }
    }
    attr(n, "nodePar") <- c(a$nodePar, lab.col=labCol, pch=NA)
  }
  n
}

dist.metric <- function(g.filt,method="hap") {
  if(method=="dist") {
    return(dist(t(g.filt)))
  }else if(method=="cor") {
    return(as.dist(1-cor(g.filt,use="pairwise.complete.obs")^2))
  }else if(method=="hap") {
    hap <- t(g.filt)	
    hap[hap==0] <- "c"
    hap[hap==1] <- NA
    hap[hap==2] <- "t"
    return(dist.dna(as.DNAbin(hap), model="N", pairwise.deletion=TRUE))
  }else if(method=="binary") {
    return(dist(t(g.filt),method=method))
  }else {
    stop(paste0("unsupported dist.metric method=",method))
  }
}

processHmb <- function(dataset,sample.snp.cutoff,site.thresh.low,dist.method="dist",hmbFile="fin.SAf2.HapMap.hmb.txt",hmnFile="fin.SAf2.HapMap.hmn.txt",label.fields=6,drop.samples=TRUE,remove.merged=TRUE,output.snps=FALSE,pdfFile="",groups=list(),remove=list(),output.nexus=1) {
  assign('sample.snp.cutoff',sample.snp.cutoff,meta)
  assign('site.thresh.low',site.thresh.low,meta)
  assign('remove',paste(remove,collapse=", "),meta)
  assign('groups',paste(groups,collapse=", "),meta)
  run.name<-paste0(sample.snp.cutoff,"_",site.thresh.low)
  cat(paste0("--\n",dataset,"\n",run.name,".",dist.method,"\n",hmbFile,"\n",hmnFile,"\n"))
  hmb <- read.table(hmbFile,header=TRUE,na.strings=".",row.names=1)
  g <- hmb[,-((ncol(hmb)-4):ncol(hmb))]
  hmn <- read.table(hmnFile,header=TRUE,na.strings=".",row.names=1)
  h <- hmn[,-(1:10)]
  cat(paste0("raw: ",ncol(h)," samples, ",nrow(h)," SNPs\n"))
  threshFile<-paste0(dataset,".",run.name,".pdf")
  assign('threshFile',threshFile,meta)
  pdf(file=threshFile,paper="A4r",width=11,height=7.5)
  image(1:nrow(h),1:ncol(h),as.matrix(h),xlab="SNPs",ylab="samples",main = paste0(dataset," (",ncol(g)," samples, ",nrow(h)," snps)"),useRaster=TRUE)
  empty.cols<-c(grep("empty",tolower(colnames(g))),grep("blank",tolower(colnames(g))))
  if(length(empty.cols)>0) {
    cat("Removing blank & empty columns: ",colnames(g)[empty.cols],"\n")
    g <- g[,-empty.cols] #to handle blank / empty columns in some datasets
    h <- h[,-empty.cols]
  }
  if(remove.merged) {
    merged.cols<-grep("_merged",colnames(g))
    if(length(merged.cols)>0) {
      cat("Removing TASSEL merged columns: ",colnames(g)[merged.cols],"\n")
      g <- g[,-merged.cols] #to handle blank merged columns in some datasets
      h <- h[,-merged.cols]
    }
  }else {
    colnames(g) <- sub("_merged_X(\\d)$","_m_0_PI0_A00_X\\1",colnames(g),perl=TRUE)
  }
  colnames(g)<-gsub("^X6m","6m",colnames(g),perl=TRUE)
  print(colnames(g))
  names.mat <- matrix(unlist(strsplit(colnames(g),split="_")),nrow=label.fields)
  print(names.mat)
  #	colnames(g) <- paste(names.mat[1,],names.mat[2,],sep='_')
  sampleSnps<-apply(g,2,function (x) sum(x!=0))
  #	colnames(g) <- paste0(names.mat[1,]," (",sampleSnps,")")
  colnames(g) <- names.mat[1,]
  
  if(length(remove)>0) {
    remove.cols<-unlist(sapply(remove,grep,colnames(g)))
    cat("Removing specified samples: ",colnames(g)[remove.cols],"\n")
    g <- g[,-remove.cols]
    h <- h[,-remove.cols]
    names.mat <- names.mat[,-remove.cols]
  }
  if(length(groups)>0) {
    categories<<-groups
  } 
  #	color.pal <<- c('black',rainbow(length(categories)-1),'gray') #global
  #color.pal <<- c('black',brewer.pal(length(categories),"Paired")) #global
  #color.pal[12] <<- 'gold'
  #color.pal <<- c('black',brewer.pal(length(categories),"Dark2")) #global
  
  samples<-apply(h,2,function(x) sum(!is.na(x)));
  samp.index <- samples>=sample.snp.cutoff
  print(table(samp.index))
  
  hist(samples,breaks=ncol(h),main=paste0(dataset," (",sum(samp.index)," / ",ncol(h)," samples, ",nrow(h)," SNPs)"), xlab="SNPs")
  abline(v=sample.snp.cutoff,col="red",lty=2)
  text(sample.snp.cutoff, 0, sample.snp.cutoff, col="red") 
  
  assign('samples',ncol(h),meta)
  assign('snps',nrow(h),meta)
  assign('samples.filtered',table(samp.index)[1],meta)
  assign('samples.remaining',table(samp.index)[2],meta)
  gg <- g[,samp.index]
  h <- h[,samp.index]
  names.mat <- names.mat[,samp.index]
  
  sites <- apply(h,1, function(x) sum(!is.na(x)));
  hist(sites,breaks=ncol(h),main = "Number of SNPs per Sample", xlab = "Samples")
  abline(v=site.thresh.low,lty=2,col="blue")
  text(site.thresh.low, 0, site.thresh.low, col="blue")
  #	abline(v=site.thresh.high,lty=2,col="darkblue")
  #	text(site.thresh.high, 0, site.thresh.high, col="darkblue")
  #	site.index <- sites > site.thresh.low & sites < site.thresh.high
  site.index <- sites > site.thresh.low
  print(table(site.index))
  assign('snps.filtered',table(site.index)[1],meta)
  assign('snps.remaining',table(site.index)[2],meta)
  #	g.minor <- gg
  g.minor <- gg[site.index,]
  h.minor <- h[site.index,]
  dim(g.minor)
  filteredSnps<-apply(g.minor,2,function (x) sum(x>0))
  #	colnames(g.minor)<-paste0(colnames(g.minor),' ',filteredSnps)
  
  samples.minor<-apply(h.minor,2,function(x) sum(!is.na(x)));
  hist(samples.minor,breaks=ncol(h.minor),main=paste0(dataset," filtered (",ncol(h.minor)," / ",ncol(h)," samples, ",nrow(h.minor)," SNPs)"), xlab="SNPs")
  #	abline(v=sample.snp.cutoff,col="red",lty=2)
  #text(sample.snp.cutoff, 0, sample.snp.cutoff, col="red") 
  
  h.na <- is.na(h)
  print(table(h.na))
  assign('genotypes',nrow(h.na)*ncol(h.na),meta)
  assign('genotypes.filled',table(h.na)[1],meta)
  assign('genotypes.na',table(h.na)[2],meta)
  
  g.minor.filtered <- g.minor
  d.dist.filtered <<- dist.metric(g.minor.filtered,method=dist.method)
  snps<-nrow(g.minor.filtered)/2
  samps<-ncol(g.minor.filtered)
  samp.names<-names(g.minor.filtered)
  cat(paste0("filt: ",samps," samples, ",snps," SNPs\n\n"))
  assign('dataset',dataset,meta)
  assign('snps.filt',snps,meta)
  assign('samples.filt',samps,meta)
  outfile.pref<-paste0(dataset,".",samps,"samp.",snps,"snp")
  if(output.snps==TRUE) {
    write.csv(g.minor.filtered,file=paste0(outfile.pref,".csv"))
  }
  
  image(1:nrow(h.minor),1:ncol(h.minor),as.matrix(h.minor),xlab="SNPs",ylab="samples",main = paste0(dataset,' filtered (',samps," samples, ",snps," snps)"),useRaster=TRUE)
  dev.off()
  
  if(pdfFile=="") {
    pdfFile<-paste0(outfile.pref,".",dist.method,".pdf")
  }
  assign('pdfFile',pdfFile,meta)
  pdf(file=pdfFile,paper="A4r",width=11,height=7.5)
  par(cex=140/(140+samps),xpd=TRUE)# adjust text size
  #with "other"
  #	legend<-c(categories,"other")
  #	legend.color<-color.pal[c(2:(length(categories)+1),1)]
  #without #other
  legend<-c(categories)
  legend.color<-color.pal[c(2:(length(categories)+1))]
  
  write.csv(as.matrix(d.dist.filtered),file=paste0(outfile.pref,".",dist.method,".csv"))
  # dendrogram
  hc <- hclust(d.dist.filtered)
  hc <- dendrapply(as.dendrogram(hc,hang=0.1),colLab2)
  plot(hc,main=paste0(dataset," Cluster Dendrogram (",samps," samples, ",snps," SNPs, ",dist.method," method)"),cex=0.07,edge.width=0.3)
  legend("bottomright", legend, col=legend.color, lty=1)
  
  # multidimensional scaling.  Principal Coordinate analysis
  pcoa <- cmdscale(d.dist.filtered, k=8)
  pcoa_df <- data.frame(pcoa)
  tot.var <- sum(sapply(pcoa_df, var))
  per.var <<- round((sapply(pcoa_df, var)/tot.var)*100,2)
  
  pcoa.color <- rep(color.pal[1],ncol(g.minor.filtered))
  for(i in 1:length(categories)) {
    pcoa.color[grep(categories[i],samp.names)]<-color.pal[i+1]
  }
  cat.glyphs<-c('U','E','K')
  glyphs<-c(21,3,24) # circle, cross, square
  pcoa.glyph <- rep(23,ncol(g.minor.filtered))
  for(i in 1:length(cat.glyphs)) {
    pcoa.glyph[grep(cat.glyphs[i],samp.names)]<-glyphs[i]
  }
  
  plotPCoA(pcoa,1,2,dataset,samp.names,pcoa.color,pcoa.glyph,legend,legend.color,dist.method)
  #	plotPCoA(pcoa,1,3,dataset,samp.names,pcoa.color,pcoa.glyph,legend,legend.color,dist.method)
  plotPCoA(pcoa,2,3,dataset,samp.names,pcoa.color,pcoa.glyph,legend,legend.color,dist.method)
  plotPCoA(pcoa,3,4,dataset,samp.names,pcoa.color,pcoa.glyph,legend,legend.color,dist.method)
  plotPCoA(pcoa,4,5,dataset,samp.names,pcoa.color,pcoa.glyph,legend,legend.color,dist.method)
  dev.off() # close off the .pdf to finish assignment
  
  if(output.nexus>0) { #output NEXUS, numerotypes (fast) if 1, genotypes (slow!) if 2, both if 3
    d<-as.matrix(t(h))
    d[is.na(d)]<-'.'
    sample<-as.list(rownames(d))
    samples<-length(sample)
    snp<-as.list(colnames(d))
    snps<-length(snp)
    sample.names<-paste(sample,sep="\n",collapse="\n")
    if(bitwAnd(output.nexus,1)) { #numerotype NEXUS
      printrow<-function(x) paste0(sample[[x]],"\t",paste0(d[x,],collapse=""))
      print(system.time(numerotypes<-paste0(lapply(1:samples,printrow),collapse="\n")))
      cat(paste0('#NEXUS\n\n[!Data from: ',dataset,' GBS hmn.txt]\n\nbegin taxa;\n\tdimensions ntax=',samples,';\n\ttaxlabels\n',sample.names,"\n\t;\nend;\n\nbegin characters;\n\tdimensions nchar=",snps,';\n\tformat missing=. symbols="012";\n\tmatrix\n',numerotypes,"\n\t;\nend;\n\nbegin assumptions;\n\toptions deftype = ord;\nend;\n"),file=paste0(outfile.pref,".nex"))
    }
    if(bitwAnd(output.nexus,2)) { #genotype NEXUS
      nucs<-c("A","C","G","T")
      alleles<-sort(apply(expand.grid(nucs,nucs),1,paste,collapse="/",sep=""))
      iupac<-ConvertIUPAC(alleles)
      names(iupac)<-alleles
      hmn2hmp<-function(x) c(c(".",substr(x,1,1),iupac[x],substr(x,3,3)))
      var<-data.frame(allele=hmn$alleles,row.names=rownames(hmn))
      geno<-as.matrix(t(apply(var,1,hmn2hmp)))
      colnames(geno)<-c(".","0","1","2")
      lookup.geno<-function(y,x) geno[[snp[[y]],d[x,y]]]
      printrow2<-function(x) paste0(sample[[x]],"\t",paste0(lapply(1:snps,lookup.geno,x),collapse=""),collapse="")
      print(system.time(genotypes<-paste0(lapply(1:samples,printrow2),collapse="\n")))
      cat(paste0('#NEXUS\n\n[!Data from: ',dataset,' GBS hmp.txt]\n\nbegin taxa;\n\tdimensions ntax=',samples,';\n\ttaxlabels\n',sample.names,"\n\t;\nend;\n\nbegin characters;\n\tdimensions nchar=",snps,';\n\tformat missing=. datatype=DNA;\n\tmatrix\n',genotypes,"\n\t;\nend;\n\nbegin assumptions;\n\toptions deftype = ord;\nend;\n"),file=paste0(outfile.pref,".geno.nex")) 
    }
  }
  cat(paste0("--\n",meta[["dataset"]],":-\n\nExperiment:\nA total of ",meta[["samples"]]," ",meta[["dataset"]]," samples with ",meta[["snps"]]," SNPs called by TASSEL UNEAK were analysed using the provided R script.\n\nSamples:\nThe following samples were identified beforehand as outliers and were removed from further analysis: ",meta[["remove"]],"\nThe remaining samples were grouped into these categories: ",meta[["groups"]],"\n\nResults:\nFrom the histogram produced in ",meta[["threshFile"]],", a per-sample minimum SNP threshold of ",meta[["sample.snp.cutoff"]]," SNPs was selected, which removed ",meta[["samples.filtered"]]," samples. The remaining genotype matrix had a call rate of ",sprintf("%5.2f%%",meta[["genotypes.filled"]]/meta[["genotypes"]]*100)," (",meta[["genotypes.filled"]]," / ",meta[["genotypes"]],"). The final genotype matrix had  ",meta[["snps.filt"]]," SNPs across ",meta[["samples.filt"]]," samples.\n\nA distance matrix was then computed, and a dendrogram plotted using hclust(), and Principal Coordinate analysis was performed via cmdscale() with pairwise plots of the top 5 dimensions.\n\nConclusion:\n"))
}

colLab2 <- function(n) { #helper function to color dendrogram labels
  if (is.leaf(n)) {
    a <- attributes(n)
    labCol<-color.pal[1]
    for(i in 1:length(categories)) {
      #			family<-matrix(unlist(strsplit(a$label,split="-")),nrow=2)[2]
      #			if(!is.na(pmatch(category.vals[i],family))) {
      if(length(grep(categories[i],a$label))) {
        labCol<-color.pal[i+1]
      }
    }
    attr(n, "nodePar") <- c(a$nodePar, lab.col=labCol, pch=NA)
  }
  n
}

processHmc <- function(dataset,dist.method="binary",hmcFile="",label.fields=5,drop.samples=TRUE,remove.merged=TRUE,output.snps=FALSE,pdfFile="",groups=list(),remove=list(),output.nexus=1) {
  cat(paste0("--\n",dataset,"\n",dist.method,"\n",hmcFile,"\n"))
  hmc <- read.table(hmcFile,header=T)
  std.head <- c("rs","HetCount_allele1","HetCount_allele2","Count_allele1","Count_allele2","Frequency")
  
  #setup matrix of sample rows 1x with allele pairs 2x
  hmc.allele <- apply(hmc[,!colnames(hmc)%in%std.head],1,function (x) unlist(strsplit(x,split="|",fixed=T)))
  sampN <- nrow(hmc.allele)/2
  snpN <- ncol(hmc.allele)
  #names.list <- list(colnames(hmc[,!colnames(hmc)%in%std.head]),paste(rep(hmc$rs,each=2),1:2,sep="_")  )
  short.names <- matrix(unlist(strsplit(colnames(hmc[,!colnames(hmc)%in%std.head]),split="_")),nr=label.fields)[1,]
  names.list <- list(short.names,paste(rep(hmc$rs,each=2),1:2,sep="_")  )
  
  #split genotypes into allele groups 2x samples 1x snps
  hmc.allele2 <- matrix(ncol=2*snpN,nrow=sampN,dimnames = names.list)
  # fill allele pairs
  hmc.allele2[,seq(1,snpN*2,2)] <- as.numeric(hmc.allele[seq(1,sampN*2,2),])
  hmc.allele2[,seq(2,snpN*2,2)] <- as.numeric(hmc.allele[seq(2,sampN*2,2),])
  #call Presence/Absense
  hmc.allele01 <- hmc.allele2
  hmc.allele01[hmc.allele2 != 0] <- 1
  
  # look at coverage across samples
  reads.samp <- rowSums(hmc.allele2)
  alleles.samp <- rowSums(hmc.allele01)
  # VISUALLY INSPECT BAD SAMPLES
  threshFile<-paste0(dataset,".hmc.pdf")
  assign('threshFile',threshFile,meta)
  if(pdfFile=="") {
    pdfFile<-paste0(dataset,'.pdf')
  }
  pdf(file=pdfFile,paper="A4r",width=11,height=7.5)
  
  plot(alleles.samp,reads.samp,xlab = "Alleles Called per Sample", ylab = "Total Reads per Sample")
  s.cuts <- 2000 #need to edit this for each experiment
  abline(v=s.cuts)
  keep <- alleles.samp>s.cuts
  table(keep)
  hmc01 <- hmc.allele01[keep,]
  
  # look at read coverage
  reads.allele <- colSums(hmc.allele2[keep,])
  samps.allele <- colSums(hmc01)
  plot(samps.allele, reads.allele,xlab = "Samples with Allele Called", ylab = "Total Reads per Allele",
       pch=".",ylim=c(0,10000))
  quantile(reads.allele,1:20/20)
  quantile(samps.allele,1:20/20)
  freq.allele <- samps.allele/sum(keep)
  # remove rare and repeat alleles
  s.cuts <- c(-0.1,0.05,0.95,1.1) # must be observed in 5% of samples but not more that 95%
  s.cuts <- c(-0.1,0.1,0.9,1.1) # must be observed in 10% of samples but not more that 90%
  
  abline(v=s.cuts*sum(keep))
  keep.allele <- cut(freq.allele,s.cuts,labels = c("rare","mid","repeat"))
  table(keep.allele)
  table(is.na(keep.allele) ) #verify nothing missing because exactly on threshold
  hmc01 <- hmc01[,keep.allele=="mid"]  # strange bug, check diminsions
  dim(hmc01)
  d.dist.filtered <<- dist.metric(hmc01,method=dist.method)
  categories<<-groups
  category.vals<<-as.character(levels(as.factor(groups)))
  #	library(RColorBrewer)
  #	color.pal<<-c('black',brewer.pal(length(category.vals),"Paired")) #big palette
  hc <- hclust(as.dist(d.dist.filtered))
  hc <- dendrapply(as.dendrogram(hc,hang=0.1),colLab2)
  par(cex=0.2)
  plot(hc)
  
}

meta<-new.env()
setwd("~/Documents/GDU output/final/aaronrun3a") #comment to output to current working dir
sample.groups<-c("^U","^UD","^UP","^UDP","^M","^MB","^MG","^ME","^P","^PM","^PH","^PUD","K","^AU","^E","^HB","^ORW","^NR","^SEW","^VR")
color.pal<<-c('black','green1','green2','green3','green4','purple1','purple2','purple3','purple4','orange1','orange2','orange3','orange4','blue1','blue2','blue3','blue4','red1','red2','red3','red4')
sample.remove<-c()

processHmb("Pscabrum SAf2 (1200,3pc)",1200,floor(517*0.03),dist.method="binary",groups=sample.groups,remove=sample.remove,remove.merged=FALSE,drop.samples=FALSE,output.nexus=1,output.snps=TRUE)

