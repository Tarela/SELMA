#!/usr/bin/env python
"""
single cell ATAC-seq clustering methods
"""
import subprocess
import sys
import os
import math
import time
import random
import string
import gzip
import numpy

def sp(cmd): 
    '''
    Call shell cmd or software and return its stdout
    '''
    a=subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell='TRUE')
    ac = a.communicate()
    return ac

def sperr(cmd):
    '''
    Call shell cmd or software and return its stdout
    '''
    a=subprocess.Popen(cmd, stderr=subprocess.PIPE, shell='TRUE')
    ac = a.communicate()
    return ac

def scClustering_correction_matrix(outname,peakMaxReads):
    scRscript="""
count_nonZero <- function(inline){
    return(length(which(inline>0)))
}
set.seed(1)

outname <- "%s"
peakMaxReads <- %s

data <- read.table(paste0(outname,"_peakXcellMat.txt"),row.names=1,header=T)
regionFeature <- read.table(paste0(outname,"_peakFeatures.txt"),row.names=4,header=T)
peaknum <- nrow(regionFeature)
useCell <- colnames(data)

x <- ((1:100)-1)/100
useY <- 12*(x)*((1-(x))**2)
if(as.numeric(peakMaxReads) != 0){
  usePeak <- regionFeature[which(regionFeature[,"cutsSum"] <= as.numeric(peakMaxReads)),]
}else{
  usePeak <- regionFeature
}
rawdata <- data[rownames(usePeak)[order(usePeak[,"avebias"])],useCell]
correctdata_raw <- c()
each <- round(nrow(rawdata)/100)
for(i in 0:98){
        tmp <- rawdata[(each*i+1):(each*(i+1)),]
        outtmp <- tmp*useY[i+1]
        correctdata_raw <- rbind(correctdata_raw,outtmp)
}
i <- 99
tmp <- rawdata[(each*i+1):(nrow(rawdata)),]
outtmp <- tmp*useY[i+1]
correctdata_raw <- rbind(correctdata_raw,outtmp)
peak_count_non0 <- apply(correctdata_raw,1,count_nonZero)
usedata <- correctdata_raw[which(peak_count_non0>0),]
y1raw <- usedata
y1scale <- y1raw  * as.numeric(sum(rawdata) / sum(y1raw))
y1scaleRound <- round(y1scale)
usedata0 <- y1scale

peak_count_non0_row <- apply(usedata0,1,count_nonZero)
peak_count_non0_col <- apply(usedata0,2,count_nonZero)
usedata <- usedata0[which(peak_count_non0_row>0),which(peak_count_non0_col>0)]

write.table(usedata,file=paste0(outname,"_correctionMat.txt"),row.names=T,col.names=T,sep="\t",quote=F)

"""%(outname,peakMaxReads)
    outf = open(outname+"_correctionMatRscript.r",'w')
    outf.write(scRscript)
    outf.close()
    tmplog = sp("Rscript %s_correctionMatRscript.r"%outname)    
#    if "simpleError" in tmplog[0].decode("ascii") and "noPackage" in tmplog[0].decode("ascii"):
#        return("noPackage")
#    else:
#        return("yesPackage")
def RDreadloci(peakloci):
    chrm = peakloci[0]
    start = int(peakloci[1])
    end = int(peakloci[2])
    return [chrm] + sorted(random.sample(range(start,end),2))

def RDreads_correctionMat(outname):
    matfile =open("%s_correctionMat.txt"%outname)
    peakfile = open("%s_peakFeatures.txt"%outname)
     
    peakinfo = {}
    for line in peakfile:
        if line.startswith("chrm") and "cutsSum" in line:
            peakcolname = line
            continue
        ll = line.split()
        peakinfo[ll[3]] = ll
    peakfile.close()

    usepeakinfo = {}

    outreads = open("%s_correctionMatRDreads.bed"%outname,'w')
    cellnames = matfile.readline().split()#[1:]
    for line in matfile:
        ll = line.split()
        peakname = ll[0]
        readscount = list(map(int,map(float,ll[1:])))#list(map(int,ll[1:]))
        usepeakinfo[peakname] = peakinfo[peakname]
        peakloci = peakinfo[peakname][:3]
        #center = int((int(peakloci[1]) + int(peakloci[2]))/2)
        #readsloci = [peakloci[0], center-25,center+25 ]
        for idx in range(len(readscount)):
            if readscount[idx] > 0:
                thiscell = cellnames[idx]
                for C in range(readscount[idx]):
                    readsloci = RDreadloci(peakloci)
                    readsinfo = readsloci + [thiscell,1]
                    outreads.write("\t".join(map(str,readsinfo))+"\n")
    outreads.close()
    
    outpeaks = open("%s_correctionMatPeaks.txt"%outname,'w')
    outpeaks.write(peakcolname)
    for peak in usepeakinfo.keys():
        newll = usepeakinfo[peak]
        outpeaks.write("\t".join(map(str,newll))+"\n")
    outpeaks.close()

def scClustering_PCAkm(outname,Percent,clusterNum,topDim,peakMaxReads,makeUMAP,SCcorrection):
    scRscript="""
count_nonZero <- function(inline){
    return(length(which(inline>0)))
}
set.seed(1)

outname <- "%s"
Percent <- %s
clusterNum <- %s
topDim <- %s
makeUMAP <- %s
SCcorrection <- %s
peakMaxReads <- %s

data <- read.table(paste0(outname,"_peakXcellMat.txt"),row.names=1,header=T)
regionFeature <- read.table(paste0(outname,"_peakFeatures.txt"),row.names=4,header=T)
peaknum <- nrow(regionFeature)
useCell <- colnames(data)

if(SCcorrection == 0){
  keep_percent <- as.numeric(Percent)
  usepeak0 <- regionFeature[ order(regionFeature[,"avebias"],decreasing=T)[ round(peaknum* (1-keep_percent/100) ) : peaknum ], ]
  if(as.numeric(peakMaxReads) != 0){
    usepeak <- usepeak0[which(usepeak0[,"cutsSum"] <= as.numeric(peakMaxReads)),]
  }else{
    usepeak <- usepeak0
  }
  usedata0 <- data[rownames(usepeak),useCell]
  
}else{
  x <- ((1:100)-1)/100
  useY <- 12*(x)*((1-(x))**2)
  if(as.numeric(peakMaxReads) != 0){
    usePeak <- regionFeature[which(regionFeature[,"cutsSum"] <= as.numeric(peakMaxReads)),]
  }else{
    usePeak <- regionFeature
  }
  rawdata <- data[rownames(usePeak)[order(usePeak[,"avebias"])],useCell]
  correctdata_raw <- c()
  each <- round(nrow(rawdata)/100)
  for(i in 0:98){
          tmp <- rawdata[(each*i+1):(each*(i+1)),]
          outtmp <- tmp*useY[i+1]
          correctdata_raw <- rbind(correctdata_raw,outtmp)
  }
  i <- 99
  tmp <- rawdata[(each*i+1):(nrow(rawdata)),]
  outtmp <- tmp*useY[i+1]
  correctdata_raw <- rbind(correctdata_raw,outtmp)
  peak_count_non0 <- apply(correctdata_raw,1,count_nonZero)
  usedata <- correctdata_raw[which(peak_count_non0>0),]
  y1raw <- usedata
  y1scale <- y1raw  * as.numeric(sum(rawdata) / sum(y1raw))
  y1scaleRound <- round(y1scale)
  usedata0 <- y1scale
}

peak_count_non0_row <- apply(usedata0,1,count_nonZero)
peak_count_non0_col <- apply(usedata0,2,count_nonZero)
usedata <- usedata0[which(peak_count_non0_row>0),which(peak_count_non0_col>0)]


if(SCcorrection == 1){
  write.table(usedata,file=paste0(outname,"_correctionMat.txt"),row.names=T,col.names=T,sep="\t",quote=F)
}

PCAdata_scale <- prcomp(t(scale(usedata)))
cellname <- as.vector(rownames(PCAdata_scale$x))
CTnum <- clusterNum
PCAdata_scale_use <- PCAdata_scale$x[,1:min(topDim,ncol(usedata))]
rownames(PCAdata_scale_use) <- colnames(usedata)

set.seed(1)
km_use_scale <- kmeans(PCAdata_scale_use, CTnum,iter.max = 100)
clusterInfo <- km_use_scale$cluster

outdata <- cbind(names(clusterInfo), clusterInfo)
colnames(outdata) <- c("cellname","cluster")
write.table(outdata[sort(rownames(outdata)),],file=paste0(outname,"_scClusters.txt"),row.names=F,col.names=T,sep="\t",quote=F)


if(makeUMAP == 1){
  if(require("umap")){
    if(nrow(outdata)<1000){
      PCH <- 16
    }else{
      PCH <- "."
    }
    pdf(file=paste0(outname,"_clusteringUMAP.pdf"))
    layout(matrix(c(1,2),nrow=1),width=c(4,1))
    par(mar=c(4,4,2,1))
      umapdata <- umap(PCAdata_scale_use)$layout
      colorMat <- rep("black",nrow(umapdata))
      names(colorMat) <- rownames(umapdata)
      rain <- rainbow(length(sort(unique(clusterInfo))))
      for(i in sort(unique(clusterInfo))){
        colorMat[names(clusterInfo)[which(clusterInfo==i)]] <- rain[i]
      }
      plot(umapdata[,1],umapdata[,2],col=colorMat[rownames(umapdata)],pch=PCH,xlab="UMAP-1",ylab="UMAP-2",main=paste0(outname," UMAP"))
    par(mar=c(4,1,2,1))
    plot(1,1,type="n",xlab="",ylab="",main="",axes=F)
    legend("center",legend=sort(unique(clusterInfo)),col=rain,bty="n",pch=16)
    dev.off()        
  }else{
      simpleError("noPackage")
  }        
}

"""%(outname,Percent,clusterNum,topDim,makeUMAP,SCcorrection,peakMaxReads)
    outf = open(outname+"_scRscript.r",'w')
    outf.write(scRscript)
    outf.close()
    tmplog = sp("Rscript %s_scRscript.r"%outname)    
    if "simpleError" in tmplog[0].decode("ascii") and "noPackage" in tmplog[0].decode("ascii"):
        return("noPackage")
    else:
        return("yesPackage")

def scClustering_ArchR(outname,GENOME,Percent,topDim,peakMaxReads,SCcorrection,makeUMAP,clusterMethod):
    check_bgzip = sp("which bgzip")
    check_tabix = sp("which tabix")
    if check_bgzip[0].decode("ascii") == "" or check_tabix[0].decode("ascii") == "" :
        return("noPackage")

    if int(SCcorrection) == 0: 
        cmd1 = """awk '{OFS="\t";print $1,$2,$3,$4,1}' %s_highQcellReads.bed | sort -k 1,1 -k 2,2g -k 3,3g -k 4,4  | bgzip > %s_ArchRReads.bed.gz"""%(outname,outname)
    else:
        cmd1 = """awk '{OFS="\t";print $1,$2,$3,$4,1}' %s_correctionMatRDreads.bed | sort -k 1,1 -k 2,2g -k 3,3g -k 4,4  | bgzip > %s_ArchRReads.bed.gz"""%(outname,outname)
    
    cmd2 = """tabix -p bed %s_ArchRReads.bed.gz"""%(outname)
    tmplog = sp(cmd1)
    tmplog = sp(cmd2)
    
    if os.path.isfile("%s_ArchRReads.bed.gz"%(outname)) and os.path.isfile("%s_ArchRReads.bed.gz.tbi"%(outname)):
        pass
    else:
        return("noTabix")
    
    scRscript="""
if(require("ArchR")){
  set.seed(1)
  
  outname <- "%s"
  GENOME <- "%s"
  Percent <- %s
  topDim <- %s
  makeUMAP <- %s
  peakMaxReads <- %s
  SCcorrection <- %s
  clusterMethod <- "%s"

  dir.create(paste0(outname,"_ArchR"))
  setwd(paste0(outname,"_ArchR"))
  set.seed(1)
  addArchRGenome(GENOME)
  addArchRThreads(threads = 1)
  
  getCN <- function(inname){
      return(strsplit(inname,"#")[[1]][2])
  }
  
  inputFiles <- c(paste0("../",outname,"_ArchRReads.bed.gz"))
  names(inputFiles)<-c("combine")
  
  ArrowFiles <- createArrowFiles(
    inputFiles = inputFiles,
    sampleNames = names(inputFiles),
    minTSS = 0, #Dont set this too high because you can always increase later
    minFrags = 0, maxFrags=1e+10,
    addTileMat = F,
    force=T,
    addGeneScoreMat = F
  )
  
  proj1 <- ArchRProject(
    ArrowFiles = ArrowFiles,
    outputDirectory = outname,
    copyArrows = FALSE 
  )
  
  if(SCcorrection==0){
    regionFeature <- read.table(paste0("../",outname,"_peakFeatures.txt"),row.names=4,header=T)
    peaknum <- nrow(regionFeature)
    keep_percent <- as.numeric(Percent)
    usepeak0 <- regionFeature[ order(regionFeature[,"avebias"],decreasing=T)[ round(peaknum* (1-keep_percent/100) ) : peaknum ], ]
    usePeak <- usepeak0[which(usepeak0[,"cutsSum"] <= peakMaxReads),]
    peakdata <- usePeak
    peaknum <- nrow(peakdata)
  }else{
    regionFeature <- read.table(paste0("../",outname,"_correctionMatPeaks.txt"),row.names=4,header=T)
    peakdata <- regionFeature[which(regionFeature[,"cutsSum"] <= peakMaxReads),]
    peaknum <- nrow(peakdata)
  }
  
  peakGR <- GRanges(seqnames=peakdata[,1],ranges=IRanges(peakdata[,2],peakdata[,3]))
  proj1 <- addPeakSet(ArchRProj=proj1, peakSet=peakGR)
  proj1 <- addPeakMatrix(proj1)
  
  set.seed(1)
  proj1 <- addIterativeLSI(
      ArchRProj = proj1,
      useMatrix = "PeakMatrix", 
      name = "IterativeLSI", 
      iterations = 2, 
      clusterParams = list( #See Seurat::FindClusters
          resolution = c(0.2), 
          sampleCells = length(proj1$cellNames), 
          n.start = 10
      ), 
      varFeatures = peaknum, 
      dimsToUse = 1:(min(topDim,length(proj1$cellNames))),
      seed=1,force=T,sampleCellsPre=20000,
      sampleCellsFinal = 20000,projectCellsPre=F
  )
  
  # clustering
  if(clusterMethod == "SEURAT"){
    proj1_seurat <- addClusters(
        input = proj1,
        reducedDims = "IterativeLSI",
        method = "Seurat",
        name = "Clusters",
        resolution = 0.8,
        force=T,seed=1
    )
  }else{
    proj1_seurat <- addClusters(
        input = proj1,
        reducedDims = "IterativeLSI",
        method = "scran",
        name = "Clusters",
        k=15,
        force=T,seed=1
    )
  }
  
  clusterInfo <- proj1_seurat$Clusters#as.numeric(Ddata$seurat_clusters)
  names(clusterInfo) <- unlist(lapply(proj1_seurat$cellNames,getCN))
  outdata <- cbind(names(clusterInfo), clusterInfo)
  colnames(outdata) <- c("cellname","cluster")
  write.table(outdata,file=paste0("../",outname,"_scClusters.txt"),row.names=F,col.names=T,sep="\t",quote=F)
  
  if(makeUMAP == 1){
    proj1 <- addUMAP(
        ArchRProj = proj1, 
        reducedDims = "IterativeLSI", 
        name = "UMAP", 
        nNeighbors = 30, 
        minDist = 0.5, 
        metric = "cosine",
        force=T
    )
    umapdata <- proj1@embeddings$UMAP$df
    rownames(umapdata) <- unlist(lapply(proj1$cellNames,getCN))
    if(nrow(outdata)<1000){
      PCH <- 16
    }else{
      PCH <- "."
    }
    pdf(file=paste0("../",outname,"_clusteringUMAP.pdf"))
    layout(matrix(c(1,2),nrow=1),width=c(4,1))
    par(mar=c(4,4,2,1))
    colorMat <- rep("black",nrow(umapdata))
    names(colorMat) <- rownames(umapdata)
    rain <- rainbow(length(sort(unique(clusterInfo))))
    for(i in seq(sort(unique(clusterInfo)))){
      ii <- sort(unique(clusterInfo))[i]
      colorMat[names(clusterInfo)[which(clusterInfo==ii)]] <- rain[i]
    }
    plot(umapdata[,1],umapdata[,2],col=colorMat[rownames(umapdata)],pch=PCH,xlab="UMAP-1",ylab="UMAP-2",main=paste0(outname," ArchR UMAP"))
    par(mar=c(4,1,2,1))
    plot(1,1,type="n",xlab="",ylab="",main="",axes=F)
    legend("center",legend=sort(unique(clusterInfo)),col=rain,bty="n",pch=16)
    dev.off()
  }
}else{
  simpleError("NoInstall")
}

"""%(outname,GENOME,Percent,topDim,makeUMAP,peakMaxReads,SCcorrection,clusterMethod)
    outf = open(outname+"_scRscript.r",'w')
    outf.write(scRscript)
    outf.close()

    tmplog = sp("Rscript %s_scRscript.r"%outname)    
    if "simpleError" in tmplog[0].decode("ascii") and "noPackage" in tmplog[0].decode("ascii"):
        return("noPackage")
    else:
        return("yesPackage")


def scClustering_APEC(outname,Percent,makeUMAP,peakMaxReads,SCcorrection):
    try:
        from APEC import clustering,plot,generate
    except:
        return("noPackage")

    Percent = int(Percent)
    plotUMAP = int(makeUMAP)
    peakMaxReads = int(peakMaxReads)

    if not os.path.isdir(outname+"_APEC"):
        os.mkdir(outname+"_APEC")
    os.chdir(outname+"_APEC")
    if not os.path.isdir("matrix"):
        os.mkdir("matrix")
    if not os.path.isdir("peak"):
        os.mkdir("peak")
    
    #if SCcorrection == 0:
    regionFeature_inf = open("../%s_peakFeatures.txt"%outname)
    #else:
    #    regionFeature_inf = open("../%s_correctionMatPeaks.txt"%outname)      

    allpeak_list = []
    for line in regionFeature_inf:
        if "avebias" in line:
            continue
        allpeak_list.append( line.split() )
    regionFeature_inf.close()

    if SCcorrection == 0:
        i=6
        sorted_allpeak_list = sorted(allpeak_list,key=lambda x:float(x[i]),reverse=True)
        peaknum = len(sorted_allpeak_list)
        keep_percent = int(Percent)
        keep_num1 = int(round(peaknum* (1-keep_percent/100)))-1
        usepeak_list = sorted_allpeak_list[keep_num1:]
        usepeakname_dict = {}
        for i in usepeak_list:
            if int(i[5]) <= peakMaxReads:
                usepeakname_dict[i[3]] = i
        inf = open("../%s_peakXcellMat.txt"%outname)

    else:
        usepeakname_dict = {}
        for i in allpeak_list:
            if int(i[5]) <= peakMaxReads :
                usepeakname_dict[i[3]] = i
        inf = open("../%s_correctionMat.txt"%outname)


    colnames = inf.readline()
    cell_outf = open("matrix/filtered_cells.csv",'w')
    cell_outf.write("\tnotes\n")
    cell_count = 0
    for cellname in colnames.split():
        if not cellname == "peakname":
            cell_count += 1
            cell_outf.write("%s\tsc\n"%cellname)
    cell_outf.close()
    
    mtx_out = open('matrix/filtered_reads.mtx.material','w')
    peak_out = open("peak/top_filtered_peaks.bed",'w')
    peak_idx = 0
    total_count = 0
    for line in inf:
        ll = line.split()
        this_peakname = ll[0]
        this_values = ll[1:]
        if this_peakname in usepeakname_dict:
            peak_idx += 1
            peak_out.write("\t".join(usepeakname_dict[this_peakname][:4])+"\n")
            for cell_idx_raw in range(len(this_values)):
                this_value = int(float(this_values[cell_idx_raw]))
                if this_value > 0:
                    cell_idx = cell_idx_raw + 1
                    total_count += 1
                    newll = [peak_idx,cell_idx,this_value]
                    mtx_out.write(" ".join(map(str,newll))+"\n")
    inf.close()
    mtx_out.close()
    peak_out.close()
    
    mtx_header_out = open('matrix/filtered_reads.mtx.header','w')
    mtx_header_out.write('%%MatrixMarket matrix coordinate integer general\n')
    mtx_header_out.write('%s %s %s\n'%(peak_idx,cell_count,total_count))
    mtx_header_out.close()
    
    inf = open("matrix/filtered_reads.mtx.header")
    inf2 = open("matrix/filtered_reads.mtx.material")
    outf = open("matrix/filtered_reads.mtx",'w')
    for line in inf:
        outf.write(line)
    for line in inf2:
        outf.write(line)
    inf.close()
    inf2.close()
    outf.close()
    
    os.chdir("../")
    
    APEC_project=outname+"_APEC"
    clustering.build_accesson(APEC_project, ngroup=600)
    clustering.cluster_byAccesson(APEC_project, nc=0, norm='probability')
    
    inf = open(outname+"_APEC/result/cluster_by_APEC.csv")
    outf = open(outname+"_scClusters.txt",'w')
    outf.write("\t".join(["cellname","cluster"])+"\n")
    for line in inf:
        ll = line.strip().split("\t")
        if ll == ['cluster']:
            continue
        outf.write(line)
    outf.close()
    inf.close()
    
    if plotUMAP==1:
        plot.plot_tsne(APEC_project, rs=0)
        tmplog = sp("cp %s %s"%(outname+"_APEC/figure/TSNE_by_APEC_with_cluster_label.pdf", outname+"_clusteringUMAP.pdf"))
    
    return("yesPackage")




def scClustering_ArchR_scran(outname,GENOME,Percent,topDim,makeUMAP):
    cmd1 = """awk '{OFS="\t";print $1,$2,$3,$4,1}' %s_highQcellReads.bed | sort -k 1,1 -k 2,2g -k 3,3g -k 4,4  | bgzip > %s_ArchRReads.bed.gz"""%(outname,outname)
    cmd2 = """tabix -p bed %s_ArchRReads.bed.gz"""%(outname)
    check_bgzip = sp("which bgzip")
    check_tabix = sp("which tabix")
    if check_bgzip[0].decode("ascii") == "" or check_tabix[0].decode("ascii") == "" :
        return("noPackage")

    tmplog = sp(cmd1)
    tmplog = sp(cmd2)

    if os.path.isfile("%s_ArchRReads.bed.gz"%(outname)) and os.path.isfile("%s_ArchRReads.bed.gz.tbi"%(outname)):
        pass
    else:
        return("noTabix")

    scRscript="""
if(require("ArchR")){
  set.seed(1)
  
  outname <- "%s"
  GENOME <- "%s"
  Percent <- %s
  topDim <- %s
  makeUMAP <- %s
  
  dir.create(paste0(outname,"_ArchR"))
  setwd(paste0(outname,"_ArchR"))
  set.seed(1)
  addArchRGenome(GENOME)
  addArchRThreads(threads = 1)
  
  getCN <- function(inname){
      return(strsplit(inname,"#")[[1]][2])
  }
  
  inputFiles <- c(paste0("../",outname,"_ArchRReads.bed.gz"))
  names(inputFiles)<-c("combine")
  
  ArrowFiles <- createArrowFiles(
    inputFiles = inputFiles,
    sampleNames = names(inputFiles),
    minTSS = 0, #Dont set this too high because you can always increase later
    minFrags = 0, maxFrags=1e+10,
    addTileMat = F,
    force=T,
    addGeneScoreMat = F
  )
  
  proj1 <- ArchRProject(
    ArrowFiles = ArrowFiles,
    outputDirectory = outname,
    copyArrows = FALSE 
  )
  
  regionFeature <- read.table(paste0("../",outname,"_peakFeatures.txt"),row.names=4,header=T)
  peaknum <- nrow(regionFeature)
  keep_percent <- as.numeric(Percent)
  usepeak <- regionFeature[ order(regionFeature[,"avebias"],decreasing=T)[ round(peaknum* (1-keep_percent/100) ) : peaknum ], ]
  peakdata <- usepeak
  
  peakGR <- GRanges(seqnames=peakdata[,1],ranges=IRanges(peakdata[,2],peakdata[,3]))
  proj1 <- addPeakSet(ArchRProj=proj1, peakSet=peakGR)
  proj1 <- addPeakMatrix(proj1)
  
  proj30 <- addIterativeLSI(
      ArchRProj = proj1,
      useMatrix = "PeakMatrix", 
      name = "IterativeLSI", 
      iterations = 2, 
      clusterParams = list( #See Seurat::FindClusters
          resolution = c(0.2), 
          sampleCells = length(proj1$cellNames), 
          n.start = 10
      ), 
      varFeatures = peaknum, 
      dimsToUse = 1:(min(topDim,length(proj1$cellNames))),
      seed=1,force=T,sampleCellsPre=20000,
      sampleCellsFinal = 20000,projectCellsPre=F
  )
  
  # clustering
  proj30_seurat <- addClusters(
      input = proj30,
      reducedDims = "IterativeLSI",
      method = "scran",
      name = "Clusters",
      k=15,
      force=T,seed=1
  )
  
  clusterInfo <- proj30_seurat$Clusters#as.numeric(Ddata$seurat_clusters)
  names(clusterInfo) <- unlist(lapply(proj30_seurat$cellNames,getCN))
  outdata <- cbind(names(clusterInfo), clusterInfo)
  colnames(outdata) <- c("cellname","cluster")
  write.table(outdata,file=paste0("../",outname,"_scClusters.txt"),row.names=F,col.names=T,sep="\t",quote=F)
  
  if(makeUMAP == 1){
    proj30 <- addUMAP(
        ArchRProj = proj30, 
        reducedDims = "IterativeLSI", 
        name = "UMAP", 
        nNeighbors = 30, 
        minDist = 0.5, 
        metric = "cosine",
        force=T
    )
    umapdata <- proj30@embeddings$UMAP$df
    rownames(umapdata) <- unlist(lapply(proj30$cellNames,getCN))
    if(nrow(outdata)<1000){
      PCH <- 16
    }else{
      PCH <- "."
    }
    pdf(file=paste0("../",outname,"_clusteringUMAP.pdf"))
    layout(matrix(c(1,2),nrow=1),width=c(4,1))
    par(mar=c(4,4,2,1))
    colorMat <- rep("black",nrow(umapdata))
    names(colorMat) <- rownames(umapdata)
    rain <- rainbow(length(sort(unique(clusterInfo))))
    for(i in seq(sort(unique(clusterInfo)))){
      ii <- sort(unique(clusterInfo))[i]
      colorMat[names(clusterInfo)[which(clusterInfo==ii)]] <- rain[i]
    }
    plot(umapdata[,1],umapdata[,2],col=colorMat[rownames(umapdata)],pch=PCH,xlab="UMAP-1",ylab="UMAP-2",main=paste0(outname," ArchR UMAP"))
    par(mar=c(4,1,2,1))
    plot(1,1,type="n",xlab="",ylab="",main="",axes=F)
    legend("center",legend=sort(unique(clusterInfo)),col=rain,bty="n",pch=16)
    dev.off()
  }
}else{
  simpleError("NoInstall")
}

"""%(outname,GENOME,Percent,topDim,makeUMAP)
    outf = open(outname+"_scRscript.r",'w')
    outf.write(scRscript)
    outf.close()
    tmplog = sp("Rscript %s_scRscript.r"%outname)    
    if "simpleError" in tmplog[0].decode("ascii") and "noPackage" in tmplog[0].decode("ascii"):
        return("noPackage")
    else:
        return("yesPackage")


def scClustering_Seurat(outname,Percent,topDim,makeUMAP):
    cmd1 = """awk '{OFS="\t";print $1,$2,$3,$4,1}' %s_highQcellReads.bed | sort -k 1,1 -k 2,2g -k 3,3g -k 4,4  | bgzip > %s_SeuratReads.bed.gz"""%(outname,outname)
    cmd2 = """tabix -p bed %s_SeuratReads.bed.gz"""%(outname)
    check_bgzip = sp("which bgzip")
    check_tabix = sp("which tabix")
    if check_bgzip[0].decode("ascii") == "" or check_tabix[0].decode("ascii") == "" :
        return("noPackage")

    tmplog = sp(cmd1)
    tmplog = sp(cmd2)

    scRscript="""
if(require("Signac") & require("Seurat") & require("Matrix")){
  set.seed(1)
  
  outname <- "%s"
  Percent <- %s
  topDim <- %s
  makeUMAP <- %s
  
  data <- read.table(paste0(outname,"_peakXcellMat.txt"),row.names=1,header=T)
  regionFeature <- read.table(paste0(outname,"_peakFeatures.txt"),row.names=4,header=T)
  peaknum <- nrow(regionFeature)
  keep_percent <- as.numeric(Percent)
  usepeak <- regionFeature[ order(regionFeature[,"avebias"],decreasing=T)[ round(peaknum* (1-keep_percent/100) ) : peaknum ], ]
  
  useCell <- as.matrix(colnames(data))
  rownames(useCell)<- useCell[,1]
  useCell[,1] <- rep(1,length(useCell))
  useCell <- as.data.frame(useCell)
  colnames(useCell) <- "is__cell_barcode"
  
  usedata0 <- data[rownames(usepeak),]
  gbm <- as.matrix(usedata0)
  sparse.gbm <- Matrix(gbm , sparse = T )
  
  counts <- sparse.gbm
  oldrownames <- rownames(counts)
  newrownames <- paste0(regionFeature[oldrownames,1],":",regionFeature[oldrownames,2]+199,"-",regionFeature[oldrownames,2]+200)
  rownames(counts) <- newrownames
  METADATA <- useCell
  
  chrom_assay <- CreateChromatinAssay(
    counts = counts,
    sep = c(":", "-"),
    fragments = paste0(outname,"_SeuratReads.bed.gz"),
    min.cells = 0,
    min.features = 0
  )
  
  pbmc <- CreateSeuratObject(
    counts = chrom_assay,
    assay = "peaks",
    meta.data = useCell
  )
  
  pbmc <- RunTFIDF(pbmc)
  pbmc <- FindTopFeatures(pbmc, min.cutoff = 'q0')
  pbmc <- RunSVD(pbmc,n=min(topDim,ncol(usedata0)))
  
  Ddata <- FindNeighbors(object = pbmc, reduction = 'lsi', dims = 1:(min(topDim,ncol(usedata0))))
  Ddata <- FindClusters(object = Ddata, verbose = FALSE, algorithm = 3)
  
  clusterInfo <- as.numeric(Ddata$seurat_clusters)
  names(clusterInfo) <- names(Ddata$seurat_clusters)
  
  outdata <- cbind(names(clusterInfo), clusterInfo)
  colnames(outdata) <- c("cellname","cluster")
  write.table(outdata,file=paste0(outname,"_scClusters.txt"),row.names=F,col.names=T,sep="\t",quote=F)
  
  if(makeUMAP == 1){
    pbmc <- RunUMAP(object = pbmc, reduction = 'lsi', dims = 1:(min(topDim,ncol(usedata0))))
    umapdata <-pbmc@reductions$umap@cell.embeddings
  
    if(nrow(outdata)<1000){
      PCH <- 16
    }else{
      PCH <- "."
    }
  
    pdf(file=paste0(outname,"_clusteringUMAP.pdf"))
    layout(matrix(c(1,2),nrow=1),width=c(4,1))
    par(mar=c(4,4,2,1))
      colorMat <- rep("black",nrow(umapdata))
      names(colorMat) <- rownames(umapdata)
      rain <- rainbow(length(sort(unique(clusterInfo))))
      for(i in sort(unique(clusterInfo))){
        colorMat[names(clusterInfo)[which(clusterInfo==i)]] <- rain[i]
      }
      plot(umapdata[,1],umapdata[,2],col=colorMat[rownames(umapdata)],pch=PCH,xlab="UMAP-1",ylab="UMAP-2",main=paste0(outname," Seurat UMAP"))
    par(mar=c(4,1,2,1))
    plot(1,1,type="n",xlab="",ylab="",main="",axes=F)
    legend("center",legend=sort(unique(clusterInfo)),col=rain,bty="n",pch=16)
    dev.off()
  }
}else{
  simpleError("NoInstall")
}
"""%(outname,Percent,topDim,makeUMAP)
    outf = open(outname+"_scRscript.r",'w')
    outf.write(scRscript)
    outf.close()
    tmplog = sp("Rscript %s_scRscript.r"%outname)    
    if "simpleError" in tmplog[0].decode("ascii") and "noPackage" in tmplog[0].decode("ascii"):
        return("noPackage")
    else:
        return("yesPackage")

def scClustering_Cicero(outname,Percent,makeUMAP):

    Percent = int(Percent)

    if not os.path.isdir(outname+"_Cicero"):
        os.mkdir(outname+"_Cicero")
    os.chdir(outname+"_Cicero")
    if not os.path.isdir("matrix"):
        os.mkdir("matrix")
    if not os.path.isdir("peak"):
        os.mkdir("peak")
    
    regionFeature_inf = open("../%s_peakFeatures.txt"%outname)
    allpeak_list = []
    for line in regionFeature_inf:
        if "avebias" in line:
            continue
        allpeak_list.append( line.split() )
    regionFeature_inf.close()
    
    i=6
    sorted_allpeak_list = sorted(allpeak_list,key=lambda x:float(x[i]),reverse=True)
    peaknum = len(sorted_allpeak_list)
    keep_percent = int(Percent)
    keep_num1 = int(round(peaknum* (1-keep_percent/100)))-1
    usepeak_list = sorted_allpeak_list[keep_num1:]
    
    usepeakname_dict = {}
    for i in usepeak_list:
        usepeakname_dict[i[3]] = i
    
    inf = open("../%s_peakXcellMat.txt"%outname)
    colnames = inf.readline()
    cell_outf = open("matrix/filtered_cells.csv",'w')
    cell_outf.write("\tnotes\n")
    cell_count = 0
    for cellname in colnames.split():
        if not cellname == "peakname":
            cell_count += 1
            cell_outf.write("%s\tsc\n"%cellname)
    cell_outf.close()
    
    mtx_out = open('matrix/filtered_reads.mtx.material','w')
    peak_out = open("peak/top_filtered_peaks.bed",'w')
    peak_idx = 0
    total_count = 0
    for line in inf:
        ll = line.split()
        this_peakname = ll[0]
        this_values = ll[1:]
        if this_peakname in usepeakname_dict:
            peak_idx += 1
            peak_out.write("\t".join(usepeakname_dict[this_peakname][:4])+"\n")
            for cell_idx_raw in range(len(this_values)):
                this_value = int(this_values[cell_idx_raw])
                if this_value > 0:
                    cell_idx = cell_idx_raw + 1
                    total_count += 1
                    newll = [peak_idx,cell_idx,this_value]
                    mtx_out.write(" ".join(map(str,newll))+"\n")
    inf.close()
    mtx_out.close()
    peak_out.close()
    
    mtx_header_out = open('matrix/filtered_reads.mtx.header','w')
    mtx_header_out.write('%%MatrixMarket matrix coordinate integer general\n')
    mtx_header_out.write('%s %s %s\n'%(peak_idx,cell_count,total_count))
    mtx_header_out.close()
    
    inf = open("matrix/filtered_reads.mtx.header")
    inf2 = open("matrix/filtered_reads.mtx.material")
    outf = open("matrix/filtered_reads.mtx",'w')
    for line in inf:
        outf.write(line)
    for line in inf2:
        outf.write(line)
    inf.close()
    inf2.close()
    outf.close()
    
    os.chdir("../")
    scRscript="""
if(require("cicero") & require("Matrix")){

  set.seed(1)
  outname <- "%s"
  Percent <- %s
  makeUMAP <- %s
  
  # read in matrix data using the Matrix package
  indata <- Matrix::readMM(paste0(outname,"_Cicero/matrix/filtered_reads.mtx"))
  indata@x[indata@x > 0] <- 1
  
  # format cell info
  cellinfo <- read.table(paste0(outname,"_Cicero/matrix/filtered_cells.csv"),header=T,row.names=1)
  cellinfo[,1] <- rownames(cellinfo)
  names(cellinfo) <- "cells"
  
  # format peak info
  peakinfo <- read.table(paste0(outname,"_Cicero/peak/top_filtered_peaks.bed"))
  names(peakinfo) <- c("chr", "bp1", "bp2","peakname")
  peakinfo$site_name <- paste(peakinfo$chr, peakinfo$bp1, peakinfo$bp2, sep="_")
  row.names(peakinfo) <- peakinfo$site_name
  row.names(indata) <- row.names(peakinfo)
  colnames(indata) <- row.names(cellinfo)
  
  # make CDS
  fd <- methods::new("AnnotatedDataFrame", data = peakinfo)
  pd <- methods::new("AnnotatedDataFrame", data = cellinfo)
  input_cds <-  suppressWarnings(newCellDataSet(indata,
                              phenoData = pd,
                              featureData = fd,
                              expressionFamily=VGAM::binomialff(),
                              lowerDetectionLimit=0))
  input_cds@expressionFamily@vfamily <- "binomialff"
  input_cds <- monocle::detectGenes(input_cds)
  
  #Ensure there are no peaks included with zero reads
  input_cds <- input_cds[Matrix::rowSums(exprs(input_cds)) != 0,] 
  
  agg_cds <- aggregate_nearby_peaks(input_cds, distance = 10000)
  agg_cds <- detectGenes(agg_cds)
  agg_cds <- estimateSizeFactors(agg_cds)
  agg_cds <- estimateDispersions(agg_cds)
  
  agg_cds <- reduceDimension(agg_cds,
                              max_components = 2,
                              norm_method = 'log',
                              num_dim = 3,
                              reduction_method = 'tSNE',
                              verbose = T)
  
  agg_cds <- clusterCells(agg_cds, verbose = F)
  cellname <- as.vector(agg_cds$cells)
  outdata <- cbind(cellname,as.numeric(agg_cds$Cluster))
  colnames(outdata) <- c("cellname","cluster")
  write.table(outdata,file=paste0(outname,"_scClusters.txt"),row.names=F,col.names=T,sep="\t",quote=F)
  
  if(makeUMAP == 1){
    pdf(file=paste0(outname,"_clusteringUMAP.pdf"))
    plot_cell_clusters(agg_cds, color_by = 'as.factor(Cluster)') + theme(text = element_text(size=8))
    dev.off()
  }
}else{
  simpleError("NoInstall")
}
"""%(outname,Percent,makeUMAP)
    outf = open(outname+"_scRscript.r",'w')
    outf.write(scRscript)
    outf.close()
    tmplog = sp("Rscript %s_scRscript.r"%outname)    
    if "simpleError" in tmplog[0].decode("ascii") and "noPackage" in tmplog[0].decode("ascii"):
        return("noPackage")
    else:
        return("yesPackage")


