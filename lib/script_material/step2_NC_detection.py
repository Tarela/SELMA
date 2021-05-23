#!/usr/bin/env python

# ------------------------------------
# Python Modual
# ------------------------------------

import os
import sys
import string
import time

# --------------------------
# custom package
# --------------------------

### tool function
from HMRpipe.Utility      import (sp,
                                   pdf_name,
                                   raise_error,
                                   wlog,
                                   ewlog,
                                   rwlog,
                                   CMD,
                                   createDIR)
# --------------------------
# main 
# --------------------------
def step2_NC_detection(conf_dict,logfile):
    '''
    analysis part, mainly Rscript
    Detect 
    '''   

    # Rscript detectNonClassic.r outname signalname usePQ cutoff alpha lambdachoice topN tmpRpackgeDIR
    if conf_dict['General']['mode'] == "binary":
        Rscript = """
### read parameter
a<-commandArgs(T)

outname <- a[1]
signalname <- unlist(strsplit(a[2],","))
cutoff <- as.numeric(a[3])
Alpha <- as.numeric(a[4])
LambdaChoice <- (a[5])
topN <- a[6]
R2cutoff <- as.numeric(a[8])
R2classic <- 0.1
tmp_rpackage_dir <- a[7]

### install R packages
if("foreach" %in% installed.packages()[,"Package"]){
    library(foreach)
}else{
    install.packages("foreach",dependencies=TRUE,lib=tmp_rpackage_dir,repos="https://mirrors.tongji.edu.cn/CRAN/")
    library(foreach,lib.loc=tmp_rpackage_dir)
}

if("glmnet" %in% installed.packages()[,"Package"]){
    library(glmnet)
}else{
    install.packages("glmnet",dependencies=TRUE,lib=tmp_rpackage_dir,repos="https://mirrors.tongji.edu.cn/CRAN/")
    library(glmnet,lib.loc=tmp_rpackage_dir)
}

library(methods)

### define functions
unilinear_only <- function(TFov,HMsig){
    # univariate linear regression 
    lmresult <- summary(lm(HMsig ~ TFov))
    fgR2 <- lmresult$adj.r.squared
    lmcoeff <- lmresult$coefficients['TFov','Estimate']
    lmP <- lmresult$coefficients['TFov',c('Pr(>|t|)')]
    return(c(fgR2,lmcoeff))
}

unilinear_permute <- function(TFov,HMsig){
    # univariate linear regression + empirical pvalue
    lmresult <- summary(lm(HMsig ~ TFov))
    fgR2 <- lmresult$adj.r.squared
    bgR2s <- c()
    for(i in 1:999){
        bgTFov <- rep(0,length(TFov))
        bgTFov[sample(length(TFov),length(which(TFov > 0)))] <- 1
        bgR2 <- summary(lm(HMsig ~ bgTFov))$adj.r.squared
        bgR2s <- c(bgR2s,bgR2)
    }
    permuteP <- (length(which(fgR2 < bgR2s))+1)/(length(bgR2s)+1)
    return(permuteP)
}

maxG<-function(cutoff,usedata){
    # function of estimating G given cutoff
    # G : inter-class variance
    # method refer to "Otus' method"
    w0 <- length(which(usedata < cutoff))/length(usedata)
    w1 <- length(which(usedata >= cutoff))/length(usedata)
    u0 <- mean(usedata[which(usedata < cutoff)])
    u1 <- mean(usedata[which(usedata >= cutoff)])
    g <- w0*w1*(u0-u1)**2
    return(g)
}

signal2cutoff <- function(rawsig){
    # function for separate NC peaks considering HMsignal
    # input data: a vector of signal (HMsignal on HMR peak, linear scale)
    # method: go through all possible cutoff, find the cutoff corresponding to maximum G
    # output: the cutoff, a vector of selected cutoff candidates (Gbins) and a vector of G value for each cutoff candidates (G) 
    
    sig <- log10(rawsig[which(rawsig>0)])
    # separate the section of (log) signal to N cutoffs 
    Gbins <- seq(min(sig)+0.1,max(sig)-0.1,0.01)
    
    # estimate G for each cutoff candidates
    G <- unlist(lapply(Gbins,maxG, sig))
    
    # select cutoff at the first time G meats its maximum value
    NCcut<-seq(min(sig)+0.1,max(sig)-0.1,0.01)[which(G==max(G))][1]

    group_detail <- rep(0, length(rawsig))
    group_detail[which(rawsig >= 10**NCcut)] <- 1
    
    # output the grouping result: list contains 3 items 
    # item1: NC cutoff
    # item2: 2 column for cutoff candidates and corresponded G
    # item3: Nrow = peak number,Ncolumn = 2, c1 for signal, c2 for group number, 0 for lowHM group (solo, non-classical), 1 for highHM group (ensemble, classical)
    
    return(list(NCcut, cbind(Gbins,G), cbind(rawsig, group_detail) ))
}

# step1 data preprocess
peakov <- read.table(paste0(outname,"_peakov.bed"))
signal <- read.table(paste0(outname,"_HMsig.bed"))

candidate_list <- as.vector(read.table(paste0(outname,"_cofactor_candidate_list.txt"))[,1])
bedncol <- ncol(peakov) - length(candidate_list)
colnames(peakov) <- c(paste0("c",seq(1:bedncol)),candidate_list)
colnames(signal) <- c(paste0("c",seq(1:bedncol)),signalname)
rownames(peakov) <- paste0("r",seq(1,nrow(peakov)))
rownames(signal) <- paste0("r",seq(1,nrow(peakov)))

TFov <- peakov[,candidate_list]
TFov[TFov > 1] <- 1
HMsig <- as.matrix(signal[,signalname])
colnames(HMsig) <- signalname
rownames(HMsig) <- rownames(signal)


### step2, prepare X,Y for model selection 
lindata <- cbind(HMsig,TFov)
bind_sum <- apply(TFov,1,sum)
## sites with 90% factors overlapped are excluded
use1_lindata <- lindata[names(bind_sum[which(bind_sum < ncol(TFov)*0.9)]),]
## factors with < 100 or > 95% cobinding events are excluded
TF_sum <- apply(use1_lindata[,colnames(TFov)],2,sum)
use_lindata <- lindata[names(bind_sum[which(bind_sum < ncol(TFov)*0.9)]),c(colnames(HMsig),names(TF_sum)[which(TF_sum>=100 & TF_sum <= nrow(use1_lindata)*0.95)])]
## form X, Y
### raw Y is used in otsu' method
rawY <- as.matrix(use_lindata[,colnames(HMsig)])
colnames(rawY) <- colnames(HMsig)
Y <- as.matrix(scale(use_lindata[,colnames(HMsig)]))
colnames(Y) <- colnames(HMsig)
X <- as.matrix(use_lindata[,(ncol(HMsig)+1):ncol(use_lindata)])
peakX <- peakov[rownames(X),1:bedncol]


### elastic net model selection
set.seed(1007)
coTFusage <- matrix(rep(0,ncol(X)*ncol(Y)),nrow=ncol(X))
rownames(coTFusage) <- colnames(X)
colnames(coTFusage) <- colnames(Y)

if(ncol(Y) == 1){
    pdf(file=paste0(outname,"_elnet_lambdaSelection.pdf"),height=6,width=6)
    par(mar=c(4,4,4,2))    
}else if(ncol(Y) == 2){
    pdf(file=paste0(outname,"_elnet_lambdaSelection.pdf"),height=6,width=12)
    par(mfrow=c(1,2),mar=c(4,4,4,2))
}else if(ncol(Y) == 3){
    pdf(file=paste0(outname,"_elnet_lambdaSelection.pdf"),height=12,width=12)
    par(mfrow=c(2,2),mar=c(4,4,4,2))
}else{
    pdf(file=paste0(outname,"_elnet_lambdaSelection.pdf"),height=12,width=12)
    par(mfrow=c(2,2),mar=c(4,4,4,2))
}

for(i in 1:ncol(Y)){
    thisY <- Y[,i]
    this_name <- colnames(Y)[i]
    cv.glmmod <- cv.glmnet(x=X,y=thisY,alpha=Alpha,family="gaussian")
    coeff_this <- coef(cv.glmmod,s=paste0("lambda.",LambdaChoice)) 
    coTF_name <- rownames(coeff_this)[which(coeff_this[,1] < 0 & rownames(coeff_this)!= "(Intercept)")]

    coTFusage[coTF_name,this_name] <- 1

    ### cross-validation curve for each HM
    if (i %in% 1:4){
        plot(cv.glmmod)
        title(this_name,line=2.5)
        if (LambdaChoice == "1se"){
            abline(v=log(cv.glmmod$lambda.1se),col="blue",lwd=2)
            legend("topleft",legend="lambda.1se",lwd=3,bty="n",col='blue')
        }else{
            abline(v=log(cv.glmmod$lambda.min),col="blue",lwd=2)
            legend("topleft",legend="lambda.min",lwd=3,bty="n",col='blue')
        }
    }
}

dev.off()

### generate R2 and coeff table
R2_mat <- matrix(rep(0,ncol(X)*ncol(Y)),nrow=ncol(X))
coeff_mat <- matrix(rep(0,ncol(X)*ncol(Y)),nrow=ncol(X))
rownames(R2_mat) <- colnames(X)
colnames(R2_mat) <- colnames(Y)
rownames(coeff_mat) <- colnames(X)
colnames(coeff_mat) <- colnames(Y)
for(TFname in colnames(X)){
    for(HMname in colnames(Y)){
        unilinResult <- unilinear_only(X[,TFname],Y[,HMname])
        R2_mat[TFname,HMname] <- unilinResult[1]
        coeff_mat[TFname,HMname] <- unilinResult[2]
    }
}

summary_table <- c()
### select candidates
# for each candidates detected in el-net step, require R2>0.1, coeff<0, R2 for any positive correlated HM signal (coeff>0) < 0.1 
for(TFname in colnames(X)){
    this_usage <- coTFusage[TFname,]
    this_R2 <- R2_mat[TFname,]
    this_coeff <- coeff_mat[TFname,]

    classicFun <- 0
    nonClassicHMidx <- c()

    for(i in 1:ncol(coTFusage)){
        if(this_usage[i]==1 & this_R2[i] >= R2cutoff & this_coeff[i] < 0){
            nonClassicHMidx <- c(nonClassicHMidx, i)
        }else if(this_R2[i] >= R2classic & this_coeff[i] > 0 ){
            classicFun <- 1
        }
    }

    if(classicFun == 0){
        for(idx in nonClassicHMidx){
            nonClassicHM <- colnames(Y)[idx]
            R2 <- this_R2[idx]
            coeff <- this_coeff[idx]
            Pval <- unilinear_permute(X[,TFname],Y[,idx])
            summary_table <- rbind(summary_table, c(TFname, nonClassicHM, Pval, R2, coeff))
        }        
    }
}

### summary and output steps
if(is.null(summary_table)){
    print("no significant candidates detected")
    write.table("no non-classical function detected",file=paste0(outname,"_filterNC.txt"),quote=F,sep="\t",row.names=F,col.names=F)
    write.table("no non-classical function detected",file=paste0(outname,'_NCsummary.txt'),quote=F,sep="\t",row.names=F,col.names=F)
}else{
    colnames(summary_table) <- c("TFname","HMname","Pval","R2","coeff")
    if(nrow(summary_table) > 1){
        summary_table <- summary_table[order(as.numeric(summary_table[,"R2"]),decreasing=TRUE),]
    }

    if(topN == "all"){
        topN <- nrow(summary_table)
    }else{
        topN <- as.numeric(topN)
    }
    
    if(topN == 1){
        out_table_raw <- summary_table
    }else{
        out_table_raw <- as.matrix(summary_table[1:topN,])
    }
    
    ## separate NC/C peaks based on HMsignal using otsu's method
    peakgroup <- c()
    for(i in 1:ncol(rawY)){
        thisY <- rawY[,i]
        thisY_NCotsu <- signal2cutoff(thisY)
        thisY_peakgroup <- thisY_NCotsu[[3]][,2]
        peakgroup <- cbind(peakgroup, thisY_peakgroup)
    }
    colnames(peakgroup) <- colnames(Y)

    ## for each predicted coTFvsHM pair, output the cobinding NC sites as a bed file
    if(!file.exists("nonClassicalPeaks/")){dir.create("nonClassicalPeaks")}
    if(topN == 1){
        pdf(file=paste0(outname,"_cofactor_HMsignal.pdf"),height=6,width=6)
        par(mar=c(4,4,2,2))    
    }else if(topN == 2){
        pdf(file=paste0(outname,"_cofactor_HMsignal.pdf"),height=6,width=12)
        par(mfrow=c(1,2),mar=c(4,4,2,2))
    }else{
        pdf(file=paste0(outname,"_cofactor_HMsignal.pdf"),height=12,width=12)
        par(mfrow=c(2,2),mar=c(4,4,2,2))
    }
    num_NCsites <- c()
    for(topnum in 1:nrow(out_table_raw)){
        coTF = out_table_raw[topnum,1]
        substrateHM = out_table_raw[topnum,2]
        coTF_cobinding <- X[,coTF]
        HMnc <- peakgroup[,substrateHM]
        cobinding_NC_peak <- peakX[which(coTF_cobinding>0 & HMnc == 0),]
        write.table(cobinding_NC_peak, file=paste0("nonClassicalPeaks/",outname,"_",coTF,"_",substrateHM,"_top",topnum,"nonclassical_peaks.bed"),quote=F,sep="\t",row.names=F,col.names=F)
        num_NCsites <- c(num_NCsites, nrow(cobinding_NC_peak))
        if(topnum %in% 1:4){
            ## boxplot compare the histone modification signal between on classical and non-classical peaks
            boxplot(Y[which(coTF_cobinding>0 & HMnc == 0)],Y[which(coTF_cobinding==0 | HMnc > 0)],
                names=c("non-classical peak","classical peak"),ylab=paste0(substrateHM," signal"),main=paste0(coTF," & ",substrateHM),
                outline=F,cex.main=1)
            legend("topleft",legend=paste0("#NCpeak = ",nrow(cobinding_NC_peak)),bty="n")
        }
    }
    dev.off()
    
    out_table <- cbind(out_table_raw, num_NCsites)
    write.table(summary_table,file=paste0(outname,"_filterNC.txt"),quote=F,sep="\t",row.names=F,col.names=T)
    write.table(out_table,file=paste0(outname,'_NCsummary.txt'),quote=F,sep="\t",row.names=F,col.names=T)
} 
  
"""
# Rscript detectNonClassic.r outname signalname usePQ cutoff alpha lambdachoice topN #tmpRpackgeDIR    
    else:# conf_dict['General']['mode'] == "signal":
        Rscript = """
### read parameter
a<-commandArgs(T)

outname <- a[1]
signalname <- unlist(strsplit(a[2],","))
cutoff <- as.numeric(a[3])
Alpha <- as.numeric(a[4])
LambdaChoice <- (a[5])
topN <- a[6]
R2cutoff <- as.numeric(a[8])
R2classic <- 0.1
tmp_rpackage_dir <- a[7]

### install R packages
if("foreach" %in% installed.packages()[,"Package"]){
    library(foreach)
}else{
    install.packages("foreach",dependencies=TRUE,lib=tmp_rpackage_dir,repos="https://mirrors.tongji.edu.cn/CRAN/")
    library(foreach,lib.loc=tmp_rpackage_dir)
}

if("glmnet" %in% installed.packages()[,"Package"]){
    library(glmnet)
}else{
    install.packages("glmnet",dependencies=TRUE,lib=tmp_rpackage_dir,repos="https://mirrors.tongji.edu.cn/CRAN/")
    library(glmnet,lib.loc=tmp_rpackage_dir)
}

library(methods)

### define functions
unilinear_only <- function(TFov,HMsig){
    # univariate linear regression 
    lmresult <- summary(lm(HMsig ~ TFov))
    fgR2 <- lmresult$adj.r.squared
    lmcoeff <- lmresult$coefficients['TFov','Estimate']
    lmP <- lmresult$coefficients['TFov',c('Pr(>|t|)')]
    return(c(fgR2,lmcoeff))
}

unilinear_permute <- function(TFov,HMsig){
    # univariate linear regression + empirical pvalue
    lmresult <- summary(lm(HMsig ~ TFov))
    fgR2 <- lmresult$adj.r.squared
    bgR2s <- c()
    for(i in 1:999){
        bgTFov <- rep(0,length(TFov))
        bgTFov[sample(length(TFov),length(which(TFov > 0)))] <- 1
        bgR2 <- summary(lm(HMsig ~ bgTFov))$adj.r.squared
        bgR2s <- c(bgR2s,bgR2)
    }
    permuteP <- (length(which(fgR2 < bgR2s))+1)/(length(bgR2s)+1)
    return(permuteP)
}

maxG<-function(cutoff,usedata){
    # function of estimating G given cutoff
    # G : inter-class variance
    # method refer to "Otus' method"
    w0 <- length(which(usedata < cutoff))/length(usedata)
    w1 <- length(which(usedata >= cutoff))/length(usedata)
    u0 <- mean(usedata[which(usedata < cutoff)])
    u1 <- mean(usedata[which(usedata >= cutoff)])
    g <- w0*w1*(u0-u1)**2
    return(g)
}

signal2cutoff <- function(rawsig){
    # function for separate NC peaks considering HMsignal
    # input data: a vector of signal (HMsignal on HMR peak, linear scale)
    # method: go through all possible cutoff, find the cutoff corresponding to maximum G
    # output: the cutoff, a vector of selected cutoff candidates (Gbins) and a vector of G value for each cutoff candidates (G) 
    
    sig <- log10(rawsig[which(rawsig>0)])
    # separate the section of (log) signal to N cutoffs 
    Gbins <- seq(min(sig)+0.1,max(sig)-0.1,0.01)
    
    # estimate G for each cutoff candidates
    G <- unlist(lapply(Gbins,maxG, sig))
    
    # select cutoff at the first time G meats its maximum value
    NCcut<-seq(min(sig)+0.1,max(sig)-0.1,0.01)[which(G==max(G))][1]

    group_detail <- rep(0, length(rawsig))
    group_detail[which(rawsig >= 10**NCcut)] <- 1
    
    # output the grouping result: list contains 3 items 
    # item1: NC cutoff
    # item2: 2 column for cutoff candidates and corresponded G
    # item3: Nrow = peak number,Ncolumn = 2, c1 for signal, c2 for group number, 0 for lowHM group (solo, non-classical), 1 for highHM group (ensemble, classical)
    
    return(list(NCcut, cbind(Gbins,G), cbind(rawsig, group_detail) ))
}

trim95 <- function(INdata){
	indata <- INdata
	indata[indata>quantile(indata,0.95)] <- quantile(indata,0.95)
	return(indata)
}
# step1 data preprocess
peakov <- read.table(paste0(outname,"_peakov.bed"))
peaksig <- read.table(paste0(outname,"_TFsig.bed"))
signal <- read.table(paste0(outname,"_HMsig.bed"))

candidate_list <- as.vector(read.table(paste0(outname,"_cofactor_candidate_list.txt"))[,1])
bedncol <- ncol(peakov) - length(candidate_list)
colnames(peakov) <- c(paste0("c",seq(1:bedncol)),candidate_list)
colnames(peaksig) <- c(paste0("c",seq(1:bedncol)),candidate_list)
colnames(signal) <- c(paste0("c",seq(1:bedncol)),signalname)

rownames(peakov) <- paste0("r",seq(1,nrow(peakov)))
rownames(peaksig) <- paste0("r",seq(1,nrow(peaksig)))
rownames(signal) <- paste0("r",seq(1,nrow(peakov)))

use_candidate_list <- c()
for(i in candidate_list){
	SDsig <- sd(as.numeric(peaksig[,i]))
    if (i != outname || SDsig == 0){
        use_candidate_list <- c(use_candidate_list,i)
    }
}

TFov <- as.matrix(peakov[,use_candidate_list])
TFov[TFov > 1] <- 1
TFsig_raw <- as.matrix(peakov[,use_candidate_list])
TFsig <- apply(TFsig_raw,2,trim95)

HMsig <- as.matrix(signal[,signalname])
colnames(HMsig) <- signalname
rownames(HMsig) <- rownames(signal)


### step2, prepare X,Y for model selection 
lindata <- cbind(HMsig,TFov)
bind_sum <- apply(TFov,1,sum)
## sites with 90% factors overlapped are excluded
use1_lindata <- lindata[names(bind_sum[which(bind_sum < ncol(TFov)*0.9)]),]
## factors with < 100 or > 95% cobinding events are excluded
TF_sum <- apply(use1_lindata[,colnames(TFov)],2,sum)
use_lindata <- lindata[names(bind_sum[which(bind_sum < ncol(TFov)*0.9)]),c(colnames(HMsig),names(TF_sum)[which(TF_sum>=100 & TF_sum <= nrow(lindata)*0.95)])]
## form X, Y
### raw Y is used in otsu' method
rawY <- as.matrix(use_lindata[,colnames(HMsig)])
colnames(rawY) <- colnames(HMsig)
Y <- as.matrix(scale(use_lindata[,colnames(HMsig)]))
colnames(Y) <- colnames(HMsig)
X_peakOV <- as.matrix(use_lindata[,(ncol(HMsig)+1):ncol(use_lindata)])
peakX <- peakov[rownames(X_peakOV),1:bedncol]
X_sig <- as.matrix(TFsig[rownames(X_peakOV),colnames(X_peakOV)])
X_sig[X_peakOV == 0] <- 0
X <- X_sig

### elastic net model selection
set.seed(1007)
coTFusage <- matrix(rep(0,ncol(X)*ncol(Y)),nrow=ncol(X))
rownames(coTFusage) <- colnames(X)
colnames(coTFusage) <- colnames(Y)

if(ncol(Y) == 1){
    pdf(file=paste0(outname,"_elnet_lambdaSelection.pdf"),height=6,width=6)
    par(mar=c(4,4,4,2))    
}else if(ncol(Y) == 2){
    pdf(file=paste0(outname,"_elnet_lambdaSelection.pdf"),height=6,width=12)
    par(mfrow=c(1,2),mar=c(4,4,4,2))
}else if(ncol(Y) == 3){
    pdf(file=paste0(outname,"_elnet_lambdaSelection.pdf"),height=12,width=12)
    par(mfrow=c(2,2),mar=c(4,4,4,2))
}else{
    pdf(file=paste0(outname,"_elnet_lambdaSelection.pdf"),height=12,width=12)
    par(mfrow=c(2,2),mar=c(4,4,4,2))
}

for(i in 1:ncol(Y)){
    thisY <- Y[,i]
    this_name <- colnames(Y)[i]
    cv.glmmod <- cv.glmnet(x=X,y=thisY,alpha=Alpha,family="gaussian")
    coeff_this <- coef(cv.glmmod,s=paste0("lambda.",LambdaChoice)) 
    coTF_name <- rownames(coeff_this)[which(coeff_this[,1] < 0 & rownames(coeff_this)!= "(Intercept)")]

    coTFusage[coTF_name,this_name] <- 1

    ### cross-validation curve for each HM
    if (i %in% 1:4){
        plot(cv.glmmod)
        title(this_name,line=2.5)
        if (LambdaChoice == "1se"){
            abline(v=log(cv.glmmod$lambda.1se),col="blue",lwd=2)
            legend("topleft",legend="lambda.1se",lwd=3,bty="n",col='blue')
        }else{
            abline(v=log(cv.glmmod$lambda.min),col="blue",lwd=2)
            legend("topleft",legend="lambda.min",lwd=3,bty="n",col='blue')
        }
    }
}

dev.off()

### generate R2 and coeff table
R2_mat <- matrix(rep(0,ncol(X)*ncol(Y)),nrow=ncol(X))
coeff_mat <- matrix(rep(0,ncol(X)*ncol(Y)),nrow=ncol(X))
rownames(R2_mat) <- colnames(X)
colnames(R2_mat) <- colnames(Y)
rownames(coeff_mat) <- colnames(X)
colnames(coeff_mat) <- colnames(Y)
for(TFname in colnames(X)){
    for(HMname in colnames(Y)){
        unilinResult <- unilinear_only(X[,TFname],Y[,HMname])
        R2_mat[TFname,HMname] <- unilinResult[1]
        coeff_mat[TFname,HMname] <- unilinResult[2]
    }
}

summary_table <- c()
### select candidates
# for each candidates detected in el-net step, require R2>0.1, coeff<0, R2 for any positive correlated HM signal (coeff>0) < 0.1 
for(TFname in colnames(X)){
    this_usage <- coTFusage[TFname,]
    this_R2 <- R2_mat[TFname,]
    this_coeff <- coeff_mat[TFname,]

    classicFun <- 0
    nonClassicHMidx <- c()

    for(i in 1:ncol(coTFusage)){
        if(this_usage[i]==1 & this_R2[i] >= R2cutoff & this_coeff[i] < 0){
            nonClassicHMidx <- c(nonClassicHMidx, i)
        }else if(this_R2[i] >= R2classic & this_coeff[i] > 0 ){
            classicFun <- 1
        }
    }

    if(classicFun == 0){
        for(idx in nonClassicHMidx){
            nonClassicHM <- colnames(Y)[idx]
            R2 <- this_R2[idx]
            coeff <- this_coeff[idx]
            Pval <- unilinear_permute(X[,TFname],Y[,idx])
            summary_table <- rbind(summary_table, c(TFname, nonClassicHM, Pval, R2, coeff))
        }        
    }
}

### summary and output steps
if(is.null(summary_table)){
    print("no significant candidates detected")
    write.table("no non-classical function detected",file=paste0(outname,"_filterNC.txt"),quote=F,sep="\t",row.names=F,col.names=F)
    write.table("no non-classical function detected",file=paste0(outname,'_NCsummary.txt'),quote=F,sep="\t",row.names=F,col.names=F)
}else{
    colnames(summary_table) <- c("TFname","HMname","Pval","R2","coeff")
    if(nrow(summary_table) > 1){
        summary_table <- summary_table[order(as.numeric(summary_table[,"R2"]),decreasing=TRUE),]
    }

    if(topN == "all"){
        topN <- nrow(summary_table)
    }else{
        topN <- as.numeric(topN)
    }
    
    if(topN == 1){
        out_table_raw <- summary_table
    }else{
        out_table_raw <- as.matrix(summary_table[1:topN,])
    }
    
    ## separate NC/C peaks based on HMsignal using otsu's method
    peakgroup <- c()
    for(i in 1:ncol(rawY)){
        thisY <- rawY[,i]
        thisY_NCotsu <- signal2cutoff(thisY)
        thisY_peakgroup <- thisY_NCotsu[[3]][,2]
        peakgroup <- cbind(peakgroup, thisY_peakgroup)
    }
    colnames(peakgroup) <- colnames(Y)

    ## for each predicted coTFvsHM pair, output the cobinding NC sites as a bed file
    if(!file.exists("nonClassicalPeaks/")){dir.create("nonClassicalPeaks")}
    if(topN == 1){
        pdf(file=paste0(outname,"_cofactor_HMsignal.pdf"),height=6,width=6)
        par(mar=c(4,4,2,2))    
    }else if(topN == 2){
        pdf(file=paste0(outname,"_cofactor_HMsignal.pdf"),height=6,width=12)
        par(mfrow=c(1,2),mar=c(4,4,2,2))
    }else{
        pdf(file=paste0(outname,"_cofactor_HMsignal.pdf"),height=12,width=12)
        par(mfrow=c(2,2),mar=c(4,4,2,2))
    }
    num_NCsites <- c()
    for(topnum in 1:nrow(out_table_raw)){
        coTF = out_table_raw[topnum,1]
        substrateHM = out_table_raw[topnum,2]
        coTF_cobinding <- X[,coTF]
        HMnc <- peakgroup[,substrateHM]
        cobinding_NC_peak <- peakX[which(coTF_cobinding>0 & HMnc == 0),]
        write.table(cobinding_NC_peak, file=paste0("nonClassicalPeaks/",outname,"_",coTF,"_",substrateHM,"_top",topnum,"nonclassical_peaks.bed"),quote=F,sep="\t",row.names=F,col.names=F)
        num_NCsites <- c(num_NCsites, nrow(cobinding_NC_peak))
        if(topnum %in% 1:4){
            ## boxplot compare the histone modification signal between on classical and non-classical peaks
            boxplot(Y[which(coTF_cobinding>0 & HMnc == 0)],Y[which(coTF_cobinding==0 | HMnc > 0)],
                names=c("non-classical peak","classical peak"),ylab=paste0(substrateHM," signal"),main=paste0(coTF," & ",substrateHM),
                outline=F,cex.main=1)
            legend("topleft",legend=paste0("#NCpeak = ",nrow(cobinding_NC_peak)),bty="n")
        }
    }
    dev.off()
    
    out_table <- cbind(out_table_raw, num_NCsites)
    write.table(summary_table,file=paste0(outname,"_filterNC.txt"),quote=F,sep="\t",row.names=F,col.names=T)
    write.table(out_table,file=paste0(outname,'_NCsummary.txt'),quote=F,sep="\t",row.names=F,col.names=T)
} 
 
"""
    createDIR("tmpPackage/")
    outf = open("tmpPackage/detectNonClassic.r",'w')
    outf.write(Rscript)
    outf.close()
    cmd = "Rscript %s %s %s %s %s %s %s %s %s"%("tmpPackage/detectNonClassic.r",
                         conf_dict['General']['outname'],
                         ",".join(conf_dict['General']['signalname']),
                         conf_dict['options']['Pvalue'],
                         conf_dict['options']['Alpha'],
                         conf_dict['options']['Lambda'],
                         conf_dict['options']['TopNcofactors'],
                         conf_dict['General']['startdir']+"tmpPackage/",
                         conf_dict['options']['Rcutoff'])
    #rwlog(cmd,logfile)
    os.system('echo "[CMD] %s " >> %s'%(cmd,logfile))
    tmpobj = sp(cmd)

    return conf_dict

