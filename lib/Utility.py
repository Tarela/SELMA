#!/usr/bin/env python
"""
Utilities used in SELMA pipeline
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

def CMD(cmd):
    os.system(cmd) 

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
   
def raise_error():
    '''
    Raise an error messgae and exit
    '''
    print('error occurs, check log file~!')
    sys.exit(1)


def detect_memory():
    meminfo={}#OrderedDict()
    try:
        with open('/proc/meminfo') as f:
            for line in f:
                meminfo[line.split(':')[0].strip()] = line.split(':')[1].strip()
        totalM = meminfo['MemTotal'].split()
        #freeM = meminfo['MemFree'].split()
        if totalM[1].lower() == "kb":
            try:
                totalM_G = int(totalM[0])/1e6
                return totalM_G
            except:
                return 'NA'
        else:
            return 'NA'    
    except:
        return 'NA'


def pdf_name(input_name):
    '''
    Change filename to pdf file name
    '''
    outputname = "_".join(input_name.split('.')[:-1])+".pdf"
    return outputname
    
def wlog(message,logfile):
    '''
    print a message and write the message to logfile
    '''
    message = "### "+message
    print(message)
    os.system('echo "%s " >> %s'%(message,logfile))
    
def ewlog(message,logfile):
    '''
    print an error message and write the error message to logfile
    then exit Dr.seq
    error messages start with [ERROR]
    '''
    print("[ERROR] %s "%(message))
    os.system('echo "[ERROR] %s " >> %s'%(message,logfile))
    raise_error()
    
def rwlog(cmd,logfile) :
    '''
    print an (shell) command line and write the command line to logfile
    then conduct the command line
    command lines start with [CMD]
    '''
    print("[CMD] %s "%(cmd))
    os.system('echo "[CMD] %s " >> %s'%(cmd,logfile))
    CMD(cmd)

def rlogonly(cmd,logfile) :
    '''
    print an (shell) command line and write the command line to logfile
    then conduct the command line
    command lines start with [CMD]
    '''
    #print "[CMD] %s "%(cmd)
    os.system('echo "[CMD] %s " >> %s'%(cmd,logfile))
    CMD(cmd)

def checkbedformat(bedfile):
    if bedfile.endswith(".bed"):
        inf = open(bedfile)
        line = inf.readline()
        inf.close()
    else:
        inf = gzip.open(bedfile,'rb')
        line = inf.readline().decode("ascii")
        inf.close()

    ll = line.strip().split("\t")
    if len(ll) < 3:
        return "fail"
    try:
        a=int(ll[1])
        b=int(ll[2])
    except:
        return "fail"
    if len(ll) < 6 and len(ll)>=3:
        return "PE"
    elif len(ll) >= 6 and ll[5] in ["+","-"]:
        return "SE"
    else:
        return "fail"
    #peaknum = int(sp("wc -l %s"%(bedfile))[0].split()[0])
    #if peaknum < cutoff:
    #    return "lesspeak"
    #return peaknum#"pass"

def split_chromosome_reads(bedfile,outname,scATAC10x,usechrom):
    if bedfile.endswith(".bed"):
        inf = open(bedfile)
    else:
        inf = gzip.open(bedfile,'rb')
    outf_chrM = open(outname + "_chrM.bed",'w')
    outf_chromatin = open(outname + "_chromatin.bed",'w')

    chrom_reads = {}
    for lineRaw in inf:
        if bedfile.endswith(".bed.gz"):
            line = lineRaw.decode("ascii")
        else:
            line = lineRaw
        ll = line.strip().split("\t")
        chrom=ll[0]
        if not chrom in usechrom:
            continue
        if not chrom in chrom_reads:
            chrom_reads[chrom] = 0
        chrom_reads[chrom] += 1
        if scATAC10x:
            if int(ll[1]) >= 4:
                newll = [ll[0], max(0,int(ll[1])-4), int(ll[2])+5, ll[3]  ]
                newline = "\t".join(map(str,newll))+"\n"
                if chrom == "chrM":
                    outf_chrM.write(newline)
                else:
                    outf_chromatin.write(newline)
        else:
            if chrom == "chrM":
                outf_chrM.write(line)
            else:
                outf_chromatin.write(line)
    inf.close()
    outf_chrM.close()
    outf_chromatin.close()
    return chrom_reads

def filter_highQcell_reads(outname,cutoff,usecells):

    cell_reads = {}
    allcells = []
    inf =open(outname + "_chromatin.bed")
    for line in inf:
        ll = line.strip().split("\t")
        chrom = ll[0]
        cellname = ll[3]
        if chrom == "chrM":
            continue
        if not cellname in cell_reads:
            cell_reads[cellname]=0
            allcells.append(cellname)
        cell_reads[cellname] += 1

    highQcells = []
    for cell in allcells:#cell_reads.keys():
        if cell_reads[cell] >= cutoff:
            highQcells.append(cell)

    inf.seek(0)

    if len(highQcells) < 100:
        return "fail"

    if len(usecells) == 0:
        usehighQcells = highQcells
    else:
#        usehighQcells = [value for value in highQcells if value in usecells]
        usehighQcells = [value for value in usecells if value in highQcells]
    if len(usehighQcells) < 100:
        finalcell = highQcells
        usetag = "highQ"
    else:
        finalcell = usehighQcells
        usetag = "highQuse"

    highQcellnum = len(finalcell)
    highQreadnum = 0
    outf = open(outname + "_highQcellReads.bed",'w')
    for line in inf:
        ll = line.strip().split("\t")
        chrom = ll[0]
        cellname = ll[3]
        if chrom == "chrM":
            continue
        if cellname in  finalcell:
            outf.write(line)
            highQreadnum += 1
    outf.close()
    inf.close()
    return [finalcell, len(usehighQcells), highQreadnum,len(cell_reads.keys()),usetag]

def readBias(bgmatrix):
    pBG = {}
    if bgmatrix.endswith(".gz"):
        inf = gzip.open(bgmatrix,mode="rb")
        for line in inf:
            ll = line.decode("ascii").strip().split("\t")
            name = ll[0]
            pBG[name] = float(ll[1])
        inf.close()
    else:
        inf = open(bgmatrix)
        for line in inf:
            ll = line.strip().split("\t")
            name = ll[0]
            pBG[name] = float(ll[1])
        inf.close()
    return pBG

def rev(seq):
    revseq = ""
    for i in seq[::-1]:
        if i == 'A':
            r = 'T'
        elif i == 'T':
            r = 'A'
        elif i == 'C':
            r = 'G'
        elif i == 'G':
            r = 'C'
        else:
            r=i#print i
        revseq += r
    return revseq

def make_nmer_dict(n):
    nmer_seq = {}
    bp = ['A','C','G','T']
    allseq = [0]*n
    allseq[0] = bp
    i=1
    while i < n:
        allseq[i] = []
        for previous_seq in allseq[i-1]:
            for add_bp in bp:
                new_seq = previous_seq + add_bp
                allseq[i].append(new_seq)
        i += 1
    for seq in allseq[n-1]:
        nmer_seq[seq] = 0
    del allseq
    return nmer_seq

def naive_kmerBias_chrM(outname,seq2bit,kmer,twoBitToFaTool,dataformat):
    flank = int(int(kmer)/2)
    pcut = make_nmer_dict(int(kmer))
    bgseq = make_nmer_dict(int(kmer))
    chrMseq = sp("%s %s:chrM stdout"%(twoBitToFaTool,seq2bit))[0].decode("ascii").replace("\n","").lstrip(">chrM")

    for i in range(len(chrMseq)):
        subseq = chrMseq[(i-flank):(i+flank)]
        if subseq in bgseq:
            bgseq[subseq]+=1

    inf = open(outname +"_chrM.bed")
    for line in inf:
        ll = line.strip().split("\t")
        if dataformat == "PE" or ll[5] == '+' :
            seq = chrMseq[(int(ll[1])-flank):(int(ll[1])+flank)]#.upper()
            if seq in pcut:
                pcut[seq]+=1
            else:
                pass
        if dataformat == "PE" or ll[5] == '-' :
            seq = rev(chrMseq[(int(ll[2])-flank):(int(ll[2])+flank)])#.upper())
            if seq in pcut:
                pcut[seq]+=1
            else:
                pass
    inf.close()        

    naiveBias = {}
    for seqtype in sorted(pcut.keys()):
        if bgseq[seqtype] == 0:
            pbias = -1
        else:  
            pbias = float(pcut[seqtype])/float(bgseq[seqtype])
        naiveBias[seqtype]=pbias
    return naiveBias

def extsummit(summitfile,outfile,extsize):
    peakcount = 0
    inf = open(summitfile)
    outf = open(outfile,'w')
    for line in inf:
        ll = line.strip().split("\t")
        if ll[0] != "chrM":
            peakcount += 1
            newll = [ll[0], max(0,int(ll[1])-extsize), int(ll[1])+extsize, ll[3],ll[4]]
            outf.write("\t".join(map(str,newll))+"\n")
    inf.close()
    outf.close()
    return peakcount

def extExternal(peakfile,outfile,extsize):
    peakcount = 0
    inf = open(peakfile)
    outf = open(outfile,'w')
    for line in inf:
        ll = line.strip().split("\t")
        if ll[0] != "chrM":
            peakcount += 1
            peakcenter = int((int(ll[1]) + int(ll[2]))/2)
            newll = [ll[0], max(0,peakcenter-extsize), peakcenter+extsize, "extPeak%s"%peakcount,"1"]
            outf.write("\t".join(map(str,newll))+"\n")
    inf.close()
    outf.close()
    return peakcount

def fetchseq_2bit(twoBitToFaTool,seq2bit,chrm,start,end):
    result = sp("%s %s:%s:%s-%s stdout"%(twoBitToFaTool,seq2bit,chrm,start,end))[0].decode("ascii").strip().split("\n")
    if len(result) >=2:
        return "".join(result[1:]).upper()
    else:
        return "NA"

def fetchseq_2bit_chrom(twoBitToFaTool,seq2bit,chrm):
    result = sp("%s %s:%s stdout"%(twoBitToFaTool,seq2bit,chrm))[0].decode("ascii").strip().split("\n")
    if len(result) >=2:
        return "".join(result[1:]).upper()
    else:
        return "NA"


def fetchsignal_bw(bwsum,bwfile,chrm,start,end):
    cmd = '%s %s %s %s %s %s'%(bwsum,bwfile,chrm,start,end,end-start)
    rawsig = sp(cmd)[0].decode("ascii").strip().split()
    if len(rawsig) == (end-start):
        sig = [0 if "n" in x.lower() else int(x) for x in rawsig ]
        return sig
    else:
        return [0]*(end-start)


def bias_exp_cleavage_DNase(outname,peakfile,biasMat,kmer,bedtools,seq_dict,totalreads,dataformat):

    Cspan = 25
    kmer=int(kmer)
    flank = int(kmer/2)

    # extend peak to peak+Cspan, split to chromosome level
    chromosome_peak_dict = {}
    plus_cut_dict = {}
    minus_cut_dict = {}
    # split merge peak
    inf = open(peakfile)
    count = 0
    for line in inf:
        ll = line.split()
        chrom = ll[0]
        count += 1
        newll = [chrom, int(ll[1]) - Cspan, int(ll[2]) + Cspan, "mergePeak%s"%count]
        if not chrom in chromosome_peak_dict:
            chromosome_peak_dict[chrom] = open("%s_mergePeaks.bed"%(chrom),'w')
            plus_cut_dict[chrom] = open("%s_plusCuts.bed"%(chrom),'w')
            minus_cut_dict[chrom] = open("%s_minusCuts.bed"%(chrom),'w')
        chromosome_peak_dict[chrom].write("\t".join(map(str,newll))+"\n")
    inf.close()
    for chrom in chromosome_peak_dict.keys():
        chromosome_peak_dict[chrom].close()
    # split reads (to plus and minus)
    inf = open(totalreads)
    count = 0
    for line in inf:
        ll = line.split()
        chrom = ll[0]
        if not chrom in chromosome_peak_dict:
            continue
        if dataformat == "PE":
            count += 1
            newll = [chrom, ll[1] ,int(ll[1])+1, "c%s"%count,".","+"]
            plus_cut_dict[chrom].write("\t".join(map(str,newll))+"\n")
            count += 1
            newll = [chrom, int(ll[2])-1 ,ll[2], "c%s"%count,".","-"]
            minus_cut_dict[chrom].write("\t".join(map(str,newll))+"\n")
        else:
            count += 1
            if ll[5] == "+":
                newll = [chrom, ll[1] ,int(ll[1])+1, "c%s"%count,".","+"]
                plus_cut_dict[chrom].write("\t".join(map(str,newll))+"\n")
            else:
                newll = [chrom, int(ll[2])-1 ,ll[2], "c%s"%count,".","-"]
                minus_cut_dict[chrom].write("\t".join(map(str,newll))+"\n")
    inf.close()
    for chrom in chromosome_peak_dict.keys():
        plus_cut_dict[chrom].close()
        minus_cut_dict[chrom].close()


    outf_plus = open(outname + "_biasExpCuts_plus.bdg",'w')
    outf_minus = open(outname + "_biasExpCuts_minus.bdg",'w')
    outf_plusCuts = open(outname + "_cleavage_plus.bdg",'w')
    outf_minusCuts = open(outname + "_cleavage_minus.bdg",'w')

    ### for each chromosome, intersect and calculate cleavage pattern
    for chrom in chromosome_peak_dict.keys():
        OVcmd1 = """%s intersect -a %s -b %s -wao > %s """%(bedtools,"%s_mergePeaks.bed"%(chrom), "%s_plusCuts.bed"%(chrom), "%s_plusCutsOnPeak.bed"%(chrom))
        OVcmd2 = """%s intersect -a %s -b %s -wao > %s """%(bedtools,"%s_mergePeaks.bed"%(chrom), "%s_minusCuts.bed"%(chrom), "%s_minusCutsOnPeak.bed"%(chrom))
        os.system(OVcmd1)
        os.system(OVcmd2)

        #### readin bias vector
        mergePeaks = []
        mergePeak_dict = {}
        thisChrom_plus_bias_dict = {}
        thisChrom_minus_bias_dict = {}
        thisChrom_plus_cuts_dict = {}
        thisChrom_minus_cuts_dict = {}
        inf = open("%s_mergePeaks.bed"%(chrom))
        for line in inf:
            ll = line.strip().split("\t")
            chrm = ll[0]
            start = int(ll[1])
            end = int(ll[2])
            peakname = ll[3]
            plus_seq_all = seq_dict[chrm][(start-flank):(end+flank)]
            minus_seq_all = rev(plus_seq_all)
            plus_single_bias_enc_vector = []
            minus_single_bias_enc_vector = []
            for pos in range(end-start):
                plus_seq = plus_seq_all[pos:(pos+kmer)]
                minus_seq = minus_seq_all[pos:(pos+kmer)]
                if len(plus_seq) == kmer and not "N" in plus_seq:
                    plus_bias_enc = 2**biasMat[plus_seq]
                else:
                    plus_bias_enc = 1
    
                if len(minus_seq) == kmer and not "N" in minus_seq:
                    minus_bias_enc = 2**biasMat[minus_seq]
                else:
                    minus_bias_enc = 0
                plus_single_bias_enc_vector.append(plus_bias_enc)
                minus_single_bias_enc_vector.append(minus_bias_enc)
            Plus_Single_encBias = numpy.array(plus_single_bias_enc_vector)
            Minus_Single_encBias = numpy.array(minus_single_bias_enc_vector[::-1]) 
            thisChrom_plus_bias_dict[peakname] = Plus_Single_encBias
            thisChrom_minus_bias_dict[peakname] = Minus_Single_encBias
            thisChrom_plus_cuts_dict[peakname] = [0]*(end-start)#Plus_Single_encBias
            thisChrom_minus_cuts_dict[peakname] = [0]*(end-start)#Minus_Single_encBias
            mergePeaks.append(peakname)
            mergePeak_dict[peakname] = [chrm,start+Cspan,end-Cspan,peakname]
        inf.close()

        ### readin cuts 
        inf = open("%s_plusCutsOnPeak.bed"%(chrom))
        for line in inf:        
            ll = line.strip().split("\t")
            chrm = ll[0]
            peak_start = int(ll[1])#+Cspan
            peak_end = int(ll[2])#-Cspan
            peakname = ll[3]
            read_start = int(ll[5])
            if read_start != -1:
                relative_pos = read_start - peak_start
            thisChrom_plus_cuts_dict[peakname][relative_pos] += 1
        inf.close()

        inf = open("%s_minusCutsOnPeak.bed"%(chrom))
        for line in inf:        
            ll = line.strip().split("\t")
            chrm = ll[0]
            peak_start = int(ll[1])#+Cspan
            peak_end = int(ll[2])#-Cspan
            peakname = ll[3]
            read_start = int(ll[5])
            if read_start != -1:
                relative_pos = read_start - peak_start
            thisChrom_minus_cuts_dict[peakname][relative_pos] += 1
        inf.close()

        # calculate biasExpCuts
        for peakname in mergePeaks:
            peakinfo = mergePeak_dict[peakname]
            chrm = peakinfo[0]
            start = peakinfo[1]
            end = peakinfo[2]
            Plus_Single_encBias = thisChrom_plus_bias_dict[peakname]
            Minus_Single_encBias = thisChrom_minus_bias_dict[peakname]
            plus_vector = thisChrom_plus_cuts_dict[peakname]
            minus_vector = thisChrom_minus_cuts_dict[peakname]

            for outpos in range(Cspan,(end-start+Cspan)):
                this_plus_single_enc = Plus_Single_encBias[outpos]
                this_minus_single_enc = Minus_Single_encBias[outpos]
                this_plus_sum_enc = sum(Plus_Single_encBias[(outpos-Cspan):(outpos+Cspan)])
                this_minus_sum_enc = sum(Minus_Single_encBias[(outpos-Cspan):(outpos+Cspan)])
    
                this_plus = plus_vector[outpos]
                this_minus = minus_vector[outpos]
                this_plus_cuts_sum = sum(plus_vector[(outpos-Cspan):(outpos+Cspan)])
                this_minus_cuts_sum = sum(minus_vector[(outpos-Cspan):(outpos+Cspan)])
    
                out_chrm = chrm
                out_start = start + outpos - Cspan
                out_end = out_start+1
                
                expcut_plus_enc = this_plus_cuts_sum * (this_plus_single_enc/this_plus_sum_enc)
                expcut_minus_enc = this_minus_cuts_sum * (this_minus_single_enc/this_minus_sum_enc)
    
                outf_plus.write("\t".join( map(str, [out_chrm,out_start,out_end,round(expcut_plus_enc,6)] ))+"\n")
                outf_minus.write("\t".join( map(str, [out_chrm,out_start,out_end,round(expcut_minus_enc,6)] ))+"\n")
                outf_plusCuts.write("\t".join( map(str, [out_chrm,out_start,out_end,round(this_plus,6)] ))+"\n")
                outf_minusCuts.write("\t".join( map(str, [out_chrm,out_start,out_end,round(this_minus,6)] ))+"\n")
    outf_plus.close()
    outf_minus.close()
    outf_plusCuts.close()
    outf_minusCuts.close()
    inf.close()


def bias_exp_cleavage_ATAC(outname,peakfile,biasMat,kmer,bedtools,seq_dict,totalreads,dataformat):

    offset=9
    Cspan = 25
    kmer=int(kmer)
    flank = int(kmer/2)

    chromosome_peak_dict = {}
    plus_cut_dict = {}
    minus_cut_dict = {}
    # split merge peak
    inf = open(peakfile)
    count = 0
    for line in inf:
        ll = line.split()
        chrom = ll[0]
        count += 1
        newll = [chrom, int(ll[1]) - Cspan, int(ll[2]) + Cspan, "mergePeak%s"%count]
        if not chrom in chromosome_peak_dict:
            chromosome_peak_dict[chrom] = open("%s_mergePeaks.bed"%(chrom),'w')
            plus_cut_dict[chrom] = open("%s_plusCuts.bed"%(chrom),'w')
            minus_cut_dict[chrom] = open("%s_minusCuts.bed"%(chrom),'w')
        chromosome_peak_dict[chrom].write("\t".join(map(str,newll))+"\n")
    inf.close()
    for chrom in chromosome_peak_dict.keys():
        chromosome_peak_dict[chrom].close()
    # split reads (to plus and minus)
    inf = open(totalreads)
    count = 0
    for line in inf:
        ll = line.split()
        chrom = ll[0]
        if not chrom in chromosome_peak_dict:
            continue
        if dataformat == "PE":
            count += 1
            newll = [chrom, ll[1] ,int(ll[1])+1, "c%s"%count,".","+"]
            plus_cut_dict[chrom].write("\t".join(map(str,newll))+"\n")
            count += 1
            newll = [chrom, int(ll[2])-1 ,ll[2], "c%s"%count,".","-"]
            minus_cut_dict[chrom].write("\t".join(map(str,newll))+"\n")
        else:
            count += 1
            if ll[5] == "+":
                newll = [chrom, ll[1] ,int(ll[1])+1, "c%s"%count,".","+"]
                plus_cut_dict[chrom].write("\t".join(map(str,newll))+"\n")
            else:
                newll = [chrom, int(ll[2])-1 ,ll[2], "c%s"%count,".","-"]
                minus_cut_dict[chrom].write("\t".join(map(str,newll))+"\n")
    inf.close()
    for chrom in chromosome_peak_dict.keys():
        plus_cut_dict[chrom].close()
        minus_cut_dict[chrom].close()


    outf_plus = open(outname + "_biasExpCuts_plus.bdg",'w')
    outf_minus = open(outname + "_biasExpCuts_minus.bdg",'w')
    outf_plusCuts = open(outname + "_cleavage_plus.bdg",'w')
    outf_minusCuts = open(outname + "_cleavage_minus.bdg",'w')

    ### for each chromosome, intersect and calculate cleavage pattern
    for chrom in chromosome_peak_dict.keys():
        OVcmd1 = """%s intersect -a %s -b %s -wao > %s """%(bedtools,"%s_mergePeaks.bed"%(chrom), "%s_plusCuts.bed"%(chrom), "%s_plusCutsOnPeak.bed"%(chrom))
        OVcmd2 = """%s intersect -a %s -b %s -wao > %s """%(bedtools,"%s_mergePeaks.bed"%(chrom), "%s_minusCuts.bed"%(chrom), "%s_minusCutsOnPeak.bed"%(chrom))
        os.system(OVcmd1)
        os.system(OVcmd2)

        #### readin bias vector
        mergePeaks = []
        mergePeak_dict = {}
        thisChrom_plus_bias_dict = {}
        thisChrom_minus_bias_dict = {}
        thisChrom_plus_cuts_dict = {}
        thisChrom_minus_cuts_dict = {}
        inf = open("%s_mergePeaks.bed"%(chrom))
        for line in inf:
            ll = line.strip().split("\t")
            chrm = ll[0]
            start = int(ll[1])
            end = int(ll[2])
            peakname = ll[3]
            plus_single_bias_enc_vector = []
            minus_single_bias_enc_vector = []
            seqall = seq_dict[chrm][(start-flank-offset):(end+flank+offset)]#fetchseq_2bit(twoBitFa,seq2bit,chrm,start-Cspan-flank,end+Cspan+flank)
            for pos in range(offset,end-start+offset):
                plus_forward_seq = seqall[pos:(pos+kmer)]
                plus_reverse_seq = rev(seqall[(pos+offset):(pos+kmer+offset)])
                minus_forward_seq = rev(seqall[(pos+1):(pos+1+kmer)])
                minus_reverse_seq = seqall[(pos+1-offset):(pos+1+kmer-offset)]
                if len(plus_forward_seq) == kmer and not "N" in plus_forward_seq:
                    plus_bias = 2**biasMat[plus_forward_seq]
                else:
                    plus_bias = 1
                if len(minus_forward_seq) == kmer and not "N" in minus_forward_seq:
                    minus_bias = 2**biasMat[minus_forward_seq]
                else:
                    minus_bias = 1
                if len(plus_reverse_seq) == kmer and not "N" in plus_reverse_seq:
                    plus_reverse_bias = 2**biasMat[plus_reverse_seq]
                else:
                    plus_reverse_bias = 1
                if len(minus_reverse_seq) == kmer and not "N" in minus_reverse_seq:
                    minus_reverse_bias = 2**biasMat[minus_reverse_seq]
                else:
                    minus_reverse_bias = 1    
                plus_cb_bias = numpy.sqrt(plus_bias * plus_reverse_bias ) 
                minus_cb_bias = numpy.sqrt(minus_bias * minus_reverse_bias)
                plus_single_bias_enc_vector.append(plus_cb_bias)
                minus_single_bias_enc_vector.append(minus_cb_bias)
                
            Plus_Single_encBias = numpy.array(plus_single_bias_enc_vector)
            Minus_Single_encBias = numpy.array(minus_single_bias_enc_vector)
            thisChrom_plus_bias_dict[peakname] = Plus_Single_encBias
            thisChrom_minus_bias_dict[peakname] = Minus_Single_encBias
            thisChrom_plus_cuts_dict[peakname] = [0]*(end-start)#Plus_Single_encBias
            thisChrom_minus_cuts_dict[peakname] = [0]*(end-start)#Minus_Single_encBias
            mergePeaks.append(peakname)
            mergePeak_dict[peakname] = [chrm,start+Cspan,end-Cspan,peakname]
        inf.close()
        ### readin cuts 
        inf = open("%s_plusCutsOnPeak.bed"%(chrom))
        for line in inf:        
            ll = line.strip().split("\t")
            chrm = ll[0]
            peak_start = int(ll[1])#+Cspan
            peak_end = int(ll[2])#-Cspan
            peakname = ll[3]
            read_start = int(ll[5])
            if read_start != -1:
                relative_pos = read_start - peak_start
            thisChrom_plus_cuts_dict[peakname][relative_pos] += 1
        inf.close()

        inf = open("%s_minusCutsOnPeak.bed"%(chrom))
        for line in inf:        
            ll = line.strip().split("\t")
            chrm = ll[0]
            peak_start = int(ll[1])#+Cspan
            peak_end = int(ll[2])#-Cspan
            peakname = ll[3]
            read_start = int(ll[5])
            if read_start != -1:
                relative_pos = read_start - peak_start
            thisChrom_minus_cuts_dict[peakname][relative_pos] += 1
        inf.close()

        # calculate biasExpCuts
        for peakname in mergePeaks:
            peakinfo = mergePeak_dict[peakname]
            chrm = peakinfo[0]
            start = peakinfo[1]
            end = peakinfo[2]
            Plus_Single_encBias = thisChrom_plus_bias_dict[peakname]
            Minus_Single_encBias = thisChrom_minus_bias_dict[peakname]
            plus_vector = thisChrom_plus_cuts_dict[peakname]
            minus_vector = thisChrom_minus_cuts_dict[peakname]

            for outpos in range(Cspan,(end-start+Cspan)):
                this_plus_single_enc = Plus_Single_encBias[outpos]
                this_minus_single_enc = Minus_Single_encBias[outpos]
                this_plus_sum_enc = sum(Plus_Single_encBias[(outpos-Cspan):(outpos+Cspan)])
                this_minus_sum_enc = sum(Minus_Single_encBias[(outpos-Cspan):(outpos+Cspan)])
    
                this_plus = plus_vector[outpos]
                this_minus = minus_vector[outpos]
                this_plus_cuts_sum = sum(plus_vector[(outpos-Cspan):(outpos+Cspan)])
                this_minus_cuts_sum = sum(minus_vector[(outpos-Cspan):(outpos+Cspan)])
    
                out_chrm = chrm
                out_start = start + outpos - Cspan
                out_end = out_start+1
                
                expcut_plus_enc = this_plus_cuts_sum * (this_plus_single_enc/this_plus_sum_enc)
                expcut_minus_enc = this_minus_cuts_sum * (this_minus_single_enc/this_minus_sum_enc)
    
                outf_plus.write("\t".join( map(str, [out_chrm,out_start,out_end,round(expcut_plus_enc,6)] ))+"\n")
                outf_minus.write("\t".join( map(str, [out_chrm,out_start,out_end,round(expcut_minus_enc,6)] ))+"\n")
                outf_plusCuts.write("\t".join( map(str, [out_chrm,out_start,out_end,round(this_plus,6)] ))+"\n")
                outf_minusCuts.write("\t".join( map(str, [out_chrm,out_start,out_end,round(this_minus,6)] ))+"\n")
    outf_plus.close()
    outf_minus.close()
    outf_plusCuts.close()
    outf_minusCuts.close()
    inf.close()
  
  

def bias_peakXcell_mat(outname,bedtools,chrom_list, kmer, biasDict, seqDict, usecells, datatype,peakminreads):

    flank = int(int(kmer)/2)
    offset=9
    peakminreads = int(peakminreads)
    #if peakmaxreads == "NA":
    #    peakmaxreads = int(1e10)
    #else:
    #    peakmaxreads = int(peakmaxreads)

    peakfile = outname + "_summitEXT.bed"
    readfile = outname + "_highQcellReads.bed"
    peakFeatures = open(outname + "_peakFeatures.txt","w")
    peakXcellMat = open(outname + "_peakXcellMat.txt","w")

    newll = ['chrm','start','end','peakname','score','cutsSum','avebias']
    peakFeatures.write("\t".join(newll)+"\n")
    
    newll = ['peakname'] + usecells
    peakXcellMat.write("\t".join(newll)+"\n")

    for chrom in chrom_list:
        if chrom == "chrM":
            continue
        cmdpeak = """awk '{OFS="\\t";if($1=="%s") print $0}' %s > %s"""%(chrom,peakfile, outname+"_tmpSCpeaks.bed")
        tmplog = sp(cmdpeak)
        try:
            chrom_peaknum = int(sp('wc -l %s'%(outname+"_tmpSCpeaks.bed"))[0].decode("ascii").split()[0])
        except:
            #wlog('no peak detected for %s, skip %s'%(chrom,chrom),logfile)
            continue 
        if chrom_peaknum == 0:
            #wlog('no peak detected for %s, skip %s'%(chrom,chrom),logfile)
            continue
        cmdread = """awk '{OFS="\\t";if($1=="%s") print $1,$2,$2+1,$4,".","+\\n"$1,$3-1,$3,$4,".","-"}' %s > %s"""%(chrom,outname+"_highQcellReads.bed", outname+"_tmpSCreads.bed")
        tmplog = sp(cmdread)

        cmdassign = "%s intersect -a %s_tmpSCpeaks.bed -b %s_tmpSCreads.bed -wao | awk '{if($NF > 0 ) print $0}'|sort -k 4,4 > %s_scOVcleavage.bed"%(bedtools,outname,outname,outname)
        tmplog = sp(cmdassign)

        inf = open("%s_scOVcleavage.bed"%outname)
        this_peak = "NA"
        for line in inf:
            ll = line.strip().split("\t")
            if this_peak == "NA":
                this_peak = ll[3]
                this_loci = ll[:5]
                cell_count = [0]*len(usecells)
                cutsSum = 0
                biasSum = 0
                read_info = ll[5:11]
                if  read_info[3] in usecells:
                    this_bias = reads_level_bias(read_info,datatype,seqDict,biasDict,flank)
                    if this_bias != "NA":
                        cutsSum += 1
                        biasSum += this_bias
                        cell_count[usecells.index(read_info[3])] += 1
            elif this_peak != ll[3]:
                if cutsSum >= peakminreads :#and cutsSum <= peakmaxreads:
                    avebias = biasSum / cutsSum#,6)
                    newll =  this_loci + [cutsSum,avebias] 
                    peakFeatures.write("\t".join(map(str,newll))+"\n")
                    newll = [this_peak] + cell_count
                    peakXcellMat.write("\t".join(map(str,newll))+"\n")
                this_peak = ll[3]
                this_loci = ll[:5]
                cell_count = [0]*len(usecells)
                cutsSum = 0
                biasSum = 0
                read_info = ll[5:11]
                if  read_info[3] in usecells:
                    this_bias = reads_level_bias(read_info,datatype,seqDict,biasDict,flank)
                    if this_bias != "NA":
                        cutsSum += 1
                        biasSum += this_bias
                        cell_count[usecells.index(read_info[3])] += 1
            else:
                read_info = ll[5:11]
                if  read_info[3] in usecells:
                    this_bias = reads_level_bias(read_info,datatype,seqDict,biasDict,flank)
                    if this_bias != "NA":
                        cutsSum += 1
                        biasSum += this_bias
                        cell_count[usecells.index(read_info[3])] += 1
        
        if cutsSum >= peakminreads :#and cutsSum <= peakmaxreads:
            avebias = biasSum / cutsSum#,6)
            newll =  this_loci + [cutsSum,avebias] 
            peakFeatures.write("\t".join(map(str,newll))+"\n")
            newll = [this_peak] + cell_count
            peakXcellMat.write("\t".join(map(str,newll))+"\n")
        
        inf.close()

    peakFeatures.close()
    peakXcellMat.close()

def reads_level_bias(region, datatype, seqdict, biasDict, flank):
    flank = int(flank)
    offset=9
    chrm = region[0]
    start = int(region[1])#-100#upstream_ext
    end = int(region[2])#center+100#start + 200#fulllen
    strand = region[5]
    if datatype == "ATAC":
        if strand == "+":
            forward_seq = seqdict[chrm][(start-flank):(start+flank)].upper()
            reverse_seq = rev(seqdict[chrm][(start+offset-flank):(start+offset+flank)].upper())
        else:
            forward_seq = rev(seqdict[chrm][(end-flank):(end+flank)].upper())
            reverse_seq = seqdict[chrm][(end-offset-flank):(end-offset+flank)].upper()
    
        if forward_seq in biasDict and reverse_seq in biasDict:
            bias_score = (biasDict[forward_seq] + biasDict[reverse_seq])/2
        else:
            bias_score = "NA"
    else:

        if strand == "+":
            forward_seq = seqdict[chrm][(start-flank):(start+flank)].upper()
        else:
            forward_seq = rev(seqdict[chrm][(end-flank):(end+flank)].upper())
    
        if forward_seq in biasDict: #and biasDict.has_key(reverse_seq):
            bias_score = biasDict[forward_seq]# + biasDict[reverse_seq])/2
        else:
            bias_score = "NA"
            
    return bias_score


def bwsigAve(bwfile,chrm,start,end,software):
    cmd = '%s %s %s %s %s 1'%(software,bwfile,chrm,start,end)
#    print sp(cmd)
    bwS = sp(cmd)[0].strip()
    if len(bwS) == 0:#bwS == "":
        return 0
    else:
        sigAve = float( bwS)#numpy.array(map(float,sp(cmd)[0].strip().split("\t")))
        return sigAve#[CpGcount,aveME]

def readAnnotation(annotation):
    '''
    read full annotation file and output as a dictionary 
    file format is fixed to UCSC full annotation format
    '''
    inf = open(annotation)
    outdict = {}
    for line in inf:
        ll = line.split()
        outdict[ll[1]] = ll[12]
    return outdict
    
     
def textformat(inp):
    '''
    transfer 1000000 to  1,000,000 for better visualization in output report
    '''
    o = ''
    comma = 0
    for i in (inp[::-1]):
        comma += 1
        o += i
        if comma%3 == 0:
            o += ','
    return o[::-1].strip(',')

def createDIR(dirname):
    '''
    check dir name and create new dir
    '''
    if not os.path.isdir(dirname):
        os.system('mkdir %s'%(dirname))

def strlatexformat(instr):
    instr = str(instr)
    outstr = instr.replace('_','\_')
    return(outstr)
    

            
       
