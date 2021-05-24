#!/usr/bin/env python
"""

Function declare:
 
def CMD                  (cmd)
def sp                   (cmd)
def sperr                (cmd)
def raise_error          ()
def detect_memory        ()
def pdf_name             (input_name)
def wlog                 (message,logfile)
def ewlog                (message,logfile)
def rwlog                (message,logfile)
def readAnnotation       (annotation)
def textformat           (inp)
def createDIR            (dirname)
def strlatexformat       (instr)
def strdis               (str1,str2)
def sample_down_transform_sam (samfile,outbed,sampledownsam,sampledownbed,sample_reads)
def transform_refgene    (refgene,ttsdis,outname)
def reform_barcode_fastq (fq,reformtxt,cbL,umiL)
def combine_reads        (barcodeF,cdsF,utr3F,utr5F,symbolF,ttsdisF,outF,dup_measure)
def generate_matrix      (refgene,inputbed,ttsdis,qcmatfull,qcmat,expmat,coverGNcutoff,umidis1)
"""
import subprocess
import sys
import os
import math
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

def checkbedformat(bedfile,cutoff):
    inf = open(bedfile)
    line = inf.readline()
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
    inf = open(bedfile)
    outf_chrM = open(outname + "_chrM.bed",'w')
    outf_chromatin = open(outname + "_chromatin.bed",'w')

    chrom_reads = {}
    for line in inf:
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
    inf =open(outname + "_chromatin.bed")
    for line in inf:
        ll = line.strip().split("\t")
        chrom = ll[0]
        cellname = ll[3]
        if chrom == "chrM":
            continue
        if not cellname in cell_reads:
            cell_reads[cellname]=0
        cell_reads[cellname] += 1

    highQcells = []
    for cell in cell_reads.keys():
        if cell_reads[cell] >= cutoff:
            highQcells.append(cell)

    inf.seek(0)

    if len(highQcells) < 100:
        return "fail"

    if len(usecells) == 0:
        usehighQcells = highQcells
    else:
        usehighQcells = [value for value in highQcells if value in usecells]


    if len(usehighQcells) < 100:
        finalcell = highQcells
    else:
        finalcell = usehighQcells

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
    return [len(usehighQcells), len(finalcell), highQreadnum,len(cell_reads.keys())]

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

def pileup_cleavage(outname,bdg2bw,csize):
    pluslog1 = sp("macs3 pileup -i %s -f BED --extsize 1 -o %s "%(outname + "_cleavage_plus.bed", outname+ "_cleavage_plus.bdg"))
    pluslog2 = sp("sort -k1,1 -k2,2n %s > %s"%(outname+ "_cleavage_plus.bdg",outname+ "_cleavage_plus_sorted.bdg" ))
    pluslog3 = sp("%s %s %s %s"%(bdg2bw,outname+ "_cleavage_plus_sorted.bdg",csize,outname+ "_cleavage_plus.bw" ))
    minuslog1 = sp("macs3 pileup -i %s -f BED --extsize 1 -o %s "%(outname + "_cleavage_minus.bed", outname+ "_cleavage_minus.bdg"))
    minuslog2 = sp("sort -k1,1 -k2,2n %s > %s"%(outname+ "_cleavage_minus.bdg",outname+ "_cleavage_minus_sorted.bdg" ))
    minuslog3 = sp("%s %s %s %s"%(bdg2bw,outname+ "_cleavage_minus_sorted.bdg",csize,outname+ "_cleavage_minus.bw" ))


def fetchseq_2bit(twoBitToFaTool,seq2bit,chrm,start,end):
    result = sp("%s %s:%s:%s-%s stdout"%(twoBitToFaTool,seq2bit,chrm,start,end))[0].decode("ascii").strip().split("\n")
    if len(result) >=2:
        return result[1]
    else:
        return "NA"
def fetchsignal_bw(bdg2bw,bwfile,chrm,start,end):
    cmd = '%s %s %s %s %s '%(bdg2bw,bwfile,chrm,start,end,end-start)
    rawsig = sp(cmd)[0].decode("ascii").strip().split()
    if len(rawsig) == (end-start):
        sig = [0 if "n" in x.lower() else int(x) for x in rawsig ]
        return sig
    else:
        return [0]*(end-start)

def bias_exp_cleavage(outname,peakfile,biasMat,kmer,bwsum,bdg2bw,twoBitFa,seq2bit):
    Cspan = 25
    flank = int(int(kmer)/2)
    inf = open(peakfile)
    for line in inf:
        ll = line.strip().split("\t")
        chrm = ll[0]
        start = int(ll[1])
        end = int(ll[2])

        if start-Cspan < 0:
            continue
        plus_vector = fetchsignal_bw(bdg2bw, outname+ "_cleavage_plus.bw", chrm,start-Cspan,end+Cspan)
        minus_vector = fetchsignal_bw(bdg2bw, outname+ "_cleavage_minus.bw", chrm,start-Cspan,end+Cspan)
        plus_single_bias_enc_vector = []
        minus_single_bias_enc_vector = []
        for pos in range(start-Cspan,end+Cspan):
            plus_seq = fetchseq_2bit(twoBitFa,seq2bit,chrm,pos-flank,pos+flank).upper()#genome[chrm][(pos-flank):(pos+flank)].upper()
            minus_seq = rev(fetchseq_2bit(twoBitFa,seq2bit,chrm,pos-flank,pos+flank)).upper()
            if len(plus_seq) == kmer and not "N" in plus_seq:
                plus_bias_enc = biasMat[plus_seq]
            else:
                plus_bias_enc = 0

            if len(minus_seq) == kmer and not "N" in minus_seq:
                minus_bias_enc = biasMat[minus_seq]
            else:
                minus_bias_enc = 0
            plus_single_bias_enc_vector.append(plus_bias_enc)
            minus_single_bias_enc_vector.append(minus_bias_enc)
        Plus_Single_encBias = numpy.array(plus_single_bias_enc_vector)
        Minus_Single_encBias = numpy.array(minus_single_bias_enc_vector)

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
            
            if this_plus > 0:
                expcut_plus_raw = this_plus_cuts_sum * (this_plus_single_raw/this_plus_sum_raw)
                expcut_plus_enc = this_plus_cuts_sum * (this_plus_single_enc/this_plus_sum_enc)
                outf.write("\t".join(map(str,[out_chrm+"_"+str(out_start)+"_"+str(out_end)+"_+",this_plus,format(expcut_plus_raw,".3e"),format(expcut_plus_enc,".3e") ]))+"\n")

            if this_minus > 0:
                expcut_minus_raw = this_minus_cuts_sum * (this_minus_single_raw/this_minus_sum_raw)
                expcut_minus_enc = this_minus_cuts_sum * (this_minus_single_enc/this_minus_sum_enc)
                outf.write("\t".join(map(str,[out_chrm+"_"+str(out_start)+"_"+str(out_end)+"_-",this_minus,format(expcut_minus_raw,".3e"),format(expcut_minus_enc,".3e") ]))+"\n")

            #outf_rawPlus.write("\t".join(map(str,[out_chrm,out_start,out_end,this_plus_single]))+"\n")
            #outf_rawMinus.write("\t".join(map(str,[out_chrm,out_start,out_end,this_minus_single]))+"\n")
            #outf_cbPlus.write("\t".join(map(str,[out_chrm,out_start,out_end,this_plus_cb]))+"\n")
            #outf_cbMinus.write("\t".join(map(str,[out_chrm,out_start,out_end,this_minus_cb]))+"\n")
            #
            #outf_rawPlus_prop.write("\t".join(map(str,[out_chrm,out_start,out_end,this_plus_single_prop]))+"\n")
            #outf_rawMinus_prop.write("\t".join(map(str,[out_chrm,out_start,out_end,this_minus_single_prop]))+"\n")
            #outf_cbPlus_prop.write("\t".join(map(str,[out_chrm,out_start,out_end,this_plus_cb_prop]))+"\n")
            #outf_cbMinus_prop.write("\t".join(map(str,[out_chrm,out_start,out_end,this_minus_cb_prop]))+"\n")
    outf.close()
    #outf_rawPlus.close()
    #outf_rawMinus.close()
    #outf_cbPlus.close()
    #outf_cbMinus.close()
    #outf_rawPlus_prop.close()
    #outf_rawMinus_prop.close()
    #outf_cbPlus_prop.close()
    #outf_cbMinus_prop.close()

    inf.close()


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
    outstr = instr.replace('_','\_')
    return(outstr)
    

            
       
