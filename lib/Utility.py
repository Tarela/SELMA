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
    peaknum = int(sp("wc -l %s"%(bedfile))[0].split()[0])
    if peaknum < cutoff:
        return "lesspeak"
    return "pass"

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
    

            
       
