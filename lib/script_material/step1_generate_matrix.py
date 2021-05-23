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
                                   rlogonly,
                                   CMD,
                                   bwsigAve,
                                   createDIR
                                   )
# --------------------------
# main 
# --------------------------

def step1_generate_matrix(conf_dict,logfile):
    '''
    generate expression matrix file 
    main data processing step, including mapping, generate expression matrix and QC matrix which is used in next step
    for fastq format : 
        STAR/bowtie2 mapping
        q30 filter, 
    for sam format:
        q30 filter     
    ''' 
    #t= time.time()
    ### generate TF overlap matrix using bed tools
    wlog("generate TF overlap matrix",logfile)
    init = 0
    for f in conf_dict['General']['usefilename']:
        if init == 0:
            cmd = '%s intersect -a %s -b %s -c '%(conf_dict['General']['bedtools'],conf_dict['General']['HMRpeak'],conf_dict['General']['peakFolder'] + f + ".bed")
            init = 1
        else:
            cmd += '| %s intersect -a - -b %s -c '%(conf_dict['General']['bedtools'],conf_dict['General']['peakFolder'] + f + ".bed")
    
    cmd += '> %s'%(conf_dict['General']['outname']+"_peakov.bed")
    rlogonly(cmd,logfile)
    ### generate TF signal matrix using bwsummary
    if conf_dict['General']['mode'] == "signal":
        wlog("generate TF signal matrix",logfile)
        inf = open(conf_dict['General']['HMRpeak'])
        outf = open(conf_dict['General']['outname']+"_TFsig.bed",'w')
        for line in inf:
            ll = line.split()
            center = int((int(ll[1]) + int(ll[2]))/2)
            start = max(0,center-1000)
            end = center + conf_dict['options']['ext']
            addsigALL = []
            for BWname in conf_dict['General']['usefilename']:
                bwsigfile = conf_dict['General']['bwfolder'] + BWname + ".bw"
                addsigALL.append(bwsigAve(bwsigfile,ll[0],start,end,conf_dict['General']['bwsummary']))
            newll = ll + addsigALL
            outf.write("\t".join(map(str,newll))+"\n")
        inf.close()
        outf.close()
    ### generate HMsig

    wlog("generate HMsignal matrix",logfile)
    inf = open(conf_dict['General']['HMRpeak'])
    outf = open(conf_dict['General']['outname']+"_HMsig.bed",'w')
    for line in inf:
        ll = line.split()
        center = int((int(ll[1]) + int(ll[2]))/2)
        start = max(0,center-conf_dict['options']['ext'])
        end = center + conf_dict['options']['ext']
        addsigALL = []
        for bwsigfile in conf_dict['General']['signalfile']:
            addsigALL.append(bwsigAve(bwsigfile,ll[0],start,end,conf_dict['General']['bwsummary']))
        newll = ll + addsigALL
        outf.write("\t".join(map(str,newll))+"\n")
    inf.close()
    outf.close()

    #s1time = time.time() -t
    #wlog("time for Step1: %s"%(s1time),logfile)
    #conf_dict['results'] = {}
    #conf_dict['results']['expmat'] = conf_dict['Step2_ExpMat']['expmat']
    #conf_dict['results']['qcmat'] = conf_dict['Step2_ExpMat']['qcmat']
    
    return conf_dict









