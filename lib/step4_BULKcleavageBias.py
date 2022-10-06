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
from SELMApipe.Utility      import (sp, 
                                   wlog,
                                   ewlog,
                                   rwlog,
                                   CMD,
                                   fetchseq_2bit_chrom,
                                   bias_exp_cleavage_DNase,
                                   bias_exp_cleavage_ATAC)

# -------------------------- 
# main 
# --------------------------
def step4_BULKcleavageBias(conf_dict,logfile):

    ### preparing mapping state dict 
    wlog('split fragments to strand specific cleavage sites',logfile)
    if conf_dict['General']['format'] == "PE":
        cmdplus = """awk '{OFS="\t";print $1,$2,$2+1,".",".","+"}' %s > %s"""%(conf_dict['General']['outname']+"_chromatin.bed", conf_dict['General']['outname']+"_cleavage_plus.bed")
        cmdminus = """awk '{OFS="\t";print $1,$3-1,$3,".",".","-"}' %s > %s"""%(conf_dict['General']['outname']+"_chromatin.bed", conf_dict['General']['outname']+"_cleavage_minus.bed")
    else:
        cmdplus = """awk '{OFS="\t";if($6=="+") print $1,$2,$2+1,".",".","+"}' %s > %s"""%(conf_dict['General']['outname']+"_chromatin.bed", conf_dict['General']['outname']+"_cleavage_plus.bed")
        cmdminus = """awk '{OFS="\t";if($6=="-") print $1,$3-1,$3,".",".","-"}' %s > %s"""%(conf_dict['General']['outname']+"_chromatin.bed", conf_dict['General']['outname']+"_cleavage_minus.bed")

    tmplog = sp(cmdplus)
    tmplog = sp(cmdminus)

    wlog("remove redundant position from the extended peak file",logfile)
    cmduni = """sort -k 1,1 -k 2,2g -k 3,3g %s | %s merge -i - > %s"""%(conf_dict['results']['peakfile'],conf_dict['General']['bedtools'],conf_dict['General']['outname']+"_summitEXTmerge.bed")
    tmplog = sp(cmduni)

    used_chrm_list= []
    inf = open(conf_dict['results']['peakfile'])
    for line in inf:
        ll = line.split()
        if not ll[0] in used_chrm_list:
            used_chrm_list.append(ll[0])
    inf.close()

    wlog('readin sequence from 2bit',logfile)
    seq_dict = {}
    inf = open(conf_dict['options']['csize']) 
    for line in inf:
        chrm = line.split()[0]
        if chrm in used_chrm_list:
            seq_dict[chrm] = fetchseq_2bit_chrom(conf_dict['General']['twoBitToFa'],conf_dict['General']['sequence'],chrm)
    inf.close()
    conf_dict['results']['seqdict'] = seq_dict

    wlog('calculate bias expected cleavages',logfile)
    if conf_dict['General']['datatype'] == "DNase":
        tmplog = bias_exp_cleavage_DNase(conf_dict['General']['outname'],
                                conf_dict['General']['outname']+"_summitEXTmerge.bed",
                                conf_dict['results']['biasMat'],
                                conf_dict['options']['kmer'],
                                conf_dict['General']['bedtools'],
                                conf_dict['results']['seqdict'],
                                conf_dict['General']['outname']+"_chromatin.bed",
                                conf_dict['General']['format'])
    else:
        tmplog = bias_exp_cleavage_ATAC(conf_dict['General']['outname'],
                                conf_dict['General']['outname']+"_summitEXTmerge.bed",
                                conf_dict['results']['biasMat'],
                                conf_dict['options']['kmer'],
                                conf_dict['General']['bedtools'],
                                conf_dict['results']['seqdict'],
                                conf_dict['General']['outname']+"_chromatin.bed",
                                conf_dict['General']['format'])
#

    wlog('pile up cleavage sites',logfile)
    #pluslog1 = sp("macs3 pileup -i %s -f BED --extsize 1 -o %s "%(conf_dict['General']['outname'] + "_cleavage_plus.bed", conf_dict['General']['outname']+ "_cleavage_plus.bdg"))
    #pluslog = sp("%s genomecov -i %s -g %s -bga -scale 1 > %s "%(conf_dict['General']['bedtools'],conf_dict['General']['outname'] + "_cleavage_plus.bed",conf_dict['options']['csize'], conf_dict['General']['outname']+ "_cleavage_plus.bdg"))
    #pluslog = sp("sort -k1,1 -k2,2n %s > %s"%(conf_dict['General']['outname']+ "_cleavage_plus.bdg",conf_dict['General']['outname']+ "_cleavage_plus_sorted.bdg" ))
    #pluslog = sp("%s %s %s %s"%(conf_dict['General']['bedGraphToBigWig'],conf_dict['General']['outname']+ "_cleavage_plus_sorted.bdg",conf_dict['options']['csize'],conf_dict['General']['outname']+ "_cleavage_plus.bw" ))
    #minuslog1 = sp("macs3 pileup -i %s -f BED --extsize 1 -o %s "%(conf_dict['General']['outname'] + "_cleavage_minus.bed", conf_dict['General']['outname']+ "_cleavage_minus.bdg"))
    #minuslog = sp("%s genomecov -i %s -g %s -bga -scale 1 > %s "%(conf_dict['General']['bedtools'],conf_dict['General']['outname'] + "_cleavage_minus.bed",conf_dict['options']['csize'], conf_dict['General']['outname']+ "_cleavage_minus.bdg"))
    #minuslog = sp("sort -k1,1 -k2,2n %s > %s"%(conf_dict['General']['outname']+ "_cleavage_minus.bdg",conf_dict['General']['outname']+ "_cleavage_minus_sorted.bdg" ))
    #minuslog = sp("%s %s %s %s"%(conf_dict['General']['bedGraphToBigWig'],conf_dict['General']['outname']+ "_cleavage_minus_sorted.bdg",conf_dict['options']['csize'],conf_dict['General']['outname']+ "_cleavage_minus.bw" ))
    pluslog = sp("sort -k1,1 -k2,2n %s > %s"%(conf_dict['General']['outname']+ "_cleavage_plus.bdg",conf_dict['General']['outname']+ "_cleavage_plus_sorted.bdg" ))
    pluslog = sp("%s %s %s %s"%(conf_dict['General']['bedGraphToBigWig'],conf_dict['General']['outname']+ "_cleavage_plus_sorted.bdg",conf_dict['options']['csize'],conf_dict['General']['outname']+ "_cleavage_plus.bw" ))
    minuslog = sp("sort -k1,1 -k2,2n %s > %s"%(conf_dict['General']['outname']+ "_cleavage_minus.bdg",conf_dict['General']['outname']+ "_cleavage_minus_sorted.bdg" ))
    minuslog = sp("%s %s %s %s"%(conf_dict['General']['bedGraphToBigWig'],conf_dict['General']['outname']+ "_cleavage_minus_sorted.bdg",conf_dict['options']['csize'],conf_dict['General']['outname']+ "_cleavage_minus.bw" ))

    pluslog = sp("sort -k1,1 -k2,2n %s > %s"%(conf_dict['General']['outname']+ "_biasExpCuts_plus.bdg",conf_dict['General']['outname']+ "_biasExpCuts_plus_sorted.bdg" ))
    pluslog = sp("%s %s %s %s"%(conf_dict['General']['bedGraphToBigWig'],conf_dict['General']['outname']+ "_biasExpCuts_plus_sorted.bdg",conf_dict['options']['csize'],conf_dict['General']['outname']+ "_biasExpCuts_plus.bw" ))
    minuslog = sp("sort -k1,1 -k2,2n %s > %s"%(conf_dict['General']['outname']+ "_biasExpCuts_minus.bdg",conf_dict['General']['outname']+ "_biasExpCuts_minus_sorted.bdg" ))
    minuslog = sp("%s %s %s %s"%(conf_dict['General']['bedGraphToBigWig'],conf_dict['General']['outname']+ "_biasExpCuts_minus_sorted.bdg",conf_dict['options']['csize'],conf_dict['General']['outname']+ "_biasExpCuts_minus.bw" ))
    return conf_dict

