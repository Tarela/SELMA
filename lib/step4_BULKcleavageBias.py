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
from Utility      import (sp,
                                   wlog,
                                   ewlog,
                                   rwlog,
                                   CMD,
                                   pileup_cleavage,
                                   bias_exp_cleavage)

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
        cmdplus = """awk '{if($6=="+") print $0}' %s > %s"""%(conf_dict['General']['outname']+"_chromatin.bed", conf_dict['General']['outname']+"_cleavage_plus.bed")
        cmdminus = """awk '{if($6=="-") print $0}' %s > %s"""%(conf_dict['General']['outname']+"_chromatin.bed", conf_dict['General']['outname']+"_cleavage_minus.bed")

    splitlog = sp(cmdplus)
    splitlog = sp(cmdminus)

    wlog('pile up cleavage sites',logfile)
    pileup_cleavage(conf_dict['General']['outname'],conf_dict['General']['bedGraphToBigWig'],conf_dict['options']['csize'])

    wlog('calculate bias expected cleavages',logfile)
    bias_exp_cleavage(conf_dict['General']['outname'],
                      conf_dict['results']['peakfile'],
                      conf_dict['results']['biasMat'],
                      conf_dict['options']['kmer'],
                      conf_dict['General']['bigWigSummary'],
                      conf_dict['General']['bedGraphToBigWig'],
                      conf_dict['General']['twoBitToFa'],
                      conf_dict['General']['sequence'])
#
#
#
#    if "chrM" in chrom_reads:
#        conf_dict['QC']["chrM_reads"] = chrom_reads["chrM"]
#    else:
#        conf_dict['QC']["chrM_reads"] = 0
#    chromatin_reads = 0
#    for chrom in chrom_reads.keys():
#        if chrom != "chrM":
#            chromatin_reads += chrom_reads[chrom]
#    conf_dict['QC']["chromatin_reads"] = chromatin_reads
#
#    if conf_dict['General']['mode'] == "sc":
#        wlog('filter high quality single cells',logfile)
#        filter_highQcell_results = filter_highQcell_reads(conf_dict['General']['outname'],int(conf_dict['options']['readcutoff']),conf_dict['options']["usecells"])
#        if filter_highQcell_results == "fail":
#            ewlog('obtain < 100 high quality cell with reads >= %s.'%(conf_dict['options']['readcutoff']),logfile)
#        if len(conf_dict['options']["usecells"]) == 0:
#            wlog('no specified cellname list inputed',logfile)
#        elif filter_highQcell_results[0] < 100 :
#            wlog('obtain < 100 cell left after highQ + cellname filtering, use highQ cell only',logfile)
#        wlog('obtain %s cells from filtering, containing %s reads'%(filter_highQcell_results[1],filter_highQcell_results[2]),logfile)
#        conf_dict['QC']['highQcellnum'] = filter_highQcell_results[1]
#        conf_dict['QC']['highQreadnum'] = filter_highQcell_results[2]
#        conf_dict['QC']['totalcellnum'] = filter_highQcell_results[3]
#
    return conf_dict

