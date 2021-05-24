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
                                   split_chromosome_reads,
                                   filter_highQcell_reads)

# --------------------------
# main 
# --------------------------
def step1_QC_format(conf_dict,logfile):

    ### preparing mapping state dict
    wlog('summarize reads count distribution',logfile)
    chrom_reads = split_chromosome_reads(conf_dict['General']['fragments'],conf_dict['General']['outname'],conf_dict['options']['scATAC10x'],conf_dict['options']['chromosome'])
    if "chrM" in chrom_reads:
        conf_dict['QC']["chrM_reads"] = chrom_reads["chrM"]
    else:
        conf_dict['QC']["chrM_reads"] = 0
    chromatin_reads = 0
    for chrom in chrom_reads.keys():
        if chrom != "chrM":
            chromatin_reads += chrom_reads[chrom]
    conf_dict['QC']["chromatin_reads"] = chromatin_reads

    if conf_dict['General']['mode'] == "sc":
        wlog('filter high quality single cells',logfile)
        filter_highQcell_results = filter_highQcell_reads(conf_dict['General']['outname'],int(conf_dict['options']['readcutoff']),conf_dict['options']["usecells"])
        if filter_highQcell_results == "fail":
            ewlog('obtain < 100 high quality cell with reads >= %s.'%(conf_dict['options']['readcutoff']),logfile)
        if len(conf_dict['options']["usecells"]) == 0:
            wlog('no specified cellname list inputed',logfile)
        elif filter_highQcell_results[0] < 100 :
            wlog('obtain < 100 cell left after highQ + cellname filtering, use highQ cell only',logfile)
        wlog('obtain %s cells from filtering, containing %s reads'%(filter_highQcell_results[1],filter_highQcell_results[2]),logfile)
        conf_dict['QC']['highQcellnum'] = filter_highQcell_results[1]
        conf_dict['QC']['highQreadnum'] = filter_highQcell_results[2]
        conf_dict['QC']['totalcellnum'] = filter_highQcell_results[3]

    return conf_dict

