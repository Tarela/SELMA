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
from SELMApipe.Utility      import (sp, 
                                   wlog,
                                   ewlog,
                                   rwlog,
                                   CMD,
                                   fetchseq_2bit_chrom,
                                   bias_peakXcell_mat)

# -------------------------- 
# main 
# --------------------------
def step4_SCpeakbias(conf_dict,logfile):

    wlog('readin sequence from 2bit',logfile)

    used_chrm_list= []
    inf = open(conf_dict['results']['peakfile'])
    for line in inf:
        ll = line.split()
        if not ll[0] in used_chrm_list:
            used_chrm_list.append(ll[0])
    inf.close()

    seq_dict = {}
    inf = open(conf_dict['options']['csize']) 
    for line in inf:
        chrm = line.split()[0]
        if chrm in used_chrm_list:
          seq_dict[chrm] = fetchseq_2bit_chrom(conf_dict['General']['twoBitToFa'],conf_dict['General']['sequence'],chrm)
    inf.close()
    conf_dict['results']['seqdict'] = seq_dict

    wlog('scan peak level bias',logfile)

    tmplog = bias_peakXcell_mat(conf_dict['General']['outname'], 
                              conf_dict['General']['bedtools'], 
                              conf_dict['options']['chromosome'],
                              conf_dict['options']['kmer'],
                              conf_dict['results']['biasMat'],
                              conf_dict['results']['seqdict'],
                              conf_dict['results']['finalcells'],
                              conf_dict['General']['datatype'],
                              conf_dict['options']['peakminreads']                              
                              )

    return conf_dict


















