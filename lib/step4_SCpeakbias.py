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
                                   fetchseq_2bit_chrom,
                                   bias_peakXcell_mat)

# -------------------------- 
# main 
# --------------------------
def step4_SCpeakbias(conf_dict,logfile):

    wlog('readin sequence from 2bit',logfile)
    seq_dict = {}
    inf = open(conf_dict['options']['csize']) 
    for line in inf:
        chrm = line.split()[0]
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
                              conf_dict['options']['peakminreads'],
                              conf_dict['options']['peakmaxreads']
                              )

    return conf_dict


















