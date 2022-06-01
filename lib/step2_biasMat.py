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
                                   readBias,
                                   naive_kmerBias_chrM)
from SELMApipe.SimplexEncoding import simplex_encoding
# --------------------------
# main 
# --------------------------
def step2_biasMat(conf_dict,logfile):

    ### obtain bias mat
    if conf_dict['options']['bias'] == "naked":
        wlog('obtain pre-processed bias matrix from naked DNA data',logfile)
        conf_dict['results']['biasMat'] = readBias(conf_dict['options']['biasfile'])
    elif conf_dict['QC']['chrM_reads'] < 500000:
        wlog('chrM reads number < 500k, obtain pre-processed bias matrix from naked DNA data',logfile)
        if not os.path.isfile(conf_dict['options']['biasfile']):
            ewlog("no naked DNA bias matrix, cannot estimate bias",logfile)
        else:
            conf_dict['results']['biasMat'] = readBias(conf_dict['options']['biasfile'])
    else:
        wlog('estimate bias matrix from mtDNA(chrM) data',logfile)
        conf_dict['results']['biasMatNaive'] = naive_kmerBias_chrM(conf_dict['General']['outname'],conf_dict['General']['sequence'],conf_dict['options']['kmer'],conf_dict['General']['twoBitToFa'],conf_dict['General']['format'])
        conf_dict['results']['biasfile'] = "%s_bias.txt"%(conf_dict['General']['outname'])
        conf_dict['results']['biasMat']= simplex_encoding(conf_dict['results']['biasMatNaive'],conf_dict['results']['biasfile'])

    return conf_dict

