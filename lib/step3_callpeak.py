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
                                   extsummit)
# --------------------------
# main 
# --------------------------
def step3_callpeak(conf_dict,logfile):

    ### callpeak
    if conf_dict['General']['genome'] == "hg38":
        gtag = "hs"
    else:
        gtag = "mm"

    if conf_dict['General']['format'] == "PE":
        macs3cmd = "macs3 callpeak -t %s -n %s -f BEDPE -g %s -q %s --keep-dup 1"%(conf_dict['General']['outname']+"_chromatin.bed",conf_dict['General']['outname'],gtag,conf_dict['options']['peakqval'])
    else:
        macs3cmd = "macs3 callpeak -t %s -n %s -f BED -g %s -q %s --keep-dup 1 --nomodel --extsize 100"%(conf_dict['General']['outname']+"_chromatin.bed",conf_dict['General']['outname'],gtag,conf_dict['options']['peakqval'])

    wlog("peak calling with macs3: %s"%macs3cmd,logfile)
    peaklog = sp(macs3cmd)

    ### ext peak from summit
    wlog("extend peak summits to +/- %sbp"%conf_dict['options']['extend'],logfile)
    conf_dict['results']['peakfile'] = conf_dict['General']['outname']+"_summitEXT.bed"
    conf_dict['QC']['peaknumTotal'] = extsummit(conf_dict['General']['outname']+"_summits.bed",conf_dict['General']['outname']+"_summitEXT.bed",int(conf_dict['options']['extend']))
    if conf_dict['QC']['peaknumTotal'] < 1000:
        ewlog("obtain < 1000 (%s) peaks, SELMA terminated"%conf_dict['QC']['peaknumTotal'],logfile)
    else:
        wlog("obtain %s peaks"%conf_dict['QC']['peaknumTotal'],logfile)

    return conf_dict

