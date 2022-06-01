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
                                   extsummit,
                                   extExternal)
# --------------------------
# main 
# --------------------------
def step3_callpeak(conf_dict,logfile):

    conf_dict['results']['peakfile'] = conf_dict['General']['outname']+"_summitEXT.bed"

    macs3callpeak = 1
    if conf_dict['options']['peak'] != "NA":
        conf_dict['QC']['peaknumTotal'] = extExternal(conf_dict['options']['peak'],conf_dict['results']['peakfile'],int(conf_dict['options']['extend']))
        if conf_dict['QC']['peaknumTotal'] < 1000:
            wlog("obtain < 1000 (%s) external inputted peaks, use macs3 to detect peaks"%conf_dict['QC']['peaknumTotal'],logfile)
            macs3callpeak = 1
        else:
            wlog("obtain %s peaks from (-p) inputted, extend peak to +/- %sbp from peak center"%(conf_dict['QC']['peaknumTotal'],conf_dict['options']['extend']),logfile)
            macs3callpeak = 0
        
    if macs3callpeak == 1:
        if conf_dict['General']['macs3'] == "NA":
            ewlog("macs3 was not installed. SELMA requires macs3 installed in the default path ($PATH) for peak calling",logfile)

        ### callpeak 
        if conf_dict['General']['genome'] == "hg38":
            gtag = "hs"
        else:
            gtag = "mm"
    
        if conf_dict['General']['format'] == "PE":
            macs3cmd = "macs3 callpeak -t %s -n %s -f BEDPE -g %s -q %s --keep-dup all"%(conf_dict['General']['outname']+"_chromatin.bed",conf_dict['General']['outname'],gtag,conf_dict['options']['peakqval'])
        else:
            macs3cmd = "macs3 callpeak -t %s -n %s -f BED -g %s -q %s --keep-dup all --nomodel --extsize 100"%(conf_dict['General']['outname']+"_chromatin.bed",conf_dict['General']['outname'],gtag,conf_dict['options']['peakqval'])
    
        wlog("peak calling with macs3: %s"%macs3cmd,logfile)
        peaklog = sp(macs3cmd)

        ### ext peak from summit
        wlog("extend peak summits to +/- %sbp"%conf_dict['options']['extend'],logfile)
        if not os.path.isfile(conf_dict['General']['outname']+"_summits.bed"):
            ewlog("no macs3 results detected, check whether macs3 was correctly installed.",logfile)
        conf_dict['QC']['peaknumTotal'] = extsummit(conf_dict['General']['outname']+"_summits.bed",conf_dict['results']['peakfile'],int(conf_dict['options']['extend']))
        if conf_dict['QC']['peaknumTotal'] < 1000:
            ewlog("obtain < 1000 (%s) peaks, SELMA terminated"%conf_dict['QC']['peaknumTotal'],logfile)
        else:
            wlog("obtain %s peaks"%conf_dict['QC']['peaknumTotal'],logfile)

    return conf_dict

