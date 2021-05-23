#!/usr/bin/env python
'''
Created on XXXX-XX-XX

@author: Tarela
'''
"""
Description:

"""

# ------------------------------------
# Python Modual
# ------------------------------------

import os,sys,re
from optparse import OptionParser
import logging,time
import string,numpy

# ------------------------------------
# error and warning
# ------------------------------------


# ------------------------------------
# Misc functions
# ------------------------------------
import subprocess
def sp(cmd):
    a=subprocess.Popen(cmd,stdout=subprocess.PIPE,shell='TRUE')
    ac = a.communicate()
    return ac
def bwsig_pattern(bwfile,chrm,start,end,points):
    cmd = 'bigWigSummary %s %s %s %s %s'%(bwfile,chrm,start,end,points)
    sigPat = sp(cmd)[0].strip().split("\t")
    return sigPat

def get_signal(inputfile,output,signalbw,extend,N,bwfolder):
    signalbw = signalbw.strip().strip(',').split(',')
    if not bwfolder:
        bwfolder = ""
    if not bwfolder.endswith('/'):
        bwfolder += '/'
    bwHs = []
    for sb in signalbw:
        if "/" in sb:#sb.startswith('/mnt/Storage'):
            bwHs.append(sb)
        else:
            bwHs.append(bwfolder + sb)
#    signalbw = signalbw.strip().strip(',').split(',')
    
    inf = open(inputfile)
    outf = open(output,'w')
    for line in inf:
        ll = line.split()
        if "_" in ll[0]:
            continue
        inputlen = len(ll)
#        print ll
        if ll[5] == "+":
            center = int(ll[1])
        else:
            center = int(ll[2])
        S = max(0,center - extend)
        E = center + extend
        for bwH in bwHs:#bwHandle:
         #   print ll[0],S,E,N
            result = bwsig_pattern(bwH,ll[0],S,E,N)
            if 'n/a' in result or len(result) != N:
                continue#print result
            if ll[5] == "+":
                ll.extend(result)
            else:
                ll.extend(result[::-1])
        if 'n/a' in ll or len(ll) != (inputlen + N*len(signalbw)):
            pass
        else: 
            outf.write("\t".join(map(str,ll))+"\n")
    inf.close()
    outf.close()


# ------------------------------------
# Main function
# ------------------------------------

def main():
    usage = "python %prog >.< "
    description = """>.<"""

    optparser = OptionParser(version="%prog 1",description=description,usage=usage,add_help_option=False)
    optparser.add_option("-h","--help",action="help",help="Show this help message and exit.")

#========major options=============
    optparser.add_option("-i","--input",dest="inputfile",type="str",default = "",
                         help="")
    optparser.add_option("-o","--output",dest="output",type="str",default = "",
                         help="")
    optparser.add_option("-w","--signalbw",dest="bw",type="str",default = "",
                         help="")
    optparser.add_option("--ext",dest="ext",type="int",default = "10000",
                         help="")
    optparser.add_option("-n","--number",dest="NUM",type="int",default = "200",
                         help="")
    optparser.add_option("--bwfolder",dest="bwfolder",type="str",
                         help="")
#========minor options=============

    (options,args) = optparser.parse_args()

    inputfile = options.inputfile
    output = options.output
    signalbw = options.bw
    bwfolder = options.bwfolder
    if not inputfile:
        optparser.print_help()
        sys.exit(1)
    
    get_signal(inputfile,output,signalbw,options.ext,options.NUM,bwfolder)


if __name__== '__main__':
    try:
        main()

    except KeyboardInterrupt:
        sys.stderr.write("User interrupt me ^_^ \n")
        sys.exit(0)

