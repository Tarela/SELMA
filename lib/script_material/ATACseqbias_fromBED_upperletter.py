#!/usr/bin/env python
#Time-stamp:<Tarela>
"""
Description:

"""

# ------------------------------------
# Python Modual
# ------------------------------------

import os,sys,re
from optparse import OptionParser
import logging
import string
import copy,time
import twobitreader
# ------------------------------------
# error and warning
# ------------------------------------


# ------------------------------------
# Misc functions
# ------------------------------------
def rev(seq):
    revseq = ""
    for i in seq[::-1]:
        if i == 'A':
            r = 'T'
        elif i == 'T':
            r = 'A'
        elif i == 'C':
            r = 'G'
        elif i == 'G':
            r = 'C'
        else:
            r=i#print i
        revseq += r
    return revseq

def make_nmer_dict(n):
    nmer_seq = {}
    bp = ['A','C','G','T']
    allseq = [0]*n
    allseq[0] = bp
    i=1
    while i < n:
        allseq[i] = []
        for previous_seq in allseq[i-1]:
            for add_bp in bp:
                new_seq = previous_seq + add_bp
                allseq[i].append(new_seq)
        i += 1
    for seq in allseq[n-1]:
  #      print seq
        nmer_seq[seq] = 0
    del allseq
    return nmer_seq
 
def SOBfetchSEQ(fullseq,seqlen):
    if seqlen == 11:
        outseq = fullseq[:3]+fullseq[4]+fullseq[8:10]+fullseq[11:15]+fullseq[16]
    elif seqlen == 10:
        outseq = fullseq[:3]+fullseq[4]+fullseq[8:10]+fullseq[11:14]+fullseq[16]
    elif seqlen == 8:
        outseq = fullseq[0]+fullseq[2]+fullseq[8:10]+fullseq[11:14]+fullseq[16]
    elif seqlen == 6:
        outseq = fullseq[0]+fullseq[2]+fullseq[8:10]+fullseq[12:14]
    elif seqlen == 4:
        outseq = fullseq[2]+fullseq[8]+fullseq[12:14]
    else:
        outseq = "NA"
    return outseq

def seqbias(peak,tag,out,sequence,kmer,biastype,shiftbp):
    flank = int(kmer/2)
    genome = twobitreader.TwoBitFile(sequence) 
    pcut = make_nmer_dict(kmer)
    #ncut = make_nmer_dict(2*flank)
    bgseq = make_nmer_dict(kmer)

    if biastype in ["flank","bagfoot"]:
        inf = open(peak)
        for line in inf:
            ll = line.strip().split("\t")
            seq = genome[ll[0]][int(ll[1]):int(ll[2])]
            for i in range(len(seq)):
                subseq = seq[(i-flank):(i+flank)]
                if bgseq.has_key(subseq):
                    bgseq[subseq]+=1
                else:
                    pass
        inf.close()
    elif biastype == "sob":
        inf = open(peak)
        for line in inf:
            ll = line.strip().split("\t")
            seq = genome[ll[0]][int(ll[1]):int(ll[2])]
            for i in range(len(seq)):
                subseq_template = seq[(i-6):(i+11)]
                if len(subseq_template) != 17:
                    continue
                subseq = SOBfetchSEQ(subseq_template,kmer)
                if bgseq.has_key(subseq):
                    bgseq[subseq]+=1
                else:
                    pass
        inf.close()
    else:
        print "biastype choose from (flank,sob,bagfoot), current biastype is %s"%(biastype)
        sys.exit()

    if shiftbp == 45:
        shiftbp_plus = 4
        shiftbp_minus = 5
    else:
        shiftbp_plus = shiftbp
        shiftbp_minus = shiftbp

    inf = open(tag)
    PEtag = 0
    for line in inf:
        ll = line.strip().split("\t")
        if "." in ll[0] or "_" in ll[0]:
            continue
        if ll[0] == "chrMT":
            chrm = 'chrM'
        else:
            chrm = ll[0]
        if len(ll) < 6:
            PEtag = 1

        if biastype == "flank":
            if PEtag == 1 or ll[5] == '+' or ll[5] == "." :
                seq = genome[chrm][(int(ll[1])-flank+shiftbp_plus):(int(ll[1])+flank+shiftbp_plus)]#.upper()
                if pcut.has_key(seq):
                    pcut[seq]+=1
                else:
                    pass
            if PEtag == 1 or ll[5] == '-' or ll[5] == ".":
                seq = rev(genome[chrm][(int(ll[2])-flank-shiftbp_minus):(int(ll[2])+flank-shiftbp_minus)])#.upper())
                if pcut.has_key(seq):
                    pcut[seq]+=1
                else:
                    pass
        elif biastype == "bagfoot":
            if PEtag == 1 or ll[5] == '+' or ll[5] == "." :
                seq = genome[chrm][(int(ll[1])-flank+shiftbp_plus):(int(ll[1])+flank+shiftbp_plus)]#.upper()
                if pcut.has_key(seq):
                    pcut[seq]+=1
                else:
                    pass
            if PEtag == 1 or ll[5] == '-' or ll[5] == ".":
                seq = genome[chrm][(int(ll[2])-flank-shiftbp_minus):(int(ll[2])+flank-shiftbp_minus)]#.upper())
                if pcut.has_key(seq):
                    pcut[seq]+=1
                else:
                    pass
        elif biastype == "sob":
            if PEtag == 1 or ll[5] == '+' or ll[5] == "." :
                rawseq = genome[ll[0]][(int(ll[1])-6):(int(ll[1])+11)]#.upper()
                if len(rawseq) != 17:
                    continue
                seq = SOBfetchSEQ(rawseq,kmer)
                if pcut.has_key(seq):
                    pcut[seq]+=1
                else:
                    pass
            if PEtag == 1 or ll[5] == '-' or ll[5] == ".":
                rawseq = rev(genome[ll[0]][(int(ll[2])-11):(int(ll[2])+6)])#.upper())
                if len(rawseq) != 17:
                    continue
                seq = SOBfetchSEQ(rawseq,kmer)
                if pcut.has_key(seq):
                    pcut[seq]+=1
                else:
                    pass            
        else:
            print "biastype choose from (flank,sob,bagfoot), current biastype is %s"%(biastype)
            sys.exit()

    inf.close()
    outf = open(out,'w')
    for seqtype in sorted(pcut.keys()):
        if bgseq[seqtype] == 0:
            pbias = -1
        else:  
            pbias = float(pcut[seqtype])/float(bgseq[seqtype])
        #nbias = float(ncut[seqtype])/float(bgseq[seqtype])
        #outf.write("\t".join(map(str,[seqtype,pcut[seqtype]]))+"\n")
        outf.write("\t".join(map(str,[seqtype,pbias,pcut[seqtype],bgseq[seqtype]]))+"\n")
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
    optparser.add_option("-p","--peak",dest="peak",type="str",
                         help="peak region for calculate background (k-mer occurancy)")
    optparser.add_option("-t","--tag",dest="tag",type="str",
                         help="6-column reads file (bed) for calculate foreground (cuts occyrancy)")              
    optparser.add_option("-o","--out",dest="out",type="str",
                         help="")              
    optparser.add_option("-s","--sequence",dest="sequence",type="str",default='/Data/Genome/hg19.2bit',
                         help="whole genome sequence in 2bit format")
    optparser.add_option("--biastype",dest="biastype",type="str",default='flank',
                         help="biastype, choose from [flank, sob, bagfooot ]")
    optparser.add_option("-k","--kmer",dest="kmer",type="int",default=6,
                         help="kmer-mer , default =6 means 6-mer bias, also used for sob method")
    optparser.add_option("--shift",dest="shift",type="int",default=0,
                         help="shift Nbp to downstream, 45 means 4for plus and 5for minus")
                         
#========minor options=============

    (options,args) = optparser.parse_args()

    peak = options.peak
    tag = options.tag
    out = options.out
    seq = options.sequence
    kmer = options.kmer
    if not tag or not peak:
        optparser.print_help()
        sys.exit(1)
    
    seqbias(peak,tag,out,seq,kmer,options.biastype, options.shift)

if __name__== '__main__':
    try:
        main()

    except KeyboardInterrupt:
        sys.stderr.write("User interrupt me ^_^ \n")
        sys.exit(0)



