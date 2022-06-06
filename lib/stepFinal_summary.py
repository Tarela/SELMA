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
                                   createDIR,
                                   textformat,
                                   strlatexformat)

# --------------------------
# main 
# --------------------------
def stepFinal_summary(conf_dict,logfile):
    wlog('Collect results',logfile)
    summarydir =  'summary/'
    createDIR(summarydir)
    if "biasfile" in conf_dict['results'] and os.path.isfile(conf_dict['results']['biasfile']):
        sp("mv %s %s"%(conf_dict['results']['biasfile'],summarydir))        
    sp("mv %s_summitEXT.bed %s"%(conf_dict['General']['outname'],summarydir))

    conf_dict['results']['umap']="NA"
    if conf_dict['General']['mode'] == "bulk":
        sp("mv %s_cleavage_plus.bw %s"%(conf_dict['General']['outname'],summarydir))
        sp("mv %s_cleavage_minus.bw %s"%(conf_dict['General']['outname'],summarydir))
        sp("mv %s_biasExpCuts_plus.bw %s"%(conf_dict['General']['outname'],summarydir))
        sp("mv %s_biasExpCuts_minus.bw %s"%(conf_dict['General']['outname'],summarydir))
    else:
        sp("gzip %s_peakXcellMat.txt"%(conf_dict['General']['outname']))
        sp("gzip %s_peakFeatures.txt"%(conf_dict['General']['outname']))
        sp("mv %s_peakXcellMat.txt.gz %s"%(conf_dict['General']['outname'],summarydir))
        sp("mv %s_peakFeatures.txt.gz %s"%(conf_dict['General']['outname'],summarydir))
        sp("mv %s_scClusters.txt %s"%(conf_dict['General']['outname'],summarydir))

        if os.path.isfile("%s_clusteringUMAP.pdf"%(conf_dict['General']['outname'])):
            if conf_dict['options']['clustermethod'] in ["APEC","Cicero"]:
                sp("mv %s_clusteringUMAP.pdf %s/%s_clusteringTSNE.pdf"%(conf_dict['General']['outname'],summarydir,conf_dict['General']['outname']))
                conf_dict['results']['umap'] = "%s_clusteringTSNE.pdf"%(conf_dict['General']['outname'])
            else:
                sp("mv %s_clusteringUMAP.pdf %s"%(conf_dict['General']['outname'],summarydir))
                conf_dict['results']['umap'] = "%s_clusteringUMAP.pdf"%(conf_dict['General']['outname'])
    tmpresult = 'tmpResults/'
    createDIR(tmpresult)
    sp("mv %s_chromatin.bed %s"%(conf_dict['General']['outname'],tmpresult))
    sp("mv %s_chrM.bed %s"%(conf_dict['General']['outname'],tmpresult))
    sp("mv %s_summits.bed %s"%(conf_dict['General']['outname'],tmpresult))
    sp("mv %s_peaks.xls %s"%(conf_dict['General']['outname'],tmpresult))
    sp("mv %s_peaks.narrowPeak %s"%(conf_dict['General']['outname'],tmpresult))
    sp("mv %s %s"%(conf_dict['options']['csize'], tmpresult))

    if conf_dict['General']['mode'] == "bulk":    
        sp("mv %s_cleavage_plus.bed %s"%(conf_dict['General']['outname'],tmpresult))
        sp("mv %s_cleavage_minus.bed %s"%(conf_dict['General']['outname'],tmpresult))
        sp("mv %s_cleavage_plus.bdg %s"%(conf_dict['General']['outname'],tmpresult))
        sp("mv %s_cleavage_minus.bdg %s"%(conf_dict['General']['outname'],tmpresult))
        sp("mv %s_cleavage_plus_sorted.bdg %s"%(conf_dict['General']['outname'],tmpresult))
        sp("mv %s_cleavage_minus_sorted.bdg %s"%(conf_dict['General']['outname'],tmpresult))
        sp("mv %s_biasExpCuts_plus.bdg %s"%(conf_dict['General']['outname'],tmpresult))
        sp("mv %s_biasExpCuts_minus.bdg %s"%(conf_dict['General']['outname'],tmpresult))
        sp("mv %s_biasExpCuts_plus_sorted.bdg %s"%(conf_dict['General']['outname'],tmpresult))
        sp("mv %s_biasExpCuts_minus_sorted.bdg %s"%(conf_dict['General']['outname'],tmpresult))
        sp("mv %s_summitEXTmerge.bed %s"%(conf_dict['General']['outname'],tmpresult))
        for chrm in conf_dict['options']['chromosome']:
            if os.path.isfile("%s_mergePeaks.bed"%chrm):
                sp("mv %s_mergePeaks.bed %s"%(chrm,tmpresult))
                sp("mv %s_plusCuts.bed %s"%(chrm,tmpresult))
                sp("mv %s_minusCuts.bed %s"%(chrm,tmpresult))
                sp("mv %s_plusCutsOnPeak.bed %s"%(chrm,tmpresult))
                sp("mv %s_minusCutsOnPeak.bed %s"%(chrm,tmpresult))
    else:
        sp("mv %s_highQcellReads.bed %s"%(conf_dict['General']['outname'],tmpresult))
        sp("mv %s_tmpSCreads.bed %s"%(conf_dict['General']['outname'],tmpresult))
        sp("mv %s_tmpSCpeaks.bed %s"%(conf_dict['General']['outname'],tmpresult))
        sp("mv %s_scOVcleavage.bed %s"%(conf_dict['General']['outname'],tmpresult))
        sp("mv %s_scRscript.r %s"%(conf_dict['General']['outname'],tmpresult))
        if conf_dict['options']['clustermethod'].upper() == "SEURAT" or conf_dict['options']['clustermethod'].upper() == "SCRAN": 
            sp("mv %s_ArchRReads.bed.gz %s"%(conf_dict['General']['outname'],tmpresult))
            sp("mv %s_ArchRReads.bed.gz.tbi %s"%(conf_dict['General']['outname'],tmpresult))
            sp("mv %s_ArchR %s"%(conf_dict['General']['outname'],tmpresult))
        if conf_dict['options']['clustermethod'].upper() == "APEC":
            sp("mv %s_APEC %s"%(conf_dict['General']['outname'],tmpresult))

        if conf_dict['options']['clustermethod'] == "ArchR":
            sp("mv %s_ArchR %s"%(conf_dict['General']['outname'],tmpresult))
        if conf_dict['options']['clustermethod'] == "APEC":
            sp("mv %s_APEC %s"%(conf_dict['General']['outname'],tmpresult))
        if conf_dict['options']['clustermethod'] == "Cicero":
            sp("mv %s_Cicero %s"%(conf_dict['General']['outname'],tmpresult))

        if conf_dict['options']['SCcorrection']:
            sp("mv %s_correctionMatRscript.r %s"%(conf_dict['General']['outname'],tmpresult))
            sp("mv %s_correctionMat.txt %s"%(conf_dict['General']['outname'],tmpresult))
            sp("mv %s_correctionMatRDreads.bed %s"%(conf_dict['General']['outname'],tmpresult))
            sp("mv %s_correctionMatPeaks.txt %s"%(conf_dict['General']['outname'],tmpresult))


    if conf_dict['options']['keeptmp']:
        wlog('--keeptmp was set, keep intermediate results',logfile)
        pass
    else:
        wlog('--keeptmp was not set, remove intermediate results',logfile)
        sp("rm -r tmpResults/")

    wlog('Generate summary reports',logfile)
    outf = open("%s_summaryReports.txt"%conf_dict['General']['outname'],'w')
    outf.write("#settings\n")
    outf.write("mode\t%s\n"%(conf_dict['General']['mode']))
    outf.write("fragments\t%s\n"%(conf_dict['General']['fragments']))
    outf.write("data format\t%s\n"%(conf_dict['General']['format']))
    outf.write("data type\t%s\n"%(conf_dict['General']['datatype']))
    outf.write("genome version\t%s\n"%(conf_dict['General']['genome']))
    outf.write("output name\t%s\n"%(conf_dict['General']['outname']))

    outf.write("\n#parameters\n")
    outf.write("peak extend size\t%s\n"%(conf_dict['options']['extend']))
    outf.write("peak qvalue\t%s\n"%(conf_dict['options']['peakqval']))
    outf.write("bias source\t%s\n"%(conf_dict['options']['bias']))
    outf.write("k-mer\t%s\n"%(conf_dict['options']['kmer']))
    if conf_dict['General']['mode'] == "sc":
        outf.write("[sc]reads cutoff\t%s\n"%(conf_dict['options']['readcutoff']))
        outf.write("[sc]%low biaspeak\t"+str(conf_dict['options']['lowbiaspeak'])+"\n")
        outf.write("[sc]peak min reads\t%s\n"%(conf_dict['options']['peakminreads']))
        outf.write("[sc]peak max reads\t%s\n"%(conf_dict['options']['peakmaxreads']))
        outf.write("[sc]topN dimensions\t%s\n"%(conf_dict['options']['topDim']))
        outf.write("[sc]clustering method\t%s\n"%(conf_dict['options']['clustermethod']))
#        if conf_dict['options']['clustermethod'] == "PCAkm":
#            outf.write("[sc]cluster number\t%s\n"%(conf_dict['options']['clusterNum']))

    outf.write("\n#QC\n")
    outf.write("total reads\t%s\n"%(conf_dict['QC']['chrM_reads']+conf_dict['QC']['chromatin_reads']))
    outf.write("chromatin reads\t%s\n"%(conf_dict['QC']['chromatin_reads']))
    outf.write("mtDNA reads\t%s\n"%(conf_dict['QC']['chrM_reads']))
    outf.write("total peaks\t%s\n"%(conf_dict['QC']['peaknumTotal']))
    if conf_dict['General']['mode'] == "sc":
        outf.write("total single-cells\t%s\n"%(conf_dict['QC']['totalcellnum']))
        outf.write("high quality single-cells\t%s\n"%(conf_dict['QC']['highQcellnum']))
        outf.write("single-cells for clustering\t%s\n"%(conf_dict['QC']['finalusecellnum']))
        outf.write("reads in single-cells for clustering\t%s\n"%(conf_dict['QC']['finalreadnum']))
        if conf_dict['QC']['scClusters'] > 0:
	        outf.write("cluster number\t%s\n"%(conf_dict['QC']['scClusters']))

    outf.write("\n#output results\n")
    outf.write("peaks (accessible regions)\t%s_summitEXT.bed\n"%(conf_dict['General']['outname']))
    if conf_dict['General']['mode'] == "bulk":
        outf.write("observed cleavage (+ strand)\t%s_cleavage_plus.bw\n"%(conf_dict['General']['outname']))
        outf.write("observed cleavage (- strand)\t%s_cleavage_minus.bw\n"%(conf_dict['General']['outname']))
        outf.write("bias expected cleavage (+ strand)\t%s_biasExpCuts_plus.bw\n"%(conf_dict['General']['outname']))
        outf.write("bias expected cleavage (- strand)\t%s_biasExpCuts_minus.bw\n"%(conf_dict['General']['outname']))
    else:
        outf.write("peak bias features\t%s_peakFeatures.txt.gz\n"%(conf_dict['General']['outname']))
        outf.write("peakXcell count\t%s_peakXcellMat.txt.gz\n"%(conf_dict['General']['outname']))
        if os.path.isfile("%s%s_scClusters.txt"%(summarydir,conf_dict['General']['outname'])):
            outf.write("single-cell cluster\t%s_scClusters.txt\n"%(conf_dict['General']['outname']))
        if os.path.isfile(summarydir + conf_dict['results']['umap']):
            if conf_dict['options']['clustermethod'] in  ["APEC","Cicero"]:
                outf.write("sc-cluster t-SNE\t%s\n"%(conf_dict['results']['umap']))
            else:
                outf.write("sc-cluster UMAP\t%s\n"%(conf_dict['results']['umap']))
    outf.close()


    ### check pdflatex
    QCdoc = """\documentclass[11pt,a4paper]{article}
\\usepackage{tabularx}
\\usepackage[english]{babel}
\\usepackage{array}
\\usepackage{graphicx}
\\usepackage{color}
\DeclareGraphicsExtensions{.eps,.png,.pdf,.ps}
\\begin{document}
\\title{SELMA summary reports for: %s}

\\vspace{-1cm}
\maketitle
\\tableofcontents
\\newpage
\\newpage
\section{Summary description}
\\begin{quotation}
Table 1 describes the input files and settings.
\end{quotation}
\\begin{table}[h]
\\small
\caption{ settings }\label{bstable}
\\begin{tabularx}{\\textwidth}{ |X|l| }

"""%(strlatexformat(conf_dict['General']['outname']))
    ### table1 prepare parameter
    QCdoc += """      
\hline
parameter & value  \\\\
\hline
mode & %s \\\\
\hline
fragment file & %s \\\\
\hline
data format & %s \\\\
\hline
data type & %s  \\\\
\hline
genome version & %s  \\\\
\hline
output name & %s  \\\\
\hline
\end{tabularx}
\end{table}
"""%(strlatexformat(conf_dict['General']['mode']),
     strlatexformat(conf_dict['General']['fragments'].split("/")[-1]),
     strlatexformat(conf_dict['General']['format']),
     strlatexformat(conf_dict['General']['datatype']),
     strlatexformat(conf_dict['General']['genome']),
     strlatexformat(conf_dict['General']['outname']))

    QCdoc += """
\\newpage
\\newpage
\section{parameters and options}
\\begin{quotation}
Table 2 describes the parameters and options.
\end{quotation}
\\begin{table}[h]
\\small
\caption{parameters and options}\label{bstable}
\\begin{tabularx}{\\textwidth}{ |X|l| }
\hline
parameter & value  \\\\
\hline
peak extend size & %s \\\\
\hline
peak qvalue & %s \\\\
\hline
bias source & %s \\\\
\hline
k-mer & %s  \\\\
\hline
"""%(strlatexformat(conf_dict['options']['extend']),
     strlatexformat(conf_dict['options']['peakqval']),
     strlatexformat(conf_dict['options']['bias']),
     strlatexformat(conf_dict['options']['kmer'])
     )
    if conf_dict['General']['mode'] == "sc":
        QCdoc += """
[sc]reads cutoff & %s  \\\\
\hline
[sc]peak min reads & %s  \\\\
\hline
[sc]peak max reads & %s  \\\\
\hline
[sc]topN dimensions & %s  \\\\
\hline
[sc]cluster methods & %s  \\\\
\hline
"""%(strlatexformat(conf_dict['options']['readcutoff']),
     strlatexformat(conf_dict['options']['peakminreads']),
     strlatexformat(conf_dict['options']['peakmaxreads']),
     strlatexformat(conf_dict['options']['topDim']),
     strlatexformat(conf_dict['options']['clustermethod']),
     )
        if conf_dict['options']['SCcorrection'] :
            QCdoc += """
[sc]bias correction & TRUE  \\\\
\hline
"""
        else:
            QCdoc += """
[sc]bias correction & FALSE  \\\\
\hline
[sc]Percent lowBias peak & %s  \\\\
\hline
"""%(conf_dict['options']['lowbiaspeak'])
       
    QCdoc += """
\end{tabularx}
\end{table}
"""

    QCdoc += """
\\newpage
\\newpage
\section{data quality}
\\begin{quotation}
Table 3 describes data Quality.
\end{quotation}
\\begin{table}[h]
\\small
\caption{data quality}\label{bstable}
\\begin{tabularx}{\\textwidth}{ |X|l| }
\hline
parameter & value  \\\\
\hline
total reads & %s \\\\
\hline
chromatin reads & %s \\\\
\hline
mtDNA reads & %s \\\\
\hline
total peaks & %s  \\\\
\hline
"""%(strlatexformat(conf_dict['QC']['chrM_reads']+conf_dict['QC']['chromatin_reads']),
     strlatexformat(conf_dict['QC']['chromatin_reads']),
     strlatexformat(conf_dict['QC']['chrM_reads']),
     strlatexformat(conf_dict['QC']['peaknumTotal'])
     )
    if conf_dict['General']['mode'] == "sc":
        QCdoc += """
[sc]total single-cells(sc) & %s  \\\\
\hline
[sc]high quality sc & %s  \\\\
\hline
[sc]sc for clustering & %s  \\\\
\hline
[sc]reads in sc for clustering & %s  \\\\
\hline
"""%(strlatexformat(conf_dict['QC']['totalcellnum']),
     strlatexformat(conf_dict['QC']['highQcellnum']),
     strlatexformat(conf_dict['QC']['finalusecellnum']),
     strlatexformat(conf_dict['QC']['finalreadnum'])
     )
        if conf_dict['QC']['scClusters'] > 0:
            QCdoc += """[sc]number of cluster & %s  \\\\
\hline           
"""%(strlatexformat(conf_dict['QC']['scClusters']))
    QCdoc += """
\end{tabularx}
\end{table}
"""

    QCdoc += """
\\newpage
\\newpage
\section{output results}
\\begin{quotation}
Table 3 describes output results (in the summary/ folder).
\end{quotation}
\\begin{table}[h]
\\small
\caption{output results}\label{bstable}
\\begin{tabularx}{\\textwidth}{ |X|l| }
\hline
parameter & value  \\\\
\hline
peaks (accessible regions) & %s \\\\
\hline
"""%(strlatexformat(conf_dict['General']['outname']+"_summitEXT.bed"))

    if conf_dict['General']['mode'] == "bulk":
        QCdoc += """
observed cleavage(+) & %s \\\\
\hline
observed cleavage(-) & %s  \\\\
\hline
bias expected cleavage(+) & %s \\\\
\hline
bias expected cleavage(-) & %s  \\\\
\hline
"""%(strlatexformat(conf_dict['General']['outname']+"_cleavage_plus.bw"),
     strlatexformat(conf_dict['General']['outname']+"_cleavage_minus.bw"),
     strlatexformat(conf_dict['General']['outname']+"_biasExpCuts_plus.bw"),
     strlatexformat(conf_dict['General']['outname']+"_biasExpCuts_minus.bw"))
    else:
        QCdoc += """
peak bias feature & %s \\\\
\hline
peakXcell count & %s  \\\\
\hline
"""%(strlatexformat(conf_dict['General']['outname']+"_peakFeatures.txt.gz"),
     strlatexformat(conf_dict['General']['outname']+"_peakXcellMat.txt.gz"))
        if os.path.isfile("%s%s_scClusters.txt"%(summarydir,conf_dict['General']['outname'])):
            QCdoc += """single-cell cluster & %s \\\\
\hline
"""%(strlatexformat(conf_dict['General']['outname']+"_scClusters.txt"))
        if os.path.isfile(summarydir + conf_dict['results']['umap']):
            if conf_dict['options']['clustermethod'] in  ["APEC","Cicero"]:
                dimRedTerm = "t-SNE"
            else:
                dimRedTerm = "UMAP"
            QCdoc += """sc-cluster %s & %s \\\\
\hline
"""%(dimRedTerm,strlatexformat(conf_dict['General']['outname']+"_clusteringUMAP.pdf"))
    QCdoc += """
\end{tabularx}
\end{table}
"""

    if os.path.isfile(summarydir + conf_dict['results']['umap']):
        if conf_dict['options']['clustermethod'] in  ["APEC","Cicero"]:
            dimRedTerm = "t-SNE"
        else:
            dimRedTerm = "UMAP"
        QCdoc += """
\\newpage
\\newpage
\section{%s scatter plot}
The 2-dim scatter plot represent the %s results. Each dot represents an individual cell and the color represents cluster labels  
\\begin{figure}[h]
        \caption{%s visualization of cells colored by clustering} \label{fig:profileunion}
        \setlength{\\abovecaptionskip}{0pt}
        \setlength{\\belowcaptionskip}{10pt}
        \centering
        {\includegraphics[width=0.8\\textwidth]{%s}}
\end{figure}
"""%(dimRedTerm, dimRedTerm,dimRedTerm, summarydir + conf_dict['results']['umap'])

    QCdoc += """
\end{document} 
"""
    latexfile = conf_dict['General']['outname'] + '_summaryReports.tex'

    outf = open(latexfile,'w')
    outf.write(QCdoc)
    outf.close()

    check_latex = sp('which pdflatex')
    if check_latex[0].decode("ascii") == "" :
        wlog('pdflatex was not installed, SELMA will not generate pdf version of summary report. Please copy the %s to an environment with pdflatex installed and complie the pdf file'%(conf_dict['General']['outname'] + '_summaryReports.tex'),logfile)
    else:
        cmd = "pdflatex %s"%(latexfile)
        tmpobj = sp(cmd)
        tmpobj = sp(cmd)
        #tmpobj = sp(cmd2)
        tmpobj = sp("rm %s_summaryReports.aux"%conf_dict['General']['outname'])
        tmpobj = sp("rm %s_summaryReports.log"%conf_dict['General']['outname'])
        tmpobj = sp("rm %s_summaryReports.toc"%conf_dict['General']['outname'])

    return conf_dict






















