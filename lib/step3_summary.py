#!/usr/bin/env python

# ------------------------------------
# Python Modual
# ------------------------------------

import os
import sys
import time
import string

# --------------------------
# custom package
# --------------------------

### tool function
from HMRpipe.Utility      import (sp,
                                   pdf_name,
                                   raise_error,
                                   wlog,
                                   ewlog,
                                   rwlog,
                                   rlogonly,
                                   CMD,
                                   createDIR,
                                   textformat,
                                   strlatexformat)
# --------------------------
# main 
# --------------------------
def step3_summary(conf_dict,logfile):
    '''
    analysis part
    mainly Rscript
    dimentional reduction + clustering
    '''
    # start
    # create section for 
    
    wlog('collect results',logfile)
# Rscript analysis.r expmat outname coverGN highvarZ selectPCcutoff rdnumber maxKnum
    summarydir =  'summary/'
    createDIR(summarydir)
    sp("mv %s_NCsummary.txt %s"%(conf_dict['General']['outname'],summarydir))
    sp("mv %s_elnet_lambdaSelection.pdf %s"%(conf_dict['General']['outname'],summarydir))
    if os.path.isfile("%s_cofactor_HMsignal.pdf"%conf_dict['General']['outname']):
        sp("mv %s_cofactor_HMsignal.pdf %s"%(conf_dict['General']['outname'],summarydir))

    tmpresult = 'tmpResults/'
    createDIR(tmpresult)
    sp("mv %s_HMsig.bed %s"%(conf_dict['General']['outname'],tmpresult))
    sp("mv %s_peakov.bed %s"%(conf_dict['General']['outname'],tmpresult))
    sp("mv %s_cofactor_candidate_list.txt %s"%(conf_dict['General']['outname'],tmpresult))
    sp("mv %s_filterNC.txt %s"%(conf_dict['General']['outname'],tmpresult))

    if conf_dict['General']['mode'] != "binary":
        sp("mv %s_TFsig.bed %s"%(conf_dict['General']['outname'],tmpresult))

    wlog('generate summary documents',logfile)
    ### initiate 
    QCdoc = """\documentclass[11pt,a4paper]{article}
\\usepackage{tabularx}
\\usepackage[english]{babel}
\\usepackage{array}
\\usepackage{graphicx}
\\usepackage{color}
\DeclareGraphicsExtensions{.eps,.png,.pdf,.ps}
\\begin{document}
\\title{Summary reports of non-classical function detection of : %s}

\\vspace{-1cm}
\maketitle
\\tableofcontents
\\newpage
\\newpage
\section{Data description}
\\begin{quotation}
Table 1 mainly describes the input files, parameters and options.
\end{quotation}
\\begin{table}[h]
\\small
\caption{parameter description}\label{bstable}
\\begin{tabularx}{\\textwidth}{ |X|l| }

"""%(strlatexformat(conf_dict['General']['outname']))
    ### table1 prepare parameter
    NcoTF = len(conf_dict['General']['usefilename'])          
    QCdoc += """      
\hline
parameter & value  \\\\
\hline
output name & %s \\\\
\hline
HMRpeak(peak filename) & %s \\\\
\hline
mode & %s \\\\
\hline
HM signal(bw filename) & \\begin{tabular}[c]{@{}l@{}}%s\end{tabular}  \\\\
\hline
\#cofactor candidates & %s \\\\
\hline
options & value \\\\
\hline
extend size & %sbp \\\\
\hline
Alpha (Elastic net) & %s \\\\
\hline
Pvalue cutoff & %s \\\\
\hline
topN cofactors & %s \\\\
\hline
"""%(strlatexformat(conf_dict['General']['outname']),
     strlatexformat(conf_dict['General']['HMRpeak'].split("/")[-1]),
     conf_dict['General']['mode'],
     strlatexformat("\\\\ ".join(conf_dict['General']['signalname'])),
     str(NcoTF),
     str(conf_dict['options']['ext']),
     str(conf_dict['options']['Alpha']),
     str(conf_dict['options']['Pvalue']),
     str(conf_dict['options']['TopNcofactors'])
     )
    QCdoc += """
\end{tabularx}
\end{table}
"""
    ### cross validation in elastic net
    QCdoc += """
\\newpage
\\newpage
\section{ElasticNet co-factor selection}
In this step we use a feature selection (elastic-net. Zou, H. and Hastie T. (2005) to select potential co-factors which corresponded to the non-classical function. Below shows the cross-validation curve for the decison of lambda in elastic-net for each histone modification substrate.  
\\begin{figure}[h]
        \caption{cross-validation curve for lambda decision} \label{fig:profileunion}
        \setlength{\\abovecaptionskip}{0pt}
        \setlength{\\belowcaptionskip}{10pt}
        \centering
        {\includegraphics[width=0.8\\textwidth]{%s}}
\end{figure}
"""%(conf_dict['General']['outname']+"_elnet_lambdaSelection.pdf")

    inf_ncsummary = open("summary/"+conf_dict['General']['outname']+"_NCsummary.txt")
    line = inf_ncsummary.readline()
    if line.startswith("no non-classical function detected"):
        QCdoc += """
\\newpage
\\newpage
\section{potential co-factors corresponded to non-classical function}
No significant co-factor was detected, indicating that the non-classical function of the HMR was not exist or none of the existing factor candidates act as a co-factor of the non-classical function.
"""
    else:
        QCdoc += """
\\newpage
\\newpage
\section{potential co-factors corresponded to non-classical function}
In summary, %s factors were predicted to potentially act as a co-factor of the non-classical function. The top%s co-factors were listed.
\subsection{summary of co-factors}
\\begin{quotation}
The corresponded histone modification substrate (HMsubstrate), empirical P-value, R-square (ordered) and the number of non-classical (NC) sites for each potential co-factor were listed below. The empirical P-value was calculated based on the comparison of foreground (observed) R-square and background R-square (distribution of random R-square generated from the 1,000 permutations of co-binding events) for each potential co-factor. The non-classical (NC) sites were defined by lower HMsubstrate signal (using Otus' method) and co-binding events of each potential co-factor.
\end{quotation}
\\begin{table}[h]
\\small
\caption{cofactor summary}\label{bstable}
\\begin{tabular}{ |l|l|l|l|l| }
    
\hline
co-factor & HMsubstrate & Pval & Rsquare & NCsites \\\\
"""%( int(sp("wc -l tmpResults/%s_filterNC.txt"%(conf_dict['General']['outname']))[0].split()[0])-1,
      int(sp("wc -l summary/%s_NCsummary.txt"%(conf_dict['General']['outname']))[0].split()[0])-1)

        for line in inf_ncsummary:
            if line.startswith("TFname"):
                continue
            ll = line.split()
            this_doc = """\hline
%s & %s & %s & %s & %s \\\\
"""%(strlatexformat(ll[0]), strlatexformat(ll[1]) ,ll[2],round(float(ll[3]),3),ll[5])
            QCdoc += this_doc
        inf_ncsummary.close()
        QCdoc += """
\hline
\end{tabular}
\end{table}
\\newpage
\\newpage
\subsection{Boxplot of HM on non-classical and classic sites}
\\begin{quotation}
Boxplot was generated to compare the difference of the histone mark (HM) signal on either non-classical or classic sites(peak). The non-classical sites were defined by lower HM signal (using Otus' method) and co-binding events of each potential co-factor. The boxplot corresponded to top co-factors were displayed.  
\end{quotation}
\\begin{figure}[h]
        \caption{boxplot cofactor HMsignal} \label{fig:profileunion}
        \setlength{\\abovecaptionskip}{0pt}
        \setlength{\\belowcaptionskip}{10pt}
        \centering
        {\includegraphics[width=0.8\\textwidth]{%s}}
\end{figure}
"""%((conf_dict['General']['outname']+"_cofactor_HMsignal.pdf"))

    QCdoc += """
\\newpage
\\newpage
\section{Output list}
\\begin{quotation}
All the main output files were described in the following table
\end{quotation}
\\begin{table}[h]
\\small
\caption{output list}\label{bstable}
\\begin{tabular}{ |l|l| }
    
\hline
description & filename \\\\
\hline
summary table of non-classical (NC) function & summary/%s \\\\
\hline
summary report (this doc) & summary/%s \\\\
\hline
cobinding matrix on HMR peaks & tmpResults/%s \\\\
\hline
histone mark signal on HMR peaks & tmpResults/%s \\\\
\hline

\end{tabular}
\end{table} 
\end{document} 

"""%(strlatexformat(conf_dict['General']['outname']+"_NCsummary.txt"),
     strlatexformat(conf_dict['General']['outname']+"_summary.pdf"),
     strlatexformat(conf_dict['General']['outname']+"_peakov.bed"),
     strlatexformat(conf_dict['General']['outname']+"_HMsig.bed")
    )

    latexfile = conf_dict['General']['outname'] + '_summary.tex'

    outf = open(summarydir+latexfile,'w')
    outf.write(QCdoc)
    outf.close()
    cmd = "pdflatex %s"%(latexfile)
    cmd2 = 'cp %s ../'%(conf_dict['General']['outname'] + '_summary.pdf')
    if conf_dict['General']['latex'] == 1:
        wlog('pdflatex was detected in default PATH, generate summary report %s'%(conf_dict['General']['outname'] + '_summary.pdf'),logfile)
        os.chdir(summarydir)
        tmpobj = sp(cmd)
        tmpobj = sp(cmd)
        tmpobj = sp(cmd2)
        tmpobj = sp("rm %s_summary.aux"%conf_dict['General']['outname'])
        tmpobj = sp("rm %s_summary.log"%conf_dict['General']['outname'])
        tmpobj = sp("rm %s_summary.toc"%conf_dict['General']['outname'])

#        for files in os.listdir(plot_folder):
#            if os.path.isfile(files) and files[-12:-4] == "_summary":
#                if not files[-4:] in ['.tex','.pdf',',png','.txt']:
#                    cmd = "rm %s"%(files)
#                    rwlog(cmd,logfile)
    else:
        wlog('pdflatex was not detected in default PATH, generate summary report .tex file in summary/ folder, you can move the whole summary/ folder to the environment with pdflatex installed and run cmd in the summary/ folder: "pdflatex %s"'%(conf_dict['General']['outname'] + '_summary.tex'),logfile)
   

    #if conf_dict['clean']:
    #    wlog('--clean pararmeter was turned on, remove internal files with large size',logfile)
    #    rwlog("rm %s "%(conf_dict['General']['outputdirectory'] + 'expmatrix/' + conf_dict['General']['outname']+'_on_symbol.bed'),logfile)
    #    rwlog("rm %s "%(conf_dict['General']['outputdirectory'] + 'expmatrix/' + conf_dict['General']['outname']+'_on_cds.bed'),logfile)
    #    rwlog("rm %s "%(conf_dict['General']['outputdirectory'] + 'expmatrix/' + conf_dict['General']['outname']+'_on_3utr.bed'),logfile)
    #    rwlog("rm %s "%(conf_dict['General']['outputdirectory'] + 'expmatrix/' + conf_dict['General']['outname']+'_on_5utr.bed'),logfile)
    #    rwlog("rm %s "%(conf_dict['General']['outputdirectory'] + 'expmatrix/' + conf_dict['General']['outname']+'_on_TTSdis.bed'),logfile)
    #    rwlog("rm %s "%(conf_dict['General']['outputdirectory'] + 'expmatrix/' + conf_dict['General']['outname']+'_combined.bed'),logfile)
    #    rwlog("rm %s "%(conf_dict['General']['outputdirectory'] + 'expmatrix/' + conf_dict['General']['outname']+'_barcode_reform.txt'),logfile)
#
    os.chdir("../")
    wlog('Step3 summary DONE, check %s for final outputs'%(summarydir),logfile)

    return conf_dict









