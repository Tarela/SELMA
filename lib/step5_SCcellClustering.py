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
                                   CMD)
from SELMApipe.scClustering import (
                          #scClustering_Seurat,
                          scClustering_correction_matrix,
                          RDreadloci,
                          RDreads_correctionMat,
                          scClustering_PCAkm,
                          scClustering_ArchR,
                          scClustering_APEC
                          )
                          #scClustering_Cicero)
# --------------------------
# main 
# --------------------------
def step5_SCcellClustering(conf_dict,logfile):

    wlog('single-cell clustering analysis',logfile)
    if conf_dict['options']['SCcorrection']:
        SCcorrection = 1
    else:
        SCcorrection = 0  
    if conf_dict['options']['clustermethod'] == "Kmeans":
        conf_dict['General']['scPackage'] = scClustering_PCAkm(conf_dict['General']['outname'],
                                       conf_dict['options']['lowbiaspeak'],
                                       conf_dict['options']['clusterNum'],
                                       conf_dict['options']['topDim'],
                                       conf_dict['options']['peakmaxreads'],
                                       int(conf_dict['options']['UMAP']),
                                       SCcorrection
                                       )
        if conf_dict['General']['scPackage']  == "noPackage":
            wlog("umap was not installed, UMAP scatter plot will not be generated",logfile)

    #elif conf_dict['options']['clustermethod'] == "Seurat": 
    #    conf_dict['General']['scPackage'] = scClustering_Seurat(conf_dict['General']['outname'],
    #                       conf_dict['options']['lowbiaspeak'],
    #                       conf_dict['options']['topDim'],
    #                       int(conf_dict['options']['UMAP']))
    #    if conf_dict['General']['scPackage']  == "noPackage":
    #        wlog("Seurat related packages were not installed, skip single-cell clustering step",logfile)

    elif conf_dict['options']['clustermethod'].upper() in ["SEURAT","SCRAN"]:
        if SCcorrection == 1:
            scClustering_correction_matrix(conf_dict['General']['outname'],
                                           conf_dict['options']['peakmaxreads'])
            RDreads_correctionMat(conf_dict['General']['outname'])

        conf_dict['General']['scPackage'] = scClustering_ArchR(conf_dict['General']['outname'],
                                                               conf_dict['General']['genome'],
                                                               conf_dict['options']['lowbiaspeak'],
                                                               conf_dict['options']['topDim'],
                                                               conf_dict['options']['peakmaxreads'],
                                                               SCcorrection,
                                                               int(conf_dict['options']['UMAP']),
                                                               conf_dict['options']['clustermethod'].upper())

        if conf_dict['General']['scPackage']  == "noPackage":
            wlog("related package (ArchR) were not installed, skip single-cell clustering step",logfile)
        if conf_dict['General']['scPackage']  == "noTabix":
            wlog("related packages (Tabix/bgzip) were not installed, skip single-cell clustering step",logfile)

    elif conf_dict['options']['clustermethod'].upper() == "APEC": 
        if SCcorrection == 1:
            scClustering_correction_matrix(conf_dict['General']['outname'],
                                           conf_dict['options']['peakmaxreads'])

        conf_dict['General']['scPackage'] = scClustering_APEC(conf_dict['General']['outname'],
                                                              conf_dict['options']['lowbiaspeak'],
                                                              int(conf_dict['options']['UMAP']),
                                                              int(conf_dict['options']['peakmaxreads']),
                                                              SCcorrection)
        if conf_dict['General']['scPackage']  == "noPackage":
            wlog("related package (APEC) were not installed, skip single-cell clustering step",logfile)

    #elif conf_dict['options']['clustermethod'] == "Cicero": 
    #    conf_dict['General']['scPackage'] = scClustering_Cicero(conf_dict['General']['outname'],
    #                       conf_dict['options']['lowbiaspeak'],
    #                       int(conf_dict['options']['UMAP']))
    #    if conf_dict['General']['scPackage']  == "noPackage":
    #        wlog("Cicero related packages were not installed, skip single-cell clustering step",logfile)

    if os.path.isfile("%s_scClusters.txt"%(conf_dict['General']['outname'])):
        tmplist = []
        inf = open("%s_scClusters.txt"%(conf_dict['General']['outname']))
        for line in inf:
            ll = line.strip().split("\t")
            if ll[1] != "cluster":
                tmplist.append(ll[1])
        inf.close()
        conf_dict['QC']['scClusters'] = len(set(tmplist))
    else:
        conf_dict['QC']['scClusters'] = 0
    return conf_dict













