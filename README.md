# SELMA: a computational framework for modeling intrinsic biases in chromatin accessibility sequencing data

Genome-wide profiling of chromatin accessibility with the assay for transposase-accessible chromatin using sequencing (ATAC-seq) or DNaseI hypersensitivity sequencing (DNase-seq) has been widely used for studying regulatory DNA elements and transcriptional regulation in many cellular systems. Efficient and thorough computational analysis is essential for extracting biological information from such high-throughput sequencing data. It has been reported that DNase cleavage of DNA has sequence preferences that can significantly affect the footprint patterns at transcription factor binding sites in genomic profiles. We found that enzymatic sequence biases commonly exist in both bulk and single-cell chromatin accessibility profiling data. Using a regular simplex encoding model, we developed a quantitative approach for accurate characterization and systematic correction of intrinsic sequence biases contained in ATAC-seq and DNase-seq data. This approach can be applied in bioinformatics for improved analysis of high-throughput chromatin accessibility sequencing.

## 0. Introduction of SELMA package
SELMA performs estimation and correction of intrinsic cleavage bias of DNaseI(DNase-seq) and Tn5(ATAC-seq) data in both bulk level and single cell level. SELMA used DNase/ATAC-seq data from either naked DNA or mtDNA to estimate the intrinsic cleavage bias, and improve the bias estimation using a simplex encoding model. SELMA provides a series of bias free analysis for the bulk/sc DNase/ATAC-seq data. For bulk data, SELMA estimates the bias expected cleavages on chromatin accessibility regions (peaks) and compares with observed cleavages. For single cell data, SELMA estimates the summarized bias score on each potential chromatin accessibility regions (peaks) and only use those peaks with less bias effect for single-cell clustering analysis. 

- Changelog<br>
v1.0.1 update (-p) option for customized peak files.<br>
v1.0.0 First version of SELMA with both sc and bulk mode.

## 1. Installation
- Package requirements<br>
SELMA requires [python](https://www.python.org) 3.6+ and [Rscript](https://www.r-project.org) v3+ to run.<br>
SELMA requires python3 packages [numpy](https://numpy.org) pre-installed.

\# for root user
```sh
$ cd SELMA
$ sudo python setup.py install  
```
\# if you are not root user, you can install SELMA at a specific location which you have write permission
```sh
$ python setup.py install --prefix /home/SELMA  # here you can replace “/home/SELMA” with any location 
$ export PATH=/home/SELMA/bin:$PATH    # setup PATH for the software
$ export PYTHONPATH=/home/SELMA/lib/python3.6/site-packages:$PYTHONPATH    # setup PYTHONPATH for module import
```
\# To check the SELMA package, just type:
```sh
$ SELMA --help  # If you see help manual, you have successfully installed SELMA
```

\# NOTE: 
- To install SELMA on MacOS, the users need to download and install Command Line Tools beforehand
- SELMA requires python3 packages [numpy](https://numpy.org) and [macs3](https://pypi.org/project/MACS3/) pre-installed. See section 4 if the users don't have macs3 or prefer other peak calling methods or provide customized peak file. 
- SELMA suggests to have bedtools (Quinlan et al., Bioinformatics. 2010) and UCSC tools (Kuhn et al., Brief Bioinform. 2013) pre-installed for data pre-processing. The SELMA package also wrapped up both tools and would install automatically if the users did not have them pre-installed in the default PATH. 
- SELMA suggests to have pdflatex installed for the summary pdf document. To install pdflatex on macOS, you can download “MacTex” from http://tug.org/cgi-bin/mactex-download/MacTeX.pkg. To install pdflatex on linux, you need to install the texlive package from https://www.tug.org/texlive/. SELMA will generate a .txt file as a simple version of summary report without pdflatex installed. A .tex file will also be generated in case the users want to make the pdf doucment later. 
- Some function (single cell clustering) of SELMA requires the related packages pre-installed (see seciton 4)

## 2. Download genome sequence database
SELMA requires genome sequence (in .2bit format) prepared for running. You can download them from the UCSC genome browser or other public domains. 
- [hg38.2bit](https://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/hg38.2bit)
- [mm10.2bit](https://hgdownload.cse.ucsc.edu/goldenpath/mm10/bigZips/mm10.2bit)


## 3. Run SELMA (usage)
#### Essential paramters
To run SELMA with default parameters, you only need to give SELMA:
-   -m MODE, --mode=MODE
Mode of SELMA, choose from sc(single-cell) or bulk
-   -i FRAGMENTs, --input_fragments=FRAGMENTs
Input fragments file in bed format, with .bed extension, for sc mode, the 4th(name) column of bed file represents the name of the corresponded individual cell
-   -f FORMAT, --format=FORMAT
Format of the fragments.bed file. choose from PE(paired-end, default) or SE(single-end)
-   -t DATATYPE, --datatype=DATATYPE
Type of sequencing data (experiments), choose from ATAC or DNase
-   -g GENOME, --genome=GENOME
genome version of the input data, default is hg38(default). Currently SELMA only support hg38 and mm10 for de novo peak detection. External peakset input (-p) is required for other genome version/speces. 
-   -s SEQUENCE, --sequence=SEQUENCE
genome sequence file in 2bit format (e.g., hg38.2bit).
-   -o OUTNAME, --outname=OUTNAME
Name of output results
-   -p PEAK, --peak=PEAK  
[optional] external peak file for the candidate chromatin accessibility regions. 

Example for run SELMA with all default parameters:

\# sc mode 
```sh
$ SELMA -m sc -i ${path}/testdata.bed.gz -g hg38 -f PE -o testsc -t ATAC --clusterMethod PCAkm -s ${path}/hg38.2bit --readCutoff 1000 --bias naked --kmer 10 --UMAP 
```

\# bulk mode 
```sh
$ SELMA -m bulk -i ${path}/testdata.bed.gz -g hg38 -f PE -o testbulk -t ATAC -s ${path}/hg38.2bit --bias naked --kmer 10 
```

## 4. Customize target chromatin accessibility (peak) regions (v1.0.1).
SELMA provides an option (-p) for the users to use their customized peak files as the target chromatin accessibility regions for the SELMA analysis. The required peak file should be in .bed format (plain text), have >=4 columns (chrom,start,end,name), and contain >=1000 peaks. This parameter was specifically designed for the users who don't like macs3 or have problems in installing macs3. However, we strongly suggested to use the same dataset (e.g., fragments.bed file) for the peak calling (with whatever methods the users prefered) to ensure sufficient cleavages/signal on the peak regions. Below is the example with external/customized peak file:

\# sc mode 
```sh
$ SELMA -m sc -i ${path}/testdata.bed.gz -p ${path}/testpeak.bed -g hg38 -f PE -o testsc -t ATAC --clusterMethod PCAkm -s ${path}/hg38.2bit --readCutoff 1000 --bias naked --kmer 10 --UMAP 
```

\# bulk mode 
```sh
$ SELMA -m bulk -i ${path}/testdata.bed.gz -p ${path}/testpeak.bed -g hg38 -f PE -o testbulk -t ATAC -s ${path}/hg38.2bit --bias naked --kmer 10 
```

## 5. Pre-processing steps for generating the input fragments file.
SELMA takes aligned fragment files (in .bed format) as input. Thus users can perform any pre-processing steps to customize the fragments files. For example, only use high quality reads with perfect alignemnt (e.g., MAPQ > 30) to run SELMA. For bulk data, only using unique paired-end fragments (unique loci) would reduce the influence from the over amplificaiton. For single-cell data, keep unique fragments in each individual cell level is also a considerable option. 


## 6. Install and use published single cell clustering methods based on SELMA bias correction. 
SELMA sc mode implements several well acknowledged cell clustering methods in the single cell clustering analysis in additional to the default PCA+Kmeans analysis. To activate these alternative methods (name, version and link listed below), users need to install the related package, and specify the method by the --clusterMethod parameter. If any methods were specified by the --clusterMethod parameter but with no related package installed, SELMA will skip the single-cell clustering analysis. 
- [Seurat v4.0.0](https://satijalab.org/seurat/)
- [scran v1.20.1](https://bioconductor.org/packages/release/bioc/html/scran.html)
- [APEC v1.2.2](https://github.com/QuKunLab/APEC)

Since Seurat and scran were developed for scRNA-seq. We used the [ArchR v1.0.1](https://www.archrproject.com) version of the two methods for scATAC-seq data. 

SELMA also provides UMAP/t-SNE visualization for the single-cell clustering analysis. You can activate this function by the --UMAP parameter. For the PCAkm method, the [umap](https://cran.r-project.org/web/packages/umap/index.html) package in R is required. 

## 7. Processed data generated for this study
The following data were generated and used in SELMA study, but they would not be generated with SELMA pipeline automatically.
- SELMA bias for [DNaseI(DNase-seq)](https://www.dropbox.com/s/ncemdhp0cee3cic/DNase_SELMAbias_10mer.txt.gz) ([backupLink](https://data.cyverse.org/dav-anon/iplant/home/tarela/SELMA/DNase_SELMAbias_10mer.txt.gz)) and [Tn5(ATAC-seq)](https://www.dropbox.com/s/x5iiy27ef80fl19/ATAC_SELMAbias_10mer.txt.gz) ([backupLink](https://data.cyverse.org/dav-anon/iplant/home/tarela/SELMA/ATAC_SELMAbias_10mer.txt.gz)). SELMA recommended and generated 10-mer simplex encoded intrinsic cleavage bias for DNaseI and Tn5. Both bias matrix were genearted by naked DNA data DNase/ATAC-seq data. Note that the bias matrix was also built-in in the SELMA package and would also be used as reference data in SELMA pipeline.  
- [Footprint bias score for consensus footprint regions](https://www.dropbox.com/s/f3m9q0fhlq4e9vc/consensusFP_biasScore.txt.gz) ([backupLink](https://data.cyverse.org/dav-anon/iplant/home/tarela/SELMA/consensusFP_biasScore.txt.gz)). The genome-wide consensus footprint regions were generated in the previous study (Vierstra et al., Nature. 2020) and were downloaded from the [public domain](http://vierstra.org/resources/dgf). We estimated a SELMA bias score for each of the consensus footprint region (last column of the file). 

## 8. Output files
1. `NAME_summaryReports.pdf` is the summary pdf file which contains information of:
     - Input file and parameter description
     - basic QC of the data
     - Summary of the SELMA bias estimation/correction results

    \#Note: This pdf file is generated only if pdflatex is pre-installed. A NAME_summaryReports.txt file is generated as well for the simple version reports. A .tex file will also be generated in case the users want to make the pdf doucment later.

2. `NAME_peaks.bed` is the peaks detected from the fragment files (with macs3). Each peak was extended to 400bp centered on the peak summit. 

3. `NAME_cleavage.bw` (bulk mode only) is the genome-wide profile of the 1bp cleavage of DNaseI/Tn5. Plus and minus strand cleavages will be separated to two files (cleavage_plus.bw, cleavage_minus.bw)

4. `NAME_biasExpCuts.bw` is the profile of the bias expected cleavage on the peak regions. Plus and minus strand cleavages will be separated to two files (biasExpCuts_plus.bw, biasExpCuts_minus.bw)

5. `NAME_peakXcell.txt.gz` (sc mode only) is the peak X cell count matrix generated from the single cell analysis. The cells were filtered by the total reads count (default >=10k reads) and the peaks were filtered based on the intrinsic cleavage bias (default top2/3 peaks with lowest bias effect). 

6. `NAME_scClustering.txt.gz` (sc mode only) is the cell clustering results using SELMA debiased peakset.

7. `NAME_bias.txt` is the bias matrix estimated with SELMA methods. This file will only be generated if the users don't use the default (--bias naked) parameter (i.e. set --bias chrM to use mtDNA reads to estimate bias instead)

## 9. Testing data and example of output files
We provided the testing data for users to test the flexibility and the power of the SELMA. The sc/bulk output can also generated with the cmd line in section 3/4 using the testing data as input. Click the file names to download (copy the backupLink ofr cmdline download). 
- testing data: [`Dropbox`](https://www.dropbox.com/s/dcgtsgww7jbrpyl/testdata.bed.gz) ([backupLink](https://data.cyverse.org/dav-anon/iplant/home/tarela/SELMA/testdata.bed.gz))
- testing peak file(optional for -p): [`Dropbox`](https://www.dropbox.com/s/a4r3gzux7v72rr9/testpeak.bed) ([backupLink](https://data.cyverse.org/dav-anon/iplant/home/tarela/SELMA/testpeak.bed))
- output for SELMA **bulk** mode with testing data input: [`Dropbox`](https://www.dropbox.com/sh/x8f29ao73t5ka8a/AADPjRgtgmW0DXJTiPMYWIS-a?dl=0)
- output for SELMA **sc** mode with testing data input: [`Dropbox`](https://www.dropbox.com/sh/a292670gqfw2uaa/AABijJfJCwIqNm1kW3tak4-da?dl=0) 


## 10. Other options of SELMA pipeline
You can also specify the following options for more accurate bias estimation and correction:
-   -\-extend=EXTEND    
[optional]Extension size from the peak summits, default is +/- 200bp from each peak summit. The macs3 detected peaks will be extended to 400bp centered on the summit of the peak for followiing analysis (e.g., biasExpected clevages in bulk mode; cell clustering in sc mode). For the process with external peak file inputted (-p), peaks will be extended from the peak center. 
-  -\-peakQval=PEAKQVAL
[optional]Qvalue cutoff in macs3 peak calling, default is 0.01 (-q 0.01 in macs3). This parameter is ignored if (-p) is activaged.
-  -\-bias=BIAS
[optional]Methods of intrinsic cleavage bias estimation, choose from naked (default, use SELMA pre-estimated bias model from naked DNA data) or chrM (use cleavages on mtDNA to estimate bias). Our pre-processed bias model from naked DNA data was proved to work perfectly in human and mouse. For other species, users can consider using the chrM option. 
-  -\-scATAC10x
[sc optional]Turn on this parameter to use 10x scATAC mode, in which the data format is assume to be PE and the reads 5'end will be shift back to represent the cleavage sites. If the fragments bed file was directly generated from the 10x cellranger-atac pipeline, users should set this parameter to ensure that the Tn5 cleavage site can be correctly captured. 
-  -\-cellnames=CELLNAMES
[sc optional]Single column file for name list of used individual cells, each line contain the name of the individual cell. This parameter is only used for sc mode. This parameter is usually ignored for common usage. 
-  -\-readCutoff=READCUTOFF
[sc optional]Reads number cutoff for high quality cells. Cells with < 10000(default) reads will be discarded in the analysis. For some low sequencing depth data, users can change this parameter to include more cells in the analysis. Setting lower number for this parameter will possibly lead to the inaccurate clustering results due to the low quality cells. 
-  -\-lowBiasPeak=LOWBIASPEAK
[sc optional]use top% peaks with less bias effect. Default is 66 (66%, use top 2/3 peaks with lowest bias for single-cell analysis). In our analysis, we showed that different % of top peaks will always improve the clustering analysis. 
-  -\-peakMinReads=PEAKMINREADS
[sc optional]peaks with < 10(default) cleavages covered in the whole fragment files will be discarded in the analysis.
-  -\-peakMaxReads=PEAKMAXREADS
[sc optional]peaks with > X cleavages covered in the whole fragment files will be discarded in the analysis. Set NA to close this function (default)
-  -\-clusterMethod=CLUSTERMETHOD
[sc optional]Method used for single cell clustering analysis. Default is PCAkm(PCA dim reduction + K+means clustering. Optional choices (Seurat,APEC,Cicero,ArchR) require related packages installed (described in section 6)
-  -\-clusterNum=CLUSTERNUM
[sc optional] number of clusters specified for clustering. Only used for PCAkm[clusterMethod] method. If users have no idea approximately how many clusters the sample should have, we suggest to set 7 as default or simply try other optional methods (described in section 6). 
-  -\-topDim=TOPDIM
[sc optional] number of dimensions (with highest Variance) used for clustering. Only used for PCAkm(PC) and Seurat (Latent variable). This number is suggested to be >=30 (deafult=30)
-  -\-UMAP
[sc optional]Turn on this parameter to generate a UMAP plot for the clustering results
-  -\-overwrite
[optional]Force overwrite, this cmd will rm existing result if set !! SELMA will terminated if there is a folder with the same name as -o in the current directory. Turn on this parameter to force SELMA running. 
-  -\-keeptmp
[optional]whether or not keep the intermediate results (tmpResults/)

# Intermediate data and scripts in the manuscript
The intermediate data and scripts in the manuscript can be download [here](https://www.dropbox.com/sh/rc0sd0x40e0dmqg/AAAgafYjM6HNYhlU185bGjjaa?dl=0). Check the README file in the folder for detailed annotation.

