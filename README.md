# SELMA: a computational method for modeling intrinsic biases in chromatin accessibility sequencing data
Genome-wide profiling of chromatin accessibility with the assay for transposase-accessible chromatin using sequencing (ATAC-seq) or DNaseI hypersensitivity sequencing (DNase-seq) has been widely used to identify regulatory DNA elements and transcription factor binding sites. However, enzymatic DNA cleavage exhibits intrinsic sequence biases that confound chromatin accessibility profiling data analysis. Simplex Encoded Linear Model for Accessible Chromatin (SELMA) is a computational method for systematic estimation of intrinsic cleavage biases from ATAC-seq and DNase-seq data. This method can be applied to both bulk and single-cell chromatin accessibility profiles for improved data analysis.

![GitHub](https://img.shields.io/github/license/Tarela/SELMA)
![GitHub](https://img.shields.io/github/v/release/tarela/SELMA)
![GitHub](https://img.shields.io/github/commit-activity/m/tarela/SELMA)
![GitHub](https://img.shields.io/github/last-commit/tarela/SELMA)
![GitHub](https://img.shields.io/github/repo-size/tarela/SELMA)

[![Python 3.6](https://img.shields.io/badge/python-3.6-blue.svg)](https://www.python.org/downloads/release/python-360/)
[![numpy 1.9](https://img.shields.io/badge/numpy-1.9-blue.svg)](https://pypi.org/project/numpy/)
[![R 3.0](https://img.shields.io/badge/R-3.0-blue.svg)](https://www.r-project.org/)

## 0. Introduction of SELMA package
# Our SELMA paper was published on Nature Communications. [link](https://www.nature.com/articles/s41467-022-33194-z) <br>
SELMA performs estimation and correction of intrinsic cleavage bias of DNaseI (DNase-seq) and Tn5 (ATAC-seq) data at both bulk and single-cell levels. SELMA uses DNase/ATAC-seq data from either naked DNA or mitochondrial DNA (mtDNA) to estimate the intrinsic cleavage bias. SELMA provides a series of bias free analysis for the bulk/sc DNase/ATAC-seq data. For bulk data, SELMA estimates the bias expected cleavages on chromatin accessibility regions (peaks) and compares with observed cleavages. For single-cell data, SELMA estimates the summarized bias score on each candidate chromatin accessibility region (peak bias score, PBS) and uses the peaks with low PBS for single-cell clustering analysis.

- Changelog<br>
v1.1.0 add single cell bias correction model. <br>
v1.0.2 optimize signal scanning step in bulk mode. <br>
v1.0.1 update (-p) option for customized peak files.<br>
v1.0.0 First version of SELMA with both single-cell(sc) and bulk mode.

## 1. Installation
- Package requirements<br>
SELMA requires [python](https://www.python.org) 3.6+ and [Rscript](https://www.r-project.org) v3+ to run.<br>
SELMA requires python3 packages [numpy](https://numpy.org) pre-installed.

\# for root user
```sh
$ cd SELMA
$ sudo python setup.py install  
```
\# if you are not root user, you can install SELMA at a specific location where you have write permission
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
- SELMA requires python3 packages [numpy](https://numpy.org) and [MACS3(v3.0.0a6)](https://pypi.org/project/MACS3/) pre-installed. See Section 4 if the user does not have MACS3 or uses customized peak file. 
- Bedtools (Quinlan et al., Bioinformatics. 2010) and UCSC tools (Kuhn et al., Brief Bioinform. 2013) are recommended for data pre-processing. The SELMA package will install both tools automatically if the users does not have them pre-installed in the default PATH. 
- pdflatex is recommended for generating the summary pdf document. To install pdflatex on macOS, you can download “MacTex” from http://tug.org/cgi-bin/mactex-download/MacTeX.pkg. To install pdflatex on linux, you need to install the texlive package from https://www.tug.org/texlive/. SELMA will generate a .txt version of summary report without pdflatex installed. A .tex file will also be generated for users to make the pdf document later. 
- Some functions (single-cell clustering) of SELMA requires the related packages pre-installed (see Seciton 4)
- The installation should be finished in about one minute.

## 2. Download the genome sequence database
SELMA requires the genome sequence (in .2bit format) prepared for running. You can download them from the UCSC genome browser or other public domains. 
- [hg38.2bit](https://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/hg38.2bit)
- [mm10.2bit](https://hgdownload.cse.ucsc.edu/goldenpath/mm10/bigZips/mm10.2bit)


## 3. Run SELMA (usage)
#### Essential paramters
To run SELMA by default parameters, you can set the following parameters:
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
-   --SCcorrection
[sc optional] Apply SELMA bias correction model to the scATAC-seq data. SELMA will use only low bias peaks for single-cell analysis if this parameter is not activated. 
-   -p PEAK, --peak=PEAK  
[optional] external peak file for the candidate chromatin accessibility regions. 

Example of running SELMA with default parameters:

\# sc mode 
```sh
$ SELMA -m sc -i ${path}/testdata_reads.bed.gz -g hg38 -f PE -o testsc -t ATAC -s ${path}/hg38.2bit --SCcorrection 
```

\# bulk mode 
```sh
$ SELMA -m bulk -i ${path}/testdata_reads.bed.gz -g hg38 -f PE -o testbulk -t ATAC -s ${path}/hg38.2bit
```

## 4. Customize target chromatin accessibility (peak) regions (v1.0.1).
SELMA provides an option (-p) to take user supplied customized peak files as the target chromatin accessibility regions for the SELMA analysis. The required peak file should be in BED format (plain text), have >=4 columns (chrom, start, end, name), and contain >=1000 peaks. This parameter was specifically designed for those who don't want to use MACS3. However, we recommend to use the same dataset (e.g., fragments.bed file) for the peak calling (with any method) to ensure sufficient cleavages/signal on the peak regions. Below is the example with an external/customized peak file: <br />
The test files (testdata_reads.bed.gz and testpeak.bed) in the following cmd lines can be downloaded via the link in section 8

\# sc mode 
```sh
SELMA -m sc -i ${path}/testdata_reads.bed.gz -p ${path}/testpeak.bed -g hg38 -f PE -o testsc -t ATAC -s ${path}/hg38.2bit
```

\# bulk mode 
```sh
SELMA -m bulk -i ${path}/testdata_reads.bed.gz -p ${path}/testpeak.bed -g hg38 -f PE -o testbulk -t ATAC -s ${path}/hg38.2bit
```

## 5. Pre-processing steps for generating the input fragments file.
SELMA takes aligned fragment files (in .bed format) as input. Users can perform any pre-processing steps to customize the fragments files. For example, keep only high-quality reads with perfect alignment (e.g., MAPQ > 30) to run SELMA. For bulk data, using unique paired-end fragments (unique loci) only can reduce the potential influence from PCR over-amplification. For single-cell data, users can keep only unique fragments in each individual cell.

## 6. Install and use published single cell clustering methods based on SELMA bias correction. 
SELMA sc mode implements several cell clustering methods in the single-cell clustering analysis in addition to the default Kmeans analysis. To activate these methods (name, version and link listed below), users need to install the related package, and specify the method by the --clusterMethod parameter. If a methods is declared by the --clusterMethod parameter but is not installed, SELMA will skip the single-cell clustering analysis.
- Seurat (required packages: [ArchR v1.0.1](https://satijalab.org/seurat/), [tabix](http://www.htslib.org/doc/tabix.html), and [bgzip](http://www.htslib.org/doc/bgzip.html))
- scran (required packages: [ArchR v1.0.1](https://satijalab.org/seurat/), [tabix](http://www.htslib.org/doc/tabix.html), and [bgzip](http://www.htslib.org/doc/bgzip.html))
- APEC (required package: [APEC v1.2.2](https://github.com/QuKunLab/APEC))

Note that [ArchR v1.0.1](https://www.archrproject.com) is used for Seurat and scran for scATAC-seq data analysis. 

SELMA also provides UMAP/t-SNE visualization for the single-cell clustering analysis. Users can activate this function by the --UMAP parameter. For the PCAkm method, the [umap](https://cran.r-project.org/web/packages/umap/index.html) package in R is required. 

## 7. Output files
1. `NAME_summaryReports.pdf` is the summary pdf file which contains information of:
     - Input file and parameter description
     - basic QC of the data
     - Summary of the SELMA bias estimation/correction results

    \#Note: This pdf file is only generated if pdflatex is pre-installed. A NAME_summaryReports.txt file is generated as well. A .tex file will also be generated in case users want to make the pdf document later.

2. `NAME_peaks.bed` is the peaks detected from the fragment files (using MACS3). Each peak was extended to a 400bp region centered at the peak summit. 

3. `NAME_cleavage.bw` (bulk mode only) is the profile of the 1bp cleavages of DNaseI/Tn5 in the peak regions. Plus and minus strand cleavages are separated to two files (cleavage_plus.bw, cleavage_minus.bw)

4. `NAME_biasExpCuts.bw` is the profile of the bias expected cleavages in the peak regions. Plus and minus strand cleavages are separated to two files (biasExpCuts_plus.bw, biasExpCuts_minus.bw)

5. `NAME_peakXcell.txt.gz` (sc mode only) is the peak by cell read count matrix generated from the single-cell analysis. Cells are filtered by the total reads count per cell (default >=10,000 reads). 

6. `NAME_scClustering.txt.gz` (sc mode only) is the cell clustering result using SELMA bias correction or debiased peakset.

7. `NAME_bias.txt` is the SELMA-estimated bias score matrix. This file will only be generated if users don't use the default parameter (--bias naked) and set --bias chrM to use mtDNA reads for bias estimatation.

## 8. Testing data and example of output files
We provided the test data for users to test SELMA. The sc/bulk output can also be generated with the cmd lines in Section 3/4 using the testing data as input. Click the file names to download (copy the backupLink for cmdline download). 
- testing data: [`Dropbox`](https://www.dropbox.com/s/9cqcrjh17cae2d4/testdata_reads.bed.gz)
- testing peak file(optional for -p): [`Dropbox`](https://www.dropbox.com/s/a4r3gzux7v72rr9/testpeak.bed)
- testing cellnames (optional for --cellnames in sc mode): [`Dropbox`](https://www.dropbox.com/s/9eb60r9xx7gbh13/testsc_cellnames.txt)
- output for SELMA **bulk** mode with testing data input: [`Dropbox`](https://www.dropbox.com/sh/x8f29ao73t5ka8a/AADPjRgtgmW0DXJTiPMYWIS-a?dl=0)
- output for SELMA **sc** mode with testing data input: [`Dropbox`](https://www.dropbox.com/sh/jl7a28w9a984tpz/AACoGYzLRnBwZLgmbGlu6bwSa?dl=0) 
- The SELMA with testing data (e.g., using sc mode) will be finished within 30 minutes.


## 9. Other parameters in the SELMA pipeline
You can also set the following parameters for more accurate bias estimation and correction:
- -\-extend=EXTEND  
[optional] Extension size from the peak summits, default is +/- 200bp from each peak summit. The MACS3 peaks will be extended to 400bp centered on the summit of the peak for the analysis (e.g., biasExpected clevages in bulk mode; cell clustering in sc mode). For the process with external peak file inputted (-p), peaks will be extended from the peak center coordinate. 
- -\-peakQval=PEAKQVAL  
[optional] Qvalue cutoff in MACS3 peak calling, default is 0.01 (-q 0.01 in MACS3). This parameter is ignored if (-p) is set.
- -\-bias=BIAS  
[optional] Bias estimation method, to be selected from naked (default, use SELMA pre-estimated bias score from naked DNA data) or chrM (use mtDNA reads to estimate bias). Naked DNA-generated bias model works fine for human and mouse. Users can consider using the chrM option for other species. 
- -\-kmer=KMER  
[optional] Length of K (K-mer length), choose from 6,8,and 10(default).
- -\-SCcorrection  
[sc optional] Apply SELMA bias correction model to the scATAC-seq data. 
- -\-scATAC10x  
[sc optional] Turn on this parameter to use 10X scATAC mode, in which the data format is assumed to be PE and the 5'end coordinate of each read will be shifted back to represent the actual cleavage sites. If the fragments bed file is directly generated from the 10X Cellranger-atac pipeline, users should set this parameter to ensure that the Tn5 cleavage site are correctly captured. 
- -\-cellnames=CELLNAMES  
[sc optional] Single column file for name list of used individual cells, each line contains the name of an individual cell. This parameter is only used for sc mode. This parameter is not used very common. 
- -\-readCutoff=READCUTOFF  
[sc optional] Reads number cutoff for high-quality cells. Cells with < 10000(default) reads will be discarded in the analysis. For samples with low sequencing depth, users can change this parameter to include more cells in the analysis. Setting a lower number for this parameter will possibly decrease the accuracy of clustering results due to the low-quality cells. 
- -\-lowBiasPeak=LOWBIASPEAK  
[sc optional] Filter peaks based on PBS and keep top% peaks with the lowest PBS for single-cell analysis. Default is 80 (80%, use top 80% peaks with lowest PBS). Note that different percentage of lowest PBS peaks will always improve the clustering analysis. 
- -\-peakMinReads=PEAKMINREADS  
[sc optional] Peaks with < 10(default) cleavages covered (across all high-quality cells) will be discarded in the analysis.
- -\-peakMaxReads=PEAKMAXREADS  
[sc optional] Peaks with > X cleavages covered (across all high-quality cells) will be discarded in the analysis. Set 0 to close this function (default)
- -\-clusterMethod=CLUSTERMETHOD  
[sc optional] Method used for single cell clustering analysis. Default is Kmeans(PCA dim reduction + K-means clustering). Optional choices (Seurat,scran, and APEC) require related packages installed (described in section 5)
- -\-clusterNum=CLUSTERNUM  
[sc optional] Number of clusters specified for K-means clustering. Only used for the PCAkm (setting by --clusterMethod) method. Default is 10. 
- -\-topDim=TOPDIM  
[sc optional] Number of dimensions (with highest Variance) used for clustering. Only used for PCAkm(PC) and ArchR (Latent variable). This number is suggested to be >=30 (deafult=60)
- -\-UMAP  
[sc optional] Turn on this parameter to generate a UMAP plot for the clustering results
- -\-overwrite  
[optional] Force overwrite, setting this parameter will remove existing result! SELMA will terminate if there is a folder with the same name as -o in the working directory. Set this parameter to force SELMA running. 
- -\-keeptmp  
[optional] Whether or not keep the intermediate results (tmpResults/)

# Reproduce cell clustering results using SELMA package
Users can reproduce one of the clustering results in the manuscript (Figure 6, Human hematopoietic cells, K-means clustering) by running SELMA with the following cmd line:
```sh
$ SELMA -m sc -i ${path}/testdata_reads.bed.gz -g hg38 -f PE -o testsc -t ATAC -s ${path}/hg38.2bit --SCcorrection --peakQval 0.1 --UMAP --overwrite --keeptmp --peakMaxReads 4000 --SCcorrection --cellnames ${path}/testsc_cellnames.txt
```
The test files (testdata_reads.bed.gz and testsc_cellnames.txt) in the cmd line can be downloaded via the link in section 8.


# Supplementary data and scripts
- SELMA bias for [DNaseI(DNase-seq)](https://www.dropbox.com/s/ncemdhp0cee3cic/DNase_SELMAbias_10mer.txt.gz) ([backupLink](https://data.cyverse.org/dav-anon/iplant/home/tarela/SELMA/DNase_SELMAbias_10mer.txt.gz)) and [Tn5(ATAC-seq)](https://www.dropbox.com/s/x5iiy27ef80fl19/ATAC_SELMAbias_10mer.txt.gz) ([backupLink](https://data.cyverse.org/dav-anon/iplant/home/tarela/SELMA/ATAC_SELMAbias_10mer.txt.gz)). SELMA-estimated 10-mer intrinsic cleavage bias scores for DNaseI and Tn5. Both bias score matrices were generated from naked DNA data DNase/ATAC-seq data. Note that these bias matrices are already built-in and used in the SELMA package.  
- [Footprint bias scores (FBSs) for ENCODE human consensus footprint regions](https://www.dropbox.com/s/f3m9q0fhlq4e9vc/consensusFP_biasScore.txt.gz) ([backupLink](https://data.cyverse.org/dav-anon/iplant/home/tarela/SELMA/consensusFP_biasScore.txt.gz)). The human genome-wide consensus footprint regions are from Vierstra et al., Nature. 2020 and were downloaded from [this link](http://vierstra.org/resources/dgf). SELMA footprint bias score (FBS) for each footprint region is in the last column of the file). 
- [Other intermediate data generated and scripts used in the manuscript](https://www.dropbox.com/sh/rc0sd0x40e0dmqg/AAAgafYjM6HNYhlU185bGjjaa?dl=0). Check the README file in the folder for detailed annotation.
- [Single cell peakXcell matrices before and after bias correction](https://www.dropbox.com/sh/v2lfv5qv6x3f4zy/AADwxqFHoH2vVCeU6ed5zDVSa?dl=0). Check the README file in the folder for detailed annotation.
- Detailed data analysis protocol for the bias estimation and evaluation (Figure1-2) of this study can be found in https://github.com/Tarela/SELMA_data_analysis

