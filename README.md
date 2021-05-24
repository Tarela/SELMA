# SELMA: a computational framework for modeling intrinsic biases in chromatin accessibility sequencing data

Genome-wide profiling of chromatin accessibility with the assay for transposase-accessible chromatin using sequencing (ATAC-seq) or DNaseI hypersensitivity sequencing (DNase-seq) has been widely used for studying regulatory DNA elements and transcriptional regulation in many cellular systems. Efficient and thorough computational analysis is essential for extracting biological information from such high-throughput sequencing data. It has been reported that DNase cleavage of DNA has sequence preferences that can significantly affect the footprint patterns at transcription factor binding sites in genomic profiles. We found that enzymatic sequence biases commonly exist in both bulk and single-cell chromatin accessibility profiling data. Using a regular simplex encoding model, we developed a quantitative approach for accurate characterization and systematic correction of intrinsic sequence biases contained in ATAC-seq and DNase-seq data. This approach can be applied in bioinformatics for improved analysis of high-throughput chromatin accessibility sequencing.

## 1. Installation
SELMA requires [python](https://www.python.org) 3.6+ and [R](https://www.r-project.org) v3+ to run.

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
\# Install pdflatex 
Before you start running, SELMA will check your computer for pdflatex. If you have already installed pdflatex, SELMA will generate a summary report in addition to regular outputs and analysis results.
- To install pdflatex on macOS, you can download “MacTex” from http://www.tug.org/mactex/. After downloading the package MacTex.pkg from http://tug.org/cgi-bin/mactex-download/MacTeX.pkg, you just double click to install MacTex and get the pdflatex.
- For linux user, you can type the following cmd line to install pdflatex on your server/computer.
```sh
$ apt-get install texlive-all
```
- The installation of pdflatex on both mac and linux requires root privilege.


\# NOTE: 
- To install SELMA on MacOS, user needs to download and install Command Line Tools beforehand
- SELMA requires a python3 package [macs3](https://pypi.org/project/MACS3/) pre-installed
- SELMA suggests to have bedtools (Quinlan et al., Bioinformatics. 2010) and UCSC tools (Kuhn et al., Brief Bioinform. 2013) pre-installed for data pre-processing. The SELMA package also wrapped up both tools and would install automatically if the users did not have them pre-installed. 
- Some function of SELMA requires the related packages pre-installed (see seciton 5)

## 2. Download genome sequence database
SELMA require genome sequence (in .2bit format) prepared for running. You can download them from UCSC genome browser. 
- [hg38.2bit](https://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/hg38.2bit)
- [mm10.2bit](https://hgdownload.cse.ucsc.edu/goldenpath/mm10/bigZips/mm10.2bit)

## 3. Download pre-estimated DNaseI/Tn5 bias matrix
With the help of SELMA, we've estiamted DNaseI/Tn5 bias using DNase-seq/ATAC-seq generated from naked DNA data. You can download them from the following link. 
- [DNaseI(DNase-seq)](https://www.dropbox.com/s/ncemdhp0cee3cic/DNase_SELMAbias_10mer.txt.gz?dl=0)
- [Tn5(ATAC-seq)](https://www.dropbox.com/s/x5iiy27ef80fl19/ATAC_SELMAbias_10mer.txt.gz?dl=0)

## 4. Run SELMA (usage)
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
genome version of the input data, choose from hg38(default) and mm10
-   -s SEQUENCE, --sequence=SEQUENCE
genome sequence file in 2bit format
-   -o OUTNAME, --outname=OUTNAME
Name of output results


Example for run SELMA with all default parameters:


sc mode 
```sh
$ SELMA -m sc -i ${path}/fragments.bed -g hg38 -f PE -o outputname -t ATAC -s ${path}/hg38.2bit --kmer 10 --bias naked --cellnames usecells.txt --scATAC10x
```

bulk mode 
```sh
$ SELMA -m bulk -i ${path}/fragments.bed -g hg38 -f PE -o outputname -t ATAC -s ${path}/hg38.2bit --kmer 10 --bias naked
```


#### options
You can also specify the following options for more accurate bias estimation and correction:
-   -\-extend=EXTEND    
[optional]Extension size from the peak summits, default is +/- 200bp from each peak summit.
-  -\-peakQval=PEAKQVAL
[optional]Qvalue cutoff in macs3 peak calling, default is 0.01 (-q 0.01)
-  -\-bias=BIAS
[optional]Methods of intrinsic cleavage bias estimation, choose from naked (default, use SELMA pre-estimated bias model from naked DNA data) or chrM (use cleavages on mtDNA to estimate bias)
-  -\-scATAC10x
[sc optional]Turn on this parameter to use 10x scATAC mode, in which the data format is assume to be PE and the reads 5'end will be shift back to represent the cleavage sites
-  -\-cellnames=CELLNAMES
[sc optional]Single column file for name list of used individual cells, each line contain the name of the individual cell. This parameter is only used for sc mode
-  -\-readCutoff=READCUTOFF
[sc optional]Reads number cutoff for high quality cells. Cells with < 10000(default) reads will be discarded in the analysis
-  -\-lowBiasPeak=LOWBIASPEAK
[sc optional]use top% peaks with less bias effect. Default is 66 (66%, use top 2/3 peaks with lowest bias for single-cell analysis).
-  -\-peakMinReads=PEAKMINREADS
[sc optional]peaks with < 10(default) covered in the whole fragment files will be discarded in the analysis.
-  -\-clusterMethod=CLUSTERMETHOD
[sc optional]Method used for single cell clustering analysis. Default is PCAkm(PCA dim reduction + K+means clustering. Optional choices (Seurat,APEC,Cicero) require related packages installed
-  -\-UMAP
[sc optional]Turn on this parameter to generate a UMAP plot for the clustering results
-  -\-h5
[sc optional]Turn on this parameter to generate the peakXcell matrix in .h5 format
-  -\-overwrite
[optional]Force overwrite, this cmd will rm existing result if set !!

## 5. Install and use published single cell clustering methods based on SELMA bias correction. 
SELMA sc mode implements several well acknowledged cell clustering methods in the singlne cell clustering analysis. To activate these alternative methods (name, version and link listed below), users need to install the related package, and specify the method by the --clusterMethod parameter
- [Seurat v4.0.0](https://satijalab.org/seurat/)
- [APEC v1.2.2](https://github.com/QuKunLab/APEC)
- [Cicero v1.8.1](https://www.bioconductor.org/packages/release/bioc/html/cicero.html)


## 6. Output files
1. `NAME_summary.pdf` is the summary pdf file which contains information of:
     - Input file and parameter description
     - basic QC of the data
     - Summary of the SELMA bias estimation/correction results

    \#Note: This pdf file is generated only if pdflatex is pre-installed. 

2. `NAME_bias.txt` is the bias matrix estimated with SELMA methods. This file will only be generated if the users don't use the --naked parameter (i.e. use mtDNA reads to estimate bias instead)

3. `NAME_peaks.bed` is the peaks detected from the fragment files (with macs3). Each peak was extended to 400bp centered on the peak summit. 

4. `NAME_cleavage.bw` (bulk mode only) is the genome-wide profile of the 1bp cleavage of DNaseI/Tn5. Plus and minus strand cleavages will be separated to two files (cleavage_plus.bw, cleavage_minus.bw)

5. `NAME_biasExpected.bw` is the profile of the bias expected cleavage on the peak regions. Plus and minus strand cleavages will be separated to two files (biasExpected_plus.bw, biasExpected_minus.bw)

6. `NAME_peakXcell.txt` (sc mode only) is the peak X cell count matrix generated from the single cell analysis. The cells were filtered by the total reads count (default >=10k reads) and the peaks were filtered based on the intrinsic cleavage bias (default top2/3 peaks with lowest bias effect). 

7. `NAME_scClustering.txt` (sc mode only) is the cell clustering results using SELMA debiased free peakset. )


## 7. Testing data and example of output files
We provided the testing data for users to test the flexibility and the power of the SELMA and the example of `scATAC_summary.pdf` and `bulkATAC_summary.pdf` which generated from 10x scATAC data and ENCODE bulkATAC data, respectively. Click the file names to download. 
- test_fragments (small input data): [`Dropbox`](https://www.google.com)
- scATAC_summary (output): [`Dropbox`](https://www.google.com)
- bulkATAC_summary (output): [`Dropbox`](https://www.google.com)

