# CDSeq: A complete deconvolution method for dissecting tissue heterogeneity
CDSeq is able to estimate both cell-type-specific gene expression profiles and sample-specific cell-type proportions simultaneously using bulk RNA-seq data.  
----------------------------------------------
# CDSeq MATLAB code version 0.1.1

Free for academic use and educational use. Please contact me at kangkai0714@gmail.com for commercial use. 

This is the MATLAB code for computational deconvolution method CDSeq. 

Reference: A novel computational complete deconvolution method using RNA-seq data (submitted)
Authors: Kai Kang, Qian Meng, Igor Shats, David Umbach, Melissa Li, Yuanyuan Li, Xiaoling Li, Leping Li.

Affiliation: National Institute of Environmental Health Sciences, Biostatistics and Computational Biology Branch and Signal Transduciton Lab.

-----------------------------------------------
# CDSeq MATLAB code version 0.1.2
Added data-dilution option to speed up the algorithm. We provided explaination in the manuscript which you can find on bioRxiv. Check out the demo.m in CDSeq_012 for details. 
-----------------------------------------------
# CDSeq MATLAB code version 0.1.4
updates:
1. I used unsigned short integer types for large vectors to save memory cost.
2. I added binary search for cell type assignment. In early version, it was linear search. 

note:
This version of CDSeq takes advantages of new MATLAB features (R2018a and later versions) 
which allows mex C++ function hanles unsigned short int type (2 bytes). 

I only compiled the code on Mac and Linux but not Windows yet. Will find a windows machine to compile later. You can also compile by yourself following the instructions below.

To compile CDSeq's C++ function, run the following two lines of commands in MATLAB: 

mex -setup C++ 

mex -R2018a CDSeqGibbsSampler.cpp

-----------------------------------------------
# Data availability
Experimental data GEO code: GSE123604

The synthetic data (6 pure cell line GEPs and 40 mixtures) are available upon request since it is too large to put them on github.
