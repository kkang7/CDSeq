This folder contains the instruction on processing fastq files of the deconvolution project
and the instruction on how to use CDSeq software
contact: Kai Kang 
email: kangkai0714@gmail.com

 
This part is not useful for users, it is a note for developer. 
================= FOLDERS =============================
Folders: 
NS50593_94_NS50596_97 foler contains the combined fastq files. Specifically, the first 20 fastq files are combined from NS50593 and NS50594, and files 21-40 are combined from NS50596 and NS50597 

The folders without name *_trimmed are the outputs without running cutadpt and only contain the first 20 samples without combining the two fastq files. So do NOT use them. 
The folders with name *_trimmed are the outputs after running cutadpt. 



The following part contains the instruction on processing fastq data
================= DATA PROCESSING =====================
First, run Combine_fast.sh to combine the two fastq files which come from the same sequencing library
Second, run Fastq_process.sh to process the data with cutadapt, STAR, and featureCounts. 




The following part contains the instructions for running CDSeq
================= Runing CDSeq ========================


----------------- files ------------------------------
1. CDSeq.m                     --- the CDSeq main function 
2. CDSeqGibbsSampler.cpp       --- the Gibbs sampler
3. CDSeqGibbsSampler.m         --- the matlab function declaration
4. CDSeqGibbsSampler.mexa64    --- the linux-compiled Gibbs sampler function 
5. CDSeqGibbsSampler.mexmaci64 --- the Mac-compiled Gibbs sampler function
6. Combine_fastq.sh            --- the script for combining the fastq files(see DATA PROCESSING section)
6. Demo.m                      --- the demo function that shows how to run the CDSeq software
7. Fastq_process.sh            --- the script for processing fastq data(see DATA PROCESSING section)
8. Mat2Vec.m                   --- convert matrix input to vector input
9. RNA2Cell.m                  --- convert RNA proportions to cell proportions
10.SyntheticMixtureData.mat    --- sample data to show the usage of CDSeq (see Demo.m)
11.cokus.cpp                   --- random number generator
12.ctrlcDetech.h               --- ctrl+c detector heaher file, but not working 
13.gene2rpkm.m                 --- convert gene rate to RPKM normalization  
14.logpost.m                   --- log posterior functin
15.munkres.m                   --- munkres algorithm for cell type association
16.read2gene.m                 --- convert read rate to gene rate


----------------- how to run CDSeq ---------------------
open Demo.m and follow the options there
if you have trouble with the mex function, try to compile the mex function on your machine
by running: 
>>mex CDSeqGibbsSampler.cpp

or email me if you have any questions. 


---------------- Windows user -------------------------
I have compiled the mex function on Mac and Linux machines, 
so if you are using a Mac or Linux, you should be able to run the code without a problem.
If you are using a Windows machine, you will need to compile the mex function before running CDSeq.
Run the following to compile the c++ function in your MATLAB 
>>mex CDSeqGibbsSampler.cpp  

If error pops up, mostly because MATLAB cannot find the c/c++ compiler. Follow the instruction in the error message to install
the compiler, then run the mex function again.

P.S. I will try to find a windows machine to compile and upload the compiled file later.

PPS: I have compiled the c++ function on Windows machine. so there should be no problem with Windows users now.

------------------ version 0.1.2---------------------
updates:
1. I added data-dilution option in the function as described in our paper. Please see demo.m for testing run. 




