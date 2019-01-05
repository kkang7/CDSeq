# Octave version
The code is modified for Octave users. 
The only part that was modified was the parallel computing using multicore. Other parts remain the same.
I did this pretty quickly, so it may have some bugs that I didn't realize. 
If you are experiencing some troubles or having some errors, please shoot me an email at kangkai0714@gmail.com

--------------------------------

# How to use:
In Octave command window 
1. Compile cpp file: (you need to compile the cpp file, the provided MATLAB compiled file may not work) 

mkoctfile --mex -DMATLAB_MEX_FILE CDSeqGibbsSampler.cpp

2. Install and load parallel package

pkg install -forge struct

pkg install -forge parallel

pkg load parallel

3. Follow the examples in demo.m 
