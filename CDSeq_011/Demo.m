% Demo code for computational deconvolution method CDSeq
% coder: Kai Kang
% CDSeq version: 0.1.1
% last updated: 11/29/2018
%-------------------------------------------------------------------------
% Reference: A novel computational complete deconvolution method using
% RNA-seq data (submitted)
% Authors: Kai Kang, Qian Meng, Igor Shats, David Umbach, Melissa Li,
% Yuanyuan Li, Xiaoling Li, Leping Li. 
% Affiliation: National Institute of Environmental Healthe Sciences
% email: kai.kang@nih.gov or kangkai0714@gmail.com
%-------------------------------------------------------------------------
%% CDSeq function
%-------------------------------------------------------------------------
% inputs:
%-------------------------------------------------------------------------
% mydata     -- RNA-seq raw read count data, genes by samples
% beta       -- hyperparameter for cell type-specific GEPs
% alpha      -- hyperparameter for cell type proportions
% T          -- number of cell types
% N          -- number of MCMC iterations
% gene_length (optional) -- effective length of the genes 
% referenceGEP(optional) -- gene expression profiles of pure cell lines, this
% is used for identify the CDSeq-estimated cell types
% poolsize    (optional) -- number of workers used for parallel computing.
% this is only used when T is a vector instead of a scalar. 
%-------------------------------------------------------------------------
% outputs:
%-------------------------------------------------------------------------
% estProp -- estimated sample specific proportions of the cell types
% estGEP  -- estimated cell-type-specific gene expression profiles
% estT(optional) -- estimated number of cell types when input T is given as a vector
% logPosterior(optional) -- log posterior of CDSeq estimates 
% cell_type_assignment(optional) -- cell type assignment for the estimated cell
% types by comparing to the given referenceGEP
% estProp_all(optional) -- all the estimated proportions for different T values, 
% only available when T is a vector
% estGEP_all(optional) -- all the estimated cell-type-specific GEPs for
% different T values, only available when T is a vector
%-------------------------------------------------------------------------
%% load data
% SyntheticMixtureData.mat contains the following variables
% A. mixture_samples  - 40 mixture samples generated using six different pure
% cell lines and 100 randomly chosen genes. they are 1)endothelial blood vessel, 
% 2)breast epithelial carcinoma,  3)B lymphocyte, 4)CD14+ leukapheresis, 
% 5)lung fibroblast, 6)normal breast.
% B. gene_length      - the effective gene length of the 100 genes.
% C. refGEP_readCount - the read counts data of the six pure cell lines and two randomly generated pseudo-cell-line expressions.
% D. true_GEP_RPKM    - true RPKM gene expressions for the six pure cell lines
% E. true_SSP         - true mixing proportions of the 40 mixture samples
% F. true_GEP_read    - true read rate for six pure cell lines
% G. true_GEP_gene    - true gene rate for six pure cell lines
% H. theta            - true mixing proportions for the 40 samples
load SyntheticMixtureData.mat 

%% ran CDSeq
mydata = mixture_samples;
beta = .5;
alpha = 5;
T = 6; 
%T = 2:8; % CDSeq will choose the most likely number of cell types
N = 700;
%% case 1
[estprop,estGEP] = CDSeq(mydata, beta, alpha, T, N);% without given the gene_length the estGEP will be reported as read rate
%%
% GEP estimation (black circles)
figure
k=1;
for i=1:6
    for j=1:6
        subplot(6,6,k)
        scatter(estGEP(:,i),true_GEP_read(:,j),5,'k'); hold on;
        k=k+1;
    end
end
% proportion estimation (red circles)
figure
k=1;
for i=1:6
    for j=1:6
        subplot(6,6,k)
        scatter(estprop(i,:),theta(j,:),5,'r'); hold on;
        k=k+1;
    end
end
%% case 2
[estprop,estGEP,estT,logpost] = CDSeq(mydata, beta, alpha, T, N, gene_length);% estT will be equal to T when T is a scalar
%%
% GEP estimation (black circles)
figure
k=1;
for i=1:6
    for j=1:6
        subplot(6,6,k)
        scatter(estGEP(:,i),true_GEP_gene(:,j),5,'k'); hold on;
        k=k+1;
    end
end
% proportion estimation (red circles)
figure
k=1;
for i=1:6
    for j=1:6
        subplot(6,6,k)
        scatter(estprop(i,:),theta(j,:),5,'r'); hold on;
        k=k+1;
    end
end
%% case 3
[estprop,estGEP,estT,logpost,cell_type_assignment] = CDSeq(mydata, beta, alpha, T, N,gene_length,refGEP_readCount);
%%
% GEP estimation (black circles)
figure
k=1;
for i=1:6
    for j=1:6
        subplot(6,6,k)
        scatter(estGEP(:,i),true_GEP_RPKM(:,j),5,'k'); hold on;
        k=k+1;
    end
end
% proportion estimation (red circles)
figure
k=1;
for i=1:6
    for j=1:6
        subplot(6,6,k)
        scatter(estprop(i,:),theta(j,:),5,'r'); hold on;
        k=k+1;
    end
end
%% case 4
TT = 2:20;
[estprop,estGEP,estT,logpost,cell_type_assignment,estprop_all,estGEP_all] = CDSeq(mydata, beta, alpha, TT, N,gene_length,refGEP_readCount);
%%
% GEP estimation (black circles)
figure
k=1;
for i=1:6
    for j=1:6
        subplot(6,6,k)
        scatter(estGEP(:,i),true_GEP_RPKM(:,j),5,'k'); hold on;
        k=k+1;
    end
end
% proportion estimation (red circles)
figure
k=1;
for i=1:6
    for j=1:6
        subplot(6,6,k)
        scatter(estprop(i,:),theta(j,:),5,'r'); hold on;
        k=k+1;
    end
end
%% case 5
TT = 2:15;
poolsize = 5;
[estprop,estGEP,estT,logpost,cell_type_assignment,estprop_all,estGEP_all] = CDSeq(mydata, beta, alpha, TT, N,gene_length,refGEP_readCount,poolsize);
%%
% GEP estimation (black circles)
figure
k=1;
for i=1:6
    for j=1:6
        subplot(6,6,k)
        scatter(estGEP(:,i),true_GEP_RPKM(:,j),5,'k'); hold on;
        k=k+1;
    end
end
% proportion estimation (red circles)
figure
k=1;
for i=1:6
    for j=1:6
        subplot(6,6,k)
        scatter(estprop(i,:),theta(j,:),5,'r'); hold on;
        k=k+1;
    end
end




