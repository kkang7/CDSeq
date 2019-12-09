%%
% CDseq: A novel computational complete deconvolution method using bulk RNA-seq 
% coder: Kai Kang
% version: 0.1.4
% last updated: 08/15/2019
%========================================================
% Major changes: modified Mat2Vec.m function to use int16 allocate 
% memoery for the vector to save memory cost
% also in CDSeqGibbsSampler.cpp, I did the similar changes in allocating 
% memory for the long vector
%%
%-------------------------------------------------------------------------
% Reference: A novel computational complete deconvolution method using
% RNA-seq data (submitted)
% Authors: Kai Kang, Qian Meng, Igor Shats, David Umbach, Melissa Li,
% Yuanyuan Li, Xiaoling Li, Leping Li. 
% Affiliation: National Institute of Environmental Healthe Sciences
% email: kai.kang@nih.gov or kangkai0714@gmail.com
%-------------------------------------------------------------------------

function [estProp, estGEP, estT, logPosterior, cell_type_assignment, estProp_all,estGEP_all] = CDSeq(mydata, beta, alpha, T, N, shrinker, gene_length, referenceGEP, poolsize)
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
fprintf('\n ---------------\n CDSeq version 0.1.4\n ---------------\n');
tic
if nargin<5 
    error('CDSeq requires at least 5 input arguments')
end

if nargin==5
    warning('Data-Dilution option is not used. It may take long time to run.')
end

if nargin>9
    error('CDSeq requires at most 9 input arguments')
end

if nargout<2
    error('CDSeq requires at least 2 output arguments')
end

if nargout>7
    error('CDSeq requires at most 7 output arguments')
end

if length(beta)>1
    error('beta is a scalar, not a vector')
end

if length(alpha)>1
    error('alpha is a scalar, not a vector')
end

if beta<=0 || alpha<=0
    error('beta and alpha have to be positive numbers')
end

if N<1
    error('N, the number of MCMC iterations, has to be greater than 1')
end

if nargin >6
    if size(mydata,1)~=length(gene_length)
        error('The number of genes in mydata has to be the same as that of in gene_length');
    end
end

if nargin==8
    if size(mydata,1)~=size(referenceGEP,1)
        error('Genes of reference gene profile should coincide with the input data')
    end
end


% check if mydata is raw count RNA-seq or normalized
colsum = sum(mydata(:,1));
if floor(colsum)~=colsum % integer check
    warning('mydata is NOT raw count data.')
end

if nargin==5
    mydata = round(mydata);
end
if nargin>=6
    mydata = ceil(mydata/shrinker);
end
N = ceil(N);

nT = length(T);

% if nargin>=7
%     if size(referenceGEP,2)<T(nT)
%         warning('The number of cell types in reference GEPs is greater than biggest input T values.')
%         %error('reference profile need to have at least %d cell types', T(nT));
%     end
% end

% The random seed
SEED = 3;
[GeneId,SampleId] = Mat2Vec(mydata);
%==============================================================
% if T is a scalar
%==============================================================
if nT==1
    if T<2, error('T, the number of cell types, has to be greater than 2');end
    if nargout>5, error('Too many output. When T is a scalar, CDSeq outputs at most five outputs.');end
    printout = 1;
    [csGEP,ssp] = CDSeqGibbsSampler( GeneId , SampleId , T , N , alpha , beta , SEED, printout);
    tssp = transpose(ssp);
    estProp = (tssp+alpha)./sum(tssp+alpha);
    estGEP_read = (csGEP+beta)./sum(csGEP+beta);
    if nargin<=6, estGEP = estGEP_read;end
    if nargin>6
%         if size(mydata,1)~=length(gene_length)
%             error('The number of genes in mydata has to be the same as that of in gene_length');
%         end
        estGEP = read2gene(estGEP_read,gene_length);        
    end
    
    if nargin>7
        if size(estGEP,1)~=size(referenceGEP,1)
            error('Genes of reference gene profile should coincide with the input data')
        end        
        % check if referenceGEP is raw count RNA-seq or normalized
        colsum = sum(referenceGEP(:,1));
        
        if floor(colsum)~=colsum % integer check
            warning('referenceGEP is NOT raw count data, the estimated proportion will be RNA proportions instead of cell proportions.')
            RawCount = 0;
        else
            RawCount = 1;
        end
        
        if RawCount==0, dt = corr(referenceGEP,estGEP);end
        if RawCount==1, dt = corr(referenceGEP,estGEP_read);end
        
        [cell_assign,~] = munkres(1-dt);
        [ro,~] = find(cell_assign==1); 
        if RawCount==1
            pure_cell_line_raw_count = referenceGEP(:,ro); % need to check
            estProp = RNA2Cell(sum(pure_cell_line_raw_count)',estProp);
            estGEP = gene2rpkm(estGEP,gene_length,referenceGEP(:,ro));% rpkm
        end
    end
    
    if nargout>=3,estT = T;end
    
    if nargout>=5 && nargin>=8,cell_type_assignment=ro;end
    
    if nargout>=4
        Beta = beta*ones(1,size(mydata,1));
        Alpha = alpha*ones(1,T); 
        logPosterior = logpost(mydata,T,estProp,estGEP,Beta,Alpha);
    end
end


%==============================================================
% if T is a vector
%==============================================================
if nT>1
    if sum(T<2)>0, error('T, the number of cell types, has to be greater than 2');end
    if nargout>=3,estT = T(1);end
    lgpst = zeros(1,nT);
    Beta = beta*ones(1,size(mydata,1));
    printout = 0;
    
    % parfor runs on each T independently and in random order
    % the outputs from all workers need to be stored in a big matrix
    [ngenes,nsamples] = size(mydata);
    sumT = sum(T);
    estProp_slice = zeros(sumT,nsamples);
    estGEP_slice = zeros(ngenes,sumT);
    last_ids = cumsum(T);
    first_ids = [1 last_ids(1:end-1)+1];
    
    % put the slice in the cells for parfor loop
    estPropSliced = mat2cell(estProp_slice,T);
    estGEPSliced = mat2cell(estGEP_slice,ngenes,T);
    
    % set the number of workers
    if nargin==9
        core_info = evalc('feature(''numcores'')');
        core_info_split = strsplit(core_info);
        numcores = str2double(core_info_split(14));
        if numcores>poolsize
            parpool(poolsize)
        else
            warning('You are requesting %d number of workers but running evalc(''feature(''''numcores'''')'') indicates that there are %d cores assigned to MATLAB, so only %d workers are available',poolsize,numcores,numcores)
            parpool(numcores-1)
        end
%         myCluster = parcluster('local');
%         if myCluster.NumWorkers<poolsize
%             warning('The input of poolsize is larger than the number of workers in local pool. Please modify Cluster Profiles to increase the number of workers.');
%             if myCluster.NumWorkers<length(T),parpool(myCluster.NumWorkers);else,parpool(length(T)); end            
%         else
%             parpool(poolsize);
%         end
    else
        core_info = evalc('feature(''numcores'')');
        core_info_split = strsplit(core_info);
        numcores = str2double(core_info_split(14));
        %fprintf('%d cores available for computing\n',numcores)
        if numcores>nT
            parpool(nT);
        else
            parpool(numcores-1);
        end
%         myCluster = parcluster('local');
%         if myCluster.NumWorkers<length(T)
%             warning('It is possible to modify Cluster Profiles to increase the number of workers to speed up the software.');
%             parpool(myCluster.NumWorkers);
%         else
%             fprintf('parpool(length(T)) here...\n')
%             parpool(length(T));
%         end
    end   
    parfor i=1:nT
        fprintf('CDSeq: running on T = %d\n',T(i));       
        [csGEP,ssp] = CDSeqGibbsSampler( GeneId , SampleId , T(i) , N , alpha , beta , SEED, printout);
        tssp = transpose(ssp);
        proptmp = (tssp+alpha)./sum(tssp+alpha);
        GEPtmp = (csGEP+beta)./sum(csGEP+beta);
        Alpha = alpha*ones(1,T(i));       
        % compute the log posterior
        lgpst(i) = logpost(mydata,T(i),proptmp,GEPtmp,Beta,Alpha);
        fprintf('lgpst(%d) = %f\n',i,lgpst(i));
        estPropSliced{i} = proptmp;
        estGEPSliced{i} = GEPtmp;
%         if i>1 && lgpst(i)>lgpst(i-1)
%             if nargout>=3,estT = T(i);end
%             estProp = proptmp;
%             estGEP = GEPtmp;
%         end
    end
    estProp_slice = cell2mat(estPropSliced);
    estGEP_slice = cell2mat(estGEPSliced);
    
    [~,idx] = max(lgpst);
    if nargout>=3,estT = T(idx);end
    estProp = estProp_slice(first_ids(idx):last_ids(idx),:);
    estGEP = estGEP_slice(:,first_ids(idx):last_ids(idx));
    
    % close workers pool
    poolobj = gcp('nocreate');
    delete(poolobj);
    
    if nargin>6
        if size(mydata,1)~=length(gene_length)
            error('The number of genes in mydata has to be the same as that of in gene_length');
        end
        estGEP = read2gene(estGEP,gene_length);
    end
    
    if nargin>=8
        if size(estGEP,1)~=size(referenceGEP,1)
            error('Genes of reference gene profile should coincide with the input data')
        end
        % check if referenceGEP is raw count RNA-seq or normalized
        colsum = sum(referenceGEP(:,1));
        
        if floor(colsum)~=colsum % integer check
            warning('referenceGEP is NOT raw count data, the estimated proportion will be RNA proportions instead of cell proportions.\n')
            RawCount = 0;
        else
            RawCount = 1;
        end
        
        dt = corr(referenceGEP,estGEP);
        [cell_assign,cost] = munkres(1-dt);
        [ro,co] = find(cell_assign);
        if nargout>=5,cell_type_assignment=ro;end
        
%         if nargout>=4
%             logPosterior = lgpst;
%         end
        
        if nargout>=6
            estProp_all = estPropSliced;
        end
        if nargout==7
            estGEP_all = estGEPSliced;
        end
        
        if RawCount==1
            pure_cell_line_raw_count = referenceGEP(:,ro); % need to check
            %fprintf('size(pure_cell_line_raw_count)=[%d,%d]\n',size(pure_cell_line_raw_count,1),size(pure_cell_line_raw_count,2));
            estProp = RNA2Cell(sum(pure_cell_line_raw_count)',estProp);
            estGEP = gene2rpkm(estGEP,gene_length,referenceGEP(:,ro));% rpkm
        end
        
    end
    
    if nargout>=4
        logPosterior = lgpst;
    end
    
    
end

timeval = toc;
if timeval<=3600
    fprintf('CDSeq completed using %.2f seconds\n',timeval);
else
    fprintf('CDSeq completed using %.2f hours\n',timeval/3600);
end
end
