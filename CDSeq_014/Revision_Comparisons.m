% comparison with other deconvolution methods
% coder: Kai Kang
% last update: 8/5/2019

% This is made for CDSeq paper revision 1. 
% Reviewers asked for more comparisons, so here it is. 


%=========================================
% UNDO
%=========================================
%%
addpath('/Users/kangk3/Dropbox/Public/NIH/Project 1/Deconvolution_Revision')
load undo_result_matlab.mat
%%
load Deconvolution_Experimental_data_shrink_10.mat
%%

figure
k=1;
for i=1:4
    for j=1:2
        subplot(4,2,k)
        k=k+1;
        scatter(theta_adjcont(:,i),undo_result(:,j))
    end
end

%%
%=========================================
% find marker genes for the four cell types
%=========================================
marker_genes = zeros(size(phi_rpkm,1),4);
for i=1:size(phi_rpkm,1)
    cell_type_id = find(phi_rpkm(i,:));
    if length(cell_type_id)==1
        marker_genes(i,cell_type_id) = 1;
    end
end
[r,c] = find(marker_genes);
marker_gene_exp = cell(1,4);
for i=1:4
    marker_gene_exp{i} = phi_rpkm(r(c==i),i);
    
end

% get the gene names
marker_gene_names = cell(1,4);
for i=1:4
    last_id = sum(marker_genes(:,i));
    marker_gene_names{i} = gene_name_list_shrink10(r(c==i));
end

% sort marker genes expressions
r_sorted = zeros(size(r));
for i=1:4
    [~,tmpid] = sort(marker_gene_exp{i},'descend');
    tmpr = r(c==i);
    r_sorted(c==i) = tmpr(tmpid);
end

% get the gene names
marker_gene_names_sorted = cell(1,4);
for i=1:4
    last_id = sum(marker_genes(:,i));
    marker_gene_names_sorted{i} = gene_name_list_shrink10(r_sorted(c==i));
end


%% plot the marker genes expression of the four cell lines
figure
xvalues = cell_order_names;
yvalues = num2cell(1:length(r));%[marker_gene_names_sorted{1} ;marker_gene_names_sorted{2}; marker_gene_names_sorted{3}; marker_gene_names_sorted{4}];
h=heatmap(xvalues,yvalues,phi_rpkm(r_sorted,:),'GridVisible','off');
h.Colormap = parula;
h.ColorLimits = [0 5];
h.FontSize = 15;
%%
save Deconvolution_Experimental_data_shrink_10.mat marker_gene_names marker_gene_names_sorted -append
%% Data dilution on expt data 
addpath('/Users/kangk3/Dropbox/Public/NIH/Project 1/Deconvolution_CDSeq_code/CDSeq_012')
% load the data from data_to_David folder
%%

cell_lines_1 = [tumor_MCF7_gene_name_count_1 normal_breast_hMECs_hTERT_gene_name_count_1 lymphocytes_namalwa_gene_name_count_1 CAFs_Hs_343T_gene_name_count_1];
cell_lines_2 = [tumor_MCF7_gene_name_count_2 normal_breast_hMECs_hTERT_gene_name_count_2 lymphocytes_namalwa_gene_name_count_2 CAFs_Hs_343T_gene_name_count_2];
cell_lines_mean = ceil((cell_lines_1 + cell_lines_2)/2);

nz_rows = any(sample_matrix_mixture_gene_names,2);
mixGEP_nz = sample_matrix_mixture_gene_names(nz_rows,:);% remove zero rows
cell_lines_mean_nz = cell_lines_mean(nz_rows,:);
gene_length_nz = gene_name_length(nz_rows);

true_read = cell_lines_mean_nz./sum(cell_lines_mean_nz);
true_gene = read2gene(true_read,gene_length_nz);
true_rpkm = gene2rpkm(true_gene,gene_length_nz,cell_lines_mean_nz);


%% 
mydata = WD;%mixGEP_nz;%
beta = 0.5;
alpha = 5;
T=4;
N=700;
gene_length = el_wd;%gene_length_nz;%
refGEP_readCount = cells_filered_wd;%cell_lines_mean_nz;%
shrinkers = (10:100)/10;%(10:5:500)/10;% working shrink 10 data so divide the sequence by 10

estprop_all = cell(1,length(shrinkers));
estGEP_all = cell(1,length(shrinkers));
cell_type_assignment_all = cell(1,length(shrinkers));
time_DatDil = zeros(1,length(shrinkers));
for i=1:length(shrinkers)
    fprintf('-------- shrinker = %d--------\n',shrinkers(i))
    tic
    [estprop_all{i},estGEP_all{i},estT,logpost,cell_type_assignment_all{i}] = CDSeq(mydata, beta, alpha, T, N, shrinkers(i),gene_length,refGEP_readCount);
    time_DatDil(i) = toc;
end

save expt_data_data_dilution_10to100_N700.mat estprop_all estGEP_all cell_type_assignment_all time_DatDil_expt 
%%
expt_prop_corr_all = zeros(4, length(shrinkers));
expt_prop_corr_min = zeros(1,length(shrinkers));

expt_gep_corr_all = zeros(4, length(shrinkers));
expt_gep_corr_min = zeros(1,length(shrinkers));

for i=1:length(shrinkers)
    prop_corr_tmp = corr(estprop_all{i}',theta_adjcont);
    expt_prop_corr_all(:,i) = diag(prop_corr_tmp(:,cell_type_assignment_all{i}));
    expt_prop_corr_min(i) = min( expt_prop_corr_all(:,i) );

    gep_corr_tmp = corr(estGEP_all{i},true_rpkm);
    expt_gep_corr_all(:,i) = diag(gep_corr_tmp(:,cell_type_assignment_all{i}));
    expt_gep_corr_min(i) = min(expt_gep_corr_all(:,i));
end
%%
save expt_data_data_dilution_10to500_N2000.mat expt_prop_corr_all expt_prop_corr_min expt_gep_corr_all expt_gep_corr_min -append
%% GEP estimation
load expt_data_data_dilution_10to500_N700.mat
%%
fs = 22;
fig=figure;
left_color = [0 0 0];
right_color = [1 0 0];
set(fig,'defaultAxesColorOrder',[left_color; right_color]);
yyaxis left
gep_time_smooth = smooth(shrinkers*10,time_DatDil_expt/3600,0.1,'lowess');
plot(shrinkers*10,time_DatDil_expt/3600,'ko',shrinkers*10,gep_time_smooth,'k-'); hold on% in hours
xlabel('Shrinkage factor','FontSize',fs)
ylabel('Running Time (hours)','FontSize',fs)
yyaxis right
ylabel('Min correlation bwteen estimated GEPs and true GEPs','FontSize',fs)
gep_min_smooth = smooth(shrinkers*10,expt_gep_corr_min,0.1,'lowess');
plot(shrinkers*10,expt_gep_corr_min,'ro',shrinkers*10,gep_min_smooth,'r-');hold on
%scatter(shrinker,prop_corr_min,'c');
legd=legend('Original Running Time','Smoothed Running Time','Original Correlations','Smoothed Correlation');
legd.FontSize=15;
xt = get(gca, 'XTick');
set(gca, 'FontSize', 16)
hold off
%%
h=gcf;
set(h,'PaperOrientation','landscape');
set(h,'PaperUnits','normalized');
set(h,'PaperPosition', [0 0 1 1]);
pdffile = sprintf('CDSeq_speedup_data_dilution_GEP_min_ExptData.pdf');
print(gcf, '-dpdf','-fillpage','-r300', pdffile); 

%% prop estimation 
fs = 22;
fig=figure;
left_color = [0 0 0];
right_color = [1 0 0];
set(fig,'defaultAxesColorOrder',[left_color; right_color]);
yyaxis left
prop_time_smooth = smooth(shrinkers*10,time_DatDil_expt/3600,0.1,'lowess');
plot(shrinkers*10,time_DatDil_expt/3600,'ko',shrinkers*10,prop_time_smooth,'k-'); hold on% in hours
xlabel('Shrinkage factor','FontSize',fs)
ylabel('Running Time (hours)','FontSize',fs)
yyaxis right
ylabel('Min correlation bwteen estimated SSPs and true SSPs','FontSize',fs)
prop_min_smooth = smooth(shrinkers*10,expt_prop_corr_min,0.1,'lowess');
plot(shrinkers*10,expt_prop_corr_min,'ro',shrinkers*10,prop_min_smooth,'r-');hold on
%scatter(shrinker,prop_corr_min,'c');
legd=legend('Original Running Time','Smoothed Running Time','Original Correlations','Smoothed Correlation');
legd.FontSize=15;
xt = get(gca, 'XTick');
set(gca, 'FontSize', 16)
hold off
%%
h=gcf;
set(h,'PaperOrientation','landscape');
set(h,'PaperUnits','normalized');
set(h,'PaperPosition', [0 0 1 1]);
pdffile = sprintf('CDSeq_speedup_data_dilution_prop_min_ExptData.pdf');
print(gcf, '-dpdf','-fillpage','-r300', pdffile); 



%%
[tmp_estprop,tmp_estGEP] = CDSeq(mydata, beta, alpha, T, N, shrinkers(3),gene_length);%,refGEP_readCount./sum(refGEP_readCount));

%%














