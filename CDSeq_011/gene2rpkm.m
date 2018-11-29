function rpkm = gene2rpkm(phi_gene,gene_effective_length,cell_line_counts)
% this function convert the gene rate to rpkm normalization
% the phi_gene and cell_line_counts are GxC matrices 
% and their order of columns should match
% coder: Kai Kang

if nargin~=3
    error('Error: phi2rpkm takes 3 inputs.\n Usage:rpkm = phi2rpkm(phi_gene,gene_effective_length,cell_line_counts, cell_order)\n');
end
if nargout~=1
    error('Error: phi2rpkm gives 1 output.\n Usage:rpkm = phi2rpkm(phi_gene,gene_effective_length,cell_line_counts)\n');
end

[g,c] = size(phi_gene);
[g2,c2] = size(cell_line_counts);
g3 = length(gene_effective_length);


if g~=g2 || g~=g3 || g2~=g3
    error('Error:the three inputs should have the same number of rows');
end

if c~=c2
    error('Error: phi_gene and cell_line_counts should have the same number of columns');
end

% if c~=t
%     error('Error: phi_gene and cell_line_counts should have the same number of columns');
% end

tmp = sum(cell_line_counts./gene_effective_length./(sum(cell_line_counts)));
rpkm = phi_gene.*tmp*10^9;

end