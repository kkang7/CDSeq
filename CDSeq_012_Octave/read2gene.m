function phi = read2gene(read_rate,gene_effective_length)
% this function convert read rate to gene rate
% coder: Kai Kang
% last update: 2/2/2018

if nargin~=2
    error('Error: read2gene takes 2 inputs');
end
if nargout~=1
    error('Error: read2gene gives 1 output');
end

[g,c] = size(read_rate);
g1 = length(gene_effective_length);
if g~=g1
    error('Error:two inputs must have same number of rows');
end
tmp = read_rate./gene_effective_length;
phi = tmp./sum(tmp);

end