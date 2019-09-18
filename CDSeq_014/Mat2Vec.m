function [GeneId, SampleId] = Mat2Vec(sampleMatrix)
%==========================================================================
% This function takes the bulk RNA-seq matrix (gene by samples or sample by gene) 
% as input and returns GeneId and SampleId as output
% coder : Kai Kang
% last update: 08/15/2019
%==========================================================================

% initialize empty vectors
[row, col] = size(sampleMatrix);
len = sum(sampleMatrix(:));

if isinf(len)
    error('sum of all entries in the data is infinity');
end
%fprintf('len=%f\n',len);
GeneId = zeros(1,len,'uint16');
SampleId = zeros(1,len,'uint16');
% if row< 2^8-1
%     GeneId = zeros(1,len,'uint8');
% elseif row < 2^16-1
%     GeneId = zeros(1,len,'uint16');
% elseif row < 2^32-1
%     GeneId = zeros(1,len,'uint32');
% end
% 
% if col < 2^8-1
%     SampleId = zeros(1,len,'uint8');
% elseif col < 2^16-1
%     SampleId = zeros(1,len,'uint16');
% elseif col < 2^32-1
%     SampleId = zeros(1,len,'uint32');
% end

% fill out the entries
first = 1;
for i=1:col
    for j=1:row
        last = first + sampleMatrix(j,i) - 1; 
        GeneId(first:last) = ones(1,sampleMatrix(j,i))*j;
        SampleId(first:last) = ones(1,sampleMatrix(j,i))*i;
        first = last+1;
    end
end



end