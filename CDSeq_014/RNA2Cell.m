function cellProp = RNA2Cell(eta,rnaProp)
% convert RNA proportion to Cell proportion
% eta is nc by 1 vector, nc is the number of cell types, 
% the elements in eta has the same order as the ones in true proportion
% rnaProp is nc by d matrix, nc is the number of cell types, d is sample
% size
% eta and rnaProp should have the same order!
% coder: Kai Kang
% last update: 2/2/2018

if nargin~=2
    error('Error: function takes 2 inputs.');
end
if nargout~=1
    error('Error: functions gives 1 output.');
end

[nc,c]= size(eta);
[nc2,d] = size(rnaProp);

if c~=1
    error('Error: size(eta,2) should be 1, size(eta,1) should be the number of cell types.');
end

if nc~=nc2
    error('Error: size(eta,1) should equal to size(rnaProp,1) which is the number of cell types.');
end

tmp = rnaProp./eta;
cellProp = tmp./sum(tmp);
end