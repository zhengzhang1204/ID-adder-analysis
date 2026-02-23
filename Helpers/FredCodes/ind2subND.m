function idxCell = ind2subND(siz,ndx)
%IND2SUBND converts index to {subscript}, basically adapted code from
%matlab's original ind2sub but have it output in cell format so i don't
%need to know how many dimensions i give it.
% you can call data(idxCell{:})
%
% fchang@fas.harvard.edu
idxCell = cell(numel(siz),1);
nout = numel(siz);
k = [1 cumprod(siz(1:end-1))];
for i = nout:-1:1
    vi = rem(ndx-1, k(i)) + 1;
    vj = (ndx - vi)/k(i) + 1;
    idxCell{i} = double(vj);
    ndx = vi;
end


