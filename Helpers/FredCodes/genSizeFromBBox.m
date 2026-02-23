function sizeData = genSizeFromBBox(BBox)
%GENSIZEFROMBBOX generates size from BBox
%ZZ? BBox is [x,y,dx,dy] while sizeData is [dy, dx]
% fchang@fas.harvard.edu

sizeData = BBox(numel(BBox)/2+1:end);
indices = 1:numel(sizeData);
indices(2) = 1;
indices(1) = 2;
sizeData = round(sizeData(indices));
end

