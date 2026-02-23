function [Dmsk,coor4Del,ParCoor,isIso]=IAB_drawGreenMask(msk,ParidxRaw,PHMap,lkg,idxList)
sz=size(lkg);
pool=find(cellfun(@(x) ~isempty(x),lkg));
isIso=true(size(PHMap,1),1);
lkgCIDClt=[];

% 1. see which parts are in green (isolated) and which parts are in blue (lineaged)
for i=2:size(PHMap,1)%for each PHpart
    if ~isnan(PHMap(i,2)) || ~isnan(PHMap(i,3))
        isIso(i)=false;
        lkgCIDClt=[lkgCIDClt;ParidxRaw{i}];
    end
end

Dmsk=false(msk);
if isempty(lkgCIDClt)
    idxSel=0;
else
    idxSel=sub2ind(sz,lkgCIDClt(:,2),lkgCIDClt(:,1));%indices of lineaged cells
end

% create 'del' mask, Dmsk;
r8=find(~ismember(pool,idxSel));
[sub4DelRow,sub4DelCol]=ind2sub(sz,pool(r8));
for iDel=1:length(r8)
    curIdx=idxList{sub4DelCol(iDel)}{sub4DelRow(iDel)};
    tempTile=Dmsk(:,(sub4DelCol(iDel)-1)*32+1:sub4DelCol(iDel)*32);
    tempTile(curIdx)=true;
    Dmsk(:,(sub4DelCol(iDel)-1)*32+1:sub4DelCol(iDel)*32)=tempTile;
end
% output2: sub-index of the cells that are isolated:
coor4Del=[sub4DelRow,sub4DelCol];

% create div cyc linkage for visualization
ParCoor=cell(size(ParidxRaw));
for iPH=1:length(ParidxRaw)
    curPar=ParidxRaw{iPH,1};
    curCoor=nan(size(ParidxRaw{iPH,1}));
    for iCID=1:size(curPar,1)
        [row,col]=ind2sub([256,32],idxList{curPar(iCID,1),1}{curPar(iCID,2),1});
        curCoor(iCID,:)=[floor(median(row)),floor(median(col))+32*(curPar(iCID,1)-1)];
    end
    ParCoor{iPH,1}=curCoor;
end
end