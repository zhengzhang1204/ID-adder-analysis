%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%function [Spt2PHMap]=getSpt2PHMap(SLines,PHLines,CIDinParts,HorLineMap)
%Stg3Neo_V2
%Author:  Zheng Zhang
%==========================================================================
%**********Output********
%Spt2PHMap:     Maps of each spot line and its corresponding PH line
%**********Input*********
%SLines:      coordinates of the yellow lines on lineage tree (Nx4 double)
%PHLines:       coordinates of the blue lines on lineage tree (Mx4 double)
%PHparts:       CIDs of all PH lines (Mx1 cell, Ax2 double in each cell)
%PHMaps:        Map of PH lines (full map, ie., Mx3 double)
%=========================================================================
% PURPOSE:
%Find out the corresponding PH line for each spot line. Some PH lines at the
%beginning may not have corresponding spot line. Some spot lines may not have
%corresponding PH line.
%=========================================================================
% OTHER INFO:
% 1. Format of Spt2PHMap:
%   For each spot line, Column 1 is the idx of the PH line that this spot line attaches to.
%   Column 2 and 3 is the PH line's ending CID.
%   Column 4-9 are for daughter's PH index and PH initial CID. 
%   [ref cell's rowIdx+PHendingCID, upper cell's rowIdx+PHinitCID, lower cell's rowIdx+PHinitCID];
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

function [Spt2PHMap]=getSpt2PHMap(SLines,PHLines,PHparts,PHMaps)
Spt2PHMap=nan(size(SLines,1),9);
for i=1:size(SLines,1)%for each spot line
    sptY=SLines(i,2);%Y pos of the current spot line on lineage tree
    idx=find(PHLines(:,2)==sptY);%index of the corresponding PH line
    if ~isempty(idx)
        Spt2PHMap(i,1:3)=[idx,PHparts{idx,1}(end,:)];
        if ~isnan(PHMaps(idx,2))%if the corr. PH line has upper daughter
            Spt2PHMap(i,4:6)=[PHMaps(idx,2),PHparts{PHMaps(idx,2),1}(1,:)];
        end
        if ~isnan(PHMaps(idx,3))%if the corr. PH line has lower daughter
            Spt2PHMap(i,7:9)=[PHMaps(idx,3),PHparts{PHMaps(idx,3),1}(1,:)];
        end
    end
end
end