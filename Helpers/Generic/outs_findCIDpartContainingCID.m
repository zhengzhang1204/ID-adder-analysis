% for a given CID, find which PH part contains it.
function [index]=outs_findCIDpartContainingCID(CID,PHparts)
index=[];
index=cellfun(@(x) find(x(:,1)==CID(1) & x(:,2)==CID(2)),PHparts,'UniformOutput',false);
index=find(cellfun(@(x) ~isempty(x), index));
end