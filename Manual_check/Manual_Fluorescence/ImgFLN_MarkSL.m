function [selSpotLines]=ImgFL_MarkSL(sptPartsInfo,CIDinParts)
selCIDCollect_end=nan(length(sptPartsInfo),3);%PHLines that the ending spot attaches to.
for iSptRow=1:length(sptPartsInfo)%for each spot line
    refCellID_end=sptPartsInfo{iSptRow,1}(end,2:4);%in which CID the last spot is. [frm, cid, lrt]
    y1=find(cellfun(@(x) ~isempty(find(x(:,1)==refCellID_end(1) & x(:,2)==refCellID_end(2), 1)),CIDinParts));
    if isempty(y1); continue; end
    y2=find(CIDinParts{y1,1}(:,1)==refCellID_end(1),1);
    selCIDCollect_end(iSptRow,:)=[y1,y2,size(CIDinParts{y1,1},1)];%[idx of PH lines in CIDinParts, idx of frame of this PH part]
end
sptSel=~isnan(selCIDCollect_end(:,1));
selSpotLines=sptPartsInfo(sptSel);
[selSpotLines,~]=sortSLines(selSpotLines,selCIDCollect_end);
end
