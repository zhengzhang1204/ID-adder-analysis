function [hasSis,curCat,curLrt,sisIdx,errTxt]=ImgFLN_simpleIdentSisChr(SLines,iSL,PHLines,map)
% Method 1: Use existing map to identify sister.
% Find in the map if the current iSL has a sister
% Confirm sister if its sister ends up in the same PH line.
hasSis=false;
curCat=[];
curLrt=[];
sisIdx=[];
errTxt=[];

%core
[iRow,iCol]=find(map==iSL);
if length(iRow)~=1
    return;
end
if ~isempty(iRow)
    sisIdx=map(iRow,3-iCol);
    if isnan(sisIdx)
        sisIdx=[];
        return;
    end

    %idx of PH line containing the end of current SLine
    idx2=find(cellfun(@(x) ~isempty(find(x(:,1)==SLines{iSL,1}(end,2) & x(:,2)==SLines{iSL,1}(end,3),1)), PHLines));

    %idx of PH line containing the end of sister SLine
    Q=~isempty(find(SLines{sisIdx,1}(:,2)==SLines{iSL,1}(end,2) & SLines{sisIdx,1}(:,3)==SLines{iSL,1}(end,3), 1));
    %idx of PH line containing the end of sister SLine
    idx1=find(cellfun(@(x) ~isempty(find(x(:,1)==SLines{sisIdx,1}(end,2) & x(:,2)==SLines{sisIdx,1}(end,3),1)), PHLines));


    if idx1==idx2
        hasSis=true;
        Lrt=getLrtsOfOvlpRegn(PHLines{idx1,1},SLines(iSL,1));
        Lrt=Lrt{1,1};
        curLrt=mean(Lrt,'omitnan');
        if curLrt>0.5
            curCat=3;
        else
            curCat=1;
        end

    end
 
    if Q
        hasSis=true;
        Lrt=getLrtsOfOvlpRegn(PHLines{idx2,1},SLines(iSL,1));
        Lrt=Lrt{1,1};
        curLrt=mean(Lrt,'omitnan');
        if curLrt>0.5
            curCat=3;
        else
            curCat=1;
        end
    end
end
if ~isempty(sisIdx) && ~hasSis
    sisIdx=[];
end
end

function [uit]=getLrtsOfOvlpRegn(curCIDpart,SLines)
uit=cell(size(SLines));
hsOtherSL=false(size(SLines,1),1);
for i=1:length(SLines)
    Lrt=nan(size(SLines{i,1},1),1);
    for iSp=1:size(SLines{i,1},1)
        idx2=find(curCIDpart(:,1)==SLines{i,1}(iSp,2) & curCIDpart(:,2)==SLines{i,1}(iSp,3),1);
        if ~isempty(idx2)
            Lrt(iSp,1)=SLines{i,1}(iSp,4);
        end
    end
    uit{i,1}=Lrt;
end
end