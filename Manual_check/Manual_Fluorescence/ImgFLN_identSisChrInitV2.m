function [SisPairMap]=ImgFLN_identSisChrInitV2(SptLines,initIdx,CIDinParts,SisPairMap)
%% Criteria to judge sister spot line:
% Within the division cycle where current spot line ends in, check its mean LRT.
% If the mean LRT says this spot line is at mid-cell, then it does not have a sister.
% If the mean LRT says this spot line is a 1/4 or 3/4 spot, then it must have a sister.
%   NOTE: For strains not THAT asynchronous (unlike cSDR), I have the following additional criteria:
%   if the spot line is alone when it enters its final division cycle, and its avgLrt is not that deviated, assign Cat=2 for it.
% If it has sister, then try to find such a candidate in this division cycle. How:
% 1. list all spot lines that appears in this division cycle.
% 2. for the spot line part that appears, calculate the mean LRT
% 3. if there are more than one spot lines that fits the criteria, choose the one that have biggest overlapping region with the current spot line.

%% Core
if length(initIdx)==1
    SisPairMap(initIdx,2)=2;
    return;
end

%1. get Lrt of all initIdx:
initEndingLrts=zeros(length(initIdx),1);
initStLrts=zeros(length(initIdx),1);
initEndingPHidx=zeros(length(initIdx),1);
for iInitIdx=1:length(initIdx)
    idx=cellfun(@(x) ~isempty(find(x(:,1)==SptLines{initIdx(iInitIdx),1}(end,2) & x(:,2)==SptLines{initIdx(iInitIdx),1}(end,3),1)), CIDinParts);%ending PH idx of currrent init SL.
    initEndingPHidx(iInitIdx)=find(idx);
    LrtE=getLrtsOfOvlpRegn(CIDinParts{idx,1},SptLines(initIdx(iInitIdx),1));%ending LRT
    LrtE=LrtE{1,1};
    initEndingLrts(iInitIdx)=median(LrtE,'omitnan');%ending LRT
    LrtS=getLrtsOfOvlpRegn(CIDinParts{1},SptLines(initIdx(iInitIdx),1));%Starting LRT
    LrtS=LrtS{1,1};
    initStLrts(iInitIdx)=median(LrtS,'omitnan');
end
upInitSLIdx=initIdx(initStLrts>0.5);
lowInitSLIdx=initIdx(initStLrts<0.5);
% 
if length(upInitSLIdx)==1 && length(lowInitSLIdx)==1 && initEndingPHidx(initIdx==upInitSLIdx)==1 && initEndingPHidx(initIdx==lowInitSLIdx)==1
    SisPairMap(upInitSLIdx,1)=lowInitSLIdx;
    SisPairMap(upInitSLIdx,2)=3;
    SisPairMap(lowInitSLIdx,1)=upInitSLIdx;
    SisPairMap(lowInitSLIdx,2)=1;  
    return;
end

if length(upInitSLIdx)==1%arbitarily define that the length here can only be 1 or 2!
    if initEndingPHidx(initIdx==upInitSLIdx)==1
        SisPairMap(upInitSLIdx,2)=3;
    else 
        SisPairMap(upInitSLIdx,2)=2;
    end
    SisPairMap(upInitSLIdx,1)=nan;
else
    endingPHInd=zeros(2,1);
    for iInitIdxU=1:length(upInitSLIdx)
        endingPHInd(iInitIdxU)=find(cellfun(@(x) ~isempty(find(x(:,1)==SptLines{upInitSLIdx(iInitIdxU),1}(end,2) & x(:,2)==SptLines{upInitSLIdx(iInitIdxU),1}(end,3),1)), CIDinParts));
    end
    if endingPHInd(1)==endingPHInd(2)
        if initEndingLrts(upInitSLIdx(1))>initEndingLrts(upInitSLIdx(2))
            SisPairMap(upInitSLIdx(1),2)=3;
            SisPairMap(upInitSLIdx(2),2)=1;
        else
            SisPairMap(upInitSLIdx(1),2)=1;
            SisPairMap(upInitSLIdx(2),2)=3;
        end
        SisPairMap(upInitSLIdx(1),1)=upInitSLIdx(2);
        SisPairMap(upInitSLIdx(2),1)=upInitSLIdx(1);
    else
        SisPairMap(upInitSLIdx(1),1)=nan;
        SisPairMap(upInitSLIdx(1),2)=2;
        SisPairMap(upInitSLIdx(2),1)=nan;
        SisPairMap(upInitSLIdx(2),2)=2;
    end
end

if length(lowInitSLIdx)==1
    if initEndingPHidx(initIdx==lowInitSLIdx)==1
        SisPairMap(lowInitSLIdx,2)=1;
    else
        SisPairMap(lowInitSLIdx,2)=2;
    end
    SisPairMap(lowInitSLIdx,1)=nan;
else
    endingPHInd=zeros(2,1);
    for iInitIdxL=1:length(lowInitSLIdx)
        endingPHInd(iInitIdxL)=find(cellfun(@(x) ~isempty(find(x(:,1)==SptLines{lowInitSLIdx(iInitIdxL),1}(end,2) & x(:,2)==SptLines{lowInitSLIdx(iInitIdxL),1}(end,3),1)), CIDinParts));
    end
    if endingPHInd(1)==endingPHInd(2)
        if initEndingLrts(lowInitSLIdx(1))>initEndingLrts(lowInitSLIdx(2))
            SisPairMap(lowInitSLIdx(1),2)=3;
            SisPairMap(lowInitSLIdx(2),2)=1;
        else
            SisPairMap(lowInitSLIdx(1),2)=1;
            SisPairMap(lowInitSLIdx(2),2)=3;        
        end
        SisPairMap(lowInitSLIdx(1),1)=lowInitSLIdx(2);
        SisPairMap(lowInitSLIdx(2),1)=lowInitSLIdx(1);
    else
        SisPairMap(lowInitSLIdx(1),1)=nan;%<---------------------------------------------------------need to re-consider
        SisPairMap(lowInitSLIdx(1),2)=2;
        SisPairMap(lowInitSLIdx(2),1)=nan;%<---------------------------------------------------------need to re-consider
        SisPairMap(lowInitSLIdx(2),2)=2;
    end
end
end




%% nested functions
function [uit]=getLrtsOfOvlpRegn(curCIDpart,SLines)
uit=cell(size(SLines));
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
