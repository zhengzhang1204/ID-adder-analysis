function [hasSis,uit,medLrt,sisIdx,errTxt]=ImgFLN_identSisChrv2(SptLines,iSL,CIDinParts)
errTxt=[];
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
%1. get idx (idx1) of PHparts where the current spotline ends in.
idx=find(cellfun(@(x) ~isempty(find(x(:,1)==SptLines{iSL,1}(end,2) & x(:,2)==SptLines{iSL,1}(end,3),1)), CIDinParts));

%2. find the overlapping region between the current SptLine and PHpart and calculate median LRT
Lrt=getLrtsOfOvlpRegn(CIDinParts{idx,1},SptLines(iSL,1));
Lrt=Lrt{1,1};Lrt(isinf(Lrt))=[];
medLrt=median(Lrt,'omitnan');

%3. According to the avgLrt, judge if the spotline is at center, 1/4 or 3/4
[jud]=judPosPbl(medLrt,1.2);
%category of position in jud: 
% 1: 0.25
% 2: 0.5
% 3: 0.75
if length(find(jud))==1 %meaning the medLrt is definitely at certain position of a cell.
    uit=find(jud);
else                    %meaning the medLrt is at an uncertain position of a cell.
    avgLrt=mean(Lrt,'omitnan');
    [jud]=judPosPbl(avgLrt,1.1);
    medLrt=avgLrt;
    if length(find(jud))==1
        uit=find(jud);
    else
        if jud(1,1)
            uit=1;
        elseif jud(1,3)
            uit=3;
        end
    end
end

%prepare output
if uit==2
    hasSis=false;
    sisIdx=[];
    return;
elseif uit==1
    hasSis=true;
elseif uit==3
    hasSis=true;
end

% additional criteria:
isUitEqual2=isEnterOvlpRegnAln(CIDinParts{idx,1},iSL,SptLines,uit);
if isUitEqual2
    uit=2;
    hasSis=false;
    sisIdx=[];
    return;
end

%4. if it may have a sister, find it in the current PH cycle.
%4.1 list all spot lines that appeared in this PH cycle.
col=[];
for iPh=1:size(CIDinParts{idx,1},1)
    idx4=find(cellfun(@(x) ~isempty(find(x(:,2)==CIDinParts{idx,1}(iPh,1) & x(:,3)==CIDinParts{idx,1}(iPh,2), 1)),SptLines));
    col=cat(1,col,idx4);
end
col=unique(col);
col(col==iSL)=[];
[Lrts]=getLrtsOfOvlpRegn(CIDinParts{idx,1},SptLines(col,1));
[Lrts,col]=remNonOvlpSLines(SptLines{iSL},SptLines(col,1),col,Lrts);%Remove the spot lines that does not overlap with the current spot line.
if isempty(col)%return if no other spot lines appear in this PH cycle.
    sisIdx=[];
    return;
end

%4.2 Calculate their meanLRT
avgLrtS=cellfun(@(x) median(x,'omitnan'),Lrts);

%4.3 find the most possible sister spot line candidate.
if uit==1
    idx3=find(avgLrtS>0.5);
elseif uit==3
    idx3=find(avgLrtS<0.5);
end

%4.4 if there are more than one sister candidates found, look for the one with biggest overlapping rate.
if length(idx3)>1
    OvlpRates=getOvlpRate(SptLines{iSL,1},SptLines(col(idx3),1));
    [~,idx4]=max(OvlpRates);
    sisIdx=col(idx3(idx4));
else
    sisIdx=col(idx3);
end

end

%% nested functions
function [lrts,col]=remNonOvlpSLines(curSLine,candSLine,col,lrts)
curEndFrm=curSLine(end,2);
candStFrm=cellfun(@(x) x(1,2), candSLine);
col(candStFrm>=curEndFrm)=[];
lrts(candStFrm>=curEndFrm)=[];
end
function [uit]=isEnterOvlpRegnAln(curCIDpart,curISL,SLines,curCat)
uit=false;
if SLines{curISL,1}(1,2)>curCIDpart(1,1)
    targetCID=SLines{curISL,1}(1,2:3);
else
    targetCID=curCIDpart(1,1:2);
end



idx1=find(cellfun(@(x) ~isempty(find(x(:,2)==targetCID(1,1) & x(:,3)==targetCID(1,2),1)),SLines));
idx1(idx1==curISL)=[];
if ~isempty(idx1)%if there are other SLs enters this div cycle, gather their mean LRTs within the div cycle.
    Lrts=getLrtsOfOvlpRegn(curCIDpart,SLines(idx1,1));
    Lrts=cellfun(@(x) median(x,'omitnan'),Lrts);
    if curCat==1
        if isempty(find(Lrts>0.49, 1))%There has other SLs enter this div cycle, but LRT are all wrong. 
            uit=true;
        end
    elseif curCat==3
        if isempty(find(Lrts<0.51, 1))%There has other SLs enter this div cycle, but LRT are all wrong. 
            uit=true;
        end

    end
else
    uit=true;%if there has no other spot lines enters this div cycle.
end
end

function  [uit]=getOvlpRate(curLine,SptLines)
uit=nan(size(SptLines,1),1);
for i=1:size(SptLines,1)
    ref=curLine(:,2:3);
    req=SptLines{i,1}(:,2:3);
    cnt=0;
    for j=1:size(ref,1)
        idx=find(req(:,1)==ref(j,1) & req(:,2)==ref(j,2),1);
        if ~isempty(idx)
            cnt=cnt+1;
        end
    end
    uit(i)=cnt./size(ref,1);
end
end
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
function [uit]=judPosPbl(x,minThr)
%definition of position: 
% 1: 0.25
% 2: 0.5
% 3: 0.75
if x<=0.5
    y1=cos(2*pi*(x-0.25));
    y2=cos(2*pi*(x-0.5));
    y3=0;
    y1(y1<0.0001)=0.0001;
    y2(y2<0.0001)=0.0001;
    rt=max([y1,y2])./min([y1,y2]);
    if rt>minThr
        uit=([y1,y2,y3]==max([y1,y2]));
    else
        uit=[y1,y2,y3]>0;
    end

else
    y1=0;
    y2=cos(2*pi*(x-0.5));
    y2(y2<0.0001)=0.0001;
    y3=cos(2*pi*(x-0.75));
    y3(y3<0.0001)=0.0001;
    rt=max([y2,y3])./min([y2,y3]);
    if rt>minThr
        uit=([y1,y2,y3]==max([y2,y3]));
    else
        uit=[y1,y2,y3]>0;
    end
    
end
end
