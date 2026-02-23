function [partsIF]=ImgFLN_sptMap2sptParts_core(Map,List)
parts={};
partsIF={};
%0. In some cases, the 1st colomn could be NaN but 2nd has a number.

%1. sptID, 2. cellID-Frm, 3. cellID-ID, 4. LRT
Marker=false(size(List,1),1);

for iRow=1:size(List,1)% for each row
    if isnan(Map(iRow,1)) && ~isnan(Map(iRow,2))
        Map(iRow,:)=[Map(iRow,2:3),nan(1,1)];
    end
    if isnan(Map(iRow,1)) || Marker(iRow)
        continue;
    end
    selRow=Map(iRow,:);
    nxtRow=selRow(~isnan(selRow));

%1. initialization
    curPart=iRow;
    partCollect={};
    StkDiv=[];%[iRow, spotID in the colomn of Map]
    StkParts={};
    isEmptied=false;
    while ~isempty(nxtRow) || ~isempty(StkDiv)
%WHILE Part 1
        if isEmptied% remove head if looped here from WHILE Part 3
            StkDiv(1)=[];
            StkParts(1)=[];
            isEmptied=false;
        end
%WHILE Part 2
        nDes=length(nxtRow);
        if nDes==1
            curPart=cat(1,curPart,nxtRow);
            if nxtRow~=-1
                Marker(nxtRow)=true;
            end
        elseif nDes==2
            StkDiv=[StkDiv;nxtRow(2)];%at which spot does the thread splitted?
            StkParts=[StkParts,{curPart}];%store curPart for 'split' cases
            curPart=cat(1,curPart,nxtRow(1));
            Marker(nxtRow(1))=true;
        else
            error(['Spot',num2str(iRow),' has more than 3 descendants!']);
        end
%WHILE Part 3
        if nxtRow~=-1
            selRow=Map(nxtRow(1),:);
            nxtRow=selRow(~isnan(selRow));
        else
            curPart=curPart(1:end-1);
            nxtRow=[];
        end
%WHILE Part 4
        if isempty(nxtRow)%means this thread is ended
            isEmptied=true;
            partCollect=[partCollect;curPart];
            if ~isempty(StkDiv)
                curPart=StkParts{1,1};
                nxtRow=StkDiv(1);
            end
        end
    end

%2. select the longest thread for this iRow, build PartInfo
    if length(partCollect)>1
        [~,iMax]=max(cellfun(@(x) length(x),partCollect));
        partCollect=partCollect(iMax);
    end
    selRow=partCollect{1,1};
    partsInfo=[selRow,List(selRow,1:3)];

%3. if the current part is too short, omit it.
    if partsInfo(end,2)-partsInfo(1,2)>1%<---note here.
        %if the purpose is to draw replication, then 1 is suitable.
        parts=[parts;partCollect];
        partsIF=[partsIF;partsInfo];
    end
end

%4. For all parts, find those that has different starting point but large mutual 
% body.
% Round 1
C=nchoosek(1:1:length(parts),2);
repParts=[];
for iC=1:size(C,1)
    mut=intersect(parts{C(iC,1)},parts{C(iC,2)});
    if ~isempty(mut)
        repParts=cat(1,repParts,[C(iC,1),C(iC,2)]);
    end
end

% Round 2
[repParts]=combDup(repParts);
del=[];
for iRP=1:length(repParts)%each pair of duplicated rows in 'parts'.
    Lcomp=zeros(length(repParts{iRP,1}),1);
    for iRP1=1:length(repParts{iRP,1})%each element in this row.
        curL=partsIF{repParts{iRP,1}(iRP1),1}(:,2);
        Lcomp(iRP1)=curL(end)-curL(1)+1;
    end
    [~,idMax]=max(Lcomp);
    sel=false(1,length(repParts{iRP,1}));
    sel(idMax)=true;
    del=[del,repParts{iRP,1}(1,~sel)];
end
parts(del)=[];
partsIF(del)=[];
end
function [uit]=combDup(in)
[~,w]=unique(in,'stable');
iDup=setdiff(1:numel(in),w);
Dup=in(iDup);
marker=false(size(in,1),1);
uit={};
for ix=1:length(iDup)
    [iRow,~]=find(in==in(iDup(ix)));
    cellout=reshape(in(iRow,:),1,[]);
    cellout(cellout==Dup(ix))=[];
    uit=[uit;{[cellout,Dup(ix)]}];
    marker(iRow)=true;
end
%end
uit=[uit;num2cell(in(~marker,:),2)];
end
