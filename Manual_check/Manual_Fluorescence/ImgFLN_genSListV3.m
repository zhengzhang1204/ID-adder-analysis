function [SList, TSCenCoor]=ImgFLN_genSListV3(Stats,InfCell,CIDinParts,maxSptN)
%create SList, get centroid and convert coordinates to each tile.
SList=zeros(length(Stats),8);%Format of SList: [Frame, Idx, PercentageInCell, oriN, cellLength, SpotY(from image roof),SpotX(local tile), CoorInCell_Y]
TSCenCoor=round(cell2mat({Stats.Centroid}'));
Frm=floor(TSCenCoor(:,1)/32)+1;%Frame in SList
SList(:,1)=Frm';
LocalCenCoor=[TSCenCoor(:,1)-(Frm-1)*32,TSCenCoor(:,2)];%Spot_Y,Spot_X in SList
SList(:,6:7)=LocalCenCoor;

%see each spot is in which cell contours.
%Format of SList: [Frame, Idx, PercentageInCell, oriN, cellLength, SpotY(from image roof),SpotX(local tile), CoorInCell_Y]
CIDpool=cell2mat(CIDinParts);
uniqFrm=unique(Frm);
for i=1:length(uniqFrm)%for each frame. 
    curCID=sort(CIDpool(CIDpool(:,1)==uniqFrm(i),2),'ascend');% CID of selected cells on this frame
    cellIdx=cellfun(@(x) find(InfCell{uniqFrm(i)}{1}==x),num2cell(curCID),'UniformOutput',false);%index of selected cells on current mask tile
    
    curSelSptIdx=find(SList(:,1)==uniqFrm(i));
    curSptIdx=sub2ind([256,32],SList(curSelSptIdx,7),SList(curSelSptIdx,6));%values based on local tile
    for j=1:length(cellIdx)%for each selected cell on this frame
        xxx=find(ismember(curSptIdx,cellIdx{j}));
        if isempty(xxx); continue; end
        SList(curSelSptIdx(xxx),2)=j;
        [cellY,~]=ind2sub([256,32],cellIdx{j});
        YmaxCell=max(cellY);
        SList(curSelSptIdx(xxx),5)=InfCell{uniqFrm(i)}{2}(j,2);
        SList(curSelSptIdx(xxx),8)=YmaxCell-LocalCenCoor(curSelSptIdx(xxx),2);
        SList(curSelSptIdx(xxx),3)=SList(curSelSptIdx(xxx),8)./SList(curSelSptIdx(xxx),5);
    end
end
TSCenCoor(SList(:,2)==0,:)=[];
SList(SList(:,2)==0,:)=[]; %remove spots that are not on the PH lineage tree

% In short cells, if two spots are too close, they should be one.
SList=squeezeSpotsInShortCells(SList);

%SList: [Frame, Idx, PercentageInCell, oriN, cellLength, SpotY(from image roof),SpotX(local tile), CoorInCell_Y), Short cell marker]
[SList,TSCenCoor]=removeStuffedSptV2(SList,TSCenCoor,maxSptN);%wrote on 23.02.20



end

function [SList,TSCenCoor]=removeStuffedSptV2(SList,TSCenCoor,maxN)
[pool,~,ic]=unique(SList(:,1:2),"rows");
del=[];
for i=1:size(pool,1)
    curSel=find(ic==i);
    if length(curSel)>maxN
        del=[del;curSel];
    end
end
if ~isempty(del)
    SList(del,:)=[];
    TSCenCoor(del,:)=[];
end
end

function [List]=squeezeSpotsInShortCells(List)
%This program judges if the cell is too small (<30% quantile). If so, check
%if there are spot pairs that is very close (arbitarily <10). This is
%because when cell is too short, the LRT is sensitive to spot location.
%Sqeeze those pairs by forcibly changing the LRTs of these two spots by 25% towards each other.
%Spots that are processed are marked in the 8th column as 2(short cell and changed) or 0(intact).
%Spots that are determined as 'Short Cells' are marked as 1(short cell).

% 1. initiation:
reg=List(1,1:2);%stores unique [Frm, cellID]
Length=List(1,5);%stores cell length of unique cellID
regN={[1]};%stores row number in list that are in this unique cellID
Marker=zeros(size(List,1),1);

% 2. each spot is in which cell:
for i=2:size(List,1)%for each spot
    idx1=find(reg(:,1)==List(i,1));
    idx2=find(reg(:,2)==List(i,2));
    if isempty(idx1)%meaning this spot is in a new frame
        reg=cat(1,reg,List(i,1:2));%update unique cellID list
        Length=cat(1,Length,List(i,5));%update cell length for this CellID
        regN=[regN;{[i]}];
        %register the row number for this cell ID
    elseif isempty(intersect(idx1,idx2))%meaning this spot is in a old frame but a new cell
        reg=cat(1,reg,List(i,1:2));%update unique cellID list
        Length=cat(1,Length,List(i,5));%update cell length for this CellID
        regN=[regN;{[i]}];
    else%meaning this spot is in a cell that is registered in 'reg'.
        selRowInReg=intersect(idx1,idx2);
        regN{selRowInReg,1}=cat(2,regN{selRowInReg,1},i);
    end
end

% 3. retrieve 0.3, 0.6 quantile of cell length to select short cells.
trio=quantile(Length,[0.3,0.6]);

% 4. for each cell (FrmCellID), see if length is too short. If short,
% modify the LRT of spots that are very close to each other.

for i=1:length(Length)% for each FrmCellID
    if Length(i,1)<trio(1) && size(regN{i},2)>1
        Ycomp=List(regN{i},[3,7]);
        %see if the distances between any pair of spots that are <=10 in regN. The
        %index of such pair(s) in regN are in 'selPair'.
        selPair=checkDistBtwSpotPairs(Ycomp(:,2),10);%<------------arbitary number detected
        if ~isempty(selPair)
            for iSP=1:size(selPair,1)
                selRowInSlist=regN{i,1}(selPair(iSP,:));
                newLRT=moveOneForth([List(selRowInSlist(1),3),...
                List(selRowInSlist(2),3)]);

                List(selRowInSlist(1),3)=newLRT(1);
                List(selRowInSlist(2),3)=newLRT(2);
                Marker(selRowInSlist(1))=2;
                Marker(selRowInSlist(2))=2;
               % disp(['Squeezed two spots in Cell [Frm ',num2str(reg(i,1)),...
               %    ' CellID ',num2str(reg(i,2)),'].']);
            end
        else
            for iSP=1:size(regN{i},2)
                Marker(regN{i}(1,iSP),1)=1;
            end
        end
    elseif Length(i,1)<trio(1) && size(regN{i},2)==1
        Marker(regN{i})=1;
    end
end
List=[List,Marker];
end

function [out]=checkDistBtwSpotPairs(in,refD)
B=nchoosek(1:1:length(in),2);
out=[];
for i=1:size(B,1)
    if abs(in(B(i,1))-in(B(i,2)))<=refD
        out=cat(1,out,B(i,:));
    end
end
end

function [out]=moveOneForth(in)
if in(1)<=in(2)
    isShift=false;
else
    isShift=true;
    in=fliplr(in);
end
out=[in(1)+(in(2)-in(1))/4,in(2)-(in(2)-in(1))/4];
if isShift; out=fliplr(out);end
end
