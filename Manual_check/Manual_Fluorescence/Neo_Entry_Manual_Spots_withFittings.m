clear;
clc;
%v2.4b @2023.04.20: fixed a problem when saving spot map.
%v2.4a @2022.12.05: fixed a dead-loop in 'F'
%v2.2  @2022.09.23
% This version works for both fast and slow growths.
%% user-specified region
EXP='G:\20260124 RdnaA_M30_atc1_dnaNyPet2';
MappingMethod=5;%[4] fast growth
                %[5] slow growth<----M12-compatibile
deHazeAmount=0.1;
atmosphericLight=0.92;
offsetBK=100; % [100] real 16-bit
              % [3100] rescaled from 11-bit
maxSptN=8;

%% core
F=dir([EXP,'\M_res\p*_resM.mat']);
for m=14:length(F)%for each pos.<----------------------------------------------change back to 1
    posStr=F(m).name(1:3);
    load([EXP,'\M_res\',F(m).name(1:3),'_resM.mat'],'PosMask');
    load([EXP,'\Y_res\',F(m).name(1:3),'_resY.mat'],'Cand','TSstack');
    load([EXP,'\M_res\',F(m).name(1:3),'_lkgM.mat'],'Infcellv2','lxg','PHLines','PHMap','PHparts');
    nSchN=length(lxg);
    if exist([EXP,'\Y_res\',F(m).name(1:3),'_lkgY.mat'],'file')
        Sx=load([EXP,'\Y_res\',F(m).name(1:3),'_lkgY.mat'],'SptMap','SptList','SptCent');
        rawMap=Sx.SptMap;
        rawSList=Sx.SptList;
        rawTSCenCoor=Sx.SptCent;
    else
        rawMap=cell(1,length(Cand));
        rawSList=[];
        createLkgY([EXP,'\Y_res\',F(m).name(1:3),'_lkgY.mat'],nSchN);%create an empty lkg_Y file.
    end
    
    for iSchN=1:nSchN%for each side channel<----------------------------------change back to 1
        if (isempty(Cand{iSchN,1}) || isempty(lxg{1,iSchN})) && isempty(rawMap{1, iSchN})
            disp([posStr,'SchN',num2str(iSchN),' is skipped.']);
            continue;
        end
        disp([posStr,'SchN',num2str(iSchN),' is being processed.']);
    %1. Generate SList from lxg, Cand and Masks. Former:IN:cand.stats, cand.intinfo,infcell,lxg.  OUT:SList, TScenCoor
        if isempty(rawMap{1, iSchN})% use existing map for drawing
            [SList,TSCenCoor]=ImgFLN_genSListV3(Cand{iSchN}.stats,Infcellv2(:,iSchN),PHparts{iSchN},maxSptN);%
        end
        PHMask=cell2mat(cellfun(@(x) x{iSchN,1},PosMask,'UniformOutput',false));
        FLTS=im2uint8(imreducehaze(mat2gray(TSstack(:,1:size(PHMask,2),iSchN)-offsetBK),deHazeAmount,"AtmosphericLight",atmosphericLight));
        %A=prepLines(Infcellv2(:,iSchN),PHparts{iSchN},FLTS);% set up
        %initial guess.
        FLTS=imadjust(FLTS,[0.3 1]);
        
    %2. Generate map from SList. Former: imgFL_MergeToLxgv3CORE. IN:SList, lxg, cand.intinfo.  OUT:Map, SList ('E:\DLstg32B_drawFLrawSList.mat','ini','TiledCent','Map')
        if ~isempty(rawMap)% use existing map for drawing
            if ~isempty(rawMap{1,iSchN})
                Map=rawMap{1,iSchN};%replace the calculated map by formerly corrected map.
                SList=rawSList{1,iSchN};
                TSCenCoor=rawTSCenCoor{1,iSchN};
                disp('Using existing Map, SList and TSCenCoor.');
            else
                [Map,SList]=imgFL_Spt2Lxgv3_simp(SList,lxg{1,iSchN}{1,2},Cand{iSchN,1}.IntInfo);
                Map=Map(:,1:3);
                Map(Map==0)=nan;
            end
        else% draw from scratch
            [Map,SList]=imgFL_Spt2Lxgv3_simp(SList,lxg{1,iSchN}{1,2},Cand{iSchN,1}.IntInfo);
        end
    %3. Draw DivTree according to lxg. (former: DL_helper_Lxg2LL.m, Since 23.02.20, this part is in IAB.)
%         [LineCol,LinePri,CID]=imgPH_simpDecipherLineage(lxg{1,iSchN}{1,2});
%         [PHparts,CIDinParts,HorLineMap,fig]=drawTreeFromLL(LineCol,LinePri,CID);
%         clear('LineCol','LinePri','CID');
%         % save temporay divTree:
%         saveas(fig,'E:\DLstg32DivTree.fig');
        curPHLines=PHLines{iSchN};
        curPHMap=PHMap{iSchN};
        curPHparts=PHparts{iSchN};

    %4. Prepare annotated image for manual correction. former:  DL_helper_AnnotSpotOnRawFL
        disp('Drawing...');
        [out,imgLineRaw,indSpot,indLine,indLineStEnd]=drawSptMapFLvImgFLV2(TSCenCoor,Map,FLTS,PHMask);%<------change this back to 1
        %out: [uint8] sz x 3, image for background
        %imgLineRaw: [uint8] sz, thin lines
        %indSpot: [cell] each row is index of a spot square
        %indLine: [cell] each row is index of a line (w=2) starting from current spot. (row idx corresponds to SList)
        %indLineStEnd: [double] n x 4, each row is starting and ending coordinates of line. (row idx corresponds to SList)

    %5. Manually correct map on the annotated image.    
        %tool = ImageFL_N(out,imgLineRaw,indSpot,indLine,indLineStEnd,TSCenCoor,['E:\',posStr,'SchN',num2str(iSchN,'%0.2i'),'.tif']);
        load('D:\Work-Program\Neo2\stills\CrossPatternFLN.mat');
        Cell4test=[{curPHparts},{curPHLines},{curPHMap}];
        Paths4save=[{['E:\',posStr,'SchN',num2str(iSchN,'%0.2i'),'.tif']};{[EXP,'\Y_res\',F(m).name(1:3),'_lkgY.mat']};{[iSchN]}];
        tool = ImageFL_Nv2(out,SList,imgLineRaw,FLTS,indSpot,indLine,indLineStEnd,TSCenCoor,pat,PHMask,Cell4test,MappingMethod,Paths4save);
%         tool.SptMap = reshapeMap(tool.SptMap);

    %6. save results in lkgM and lkgY.    
        saveLkgMY([EXP,'\Y_res\',posStr,'_resY.mat'],[EXP,'\Y_res\',posStr,'_lkgY.mat'],[EXP,'\Y_res\',posStr],...
            tool.SptMap,tool.SList,tool.SptCent,...%curSptMap,curSptList,curSptCent
            tool.CIDinParts,tool.HorLineMapOut,PHparts,...%curPHLines,curPHMap,curPHparts,
            tool.SptLines,tool.SLinesMap,tool.SptParts,iSchN,...%curSLines,curSLinesMap,curSLinesparts,iSchN
            tool.Paths);
        
    end
end

%% nested functions:
function []=saveLkgMY(pthYres,pthY,pthTree,curSptMap,curSptList,curSptCent,curPHLines,curPHMap,curPHparts,curSLines,curSLinesMap,curSLinesparts,iSchN,Paths)
% Yres parameters:
if isempty(curSptMap)
    load(pthYres,'Cand');
    Cand{iSchN}=[];
    save(pthYres,'Cand','-append');
end
% if ~exist('PHLines','var')
%     load(pthM,'lxg');
%     PHLines=cell(size(lxg));
%     PHMap=cell(size(lxg));
%     PHparts=cell(size(lxg));
%     clear('lxg');
% end
% PHLines{1,iSchN}=curPHLines;
% PHMap{1,iSchN}=curPHMap;
% PHparts{1,iSchN}=curPHparts;
% save(pthM,'PHLines','PHMap','PHparts','-append');

% fork parameters: 2-1. spot parameters:
load(pthY,'SptMap','SptList','SptCent');
SptMap{1,iSchN}=curSptMap;
SptList{1,iSchN}=curSptList;
SptCent{1,iSchN}=curSptCent;

% fork parameters: 2-2. spot lines parameters:
load(pthY,'SLines','SLinesMap','SLinesparts');
SLines{1,iSchN}=curSLines;
SLinesMap{1,iSchN}=curSLinesMap;
SLinesparts{1,iSchN}=curSLinesparts;

% save Tree
ImgFLN_drawSPTtree(curSLinesparts,curSLinesMap,0.05,'invisible',[pthTree,'s',num2str(iSchN,'%0.2i'),'_Tree.fig'],Paths);

save(pthY,'SptMap','SptList','SptCent','SLines','SLinesMap','SLinesparts','-append');
end
function []=createLkgY(path,nSchN)
%2-1(spots): SptMap,SptList,SptCent,
%2-2(spot lines): SLines,SLinesMap,SLinesparts,Tree
%additional: Tree
SptMap=cell(1,nSchN);
SptList=cell(1,nSchN);
SptCent=cell(1,nSchN);

SLines=cell(1,nSchN);
SLinesMap=cell(1,nSchN);
SLinesparts=cell(1,nSchN);

save(path,"SptMap","SptList","SptCent","SLines","SLinesMap","SLinesparts");
end
function [SList, TSCenCoor]=genSList(Stats,IntInfo,InfCell,lxg,maxSptN)
%create SList, get centroid and convert coordinates to each tile.
SList=zeros(length(Stats),8);%Format of SList: [Frame, Idx, PercentageInCell, oriN, cellLength, SpotY(from image roof),SpotX(local tile), CoorInCell_Y]
TSCenCoor=round(cell2mat({Stats.Centroid}'));
iFrm=floor(TSCenCoor(:,1)/32)+1;%Frame in SList
SList(:,1)=iFrm';
LocalCenCoor=[TSCenCoor(:,1)-(iFrm-1)*32,TSCenCoor(:,2)];%Spot_Y,Spot_X in SList
SList(:,6:7)=LocalCenCoor;

%see each spot is in which cell contours.
%Format of SList: [Frame, Idx, PercentageInCell, oriN, cellLength, SpotY(from image roof),SpotX(local tile), CoorInCell_Y]
for i=1:length(iFrm)%for each frame. 
    rngCell=[InfCell{iFrm(i),1}{1,1}(:,3),InfCell{iFrm(i),1}{1,1}(:,4)];
    ndx=find(LocalCenCoor(i,2)*ones(size(rngCell,1),1)<rngCell(:,2) & LocalCenCoor(i,2)*ones(size(rngCell,1),1)>=rngCell(:,1), 1);
    if isempty(ndx) || length(ndx)~=1
        ndx=LookUpInPixList(LocalCenCoor(i,:),InfCell{iFrm(i),1}{1,2});
    end
    if isempty(ndx)
        continue;
    end
    SList(i,2)=ndx;

    YmaxCell=max(InfCell{iFrm(i),1}{1,2}{ndx,1}(:,2),[],'all');
    SList(i,5)=YmaxCell-InfCell{iFrm(i),1}{1,1}(ndx,3);
    SList(i,8)=YmaxCell-LocalCenCoor(i,2);
    SList(i,3)=SList(i,8)/SList(i,5);
end

% In short cells, if two spots are too close, they should be a pair.
SList=squeezeSpotsInShortCells(SList);
%SList: [Frame, Idx, PercentageInCell, oriN, cellLength, SpotY(from image roof),SpotX(local tile), CoorInCell_Y), Short cell marker]
[SList,TSCenCoor]=removeStuffedSpt(SList,TSCenCoor,maxSptN);
end
function [out,out2]=removeStuffedSpt(in,in2,maxN)
frms=unique(in(:,1));
frms(isnan(frms))=[];
out=[];
out2=[];
for i=1:length(frms)%for each frame index
    curPart=in(in(:,1)==frms(i),:);
    curIn2=in2(in(:,1)==frms(i),:);
    check=in(in(:,1)==frms(i),2);
    uCh=unique(check);
    del=[];
    for j=1:length(uCh)
        if length(find(check==uCh(j)))>maxN
            del=[del;find(check==uCh(j))];

        end
    end
    if ~isempty(del)
        curPart(del,:)=[];
        curIn2(del,:)=[];
    end
    out=cat(1,out,curPart);
    out2=cat(1,out2,curIn2);
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

function [out]=reshapeMap(in)
if isempty(in)
    out=[];
    return
end
out=nan(length(in),3);
for i=1:length(in)
    A=length(in{i,1});
    if A~=0
        out(i,1:A)=in{i,1};
    end
end
end

function [LineCollect,LinePrime,cellID]=imgPH_simpDecipherLineage(lxg)
% 1. Initialization

%fullize the lxg if the last frame is not all empty.
if ~isempty(find(~cellfun(@(x) isempty(x),lxg(:,end)), 1))
    lxg=[lxg ,cell(size(lxg,1),1)];
end

%remove entries that are not originated from Entry {1,1}.
A=lxg;%{1,nSCh}{1,2};%change here for testing. Take 5 for instance, lxg{1,5}{1,2}
Paridx={};%OUTPUT

init=A{1,1};
inTime=1;
stkDiv=[];%dynamic stack
stkCid=[];%dynamic stack
stkDrCo=[];%dynamic stack
stkDivT=[];%dynamic stack

curdrawCoor=[0,0];%dynamic stack [x,y]
curdivTime=0;%dynamic stack

isDivCyc=false;
isFirstRun=true;

LinePrime=[];% Y of previous ending point.
LineCollect=[];% [x0,y0,x1,y1].
cellID=[];
curLP=[];

figure('Visible','off');
hold on
while ~isempty(stkDiv) || isFirstRun
    %delete the first row of tempStk (that was used in previous round). Begin
    if ~isFirstRun
        stkDiv=stkDiv(2:end,:);
        stkDivT=stkDivT(2:end,:);
        stkDrCo=stkDrCo(2:end,:);
    end
    %delete the first row of tempStk. End
    for t=inTime:size(A,2)%for each time point
        if length(init)==2
            disp(['SChn',num2str(nSCh),' divided at 1st frame!']);
            error('out!')
        end

        if length(A{init,t})==1 && ~isDivCyc
            %<----stopping here indicates that the cell divided at the 1st frame. Use ManCheck to correct.
            init=A{init,t};
            LinePrime=cat(1,LinePrime,curdrawCoor(1,2));
            plot([curdrawCoor(1,1),t],[curdrawCoor(1,2),curdrawCoor(1,2)],'b');
            LineCollect=cat(1,LineCollect,[curdrawCoor(1,1),curdrawCoor(1,2),t,curdrawCoor(1,2)]);
            cellID=cat(1,cellID,[t,init]);
            curdrawCoor=[t,curdrawCoor(1,2)];
            continue;
        elseif length(A{init,t})==2 && ~isDivCyc
            isDivCyc=true;
            stkDiv=[stkDiv;[init,t]];
            stkCid=[stkCid;[init,t+1]];
            init=A{init,t}(1);
            plot([curdrawCoor(1,1),t],[curdrawCoor(1,2),curdrawCoor(1,2)],'b');
            LineCollect=cat(1,LineCollect,[curdrawCoor(1,1),curdrawCoor(1,2),t,curdrawCoor(1,2)]);
            cellID=cat(1,cellID,[t,init]);
            LinePrime=cat(1,LinePrime,curdrawCoor(1,2));
            stkDrCo=[stkDrCo;[t,curdrawCoor(1,2)]];
            stkDivT=[stkDivT;curdivTime];
            curLP=curdrawCoor(1,2);
            curdivTime=curdivTime+1;
            curdrawCoor=[t,curdrawCoor(1,2)+1/(2^curdivTime)];
        elseif length(A{init,t})==1 && isDivCyc
            init_old=init;
            init=A{init,t}(1);
            stkCid=[stkCid;[init,t+1]];
            LinePrime=cat(1,LinePrime,curLP);
            plot([curdrawCoor(1,1),t],[curdrawCoor(1,2),curdrawCoor(1,2)],'b');
            LineCollect=cat(1,LineCollect,[curdrawCoor(1,1),curdrawCoor(1,2),t,curdrawCoor(1,2)]);
            cellID=cat(1,cellID,[t,init_old]);
            curdrawCoor(1,:)=[t,curdrawCoor(1,2)];
            curLP=curdrawCoor(1,2);
            
        elseif length(A{init,t})==2 && isDivCyc
            stkDiv=[stkDiv;[init,t]];
            cellID=cat(1,cellID,[t,init]);
            init=A{init,t}(1);
            Paridx=[Paridx;{stkCid}];
            stkCid=[init,t+1];
            LinePrime=cat(1,LinePrime,curLP);
            plot([curdrawCoor(1,1),t],[curdrawCoor(1,2),curdrawCoor(1,2)],'b');
            LineCollect=cat(1,LineCollect,[curdrawCoor(1,1),curdrawCoor(1,2),t,curdrawCoor(1,2)]);
            
            stkDrCo=[stkDrCo;[t,curdrawCoor(1,2)]];
            stkDivT=[stkDivT;curdivTime];

            curdivTime=curdivTime+1;
            curdrawCoor(1,:)=[t,curdrawCoor(1,2)+1/(2^curdivTime)];
        elseif isempty(A{init,t})

            LinePrime=cat(1,LinePrime,curLP);
            plot([curdrawCoor(1,1),t],[curdrawCoor(1,2),curdrawCoor(1,2)],'b');
            LineCollect=cat(1,LineCollect,[curdrawCoor(1,1),curdrawCoor(1,2),t,curdrawCoor(1,2)]);
            cellID=cat(1,cellID,[t,init]);
            if isempty(stkDiv)
                break;
            end

            init=A{stkDiv(1,1),stkDiv(1,2)}(2);
            if init==-1% divLost case
                while init==-1 && ~isempty(stkDiv)
                    stkDiv(1,:)=[];
                    stkDrCo(1,:)=[];
                    stkDivT(1,:)=[];
                    if ~isempty(stkDiv)
                        init=A{stkDiv(1,1),stkDiv(1,2)}(2);
                        %add the follow to test BEG:
                        inTime=stkDiv(1,2)+1;
                        isDivCyc=true;
                        curLP=stkDrCo(1,2);
                        stkCid=[init,inTime];
                        curdrawCoor=[inTime-1,stkDrCo(1,2)-1/(2^(stkDivT(1)+1))];
                        curdivTime=stkDivT(1)+1;
                        %add the follow to test END.
                    end
                end
                break;
            end

            inTime=stkDiv(1,2)+1;
            isDivCyc=true;
            curLP=stkDrCo(1,2);
            stkCid=[init,inTime];
            curdrawCoor=[inTime-1,stkDrCo(1,2)-1/(2^(stkDivT(1)+1))];
            curdivTime=stkDivT(1)+1;

            break;
        end
    end
    isFirstRun=false;
end
hold off
close('all');
end
function [newParts,parts_CID,Map,F]=drawTreeFromLL(LC,LP,cellID)
marker=true(size(LC,1),1);
marker(1)=false;
parts=[];
StY=[];
parts_CID={};

%1. retrieve horizontal lines [x0,y0,x1,y1];
%acquire: parts
curPart=LC(1,:);
curStY=LP(1);
curCID=cellID(1,:);

while ~isempty(curPart)
    for iIn=1:size(LC,1)
        if isequal(LC(iIn,1:2),curPart(1,3:4))
            curPart(1,3:4)=LC(iIn,3:4);
            marker(iIn)=false;
            curCID=cat(1,curCID,cellID(iIn,:));
            %continue;
        end
    end
    parts=cat(1,parts,curPart);
    StY=cat(1,StY,curStY);
    parts_CID=[parts_CID;{curCID}];
    newStRow=find(marker,1,"first");

    if ~isempty(newStRow)
        curPart=LC(newStRow,:);
        marker(newStRow)=false;
        curStY=LP(newStRow);
        curCID=cellID(newStRow,:);
    else
        curPart=[];
        curStY=[];
        curCID=[];
    end
end

%2. retrive vertical lines.
Map=zeros(size(parts,1),1);
for iH=2:size(parts,1)
    selRow=find(parts(:,4)==StY(iH));
    if ~isempty(selRow) && length(selRow)==1
        Map(iH)=selRow;
    else
        error('not unique in Y.')
    end
end

%3. reshuffle Y
[~,I]=sort(parts(:,2));
newParts=parts;
for iD=1:size(parts,1)
    newParts(I(iD),2)=iD;
    newParts(I(iD),4)=iD;
end

%4. draw lineage tree for PH cell contours
F=figure('visible','off');
set(F,'position',[10,10,900,1200]);
grid on
hold on
% invert the tree plot upside down.
Inv=max(newParts(:,2));
newParts(:,2)=Inv-newParts(:,2);
newParts(:,4)=newParts(:,2);
% draw tree
for iD=1:size(parts,1)
    plot(newParts(iD,[1,3]),newParts(iD,[2,4]),'b','LineWidth',2);
    if Map(iD)~=0
        plot([newParts(iD,1),newParts(Map(iD),3)],[newParts(iD,2),newParts(Map(iD),4)],'b');
    end
end
end

