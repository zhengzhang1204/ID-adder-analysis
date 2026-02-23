clc;
clear;
%ONGOING TESTING
% Version 2.1 @23.02.20 PHfig and everything in lkgM are saved using this codes now.

%% user-specified region
EXP='G:\20260124 RdnaA_M30_atc1_dnaNyPet2\';%
cutEndY=249;%<-set this back to 0, can be 230.
cutEndFrame=0;%<-set this back to 0!
minArea=140;
majorL=14;
minY=6; %contours within minY will be elliminated.

%% core
F=dir([EXP,'\M_res\p*_resM.mat']);
for m=15:length(F)%for each pos.<----------------------------------------------change back to 1
    load(fullfile(F(m).folder,F(m).name),'PosMask','uitCIx','uitPH');%masks for each side-channel; X positions of each side-channel; PH image of each side-channel
    if exist(fullfile(F(m).folder,[F(m).name(1:3),'_lkgM.mat']),"file")% see if this pos is analyzed
        load(fullfile(F(m).folder,[F(m).name(1:3),'_lkgM.mat']));
    else
        lxg=cell(1,size(uitCIx,2));
        Infcellv2=cell(size(uitCIx,3),size(uitCIx,2));
    end
    if ~exist('PHLines','var')
        PHLines=cell(1,size(uitCIx,2));
        PHMap=cell(1,size(uitCIx,2));
        PHparts=cell(1,size(uitCIx,2));
    end
    nTile=min(cellfun(@(x) length(x), PosMask),[],2);%how many side channels in this pos?
    PosMaskN=PosMask;%new masks for all side-chennels in this position.
    for k=1:nTile%for each side channel<--------------------------------------change back to 1
        msk=[];
        sw=true;
        sw2=true;
        eptcnt=0;
        out2=zeros(256,32*size(uitPH,3));
        for iFrm=1:size(uitPH,3)%for each frame
            out2(:,(iFrm-1)*32+1:iFrm*32)=uitPH(:,(k-1)*32+1:k*32,iFrm);%compiled PH images for current side-channel.
        end
        
        %If the first frame is empty, omit this side channel
        F1=PosMask{1,1}{k,1};
        if cutEndY~=0
            newF1=F1;
            newF1(cutEndY:end,:)=false;
            L=regionprops(newF1,'Area');
            if isempty(find(cellfun(@(x) x>minArea, {L.Area}),1))
                sw2=false;
            end
        else
            if isempty(find(F1, 1))
                sw2=false;
            end   
        end
        if ~sw2; PosIsManChecked(1,k)=true; continue; end

        %create TS mask:
        msk=cell2mat(cellfun(@(x) x{k,1}, PosMask,'UniformOutput',false));
       
        %apply cut off for Y-end and T-end
        if cutEndY~=0% && ~PosIsManChecked(1,k)%<--------change this back to 'no %'
            msk(cutEndY:end,:)=0;
            msk=modCutOffY(msk,cutEndY);%for those mother cell that has been cut by cutEndY, imopen them.
            PosIsManChecked(1,k)=true;
        end
        if cutEndFrame~=0
            msk(:,cutEndFrame*32+1:end)=0;
        end

        %open GUI to manually draw masks:
        Pos=F(m).name(1:3);
        tex=[Pos,'(',num2str(m,'%0.2i'),')SdChn', num2str(k,'%0.2i')];
        tex4SvFig=[Pos,'S', num2str(k,'%0.2i')];
        xt=(10:32:32*(length(PosMask)-1)+10);
        yt=10*ones(1,length(PosMask));
        Label=num2cell(1:1:length(PosMask));
        TxtIn={tex,xt,yt,Label};
        Msk=imgABv3(mat2gray(out2),1,msk,{tex,xt,yt,Label,Pos(2:3),num2str(k,'%0.2i')},minArea,majorL,minY,cutEndY);
        %decipher the obtained mask:
        mskN=Msk.LabelMasks>0.1;
        nTile=size(mskN,2)/32;
        for j=1:nTile%for each time point
            TempMask=mskN(:,(j-1)*32+1:j*32);
            L=regionprops(TempMask,'Area','MajorAxisLength','PixelIdxList','PixelList');
            for nL=1:length(L)
                if L(nL).Area<minArea || L(nL).MajorAxisLength<majorL
                    TempMask(L(nL).PixelIdxList)=false;
                end
                %check if the region is at edge.
                judge=find(L(nL).PixelList(:,2)<5);
                if length(judge)>5
                    TempMask(L(nL).PixelIdxList)=false;
                end
            end
            PosMaskN{1,j}{k,1}= TempMask;
        end
        PosMask=PosMaskN;

        %3. Draw DivTree according to lxg.
        if isempty(Msk.LineageToSave)
            lxg{k}=[];
            Infcellv2(:,k)=cell(size(Infcellv2,1),1);
            PHLines{1,k}=[];
            PHMap{1,k}=[];
            PHparts{1,k}=[];
            save(fullfile(F(m).folder,[F(m).name(1:3),'_lkgM.mat']),'lxg', 'Infcellv2','PHLines','PHMap','PHparts');
            continue;
        end
        [LineCol,LinePri,CID]=imgPH_simpDecipherLineage(Msk.LineageToSave{2});
        [curPHparts,CIDinParts,HorLineMap,fig]=drawTreeFromLL(LineCol,LinePri,CID);
        clear('LineCol','LinePri','CID');
        if ~exist([EXP,'\Figs\'],"dir")
            mkdir([EXP,'\Figs\']);
        end
        saveas(fig,[EXP,'\Figs\',tex4SvFig,'.fig']);
        
        %4. save 'lxg', 'InfCell','PHLines','PHMap','PHparts' to lkgM.mat
        if ~isempty(Msk.LineageToSave)
            lxg{k}=Msk.LineageToSave(1:3);
            Infcellv2(:,k)=Msk.LineageToSave{5};
            PHLines{1,k}=curPHparts;
            PHMap{1,k}=HorLineMap;
            PHparts{1,k}=CIDinParts;
            save(fullfile(F(m).folder,[F(m).name(1:3),'_lkgM.mat']),'lxg', 'Infcellv2','PHLines','PHMap','PHparts');
        end

        % save 'PosMask','PosIsManChecked' to resM.mat
        save(fullfile(F(m).folder,F(m).name),'PosMask','-append');
        
    end
end

%% nested functions:
function [new]=modCutOffY(new,cutEndY)
patch=false(size(new));
patch(cutEndY-1:end,:)=true;
dif=new&patch;
Ldif=regionprops(dif,'Centroid');

for iDif=1:length(Ldif)
    curFrm=floor(Ldif(iDif).Centroid(1)/32)+1;
    curTile=logical(new(:,(curFrm-1)*32+1:curFrm*32));
    Lcells=regionprops(curTile,'Centroid','PixelIdxList');
    Y=cellfun(@(x) x(2),{Lcells.Centroid});
    [~,iMax]=max(Y);
    if isempty(iMax)
        new(:,(curFrm-1)*32+1:curFrm*32)=double(curTile);
        continue;
    end
    curTile(Lcells(iMax).PixelIdxList)=false;%remove mother cell on curTile
    patchTile=false(256,32);%imopen shape on patchTile
    patchTile(Lcells(iMax).PixelIdxList)=true;
    patchTile=imopen(patchTile,strel("disk",4));
    curTile(patchTile)=true;%patch the curTile
    new(:,(curFrm-1)*32+1:curFrm*32)=double(curTile);%patch the newMask;
end
end

function [uit]=genTilizedPH(pathEXP,pathPOS,SchN,nFrm,CI)
uit=[];  
disp('c');
for iFrm=1:nFrm%for each frame
    J=imread([pathEXP(1:end-5),'M_crop\',pathPOS(1:3),'\',pathPOS(1:3),'_t',num2str(iFrm,'%0.4i'),'.tif']);
    %if ~isnan(CI{1,iFrm}(SchN,1))
    if ~isnan(CI(1,SchN,iFrm))
        %uit=cat(2,uit,mat2gray(J(:,CI{1,iFrm}(SchN,1):CI{1,iFrm}(SchN,2))));
        uit=cat(2,uit,mat2gray(J(:,CI(1,SchN,iFrm):CI(2,SchN,iFrm))));
    else
        uit=cat(2,uit,mat2gray(zeros(256,32)));
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

Map=fullHLM(Map,newParts);%Format: [prevID, sucID1_high, sucID2_low] (IDs are in CIDinParts)
end

function [uit]=fullHLM(in,PHLines)
uit2=nan(size(in,1),2);
for iSL=1:size(in,1)
    xx=find(in==iSL);
    if isempty(xx)
        continue;
    end
    mom_Y=PHLines(iSL,2);
    if length(xx)==2
        if PHLines(xx(1),2)<PHLines(xx(2),2)
            xx=flipud(xx);
        end
        uit2(iSL,:)=xx';
    elseif length(xx)==1
        if PHLines(xx(1),2)>mom_Y
            uit2(iSL,1)=xx;
        else
            uit2(iSL,2)=xx;
        end

    end
end
uit=[in,uit2];
end
