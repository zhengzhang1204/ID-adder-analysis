%% This version is specifically for SLOW growing scenario. Featuring:
%1. refilling the missing part in spot lines by averaging.
%2. new drawing method: 

function [figTree]=ImgFLN_testGenTree(sptPartsInfo,HorLineMap,CIDinParts,PHparts,MappingMethod)
figTree=[];
% text='';
clc;
disp('Problems:')

selCIDCollect_end=nan(length(sptPartsInfo),3);%PHLines that the ending spot attaches to.

% For each row in horizontal line map, add a column for succeeding PH lines.
HorLineMap=fullHLM(HorLineMap,PHparts);%Format: [prevID, sucID1_high, sucID2_low] (IDs are in CIDinParts)

%% Section 1: for each spot line, find the cell line in which the ending spot is.
% (sptPartsInfo contains all possible spot lines. This section select spot lines in interested cell lines.)

for iSptRow=1:length(sptPartsInfo)%for each spot line
    refCellID_end=sptPartsInfo{iSptRow,1}(end,2:4);%in which CID the last spot is. [frm, cid, lrt]
    refCellID_beg=sptPartsInfo{iSptRow,1}(1,2:4);%in which CID the first spot is. [frm, cid, lrt]
% check the last spot
    y1=find(cellfun(@(x) ~isempty(find(x(:,1)==refCellID_end(1) & x(:,2)==refCellID_end(2), 1)),CIDinParts));
    if isempty(y1); continue; end
    y2=find(CIDinParts{y1,1}(:,1)==refCellID_end(1),1);
    selCIDCollect_end(iSptRow,:)=[y1,y2,size(CIDinParts{y1,1},1)];%[idx of PH lines in CIDinParts, idx of frame of this PH part]
% check the first spot
    x1=find(cellfun(@(x) ~isempty(find(x(:,1)==refCellID_beg(1) & x(:,2)==refCellID_beg(2), 1)),CIDinParts));
    if isempty(x1); continue; end
    x2=find(CIDinParts{x1,1}(:,1)==refCellID_beg(1),1);
end
clear('x1','x2','y1','y2','refCellID_beg','refCellID_end');

%% Section 2: For those interested spot lines, estimate its missing parts.
sptSel=find(~isnan(selCIDCollect_end(:,1)));
for i=1:length(sptSel)%for each selected spot line in sptPartsInfo:<-------------change this back to 1
    cnct=diff(sptPartsInfo{sptSel(i),1}(:,2));
    [fSLPHI,TOR]=fullSptLinePHInfo(sptPartsInfo{sptSel(i),1},CIDinParts,HorLineMap,sptPartsInfo{sptSel(i),1}(1,1));%Turn Over Row
    if isempty(fSLPHI) && isempty(TOR)
        return;
    end
    if ~isequal(cnct,ones(size(sptPartsInfo{sptSel(i),1},1)-1,1))%means this spot line has missing frame.
        missIdx=find(cnct~=1)+1;
        insertant=cell(length(missIdx)+1,1);
        for j=1:length(missIdx)
            stt=sptPartsInfo{sptSel(i),1}(missIdx(j)-1,2);
            sty=sptPartsInfo{sptSel(i),1}(missIdx(j)-1,4);
            ett=sptPartsInfo{sptSel(i),1}(missIdx(j),2);
            ety=sptPartsInfo{sptSel(i),1}(missIdx(j),4);
            if length(TOR)==3
                jud=[stt<fSLPHI(TOR(1),1) & ett>fSLPHI(TOR(1),1),...
                    stt<fSLPHI(TOR(2),1) & ett>fSLPHI(TOR(2),1),...
                    stt<fSLPHI(TOR(3),1) & ett>fSLPHI(TOR(3),1)];
            elseif length(TOR)==2
                jud=[stt<fSLPHI(TOR(1),1) & ett>fSLPHI(TOR(1),1),...
                    stt<fSLPHI(TOR(2),1) & ett>fSLPHI(TOR(2),1),false];
            elseif length(TOR)==1
                jud=[stt<fSLPHI(TOR(1),1) & ett>fSLPHI(TOR(1),1),false,false];
            end
            
            %if stt<fSLPHI(TOR,1) && ett>fSLPHI(TOR,1)
            trIdx=find(jud);
            if find(jud)
                 nx1=stt+1:1:fSLPHI(TOR(trIdx),1);
                 nx2=fSLPHI(TOR(trIdx),1)+1:1:ett-1;
                 nx=stt+1:1:ett-1;
                if fSLPHI(TOR(trIdx)+1,3)==1
                    reFillY1=interp1([stt,ett],[sty,ety/2],nx);
                    reFillY1=reFillY1(1:length(nx1));
                    reFillY2=interp1([stt,ett],[sty*2,ety],nx);
                    reFillY2=reFillY2((end-length(nx2))+1:end);
                    reFillY=[reFillY1,reFillY2];
                elseif fSLPHI(TOR(trIdx)+1,3)==2

                    reFillY1=interp1([stt,ett],[sty,(1+ety)/2],nx);
                    reFillY1=reFillY1(1:length(nx1));
                    reFillY2=interp1([stt,ett],[(1-sty)*2,ety],nx);
                    reFillY2=reFillY2((end-length(nx2))+1:end);
                    reFillY=[reFillY1,reFillY2];
                else
                    text='';%<-------------------------------------------------return 2
                end
            else
                nx=stt+1:1:ett-1;
                reFillY=interp1([stt,ett],[sty,ety],nx);
            end
            %adding the data in previous block to current insertant.
            if j==1
                prev=sptPartsInfo{sptSel(i),1}(1:missIdx(j)-1,:);
            else
                prev=sptPartsInfo{sptSel(i),1}(missIdx(j-1):missIdx(j)-1,:);
            end
            mx=zeros(length(nx),1);
            for ix=1:length(nx)
                mx(ix)=fSLPHI(fSLPHI(:,1)==nx(ix),2);
            end
            insertant{j,1}=[prev;[zeros(length(nx),1),nx',mx,reFillY']];

            %adding the remaining data block to the end of current insertant.
            if j==length(missIdx)
                insertant{j+1,1}=sptPartsInfo{sptSel(i),1}(missIdx(j):end,:);
            end
        end

        %replace the content in sptPartsInfo by the refilled version.
        out=[];
        for j=1:length(insertant)
            out=cat(1,out,insertant{j,1});
        end
        sptPartsInfo{sptSel(i),1}=out;
    end
end
clear('stt','sty','ett','ety','cnct','fSLPHI','TOR','nx','nx1','nx2','jud');

%% Section 3 (CORE)：
% 3.1 remove spot lines info that are not in interested cell lines.
selSpotLines=sptPartsInfo(sptSel);
selCIDCollect_end=selCIDCollect_end(sptSel,:);
clear('sptPartsInfo');

% 3.2 Mapping the lineage of the selected spot lines.
%(selected means the spots lines are in selected PH parts.)
%MappingMethod:
%4: slow wt.
%5: slow asynchronous.

if MappingMethod==4
    [selSptMap,iniThread]=selSptLinSlowGRWT(selSpotLines,CIDinParts,HorLineMap,selCIDCollect_end);
elseif MappingMethod==5
    [selSptMap,iniThread]=selSptLinSlowGRAS(selSpotLines,CIDinParts,HorLineMap,selCIDCollect_end); 
end

%validate the spot map: every spot line should have only one mother spot line
selA=selSptMap;
selA(isnan(selA))=[];
selA=sort(selA);
idx3=find(diff(selA)==0);
if ~isempty(idx3)
    for inx=1:length(idx3)
        disp(['Spot Line #',num2str(selA(idx3(inx))),' starting from Spot#', num2str(selSpotLines{selA(idx3(inx))}(1,1)),' has two mother lines.'])
    end
    %<-------------------------------------------------return 2
end

%selSptMap format: each row: [IDX_higher_LRT, IDX_lower_LRT] (IDX: in selSpotLines)

% 3.3 For the first 2 rounds, identify the PH indice
PHindice=init2PHidx(PHparts,HorLineMap);
%PHindice format: Index for up1,low1,upup,uplow,lowup,lowlow

% 3.4 Retrieve Raw spot line parts in the format of fram+Y_pos.
[SPTparts,selSptLineAttY_pos]=getSptTreeData(selSptMap,selSpotLines,PHindice,PHparts,HorLineMap,iniThread);%SPTparts format:
%each row: [spot line start frm, Y_pos, spot line ending frm, Y_pos].

%validate the SPTparts: every spot line should appear on this SPTparts map
for iSL=1:size(SPTparts,1)
    if isnan(SPTparts(iSL,1))%<------------------------------------------------------return 3
        disp(['Spot Line starting from Spot#', num2str(selSpotLines{iSL,1}(1,1)),' (Frame#',num2str(selSpotLines{iSL,1}(1,2)),', id#',num2str(iSL),') disappeared!']);
    end
end

%validate the SPTparts: every division event should have a corresponding spot line
for iPL=1:size(PHparts,1)
    if ~isnan(HorLineMap(iPL,2)) || ~isnan(HorLineMap(iPL,3))
        curY_PH=PHparts(iPL,2);
        idx=find(selSptLineAttY_pos==curY_PH,1);
        if isempty(idx)%<------------------------------------------------------return 4
            if iPL~=1
                disp(['Divison on Frame ',num2str(CIDinParts{iPL,1}(1,1)),' Cell #', num2str(CIDinParts{iPL,1}(1,2)),' misses a spotline!']);
            end
        else
            if selSpotLines{idx,1}(end,2)>CIDinParts{iPL,1}(end,1)
                disp(['Spot Line starting from Spot#',num2str(selSpotLines{idx,1}(1,1)),' (Frame#',num2str(selSpotLines{idx,1}(1,2)),', id#',num2str(idx),') is too long!'])
            end 
        end

    end
end


%3.5 Optimize the SPTparts for drawing lines. Prepare additional line for
%drawing.
[figTree]=ImgFLN_drawSPTtree(SPTparts,selSptMap,0.05,'visible');
%uit={CIDinParts,HorLineMap,PHparts,selCIDCollect_beg,selCIDCollect_end,selSpotLines,selSptMap,SPTparts};
end
%% nested functions:

function [SPTLines,SPTAttachPHLines]=getSptTreeData(Map,Lines,PHidComplex,PHLines,PHMap,threads)%Step 3.4

%if ~isempty(threads)% initial thread
if length(threads)~=1
    initThreadY=cellfun(@(x) mean(x(1:2,4)), Lines(threads));
    [initThreadY,idx1]=sort(initThreadY,'descend');
    threads=threads(idx1);
    %
    initUp=threads(initThreadY>0.5);
    initLow=threads(initThreadY<=0.5);
    toDivU=zeros(length(initUp),1);
    toDivL=zeros(length(initLow),1);
    threads=[initUp;initLow];%may be useless
    if length(initUp)==2
        toDivU(1)=22;
        toDivU(2)=21;
    elseif length(initUp)==1
        toDivU(1)=2;
    elseif length(initUp)==3
        %if enters here, I should judge the category by hand
        toDivU(1)=42;
        toDivU(2)=412;
        toDivU(3)=411;
    else
        error('ER0');
    end

    if length(initLow)==2
        toDivL(1)=12;
        toDivL(2)=11;
    elseif length(initLow)==1
        toDivL(1)=1;
    else
        error('ER0');
    end
    toDiv=[toDivU;toDivL];%toDiv is just a code.
    clear('idx1');

%2 identify the Y_pos for each (initial) thread.
    thd_YU=zeros(length(initUp),3);%format: each row:[Y_pos, appendix, PH_ID]
    for iUp=1:length(initUp)
        curCat=toDivU(iUp);
        if curCat==2
            thd_YU(iUp,1)=PHLines(PHidComplex(1),2);
            thd_YU(iUp,3)=PHidComplex(1);
        elseif curCat==22
            if ~isnan(PHidComplex(3))
                thd_YU(iUp,1)=PHLines(PHidComplex(3),2);
                thd_YU(iUp,3)=PHidComplex(3);
            else
                thd_YU(iUp,1)=PHLines(PHidComplex(1),2);
                thd_YU(iUp,2)=1;
                thd_YU(iUp,3)=PHidComplex(1);
            end
        elseif curCat==21
            if ~isnan(PHidComplex(4))
                thd_YU(iUp,1)=PHLines(PHidComplex(4),2);
                thd_YU(iUp,3)=PHidComplex(4);
            else
                thd_YU(iUp,1)=PHLines(PHidComplex(1),2);
                thd_YU(iUp,3)=PHidComplex(1);
                thd_YU(iUp,2)=-1;
            end
%         elseif curCat==42
%         elseif curCat==422
%         elseif curCat==421
        else
            error('ER1');
        end
    end

    thd_YL=zeros(length(initLow),2);
    for iLow=1:length(initLow)
        curCat=toDivL(iLow);
        if curCat==1
            thd_YL(iLow,1)=PHLines(PHidComplex(2),2);
            thd_YL(iLow,3)=PHidComplex(2);
        elseif curCat==12
            if ~isnan(PHidComplex(5))
                thd_YL(iLow,1)=PHLines(PHidComplex(5),2);
                thd_YL(iLow,3)=PHidComplex(5);
            else
                thd_YL(iLow,1)=PHLines(PHidComplex(2),2);
                thd_YL(iLow,2)=1;
                thd_YL(iLow,3)=PHidComplex(2);
            end
        elseif curCat==11
            if ~isnan(PHidComplex(6))
                thd_YL(iLow,1)=PHLines(PHidComplex(6),2);
                thd_YL(iLow,3)=PHidComplex(6);
            else
                thd_YL(iLow,2)=-1;
                thd_YL(iLow,1)=PHLines(PHidComplex(2),2);
                thd_YL(iLow,3)=PHidComplex(2);
            end
        else
            error('ER1');
        end
    end
    thd_Y=[thd_YU;thd_YL];
    clear('iUp','iLow','curCat','thd_YL','thd_YU','toDivU','toDivL');
else
    idx1Thread=find(PHLines(:,1)==0);
    thd_Y=[PHLines(idx1Thread,2),0,idx1Thread];
end

%testing added on 22.04.11_BEG
%for the initial threads, if there has no attached PH Y-line, modify it by
%appendix:
% for ixx=1:size(thd_Y)
%     if thd_Y(ixx,2)~=0
%         thd_Y(ixx,1)=thd_Y(ixx,1)+sign(thd_Y(ixx,2))/(2^abs(thd_Y(ixx,2)));
%         thd_Y(ixx,2)=0;
%     end
% end

%testing ~ END


%3 follow initial thread and identify all daughter's Y_pos
SPTLines=nan(size(Lines,1),4);%stores the spot line location in final tree graph.
SPTAttachPHLines=nan(size(Lines,1),1);%stores the Y_pos of PH that this spot line attaches to.
for iTh=1:length(threads)%for each (initial) thread
    curThr=threads(iTh);%i.e., mom's idx
    tempStk=[];
    isIni=true;
%fill up output info in SPTLines for initial threads:
    SPTLines(threads(iTh),:)=[1,thd_Y(iTh,1),Lines{threads(iTh),1}(end,2),thd_Y(iTh,1)];
    SPTAttachPHLines(threads(iTh),:)=thd_Y(iTh,1);
    curY=thd_Y(iTh,:);%curY format: [Y_pos, appendix, PH_ID]
    while isIni || ~isempty(tempStk) || ~isempty(curY)
    %clear up upon WHILE entry:
        isIni=false;

    %identify daughter thread:
        DaugIdx=Map(curThr,:);%curThr: mom's idx
        DaugIdxZ=DaugIdx(~isnan(DaugIdx));
        
        if isempty(DaugIdxZ)
            if isempty(tempStk)
                break;
            end
            curY=tempStk(1,1:3);
            curDaugIdx=[tempStk(1,4), -1];
            tempStk=tempStk(2:end,:);

        elseif length(DaugIdxZ) ==1
            curDaugIdx=[Map(curThr,~isnan(DaugIdx)),1];%curDaugIdx format: [row idx
            % of daugher spot line in SPTLines, +/-1]
        elseif length(DaugIdxZ) ==2
            curDaugIdx=[Map(curThr,1),1];%current daughter spot line row idx in
            % SPTLines.
            tempStk=[tempStk;[curY,Map(curThr,2)]];%i.e., stored curY
        end

    %follow one daughter thread:
        %a) identify Y_pos and PH_partID of mother thread (momY_pos,mom_PHpID)
        %momY_pos=curY(1)+1/(2^abs(curY(2)))*sign(curY(2));
        momY_pos=curY(1);
        mom_app=curY(2);%mom aligned to: 0: the PH line, +/-0.5, +/-0.25
        mom_PHpID=curY(3);

        %b) look for division event in mon_PHpID,according to 1 or -1, look for daug_PHpID
        if curDaugIdx(2)==1 || (curDaugIdx(2)==-1 && mom_app>0)
            daug_PHpID=PHMap(mom_PHpID,2);
        elseif curDaugIdx(2)==-1
            daug_PHpID=PHMap(mom_PHpID,3);
        end

        %c) if daug_PHpID exists, its Y_pos is daugY_pos
        %   if does not exist, its Y_pos is momY_pos plus appendix.
        if ~isnan(daug_PHpID)
            curDaugY_pos=PHLines(daug_PHpID,2);
            curDaug_app=0;
            curPHpID=daug_PHpID;
            curAttachY_Pos=PHLines(daug_PHpID,2);
        else
            curPHpID=mom_PHpID;
            curDaug_app=(abs(curY(2))+1)*curDaugIdx(2);
            curDaugY_pos=momY_pos+1/(2^abs(curDaug_app))*sign(curDaug_app);
            curAttachY_Pos=momY_pos;
        end
        %d)  store the spot line with Fromat:[start frm, Y_pos, end frm, Y_pos]
        SPTLines(curDaugIdx(1),:)=[Lines{curDaugIdx(1),1}(1,2),...
            curDaugY_pos,Lines{curDaugIdx(1),1}(end,2),curDaugY_pos];
        SPTAttachPHLines(curDaugIdx(1),:)=curAttachY_Pos;

    %assign new conditions for WHILE loop:
        curThr=curDaugIdx(1);
        curY=[curDaugY_pos,curDaug_app,curPHpID];
     end
end
end
function [uit]=init2PHidx(PHLines,PHMap)
%This program returns the idx of 1stDiv and 2ndDiv PH lines.
%initiation
up1=nan;
low1=nan;
upup=nan;
uplow=nan;
lowup=nan;
lowlow=nan;
%calculate round 1
idxR1=PHMap(1,2:3);
idxR1(isnan(idxR1))=[];
if ~isempty(idxR1)%normally, enters this section
    for ix=1:length(idxR1)
        if PHLines(idxR1(ix),2)>PHLines(1,2)
            up1=idxR1(ix);
        else
            low1=idxR1(ix);
        end
    end
end
%calculate round 2
if ~isempty(up1)
    idxR2u=PHMap(up1,2:3);
    idxR2u(isnan(idxR2u))=[];
    if ~isempty(idxR2u)
        for ix=1:length(idxR2u)
            if PHLines(idxR2u(ix),2)>PHLines(up1,2)
                upup=idxR2u(ix);
            else
                uplow=idxR2u(ix);
            end
        end
    end
end
if ~isempty(low1)
    idxR2l=PHMap(low1,2:3);
    idxR2l(isnan(idxR2l))=[];
    if ~isempty(idxR2l)
        for ix=1:length(idxR2l)
            if PHLines(idxR2l(ix),2)>PHLines(low1,2)
                lowup=idxR2l(ix);
            else
                lowlow=idxR2l(ix);
            end
        end
    end
end
%generate output.
uit=[up1,low1,upup,uplow,lowup,lowlow];
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
function [pX,turnOverRow]=fullSptLinePHInfo(sptPart,CIDPart,PHMap,i4disp)
%1. retrieve the initial PH idx
pX=[];
selPHidx=find(cellfun(@(x) x(1,1)<=sptPart(1,2) & x(end,1)>=sptPart(1,2),CIDPart));
for i=1:length(selPHidx)
    sel=find(CIDPart{selPHidx(i),1}(:,1)==sptPart(1,2) & CIDPart{selPHidx(i),1}(:,2)==sptPart(1,3), 1);
    if ~isempty(sel)
        initPHidx=selPHidx(i);
        break;
    end
end

%2 follow the initial PH idx until spot line ends
Frm=sptPart(1,2):1:sptPart(end,2);
uit=zeros(length(Frm),3);% uit format: each row: [Frm, ID, Portion] Portion is where the spot is in the cell.

initCIDP=CIDPart{initPHidx,1};
if size(initCIDP,1)-sel+1<length(Frm)%if this spot line spans two or more PH parts.
    p1=initCIDP(sel:end,:);
    p1X=zeros(size(p1,1),1);
    Frm1EinSP=find(sptPart(:,2)==initCIDP(end,1) & sptPart(:,3)==initCIDP(end,2), 1);
    if isempty(Frm1EinSP)
        Frm1EinSP=find(sptPart(:,2)<initCIDP(end,1), 1, 'last' );
    end

    if median(sptPart(1:Frm1EinSP,4))<0.5
        nxtPHidx=PHMap(initPHidx,3);
        V=1;%code to indicate where the spot is in the cell
    else
        nxtPHidx=PHMap(initPHidx,2);
        V=2;
    end


    rem=length(Frm)-size(p1,1);
    if rem<=size(CIDPart{nxtPHidx,1},1)
        p2=CIDPart{nxtPHidx,1}(1:rem,:);
        p2Y=V*ones(size(p2,1),1);
        turnOverRow=size(p1,1);
        pX=[[p1,p1X];[p2,p2Y]];
    else%means this spot line spans 3 PH lines.
        
        p2=CIDPart{nxtPHidx,1};
        p2Y=V*ones(size(p2,1),1);
        Frm2EinSP=find(sptPart(:,2)==p2(end,1) & sptPart(:,3)==p2(end,2));
        rem2=length(Frm)-size(p1,1)-size(p2,1);%crop the first rem2 rows from the thirt PH part.
        if median(sptPart(Frm2EinSP-4:Frm2EinSP,4))<0.5
            nxt2PHidx=PHMap(nxtPHidx,3);
            V=1;
        else
            nxt2PHidx=PHMap(nxtPHidx,2);
            V=2;
        end

        if rem2<=size(CIDPart{nxt2PHidx,1},1)
            p3=CIDPart{nxt2PHidx,1}(1:rem2,:);
            p3Y=V*ones(size(p3,1),1);
            turnOverRow=[size(p1,1),size(p1,1)+size(p2,1)];
            pX=[[p1,p1X];[p2,p2Y];[p3,p3Y]];
        else%means this spot line spans 4 PH lines.

            p3=CIDPart{nxt2PHidx,1};
            p3Y=V*ones(size(p3,1),1);
            Frm3EinSP=find(sptPart(:,2)==p3(end,1) & sptPart(:,3)==p3(end,2));
            rem3=length(Frm)-size(p1,1)-size(p2,1)-size(p3,1);%crop the first rem3 rows from the thirt PH part.
            if median(sptPart(Frm3EinSP-4:Frm3EinSP,4))<0.5
                nxt3PHidx=PHMap(nxt2PHidx,3);
                V=1;
            else
                nxt3PHidx=PHMap(nxt2PHidx,2);
                V=2;
            end
            if rem3>size(CIDPart{nxt3PHidx,1},1)
                disp(['Spot line starting from #',num2str(i4disp),' spans too long!']);%error return.
                pX=[];
                turnOverRow=[];
                return
            else
                p4=CIDPart{nxt3PHidx,1}(1:rem3,:);
            end
            p4Y=V*ones(size(p4,1),1);
            turnOverRow=[size(p1,1),size(p1,1)+size(p2,1),size(p1,1)+size(p2,1)+size(p3,1)];
            pX=[[p1,p1X];[p2,p2Y];[p3,p3Y];[p4,p4Y]];
        end


    end
    
else%if this spot line exists in the same PH part
    selE=find(initCIDP(:,1)==sptPart(end,2) & initCIDP(:,2)==sptPart(end,3));
    pA=initCIDP(sel:selE,:);
    pX=[pA,zeros(size(pA,1),1)];
    turnOverRow=1;
end
end