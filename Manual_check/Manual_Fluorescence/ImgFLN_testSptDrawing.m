% Version 2.02 @22.09.16
% Version 2.03 @23.03.10 [not finished yet] 1. Mark wrong spots in red fonts. 
function [SPTparts,selSptMap,selSptMapTemp,selSpotLines,HLMout,errSpots,sw]=ImgFLN_testSptDrawing(SLMapArch,sptPartsInfo,HorLineMap,CIDinParts,PHparts,MappingMethod)
sw=true;
SPTparts=[];
selSptMap=[];
selSpotLines=[];
selSptMapTemp=[];
clc;
disp('Problems:')
selCIDCollect_end=nan(length(sptPartsInfo),3);%PHLines that the ending spot attaches to.

% For each row in horizontal line map, add a column for succeeding PH lines.
if size(HorLineMap,2)==3
    HLMout=HorLineMap;
else
    HorLineMap=fullHLM(HorLineMap,PHparts);%Format: [prevID, sucID1_high, sucID2_low] (IDs are in CIDinParts)
    HLMout=HorLineMap;
end
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
    x1=find(cellfun(@(x) ~isempty(find(x(:,1)==refCellID_beg(1) & x(:,2)==refCellID_beg(2), 1)),CIDinParts), 1);
    if isempty(x1); continue; end
%     x2=find(CIDinParts{x1,1}(:,1)==refCellID_beg(1),1);to be deleted
end
clear('x1','y1','y2','refCellID_beg','refCellID_end');

%% Section 2: For those interested spot lines, estimate its missing parts.
sptSel=find(~isnan(selCIDCollect_end(:,1)));
for i=1:length(sptSel)%for each selected spot line in sptPartsInfo:<-------------change this back to 1
%     if i==8
%         disp('a');
%     end

    tempUit=fullizeSLines(sptPartsInfo{sptSel(i),1},CIDinParts,HorLineMap);
    if ~isempty(tempUit)

        sptPartsInfo{sptSel(i)}=tempUit;
    else

        disp(['Spot line from #', num2str(sptPartsInfo{sptSel(i)}(1,1)),' on Frame ',num2str(sptPartsInfo{sptSel(i)}(1,2)),' has fullize problem.']);
    end
end

%% Section 3 (CORE)：
% 3.1 remove spot lines info that are not in interested cell lines.
selSpotLines=sptPartsInfo(sptSel);
selCIDCollect_end=selCIDCollect_end(sptSel,:);
[selSpotLines,selCIDCollect_end]=sortSLines(selSpotLines,selCIDCollect_end);
clear('sptPartsInfo');

% 3.2 Mapping the lineage of the selected spot lines.
%(selected means the spots lines are in selected PH parts.)
%MappingMethod:
%4: fast condition.
%5: slow condition.

if MappingMethod==4
    [selSptMap,iniThread,errTxt]=selSptLinFastGRAS(selSpotLines,CIDinParts,HorLineMap,selCIDCollect_end);
elseif MappingMethod==5 %for M12 and below
    [selSptMap,iniThread,errTxt]=selSptLinSlowGRASv2(selSpotLines,CIDinParts,HorLineMap,selCIDCollect_end); 
end
if ~isempty(errTxt)
    disp(errTxt);
    sw=false;
    SPTparts=[];    
    return;
end

if ~isempty(SLMapArch)
    selSptMap=SLMapArch;
end

%validate the spot map: every spot line should have only one mother spot line
errSpots=[];
selA=selSptMap;
selA(isnan(selA))=[];
selA=sort(selA);
idx3=find(diff(selA)==0);
if ~isempty(idx3)
    for inx=1:length(idx3)
        errSpots=cat(1,errSpots,selSpotLines{selA(idx3(inx))}(1,1));
        disp(['Spot Line starting from Spot#', num2str(selSpotLines{selA(idx3(inx))}(1,1)),' (Frame#',num2str(selSpotLines{selA(idx3(inx))}(1,2)),', id#',num2str(selA(idx3(inx))),') has two mother lines.']);
        sw=false;
    end
    %<-------------------------------------------------return 1
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
    if isnan(SPTparts(iSL,1))%<------------------------------------------------------return 2
        sw=false;
        errSpots=cat(1,errSpots,selSpotLines{iSL,1}(1,1));
        disp(['Spot Line starting from Spot#', num2str(selSpotLines{iSL,1}(1,1)),' (Frame#',num2str(selSpotLines{iSL,1}(1,2)),', id#',num2str(iSL),') disappeared!']);
    end
end

%validate the SPTparts: every division event should have a corresponding spot line
for iPL=1:size(PHparts,1)
    if ~isnan(HorLineMap(iPL,2)) || ~isnan(HorLineMap(iPL,3))
        curY_PH=PHparts(iPL,2);
        idx=find(selSptLineAttY_pos==curY_PH,1);
        if isempty(idx)%<------------------------------------------------------return 3
%             if iPL~=1% && MappingMethod==5
%                 sw=false;
%                 disp(['Divison on Frame ',num2str(CIDinParts{iPL,1}(1,1)),' Cell #', num2str(CIDinParts{iPL,1}(1,2)),' misses a spotline!']);
%             end
        else
            if selSpotLines{idx,1}(end,2)>CIDinParts{iPL,1}(end,1)
                sw=false;%<------------------------------------------------------return 4
                errSpots=cat(1,errSpots,selSpotLines{idx,1}(1,1));
                disp(['Spot Line starting from Spot#',num2str(selSpotLines{idx,1}(1,1)),' (Frame#',num2str(selSpotLines{idx,1}(1,2)),', id#',num2str(idx),') is too long!'])
            end 
        end

    end
end
if ~sw
    SPTparts=[];
    selSptMapTemp=selSptMap;
    
    selSpotLines=[];
end


end
%% nested functions:
function [SPTLines,SPTAttachPHLines]=getSptTreeData(Map,Lines,PHidComplex,PHLines,PHMap,threads)%Step 3.4

%if ~isempty(threads)% initial thread
if length(threads)>1
    initThreadY=cellfun(@(x) mean(x(1:2,4)), Lines(threads));
    [initThreadY,idx1]=sort(initThreadY,'descend');
    threads=threads(idx1);
    %
    initUp=threads(initThreadY>0.5);
    initLow=threads(initThreadY<=0.5);
    toDivU=zeros(length(initUp),1);
    toDivL=zeros(length(initLow),1);
    %threads=[initUp;initLow];%may be useless
    if length(threads)==2 && (isempty(initUp) || isempty(initLow))
        [~,idxUp]=max(initThreadY);[~,idxLow]=min(initThreadY);
        initUp=threads(idxUp);initLow=threads(idxLow);
    end
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
            if find(isnan(DaugIdx))==1
                curDaugIdx=[Map(curThr,~isnan(DaugIdx)),-1];%curDaugIdx format: [row idx
                % of daugher spot line in SPTLines, +/-1]
            else
                curDaugIdx=[Map(curThr,~isnan(DaugIdx)),1];%curDaugIdx format: [row idx
                % of daugher spot line in SPTLines, +/-1]
            end
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
%         if isnan(daug_PHpID)
%             error('getSptTreeData Wrong!');
%         end

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
if ~isnan(up1)
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
if ~isnan(low1)
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
function [uit]=fullizeSLines(SLraw,PHLines,PHMap)
if (SLraw(end,2)-SLraw(1,2)+1)==size(SLraw,1)%meaning this SLine is full.
    uit=SLraw;
    return;
end

uit=zeros(SLraw(end,2)-SLraw(1,2)+1,4);
uit(1,:)=SLraw(1,:);
curPointer=2;
try
for curIdx=2:size(SLraw,1)
    if SLraw(curIdx,2)-SLraw(curIdx-1,2)==1%meaning no gap
        uit(curPointer,:)=SLraw(curIdx,:);
        curPointer=curPointer+1;
        continue;
    end

    % Now, the gap is between these frames: SLraw(curIdx-1,2) and SLraw(curIdx,2).
    % see if SLraw(curIdx-1,2:3) and SLraw(curIdx,2:3) are the same PHLine.
    gapBegPHLineIdx=find(cellfun(@(x) ~isempty(find(x(:,1)==SLraw(curIdx-1,2) & x(:,2)==SLraw(curIdx-1,3), 1)), PHLines));
    gapEndPHLineIdx=find(cellfun(@(x) ~isempty(find(x(:,1)==SLraw(curIdx,2) & x(:,2)==SLraw(curIdx,3), 1)), PHLines));
    if gapBegPHLineIdx==gapEndPHLineIdx %if the beginging and ending CIDs of the gap are in the same PHLine.
        ratRefill=interp1([SLraw(curIdx-1,2),SLraw(curIdx,2)],[SLraw(curIdx-1,4),SLraw(curIdx,4)],[SLraw(curIdx-1,2)+1:1:SLraw(curIdx,2)-1]);
        idx1=find(PHLines{gapBegPHLineIdx}(:,1)==SLraw(curIdx-1,2)+1,1);
        idx2=find(PHLines{gapBegPHLineIdx}(:,1)==SLraw(curIdx,2)-1,1);
        %insert=PHLines{gapBegPHLineIdx}(idx1:idx2,:);
%         uit(curPointer:curPointer+length(ratRefill)-1,:)=[zeros(length(ratRefill),1),(SLraw(curIdx-1,2)+1:1:SLraw(curIdx,2)-1)',...
%             ones(length(ratRefill),1)*SLraw(curIdx,3),ratRefill'];
        uit(curPointer:curPointer+length(ratRefill)-1,:)=[zeros(length(ratRefill),1),PHLines{gapBegPHLineIdx}(idx1:idx2,:),ratRefill'];
        uit(curPointer+length(ratRefill),:)=SLraw(curIdx,:);
        curPointer=curPointer+length(ratRefill)+1;
    else %if the beginging and ending CIDs of the gap are NOT in the same PHLine.
        PHIdxRoute=[gapBegPHLineIdx,Tool_getMapLinkage(gapBegPHLineIdx,gapEndPHLineIdx,PHMap),gapEndPHLineIdx];%the PH idx of the flanked (refilling) part + init and ending PH parts.
        posCode=zeros(1,length(PHIdxRoute)-1);%posCode: [1] upper daughter, [2] lower daughter
        for i=2:length(PHIdxRoute)%start from 2 intentionally
            posCode(i-1)=find(PHMap(PHIdxRoute(i-1),2:3)==PHIdxRoute(i),1);
        end

        % project the ending ratio to the begining cell:
        ratBeg=SLraw(curIdx-1,4);
        ratEnd=SLraw(curIdx,4);
        ratEndProj=ratEnd;
        for i=length(posCode):-1:1
            if posCode(i)==2
                ratEndProj=ratEndProj/2;
            else
                ratEndProj=(ratEndProj+1)/2;
            end
        end

        % calculate interpolation for rats in between.
        ratRefill=interp1([SLraw(curIdx-1,2),SLraw(curIdx,2)],[ratBeg,ratEndProj],...
            [SLraw(curIdx-1,2)+1:1:SLraw(curIdx,2)-1]);

        % prepare the rescaling map for the refilling part:
        refillFrames=[SLraw(curIdx-1,2)+1:1:SLraw(curIdx,2)-1]';
        gapBegIdxInPH=find(PHLines{gapBegPHLineIdx}(:,1)==SLraw(curIdx-1,2),1);
        gapEndIdxInPH=find(PHLines{gapEndPHLineIdx}(:,1)==SLraw(curIdx,2),1);
        refPointer=size(PHLines{gapBegPHLineIdx},1)-gapBegIdxInPH+1;
        
        % fill up the starting part.
        refBlock=zeros(length(refillFrames),4);%refilling block
        refBlock(1:refPointer-1,2:3)=PHLines{gapBegPHLineIdx}(gapBegIdxInPH+1:end,:);
        refBlock(1:refPointer-1,4)=ratRefill(1:refPointer-1)';
        % fill up the middle part.
        for iRt=2:length(PHIdxRoute)-1%start from 2 intentionally
            refBlock(refPointer:refPointer+size(PHLines{PHIdxRoute(iRt)},1)-1,2:3)=PHLines{PHIdxRoute(iRt)};
            curRoute=posCode(1:iRt-1);
            curRats=ratRefill(refPointer:refPointer+size(PHLines{PHIdxRoute(iRt)},1)-1)';
            for iPC=1:length(curRoute)
                if curRoute(iPC)==2
                    curRats=curRats*2;
                else
                    curRats=curRats*2-1;
                end
            end
            curRats(curRats>=1)=0.99;
            curRats(curRats<=0)=0.01;
            refBlock(refPointer:refPointer+size(PHLines{PHIdxRoute(iRt)},1)-1,4)=curRats;
            refPointer=refPointer+size(PHLines{PHIdxRoute(iRt)},1);
        end
        % fill up the ending part.
        refBlock(refPointer:end,2:3)=PHLines{gapEndPHLineIdx}(1:gapEndIdxInPH-1,:);
        curRats=ratRefill(refPointer:end)';
        for iPC=1:length(posCode)
            if posCode(iPC)==2
                curRats=curRats*2;
            else
                curRats=curRats*2-1;
            end
        end
        curRats(curRats>=1)=0.99;
        curRats(curRats<=0)=0.01;
        refBlock(refPointer:end,4)=curRats;

        %wrap up calculation and insert it into uit
        uit(curPointer:curPointer+length(ratRefill)-1,:)=refBlock;
        uit(curPointer+length(ratRefill),:)=SLraw(curIdx,:);
        curPointer=curPointer+length(ratRefill)+1;
    end
end
catch
    uit=[];
end
end

