%This mapping method is for SLOW growth WT.
function [uit,initIdx,errTxt]=selSptLinFastGRAS(SptLines,CIDinParts,HorLineMap,sptEndCIDidx)
errTxt=[];
%create a map for selected spot lines.
% uit=nan(size(SptLines,1),2);
% startCID=cell2mat(cellfun(@(x) x(1,2:3),SptLines,'UniformOutput',false));
initIdx=find(cellfun(@(x) x(1,2)==1, SptLines));
% marker=false(length(SptLines),1);

% in: SptLines,CIDinParts,HorLineMap. 
% out: selSptLineMap
%selSptMap=selSptLin(SptLines,CIDinParts,HorLineMap);%selSptMap format: 
%each row: [IDX_higher_LRT, IDX_lower_LRT] (IDX: in selSpotLines)

%create a map for selected spot lines.
marker=true(length(SptLines),1);
rng=cell2mat(cellfun(@(x) [x(1,2),x(end,2)], SptLines,'UniformOutput',false));
uit=nan(size(SptLines,1),2);

for iSL=1:length(SptLines)%for each selected spot line:<-------CHANGE BACK TO 1.
    rngN=rng;
    rngN(iSL,:)=[0,0];
    if ~marker(iSL)
        continue;
    end
    
%1. check the ending 3 spots to see if there are very close lines.
%1.1 extract reference info (info of 3+1+3 spots)
    if size(SptLines{iSL,1},1)>=4
        ref1=SptLines{iSL,1}(end-3:end,:);
    else
        ref1=SptLines{iSL,1};
    end
    ref2=extendedRef(SptLines{iSL,1}(end,:),CIDinParts,HorLineMap);
    ref=[ref1;ref2];
    if ~isempty(ref2)
        daugIniCellID_exp=ref2(1,2:3);%expected cellID for daughter spot line.
    else
        daugIniCellID_exp=[];
    end
    if isempty(ref2) || size(ref2,1)<3
        isNoDaugSptLine=true;
    else
        isNoDaugSptLine=false;
    end
    clear('ref1','ref2');

%1.2 Coarse screen for candidates by their existing frames.
    matchCol=[];
    for iRf=1:size(ref,1)
        matchCol=cat(1,matchCol,find(rngN(:,1)<=ref(iRf,2) & rngN(:,2)>=ref(iRf,2)));
    end
    matchCol=unique(matchCol);

    for iRf=1:length(matchCol)
        if SptLines{matchCol(iRf),1}(1,2)<=SptLines{iSL,1}(1,2)
            matchCol(iRf)=nan;
        end
    end
    matchCol(isnan(matchCol))=[];

%1.3 If any candidate exists, for each, analyze the pbb to be daughter chr.
    res=false(1,length(matchCol));
    if ~isempty(matchCol)
        pbb=zeros(size(ref,1),length(matchCol));
        if ~isNoDaugSptLine
            for iRf=1:size(ref,1)%for each spot in ref
                for iMc=1:length(matchCol)%for each candidate spot line
                    % if the initial cell of this spot line is daugIniCellID_exp, 
                    if isequal(SptLines{matchCol(iMc)}(1,2:3),daugIniCellID_exp)
                        sw=1;
                    else
                        sw=0;
                    end
                    % Compare LRTs between the reference and candidates.
                    pbb(iRf,iMc)=getPbb(ref(iRf,:),SptLines{matchCol(iMc)},sw);
                end
            end
        end
%1.4 Judge if there is any close daughter chr.
        for iMc=1:length(matchCol)
            if ~isempty(find(pbb(:,iMc)>0.3, 1)) %ARBITARY VALUE 0.8->0.6->0.58->0.3

                res(1,iMc)=true;
            end
        end
%1.5 If there are 3 candidates, delete the most unlikely one.
        isInitOnExpFrm=false(1,3);
        if length(matchCol(res))==3
            resT=find(res);
            for i=1:3 %intended to be 3
                if isequal(pbb(1:4,resT(i)),zeros(4,1)) &&...
                        isequal(pbb(5:size(pbb,1)-1,resT(i))>0.85,true(size(pbb,1)-5,1))
                    isInitOnExpFrm(i)=true;
                end
            end
        end
        if length(find(isInitOnExpFrm))==2
            res(resT(~isInitOnExpFrm))=false;
        end
    end
    %'true' in res are daughter chrs.
    FrmSt=SptLines{iSL,1}(1,2);
    if ~isempty(find(res, 1))
        selIdx=matchCol(res);
        FilUit=[];
        daug_lrt=zeros(length(selIdx),1);
        for ix=1:length(selIdx)
            daug_lrt(ix,1)=SptLines{selIdx(ix),1}(1,4);
            candFrmSt=SptLines{selIdx(ix),1}(1,2);
            if candFrmSt>FrmSt
                FilUit=cat(2,FilUit,selIdx(ix));
            end
        end
        %Judge Lrt with mom BEG:
        rowInMom=find(SptLines{iSL,1}(:,2)==candFrmSt);
        if isempty(rowInMom) && SptLines{iSL,1}(end,2)+1==candFrmSt
            momLrt=SptLines{iSL, 1}(end,4);
        elseif ~isempty(rowInMom)
            momLrt=SptLines{iSL, 1}(rowInMom,4);
        else
            momLrt=NaN;
        end
        %Judge Lrt with mom END.
        [~,iDL]=sort(daug_lrt,'descend');
        if length(FilUit)==2
            uit(iSL,:)=FilUit(iDL);
        elseif length(FilUit)==1 && ~isnan(momLrt)
            if SptLines{FilUit,1}(1,4)>momLrt
                uit(iSL,1)=FilUit;
            else
                uit(iSL,2)=FilUit;
            end
        elseif length(FilUit)==1 && isnan(momLrt)
            if SptLines{FilUit,1}(1,4)>0.5
                uit(iSL,1)=FilUit;
            else
                uit(iSL,2)=FilUit;
            end
        end
    end
end


end

%% nested functions:
function [uit]=extendedRef(endSpt,PHParts,PHLineMap)
%this nested function extends the current spot line to another 3 frames,
%assuming the LRT does not change. This considers cell division. 
uit=[];
selPHidx=find(cellfun(@(x) ~isempty(find(x(:,1)==endSpt(2) & x(:,2)==endSpt(3), 1) ),PHParts));%in which CID does this spot line end up?
if ~isempty(selPHidx)
    selPHP=PHParts{selPHidx,1};
    iRow=find(selPHP(:,1)==endSpt(2) & selPHP(:,2)==endSpt(3));
    if size(selPHP,1)-iRow>=3
        uit=[zeros(3,1),selPHP(iRow+1:iRow+3,:),endSpt(4)*ones(3,1)];
    else
        uit1=[zeros(size(selPHP,1)-iRow,1),selPHP(iRow+1:end,:),endSpt(4)*ones(size(selPHP,1)-iRow,1)];
        rem=3-size(selPHP,1)+iRow;
        %look for the descencdant of the current PHpart.
        if endSpt(4)<=0.5
            nextPHidx=PHLineMap(selPHidx,3);%
            nextLrt=endSpt(4)*2;
        else
            nextPHidx=PHLineMap(selPHidx,2);%
            nextLrt=1-(1-endSpt(4))*2;
        end
        if isnan(nextPHidx)%no PH descendant
            uit=[];
            return;
        end
        if size(PHParts{nextPHidx,1},1)>=rem
            uit2=[zeros(rem,1),PHParts{nextPHidx,1}(1:rem,:),nextLrt*ones(rem,1)];
        else
            uit2=[zeros(size(PHParts{nextPHidx,1},1),1),PHParts{nextPHidx,1},nextLrt*ones(size(PHParts{nextPHidx,1},1),1)];
        end
        uit=[uit1;uit2];
    end
end
end
function [uit]=getPbb(ref,cand,sw)
idx=find(cand(:,2)==ref(2) & cand(:,3)==ref(3), 1);
if isempty(idx)
    uit=0;
    return;
end

%Core: calc pbb for cand(idx,4) and ref(4)
if sw==0
    if cand(idx,4)>ref(4)-0.13 && cand(idx,4)<ref(4)+0.13
        uit=cos(4*pi*(cand(idx,4)-ref(4)));
    else
        uit=0;
    end
else% if the candidate spot line starts from the next frame cell
    if cand(idx,4)>ref(4)-0.25 && cand(idx,4)<ref(4)+0.25
        uit=cos(pi*(cand(idx,4)-ref(4)));
    else
        uit=0;
    end
end
end
