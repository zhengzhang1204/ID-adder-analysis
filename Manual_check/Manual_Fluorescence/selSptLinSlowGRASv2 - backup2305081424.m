%This mapping method is for medium below M12, also works for asynchronous.
%Version 2.5 @22.11.29
%Version 2.4a@22.11.02
%Version 2.3 @22.10.17
%Version 2.2 @22.09.29
%Version 2.1 @22.09.23

function [uit,initIdx,errTxt]=selSptLinSlowGRASv2(SptLines,CIDinParts,HorLineMap,sptEndCIDidx)
errTxt=[];
%create a map for selected spot lines.
uit=nan(size(SptLines,1),2);
startCID=cell2mat(cellfun(@(x) x(1,2:3),SptLines,'UniformOutput',false));
initIdx=find(cellfun(@(x) x(1,2)==1, SptLines));

if isempty(initIdx)%very slow growth / no initIdx in cell(1,1)
    for iT=1:size(CIDinParts{1},1)
        initIdx=find(cellfun(@(x) x(1,2)==CIDinParts{1}(iT,1) & x(1,3)==CIDinParts{1}(iT,2), SptLines));
        if isempty(initIdx)
            continue;
        else
            break;
        end
    end
end
marker=false(length(SptLines),1);
SisterPairMap=nan(length(SptLines),2);

[SisterPairMap]=ImgFLN_identSisChrInitV2(SptLines,initIdx,CIDinParts,SisterPairMap);

for iSL=1:length(SptLines)%for each selected spot line:<-------CHANGE BACK TO 1.
%     if iSL==42
%         disp('R');
%     end
    if marker(iSL)
        continue;
    end
    curSptLine=SptLines{iSL,1};
%0. determine if this is the initial thread : start from Cell [1,1]
%     [iR,iC]=find(uit==iSL);
%     if ~isempty(iR) && length(iR)==1
%         avoider=uit(iR,3-iC);
%     else
%         avoider=[];
%     end
%1. see in which cell does this spot line end?
    curSLendPHIdx=find(cellfun(@(x) ~isempty(find(x(:,1)==curSptLine(end,2) & x(:,2)==curSptLine(end,3),1)), CIDinParts));
    if isempty(curSLendPHIdx) || length(curSLendPHIdx)~=1
        error(['SpotLine',num2str(iSL),'ends nowhere!']);
    end
    curSLLrt=getLrtsOfOvlpRegn(CIDinParts{curSLendPHIdx},SptLines(iSL));% LRT of current SLine in its ending PH part.
    curSLLrt=mean(curSLLrt{1},'omitnan');

%2 judge if this spot line has a sister, which resides in the same cell and spans almost the same frame range.
    %get all spot lines that emerges within the interested spot line:

    desRowUp=nan;
    desRowLow=nan;
    desRowUpArch=nan;
    desRowLowArch=nan;
    %use existing map to overwrite the above values: [hasSis,curCat,sisIdx,errTxt]
    tes=find(SisterPairMap(:,1)==iSL,1);
    %if ismember(iSL,initIdx); tes=[]; end %if this is an initial thread, proceed with normal sisterFinding.
    if ~isempty(tes)
        if sptEndCIDidx(tes,1)~=sptEndCIDidx(iSL,1)
            tes=[];
        end
    end
    if isempty(tes)
        [hasSis,curCat,curLrt,sisIdx,errTxt]=ImgFLN_simpleIdentSisChr(SptLines,iSL,CIDinParts,uit);% Method 1: Use existing map to identify sister.
        if ~hasSis && ~ismember(iSL,initIdx)
            [hasSis,curCat,curLrt,sisIdx,errTxt]=ImgFLN_identSisChrv2(SptLines,iSL,CIDinParts);% Method 2: check if this spot line has a sister in the same cell?
        end
        if isempty(sisIdx) && hasSis
            hasSis=false;
            sisIdx=[];
        end
        if hasSis
            SisterPairMap(iSL,1)=sisIdx;
            SisterPairMap(iSL,2)=curCat;
        end
    else
        hasSis=1;
        SisterPairMap(iSL,1)=tes;
        sisIdx=tes;
        if SisterPairMap(tes,2)<10% cat 1 or 3
            curCat=4-SisterPairMap(tes,2);
        elseif SisterPairMap(tes,2)>10 && SisterPairMap(tes,2)<30 % cat 11 or 13
            curCat=24-SisterPairMap(tes,2);
        elseif SisterPairMap(tes,2)>30% cat 31 or 33
            curCat=64-SisterPairMap(tes,2);
        end
        SisterPairMap(iSL,2)=curCat;
        curLrt=curSLLrt;
        errTxt=[];
    end




    %check if curSL and Mother SL end up in the same PHparts:
    momSLEndPHIdx=sptEndCIDidx(uit(:,1)==iSL | uit(:,2)==iSL,1);
    curSLEndPHIdx=sptEndCIDidx(iSL,1);
    if momSLEndPHIdx==curSLEndPHIdx
        momSLIdx=find(uit(:,1)==iSL | uit(:,2)==iSL,1);
        C=uit(momSLIdx,:);
        C(C==iSL)=[];%SL idx of the current SLine's sister
        if ~isnan(SisterPairMap(momSLIdx,1))
            momSisSLEndPHIdx=sptEndCIDidx(SisterPairMap(momSLIdx,1));
            if momSisSLEndPHIdx==curSLEndPHIdx
                if ~ismember(momSLIdx,initIdx)
                    if ~isempty(find(uit(:,1)==momSLIdx,1)); code=3; else; code=1; end
                else
                    if SisterPairMap(momSLIdx,2)==3; code=3; else; code=1; end
                end
                if ~isempty(find(uit(:,1)==iSL,1)); code =code*10+3; else; code=code*10+1; end
                sisIdx=C;
                hasSis=true;
                curLrt=curSLLrt;
                curCat=code;
                SisterPairMap(iSL,2)=code;
            end
        end
    end
    %here, the current spot line may have a sister, but it could be in mother or daughter cells. So, you may have [true, [empty]].
    if ~isempty(errTxt)
        errTxt=['Spot Line starting from Spot#',num2str(SptLines{iSL,1}(1,1)), ' (Frame#', num2str(SptLines{iSL,1}(1,2)),', id#',num2str(iSL),') ends at weird location.'];
        uit=[];
        initIdx=[];
        return;
    end

    if hasSis% the current spot line has a sister:
%HASSIS-1 find spots that overlapping with current spot line, but must be in the same PH part of the cell. These are possible daughter SLines.
        col=[];%stores indice of all possible spot lines the starts within the 2nd half of the current spot line.
        for i=size(curSptLine,1)-floor(size(curSptLine,1)*0.5):size(curSptLine,1)
            jux=find(curSptLine(i,2)==startCID(:,1) & curSptLine(i,3)==startCID(:,2));
            jux=removeSisIdx(iSL,jux,uit);%remove already identified sister idx
            jux=removeDaugEndsInOtherCID(curSptLine(end,2:3),jux,SptLines);%remove candidate SLs that don't exist in the cell where mother SL ends.
            col=cat(1,col,jux');
        end
        col=unique(col);
        if ~isempty(col)
            col(col==sisIdx)=[];
        end
        if ~isempty(col)%if there are spot lines that start within the 2nd half of the current spot line.
            %Lrts=cellfun(@(x) x(1,4), SptLines(col,1));
            Lrts=getLrtsOfOvlpRegn(CIDinParts{curSLendPHIdx},SptLines(col,1));
            Lrts=cellfun(@(x) mean(x,'omitnan'),Lrts);
            if curCat==1
                selCol=find(Lrts<curSLLrt+0.25);
            elseif curCat==3
                selCol=find(Lrts>curSLLrt-0.25);
            elseif curCat==33
                selCol=find(Lrts>max([0.75,curSLLrt-0.125]));
            elseif curCat==31
                selCol=find(Lrts>max([curSLLrt-0.125, 0.5]) & Lrts<min([curSLLrt+0.125,0.75]));
            elseif curCat==13
                selCol=find(Lrts>max([curSLLrt-0.125, 0.25]) & Lrts<min([curSLLrt+0.125,0.5]));
            elseif curCat==11
                selCol=find(Lrts<min([curSLLrt+0.125,0.25]));
            else
                error('testing');
            end
            col=col(selCol);
            Lrts=Lrts(selCol);       
            
            % distribute idx of selected spot lines as daughter spot lines
            
            idxUp=find(Lrts>curSLLrt);
            if length(idxUp)>1
                [~,idxUp2]=min(Lrts(idxUp));
                desRowUp=col(idxUp(idxUp2));
            elseif length(idxUp)==1
                desRowUp=col(idxUp);
            elseif isempty(idxUp)
                desRowUp=nan;
            end

            idxLow=find(Lrts<curSLLrt);
            if length(idxLow)>1
                [~,idxLow2]=max(Lrts(idxLow));
                desRowLow=col(idxLow(idxLow2));
            elseif length(idxLow)==1
                desRowLow=col(idxLow);
            elseif isempty(idxLow)
                desRowLow=nan;
            end
        else
            desRowUp=nan;
            desRowLow=nan;
        end%<- to 2.1a
        desRowUpArch=desRowUp;
        desRowLowArch=desRowLow;

%HASSIS-2 if not all the daughter spot lines are overlapping with mother spot line (i.e., current spot line), then find them in the remaining part of this PH cycle.
        if isnan(desRowLow) || isnan(desRowUp)
            xxx=find(CIDinParts{curSLendPHIdx,1}(:,1)>=curSptLine(end,2),1);
            curCIDrem=CIDinParts{curSLendPHIdx,1}(xxx:end,:);%xxx=0 means the spot end at the PHpart's ending frame
            cand=[];
            for i=1:size(curCIDrem,1)
                jux=find(curCIDrem(i,1)==startCID(:,1) & curCIDrem(i,2)==startCID(:,2));
                if ~isempty(jux)
                    cand=cat(1,cand,jux);
                end
            end
            cand=unique(cand);% possible daughter SLs starting in the remaining part of curCIDpart.
            cand(cand==sisIdx)=[];
            %if the starting frame of sister SL does not overlap with curent SL, means less-dnaA-type
            if ~isnan(sisIdx)
                if SptLines{sisIdx}(1,2)>SptLines{iSL}(end,2)%if the starting frame of sister SL does not overlap with curent SL, means less-dnaA-type
                    %if the CID where current SL ends has only this SL:
                    foundSLIdx=find(cellfun(@(x) ~isempty(find(x(:,2)==SptLines{iSL}(end,2) & x(:,3)==SptLines{iSL}(end,3), 1)),SptLines));
                    if length(foundSLIdx)==1%This SL is the only replicating SL in this CID. It is placed at cell center, but its actual Lrt should be 1/4 or 3/4.
                        if ~isempty(find(uit(:,1)==iSL, 1))%current SL should be the upper SL
                            curCat=3;curSLLrt=0.75;
                        else
                            curCat=1;curSLLrt=0.25;
                        end
                    end
                end
            end


            %Screen out the candidates that has wrong position according to mother's category.
            if ~isempty(cand)
                %Lrts=cellfun(@(x) x(1,4),SptLines(cand,1));
                Lrts=getLrtsOfOvlpRegn(CIDinParts{curSLendPHIdx},SptLines(cand,1));
                Lrts=cellfun(@(x) mean(x,'omitnan'),Lrts);
                if curCat==1
                    selCol=find(Lrts<curSLLrt+0.25);
                elseif curCat==3  
                    selCol=find(Lrts>curSLLrt-0.25);
                elseif curCat==33
                    selCol=find(Lrts>max([0.75,curSLLrt-0.125]));
                elseif curCat==31
                    selCol=find(Lrts>max([curSLLrt-0.125, 0.5]) & Lrts<min([curSLLrt+0.125,0.75]));
                elseif curCat==13
                    selCol=find(Lrts>max([curSLLrt-0.125, 0.25]) & Lrts<min([curSLLrt+0.125,0.5]));
                elseif curCat==11
                    selCol=find(Lrts<min([curSLLrt+0.125,0.25]));
                else
                    error('testing');
                end
                cand=cand(selCol);
                Lrts=Lrts(selCol); 

                % distribute idx of selected spot lines as daughter spot lines
                if length(cand)==2 && isnan(desRowUp) && isnan(desRowLow)%if there are 2 candidates and both desendant slots are empty:
                    if Lrts(1)>Lrts(2)
                        desRowUp=cand(1); desRowLow=cand(2);
                    else
                        desRowUp=cand(2); desRowLow=cand(1);
                    end
                else%Otherwise, fill up Up and Low one by one:
                    idxUp=find(Lrts>curSLLrt);
                    if length(idxUp)>1 
                        [stFrmsUp,shef]=sort(cellfun(@(x) x(1,2),SptLines(cand(idxUp))),'ascend');
                        if stFrmsUp(1)<stFrmsUp(2)-5
                            desRowUp=cand(idxUp(shef(1)));
                        else
                            [~,idxUp2]=min(Lrts(idxUp));desRowUp=cand(idxUp(idxUp2));
                        end
                    elseif length(idxUp)==1
                        desRowUp=cand(idxUp);
                    elseif isempty(idxUp)
                        desRowUp=nan;
                    end
        
                    idxLow=find(Lrts<curSLLrt);
                    if length(idxLow)>1
                        [stFrmsLow,shef]=sort(cellfun(@(x) x(1,2),SptLines(cand(idxLow))),'ascend');
                        if stFrmsLow(1)<stFrmsLow(2)-5
                            desRowLow=cand(idxLow(shef(1)));
                        else
                            [~,idxLow2]=max(Lrts(idxLow));
                            desRowLow=cand(idxLow(idxLow2));
                        end
                    elseif length(idxLow)==1
                        desRowLow=cand(idxLow);
                    elseif isempty(idxLow)
                        desRowLow=nan;
                    end
                end

            end
        end%<-to 2.2a
    
%HASSIS-3 if one or both daughter spot lines are in the daughter cells:
        if isnan(desRowUp) || isnan(desRowLow)
            if curCat==1 || curCat==13 || curCat==11
                refPHidx1=HorLineMap(curSLendPHIdx,3);
            elseif curCat==3 || curCat==31 || curCat==33
                refPHidx1=HorLineMap(curSLendPHIdx,2);

            else
                error('testing');
            end
            if ~isnan(refPHidx1)
                %2.3a-1: list all spot lines that start in PHpart refPHidx1 
                col=[];%stores indice of all possible spot lines the starts within the given daughter div cycle.
                stFrm=[];
                for i=1:size(CIDinParts{refPHidx1,1},1)
                    jux=find(CIDinParts{refPHidx1,1}(i,1)==startCID(:,1) & CIDinParts{refPHidx1,1}(i,2)==startCID(:,2));
                    col=cat(1,col,jux);
                end
                col=unique(col);%list all SLs starting from refPHidx1
                colInf=[cellfun(@(x) x(1,4),SptLines(col,1)),startCID(col,1),col]; %[lrt, start frame,idx]
                if curCat==1 || curCat==3
                    colInfUp=colInf(colInf(:,1)>0.5,:);
                    colInfLow=colInf(colInf(:,1)<0.5,:);
                elseif curCat==13 || curCat==33
                    colInfUp=colInf(colInf(:,1)>0.75 ,:);
                    colInfLow=colInf(colInf(:,1)>0.5 & colInf(:,1)<0.75,:);
                elseif curCat==11 || curCat==31
                    colInfUp=colInf(colInf(:,1)<0.5 & colInf(:,1)>0.25,:);
                    colInfLow=colInf(colInf(:,1)<0.25,:);
                end
                [~,idxUp]=sort(colInfUp(:,2),'ascend');
                colInfUp=colInfUp(idxUp,:);
                [~,idxLow]=sort(colInfLow(:,2),'ascend');
                colInfLow=colInfLow(idxLow,:);
        
                %2.3a-2 find the missing daughter spot line in this given division cycle
                if ~isnan(desRowLow) && isempty(colInfUp)% if low is found previously but up is still empty, find up in the following upper grand-daughter cell
                    refPHidx2=HorLineMap(refPHidx1,2);
                    if ~isnan(refPHidx2)
                        col2=[];
                        for i2=1:size(CIDinParts{refPHidx2,1},1)
                            jux2=find(CIDinParts{refPHidx2,1}(i2,1)==startCID(:,1) & CIDinParts{refPHidx2,1}(i2,2)==startCID(:,2));
                            col2=cat(1,col2,jux2);
                        end 
                        col2=unique(col2);
                        if ~isempty(col2)
                            [~,idx2]=sort(cellfun(@(x) x(1,2),SptLines(col2,1)));
                            desRowUp=col2(idx2);
                        else
                            desRowUp=nan;
                        end
                    end
                elseif ~isnan(desRowUp) && isempty(colInfLow)% if up is found previously but low is still empty, find low in the following lower grand-daughter cell.
                    refPHidx2=HorLineMap(refPHidx1,3);
                    if ~isnan(refPHidx2)
                        col2=[];
                        for i2=1:size(CIDinParts{refPHidx2,1},1)
                            jux2=find(CIDinParts{refPHidx2,1}(i2,1)==startCID(:,1) & CIDinParts{refPHidx2,1}(i2,2)==startCID(:,2));
                            col2=cat(1,col2,jux2);
                        end 
                        col=unique(col);
                        if ~isempty(col)
                            [~,idx2]=sort(cellfun(@(x) x(1,2),SptLines(col,1)));
                            desRowLow=col(idx2);
                        else
                            desRowLow=nan;
                        end
                    end
                else
                    if isnan(desRowLow)
                        if ~isempty(colInfLow)
                            desRowLow=colInfLow(1,3);
                        end
                    end
                    if isnan(desRowUp)
                        if ~isempty(colInfUp)
                            desRowUp=colInfUp(1,3);
                        end
                    end
                end
            end%<-to 2.3a
        end

    else% the current spot line does not have a sister:
%NOSIS-1 find spots that starts within the last half of current spot line
        for i=size(curSptLine,1)-floor(size(curSptLine,1)*0.5):size(curSptLine,1)
            jux=find(curSptLine(i,2)==startCID(:,1) & curSptLine(i,3)==startCID(:,2));
            jux=removeSisIdx(iSL,jux,uit);%remove already identified sister idx, may not be useful here in 2.1b
            jux=removeDaugEndsInOtherCID(curSptLine(end,2:3),jux,SptLines);
            if ~isempty(jux)

                for iJux=1:length(jux)
                    %desLrt=SptLines{jux(iJux),1}(1,4);   
                    Lrts=getLrtsOfOvlpRegn(CIDinParts{curSLendPHIdx},SptLines(jux(iJux),1));
                    desLrt=cellfun(@(x) mean(x,'omitnan'),Lrts);
                    RoO=(curSptLine(end,2)-curSptLine(i,2))/size(curSptLine,1);%the Ratio of the Overlapping period over length of spot line
                    if RoO<0.5% for slow growth condition, the mother-daughter overlapping ratio is set to <0.5
                        if SptLines{jux(iJux)}(SptLines{jux(iJux)}(:,2)==SptLines{iSL}(end,2),4)>SptLines{iSL}(end,4)
                            desRowUp=jux(iJux);
                        else
                            desRowLow=jux(iJux);
                        end
                    end
                end
            end
        end
        desRowUpArch=desRowUp;
        desRowLowArch=desRowLow;

%NOSIS-2 if not all the daughter spot lines are overlapping with mother spot line (i.e., current spot line), then find them in the remaining part of this PH cycle.
        if isnan(desRowLow) || isnan(desRowUp)
            thr=[];
            % find candidate spot lines starting in the remaining part of curCIDpart:
            xxx=find(CIDinParts{curSLendPHIdx,1}(:,1)>curSptLine(end,2),1);
            curCIDrem=CIDinParts{curSLendPHIdx,1}(xxx:end,:);%xxx=0 means the spot end at the PHpart's ending frame
            cand=[];
            for i=1:size(curCIDrem,1)
                jux=find(curCIDrem(i,1)==startCID(:,1) & curCIDrem(i,2)==startCID(:,2));
                if ~isempty(jux)
                    cand=cat(1,cand,jux);
                end
            end
            cand=unique(cand);% possible daughter SLs starting in the remaining part of curCIDpart.

            % assign daughter spot line idx according to candidates' initial LRTs
            vSisIdx=find(uit(:,1)==iSL | uit(:,2)==iSL);
            if ~isempty(vSisIdx) && ~isnan(vSisIdx(1,1)) && size(vSisIdx,1)==1
                vSisIdx=uit(vSisIdx,3-find(uit(vSisIdx,:)==iSL));cand(cand==vSisIdx)=[];
                if ~isnan(vSisIdx)
                    if SptLines{vSisIdx}(1,2)>SptLines{iSL}(end,2)%if the starting frame of sister SL does not overlap with curent SL, means less-dnaA-type
                         %if the CID where current SL ends has only this SL:
                         foundSLIdx=find(cellfun(@(x) ~isempty(find(x(:,2)==SptLines{iSL}(end,2) & x(:,3)==SptLines{iSL}(end,3), 1)),SptLines));
                         if length(foundSLIdx)==1%This SL is the only replicating SL in this CID. It is placed at cell center, but its actual Lrt should be 1/4 or 3/4.
                             if ~isempty(find(uit(:,1)==iSL, 1))%current SL should be the upper SL
                                 curCat=3;thr=0.75;
                            else
                                curCat=1;thr=0.25;
                            end
                        end
                    end
                end
            end

            if ~isempty(cand)
                stFrm=cellfun(@(x) x(1,2), SptLines(cand,1));
                [stFrm,idx]=sort(stFrm,"ascend");
                cand=cand(idx);
                Lrts=getLrtsOfOvlpRegn(CIDinParts{curSLendPHIdx},SptLines(cand,1));
                Lrts=cellfun(@(x) mean(x,'omitnan'),Lrts);
                if isempty(thr)
                    if ~isnan(SisterPairMap(iSL,2))
                        if SisterPairMap(iSL,2)==3
                            thr=0.75;
                        elseif SisterPairMap(iSL,2)==1
                            thr=0.25;
                        elseif SisterPairMap(iSL,2)==2 || isnan(SisterPairMap(iSL,2))
                            thr=0.5;
                        end
                    else
                        if SisterPairMap(iSL,2)==3
                            thr=0.75;
                        elseif SisterPairMap(iSL,2)==1
                            thr=0.25;
                        elseif SisterPairMap(iSL,2)==2 || isnan(SisterPairMap(iSL,2))
                            thr=0.5;
                        end
                    end
                end
                
                if isnan(desRowLow)
                    candLow=cand(Lrts<=thr);
                    stFrmLow=stFrm(Lrts<=thr);
                    if ~isempty(candLow) && length(candLow)>1
                        [~,iLow]=min(stFrmLow);%choose the earliest one.
                        desRowLow=candLow(iLow);
                    elseif ~isempty(candLow) && length(candLow)==1
                        desRowLow=candLow;
                    end
                end
                if isnan(desRowUp)
                    candUp=cand(Lrts>thr);
                    stFrmUp=stFrm(Lrts>thr);
                    if ~isempty(candUp) && length(candUp)>1
                        [~,iUp]=min(stFrmUp);
                        desRowUp=candUp(iUp);
                    elseif ~isempty(candUp) && length(candUp)==1
                        desRowUp=candUp;
                    end
                end
                if length(stFrm)>1%under slow GR, length(cand) is always <2
                    if curSptLine(end,2)+20>=stFrm(2)
                        if Lrts(2)>thr && isnan(desRowUp)
                            desRowUp=cand(2);
                        elseif Lrts(2)<=thr && isnan(desRowLow)
                            desRowLow=cand(2);
                        end
                    end
                end
            end
        end

%NOSIS-3 if one or both daughter spot lines are in the daughter cells:
        [iRow,iCol]=find(uit==iSL,1);%look for its sister.
        if ~isempty(iRow)
            sisInUit=uit(iRow,3-iCol);
            if ~isnan(sisInUit)
                if SptLines{sisInUit}(1,2)<SptLines{iSL}(end,2)
                    iRow=[];iCol=[];sisInUit=[];
                end
            end

        else
            sisInUit=[];
        end

        if isempty(iRow)
            if isnan(desRowLow)
                lowCID=HorLineMap(curSLendPHIdx,3);%lower descendant
                if ~isnan(lowCID)
                    PHdesCID=CIDinParts{lowCID,1};
                    for i=1:size(PHdesCID,1)
                        jux=find(PHdesCID(i,1)==startCID(:,1) & PHdesCID(i,2)==startCID(:,2));
                        if ~isempty(jux)
                            [~,idx]=sort(startCID(jux,1),'ascend');
                            jux=jux(idx);
                            desRowLow=jux(1);
                            break;
                        end
                    end
                end
            end
        
            if isnan(desRowUp)
                upCID=HorLineMap(curSLendPHIdx,2);%upper descendant
                if ~isnan(upCID)
                    PHdesCID=CIDinParts{upCID,1};
                    for i=1:size(PHdesCID,1)
                        jux=find(PHdesCID(i,1)==startCID(:,1) & PHdesCID(i,2)==startCID(:,2));
                        if ~isempty(jux)
                            [~,idx]=sort(startCID(jux,1),'ascend');
                            jux=jux(idx);
                            desRowUp=jux(1);
                            break;
                        end
                    end
                end
            end
        else
            if iCol==2%meaning its descendants are in the lower cell
                lookInPHidx=HorLineMap(curSLendPHIdx,3);%lower descendant index

            elseif iCol==1%meanning its descendants are in the upper cell
                lookInPHidx=HorLineMap(curSLendPHIdx,2);%upper descendant index
            end
            if ~isnan(lookInPHidx)
                PHdesCID=mat2cell(CIDinParts{lookInPHidx,1},ones(1,size(CIDinParts{lookInPHidx,1},1)));
                out=cell2mat(cellfun(@(x) find(x(1)==startCID(:,1) & x(2)==startCID(:,2)),PHdesCID,'UniformOutput',false));
                for i=1:size(out,1)
                    Lrt=getLrtsOfOvlpRegn(CIDinParts{lookInPHidx},SptLines(out(i),1));
                    indexInf=Lrt{1}==Inf;
                    Lrt{1}(indexInf)=[];
                    Lrt=mean(Lrt{1},'all','omitnan');
                    if isnan(desRowUp) && Lrt>0.5
                        desRowUp=out(i);
                    end
                    if isnan(desRowLow) && Lrt<0.5
                        desRowLow=out(i);
                    end
                    if ~isnan(desRowUp) && ~isnan(desRowLow)
                        break;
                    end
                end
            end



        end
    end%<-to: if hasSis
    if ~isnan(desRowLowArch)
        desRowLow=desRowLowArch;
    end
    if ~isnan(desRowUpArch)
        desRowUp=desRowUpArch;
    end
    uit(iSL,1)=desRowUp(1);
    uit(iSL,2)=desRowLow(1);

end
%0-1 in case no initial thread are found by the original definition, 
% register any spotline appears in the first 3 frames or the initial PH cycle(frame<3) as initial threads.

if isempty(initIdx)
    initIdx=unique(find(startCID(:,1)==2 | startCID(:,1)==3));
end
if isempty(initIdx)
    stFrm=4;
    while isempty(initIdx)
        initIdx=find(startCID(:,1)==stFrm);
        stFrm=stFrm+1;
    end
end

%0-2 in case multiple initial threads are found, determine if there are any
%daughter spot line. If yes, remove them.
del=[];
for i=1:length(initIdx)
    jux=find(uit==initIdx(i),1);
    if ~isempty(jux)
       del=[del,i]; 
    end
end
if ~isempty(del)
    initIdx(del)=[];
end



end

%% nested functions:

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
    Lrt(isinf(Lrt))=[];
    uit{i,1}=Lrt;
end
end

function [uit]=removeSisIdx(iSL,jux,map)
uit=[];
for i=1:length(jux)
    idx1=find(map(:,1)==jux(i) & map(:,2)==iSL,1);
    idx2=find(map(:,2)==jux(i) & map(:,1)==iSL,1);
    if isempty(idx1) && isempty(idx2)
        uit=cat(2,uit,jux(i));
    end
end
end
function [uit]=removeDaugEndsInOtherCID(refCID,jux,SLines)
uit=jux;
for i=1:length(jux)
    idx=find(SLines{jux(i),1}(:,2)==refCID(1) & SLines{jux(i),1}(:,3)==refCID(2), 1 );
    if isempty(idx)
        uit(uit==jux(i))=[];
    end
end
end
