%This is for calculating ori-related information.
%This is to replace identifyUnconnectedLink.m

%in: PHLines, PHMap, SLines, SLMap
%out: oriLines, oriMap, oriNummat

function [oriLines, oriMap]=getOriInfo(PHLines, PHMap, SLines, SLMap, SptEndFrm)%in: both maps are full map.
%1. get virtual SLines:
%   When the current SLine overlaps with its mother SLine, vSLine is the current SLine.
%   When the current SLine does not overlap with its mother SLine, vSLine starts from the mother's ending frame.
oriLines=cell(length(SLines),3);

for iSL=1:length(SLines)
    oriLines{iSL,3}=SLines{iSL,1}(1,1);
    curSLStEndFrm=[SLines{iSL,1}(1,2),SLines{iSL,1}(end,2)];%current SLine start ending frame.
    curSLStPHIdx=find(cellfun(@(x) ~isempty(find(x(:,1)==SLines{iSL,1}(1,2) & x(:,2)==SLines{iSL,1}(1,3), 1)), PHLines));
    curSLEndPHIdx=find(cellfun(@(x) ~isempty(find(x(:,1)==SLines{iSL,1}(end,2) & x(:,2)==SLines{iSL,1}(end,3), 1)), PHLines));

%1.1 get the Part 2 of CIDs for upper ori: from the start of current SL to the start of upper daughter SLine.
    if ~isnan(SLMap(iSL,2))%the current SLine has an upper daughter:
        daugUpStFrm=SLines{SLMap(iSL,2),1}(1,2);
        daugUpStIdx=find(cellfun(@(x) ~isempty(find(x(:,1)==SLines{SLMap(iSL,2),1}(1,2) & x(:,2)==SLines{SLMap(iSL,2),1}(1,3), 1)), PHLines));
        if daugUpStFrm>=curSLStEndFrm(2)+1%there is a gap between current SLine and its upper daughter SLine.
            if daugUpStIdx==curSLEndPHIdx
                CIDs=[SLines{iSL,1}(:,2:3);...
                PHLines{daugUpStIdx,1}(find(PHLines{daugUpStIdx,1}(:,1)==curSLStEndFrm(2)+1,1):find(PHLines{daugUpStIdx,1}(:,1)==daugUpStFrm,1),:)]; 
%                 if daugUpStFrm==curSLStEndFrm(2)+1
%                     CIDs=CIDs(1:end-1,:);
%                 end
            else
                route=Tool_getMapLinkage(curSLEndPHIdx,daugUpStIdx,PHMap);
                CIDsUpMid=[];
                for iv=1:length(route)
                    CIDsUpMid=[CIDsUpMid;PHLines{route(iv),1}];
                end
                CIDs=[SLines{iSL,1}(:,2:3);...
                    PHLines{curSLEndPHIdx,1}(find(PHLines{curSLEndPHIdx,1}(:,1)==curSLStEndFrm(2)+1,1):end,:);...
                    CIDsUpMid;...
                    PHLines{daugUpStIdx,1}(1:find(PHLines{daugUpStIdx,1}(:,1)==daugUpStFrm,1),:)];
%                 if daugUpStFrm==curSLStEndFrm(2)+1
%                     CIDs=CIDs(1:end-1,:);
%                 end
            end
        else%the current SLine and its upper daughter SLine overlap.
            CIDs=SLines{iSL,1}(1:find(SLines{iSL,1}(:,2)==SLines{SLMap(iSL,2),1}(1,2),1),2:3);
        end
    else%the current SLine does not have an upper daughter SLine.
        if ~isnan(PHMap(curSLEndPHIdx,2))%the PHLine containing the end of current SLine has an upper daughter cell.
            CIDs=[SLines{iSL,1}(1:end-1,2:3);...
                PHLines{curSLEndPHIdx,1}(find(PHLines{curSLEndPHIdx,1}(:,1)==curSLStEndFrm(2),1):end,:);...
                PHLines{PHMap(curSLEndPHIdx,2)}];
        else%the PHLine containing the end of current SLine has an upper daughter cell.
            CIDs=[SLines{iSL,1}(1:end-1,2:3);...
                PHLines{curSLEndPHIdx,1}(find(PHLines{curSLEndPHIdx,1}(:,1)==curSLStEndFrm(2),1):end,:)];
        end
        idxU=find(CIDs(:,1)>SptEndFrm,1);
        if ~isempty(idxU)
            CIDs=CIDs(1:idxU,:);
        end
    end
    oriLines{iSL,1}=CIDs;%(1:end-1,:);
    clear('CIDs');
%1.2 get the Part 2 of CIDs for lower ori: from the start of current SL to the start of lower daughter SLine.
    if ~isnan(SLMap(iSL,3))%the current SLine has an lower daughter:
        daugLowStFrm=SLines{SLMap(iSL,3),1}(1,2);
        daugLowStIdx=find(cellfun(@(x) ~isempty(find(x(:,1)==SLines{SLMap(iSL,3),1}(1,2) & x(:,2)==SLines{SLMap(iSL,3),1}(1,3), 1)), PHLines));
        if daugLowStFrm>=curSLStEndFrm(2)+1%there is a gap between current SLine and its lower daughter SLine.
            if daugLowStIdx==curSLEndPHIdx
                CIDs=[SLines{iSL,1}(:,2:3);...
                PHLines{daugLowStIdx,1}(find(PHLines{daugLowStIdx,1}(:,1)==curSLStEndFrm(2)+1,1):find(PHLines{daugLowStIdx,1}(:,1)==daugLowStFrm,1),:)]; 
%                 if daugUpStFrm==curSLStEndFrm(2)+1
%                     CIDs=CIDs(1:end-1,:);
%                 end
            else
                route=Tool_getMapLinkage(curSLEndPHIdx,daugLowStIdx,PHMap);
                CIDsLowMid=[];
                for iv=1:length(route)
                    CIDsLowMid=[CIDsLowMid;PHLines{route(iv),1}];
                end
                CIDs=[SLines{iSL,1}(:,2:3);...
                    PHLines{curSLEndPHIdx,1}(find(PHLines{curSLEndPHIdx,1}(:,1)==curSLStEndFrm(2)+1,1):end,:);...
                    CIDsLowMid;...
                    PHLines{daugLowStIdx,1}(1:find(PHLines{daugLowStIdx,1}(:,1)==daugLowStFrm,1),:)];
%                 if daugUpStFrm==curSLStEndFrm(2)+1
%                     CIDs=CIDs(1:end-1,:);
%                 end
            end
        else%the current SLine and its lower daughter SLine overlap.
            CIDs=SLines{iSL,1}(1:find(SLines{iSL,1}(:,2)==SLines{SLMap(iSL,3),1}(1,2),1),2:3);
        end
    else%the current SLine does not have an lower daughter SLine.
        if ~isnan(PHMap(curSLEndPHIdx,3))%the PHLine containing the end of current SLine has an upper daughter cell.
            CIDs=[SLines{iSL,1}(1:end-1,2:3);...
                PHLines{curSLEndPHIdx,1}(find(PHLines{curSLEndPHIdx,1}(:,1)==curSLStEndFrm(2),1):end,:);...
                PHLines{PHMap(curSLEndPHIdx,3)}];
        else%the PHLine containing the end of current SLine has a lower daughter cell.
            CIDs=[SLines{iSL,1}(1:end-1,2:3);...
                PHLines{curSLEndPHIdx,1}(find(PHLines{curSLEndPHIdx,1}(:,1)==curSLStEndFrm(2),1):end,:)];
        end
        idxL=find(CIDs(:,1)>SptEndFrm,1);
        if ~isempty(idxL)
            CIDs=CIDs(1:idxL,:);
        end
    end
    oriLines{iSL,2}=CIDs;%(1:end-1,:);
    clear('CIDs');         
end

%2. create oriMap (FULL) from SLMap
oriMap=cell(size(SLMap,1)*2,4);
for iSL=1:size(SLMap,1)
        oriMap{(iSL-1)*2+1,1}=[iSL,1];
        oriMap{iSL*2,1}=[iSL,2];    
%filling up the ancestor info.
    if isnan(SLMap(iSL,1))
        oriMap{(iSL-1)*2+1,2}=nan;
        oriMap{iSL*2,2}=nan;
    else
        index=find(SLMap(SLMap(iSL,1),2:3)==iSL,1);
        oriMap{(iSL-1)*2+1,2}=[SLMap(iSL,1),index];
        oriMap{iSL*2,2}=[SLMap(iSL,1),index];
    end
%filling up the upper descendant info.
    if ~isnan(SLMap(iSL,2))
        oriMap{(iSL-1)*2+1,3}=[SLMap(iSL,2),1];
        oriMap{(iSL-1)*2+1,4}=[SLMap(iSL,2),2];  
    else
        oriMap{(iSL-1)*2+1,3}=nan;
        oriMap{(iSL-1)*2+1,4}=nan; 
    end
%filling up the lower descendant info.
    if ~isnan(SLMap(iSL,3))
        oriMap{iSL*2,3}=[SLMap(iSL,3),1];
        oriMap{iSL*2,4}=[SLMap(iSL,3),2];  
    else
        oriMap{iSL*2,3}=nan;
        oriMap{iSL*2,4}=nan; 
    end
end

%3. create origin MAT (using the OriMap and oriLines)
% for iOL=1:length(oriLines)
%     for k=1:2
%         curCIDs=oriLines{iOL,k};
%         for j=1:size(curCIDs,1)
%             if isnan(oriMat(curCIDs(j,2),curCIDs(j,1)))
%                 oriMat(curCIDs(j,2),curCIDs(j,1))=0;
%             end
%                 oriMat(curCIDs(j,2),curCIDs(j,1))=oriMat(curCIDs(j,2),curCIDs(j,1))+1;
%             
%         end
%     end
% end

% for i=1:size(oriLines,1)
%     for j=1:2%intentionally 1:2
%         curOL=oriLines{i,j};
%         for k=1:size(curOL,1)
%             oriMat(curOL(k,2),curOL(k,1))=oriMat(curOL(k,2),curOL(k,1))+1;
%         end
%     end
% end



%4. for the oriLines that has a daughter oriLine, crop out the last CID.
for i=1:size(oriMap,1)
    if ~isnan(oriMap{i,3})
        idx=oriMap{i,1};
        oriLines{idx(1),idx(2)}=oriLines{idx(1),idx(2)}(1:end-1,:);
    end
end







end