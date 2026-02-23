function [uit1,uit2]=outputtingL(PHMap,PHLines,CIDinit,NCRmat,win)%In: PHMap,PHLines,CIDinit,NCRmat,win. Out:[NCRDiv,NCRInit] (N x Window double)
uit1=[];
uit2=[];
% for each division:
for iPHL=1:length(PHLines)
    if ~isnan(PHMap(iPHL,2)) || ~isnan(PHMap(iPHL,3))
        if size(PHLines{iPHL,1},1)>=win
            CIDsDiv=PHLines{iPHL,1}(size(PHLines{iPHL,1},1)-win+1:end,:);
        else
            if iPHL==1
                continue;
            end
            CIDsP1=PHLines{iPHL,1};
            if size(PHLines{PHMap(iPHL,1),1},1)>(win-size(PHLines{iPHL,1},1))
                CIDsP2=PHLines{PHMap(iPHL,1),1}((size(PHLines{PHMap(iPHL,1),1},1)-win+size(PHLines{iPHL,1},1)+1):end,:);
            else
                CIDsP2=PHLines{PHMap(iPHL,1),1};
            end
            CIDsDiv=[CIDsP2;CIDsP1];
        end
        cur1=(cellfun(@(x) NCRmat(x(2), x(1)),num2cell(CIDsDiv,2)))';
        if length(cur1)<win
            cur1=cat(2,nan(1,win-length(cur1)),cur1);
        end
        uit1=cat(1,uit1,cur1);
    end
end

% for each initiation:
for iSL=1:size(CIDinit,1)
    if CIDinit(iSL,1)<=win
        continue;
    end
    PHidx=find(cellfun(@(x) ismember(CIDinit(iSL,:), x,'rows'), PHLines));
    stFrm=find(PHLines{PHidx,1}(:,1)==CIDinit(iSL,1));
    linkage=fliplr([1,Tool_getMapLinkage(1,PHidx,PHMap),PHidx]);
    cnt=stFrm-1;
    curLink=1;
    CIDsInit=PHLines{PHidx,1}(1:stFrm-1,:);
    while cnt<win && curLink<=length(linkage)
        curLink=curLink+1;
        cnt=cnt+size(PHLines{linkage(curLink),1},1);
        CIDsInit=[PHLines{linkage(curLink),1};CIDsInit];
    end
    if size(CIDsInit,1)>win
        CIDsInit=CIDsInit(end-win+1:end,:);
    end
    cur2=(cellfun(@(x) NCRmat(x(2), x(1)),num2cell(CIDsInit,2)))';
    if length(cur2)<win
        cur2=cat(2,nan(1,win-length(cur2)),cur2);
    end
    uit2=cat(1,uit2,cur2);
end
end

