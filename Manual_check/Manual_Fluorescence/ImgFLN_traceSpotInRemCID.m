function [uit]=ImgFLN_traceSpotInRemCID(sptInCID,SList,initSpt,map,mode)
uit=[];%format: [starting Spt, ending Spt]
prevRat=SList(initSpt,3);
prevSpt=initSpt;
%1. create new candidate spot lines:
if mode(1)==1
    for i=2:length(sptInCID)
        if length(sptInCID{i})==1
            uit=[uit;[prevSpt,sptInCID{i}]];
            prevSpt=sptInCID{i};
            prevRat=(prevRat+SList(sptInCID{i},3))/2;
        elseif isempty(sptInCID{i})
            continue;
        else
            [~,idx]=min(abs(SList(sptInCID{i},3)-prevRat));
            uit=[uit;[prevSpt,sptInCID{i}(idx)]];
            prevSpt=sptInCID{i}(idx);
            prevRat=0.7*prevRat+0.3*SList(sptInCID{i}(idx),3);
        end
    end

    %2. create new candidate spot lines:
    del=[];
    for i=1:size(uit,1)
        if ismember(uit(i,2),map{uit(i,1)})
            del=[del;i];
        end
    end
    uit(del,:)=[];
end
if mode(1)==2
    sptInInitCID=sptInCID{1};
    sptCandIdxUp=find(SList(sptInInitCID,3)>0.5);
    sptCandIdxLow=find(SList(sptInInitCID,3)<=0.5);
    if ~isempty(sptCandIdxUp)
        [~,kup]=min(SList(sptInInitCID(sptCandIdxUp),3)-0.75);
        curCat3=sptInInitCID(sptCandIdxUp(kup));
    else
        curCat3=[];
    end
    if ~isempty(sptCandIdxLow)
        [~,klow]=min(SList(sptInInitCID(sptCandIdxLow),3)-0.25);
        curCat1=sptInInitCID(sptCandIdxLow(klow));
    else
        curCat1=[];
    end

%     sptInInitCID(sptInInitCID==initSpt)=[];
    traceUp=[];
    traceLow=[];
%     if SList(initSpt,3)>0.5%current spot is upper
%         if length(sptInCID{1})==2
%             curCat3=sptInInitCID;
%         elseif length(sptInCID{1})>2
%             [~,idx]=min(SList(sptInInitCID,3)-0.25,[],"all");
%             curCat3=sptInInitCID(idx);
%         else
%             curCat3=initSpt;
%         end
%     else%current spot is lower
%         if length(sptInCID{1})==2
%             curCat1=sptInInitCID;
%         elseif length(sptInCID{1})>2
%             [~,idx]=min(SList(sptInInitCID,3)-0.75,[],"all");
%             curCat1=sptInInitCID(idx);
%         else
%             curCat1=initSpt;
%         end
%     end

    for i=2:length(sptInCID)
        if isempty(sptInCID{i})
            continue;
        else
            idxs=find(SList(sptInCID{i},3)>0.5);
            if ~isempty(idxs)
                [~,minUpIdx]=min(SList(sptInCID{i}(idxs),3)-0.75);
                nextCat3=sptInCID{i}(idxs(minUpIdx));
            else
                nextCat3=[];
            end

            idxs=find(SList(sptInCID{i},3)<=0.5);
            if ~isempty(idxs)
                [~,minUpIdx]=min(SList(sptInCID{i}(idxs),3)-0.25);
                nextCat1=sptInCID{i}(idxs(minUpIdx));
            else
                nextCat1=[];
            end
        end
        if ~isempty(curCat1) && ~isempty(nextCat1)
            traceLow=cat(1,traceLow,[curCat1,nextCat1]);
            curCat1=nextCat1;
            nextCat1=[];
        elseif isempty(curCat1) && ~isempty(nextCat1)
            curCat1=nextCat1;
            nextCat1=[];
        end

        if ~isempty(curCat3) && ~isempty(nextCat3)
            traceUp=cat(1,traceUp,[curCat3,nextCat3]);
            curCat3=nextCat3;
            nextCat3=[];
        elseif isempty(curCat3) && ~isempty(nextCat3)
            curCat3=nextCat3;
            nextCat3=[];
        end
    end

    %2. create new candidate spot lines:
    del=[];
    for i=1:size(traceUp,1)
        if ismember(traceUp(i,2),map{traceUp(i,1)})
            del=[del;i];
        end
    end
    traceUp(del,:)=[];

    del=[];
    for i=1:size(traceLow,1)
        if ismember(traceLow(i,2),map{traceLow(i,1)})
            del=[del;i];
        end
    end
    traceLow(del,:)=[];
    uit=[traceUp;traceLow];


end




end