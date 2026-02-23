% This function is to calculate the mean NCR in a time window ('win') and extract the distance from the window right edge to the closest initiatio and division.
function [uit]=outputtingK_edge(PHMap,PHLines,CIDinit,NCRmat,win)%In: PHMap(nx3),PHLines,SLines,NCRmat. Out:[NCR, T2nI, T2nD]

% gathering ininitiation events in each PHline.
PHLinit=cell(size(PHLines));
for iPHL=1:length(PHLines)
    ixxo=intersect(PHLines{iPHL,1},CIDinit,'rows');%CIDinit is the starting CID of each SpotLine.
    PHLinit{iPHL,1}=ixxo(:,1)-PHLines{iPHL,1}(1,1)+1;
end

% calculating CIDs in Window and its time to next initiation (T2nI)
uit=[];
for iPHL=1:length(PHLines)%for each PHLines
    if size(PHLines{iPHL,1},1)<=win
        continue;
    end
    for j=1:size(PHLines{iPHL,1},1)%for each CID in current PHLine
        if j<=size(PHLines{iPHL,1},1)-win%the window is completely inside the current PHLine
            CIDs=PHLines{iPHL,1}(j:j+win-1,:);
            %look for T2nI in current PHlines and the earliest one in any daughter
            idx1=find(PHLinit{iPHL,1}>j,1);%was:j+floor(win/2)-1
            if ~isempty(idx1)
                curT2nI=PHLinit{iPHL,1}(idx1)-j;
            else
                if ~isnan(PHMap(iPHL,2)) && ~isempty(PHLinit{PHMap(iPHL,2),1})
                    idx1u=PHLinit{PHMap(iPHL,2),1}(1)+size(PHLines{iPHL,1},1)-j;
                else
                    idx1u=[];
                end
                if ~isnan(PHMap(iPHL,3)) && ~isempty(PHLinit{PHMap(iPHL,3),1})
                    idx1l=PHLinit{PHMap(iPHL,3),1}(1)+size(PHLines{iPHL,1},1)-j;
                else
                    idx1l=[];
                end
                if ~isempty(idx1u) && ~isempty(idx1l)
                    curT2nI=min([idx1u,idx1l]);
                elseif ~isempty(idx1u)
                    curT2nI=idx1u;
                elseif ~isempty(idx1l) 
                    curT2nI=idx1l;
                else
                    curT2nI=nan;
                end
            end
            
        else%the window is not completely inside the current PHLine
            CIDsP1=PHLines{iPHL,1}(j:end,:);
            %look for T2nI in current PHlines and the earliest one in any daughter
            idx1=find(PHLinit{iPHL,1}>j+floor(win/2)-1,1);
            if ~isempty(idx1)
                curT2nI=PHLinit{iPHL,1}(idx1)-j-floor(win/2);
            else
                curT2nI=[];
            end
            CIDsP2a=[];
            choose=[];
            if ~isnan(PHMap(iPHL,2)) 
                curT2nDU=[];
                if size(PHLines{PHMap(iPHL,2),1},1)>=(win-size(CIDsP1,1)+1)
                    CIDsP2a=PHLines{PHMap(iPHL,2),1}(1:win-size(CIDsP1,1),:);%was 1:win-size(CIDsP1,1)+1
                    if isnan(PHMap(PHMap(iPHL,2),2)) && isnan(PHMap(PHMap(iPHL,2),3))
                        curT2nDU=nan;
                    else
                        curT2nDU=size(PHLines{PHMap(iPHL,2),1},1)+size(PHLines{iPHL,1},1)-j;
                    end
                end
                if ~isempty(PHLinit{PHMap(iPHL,2),1}) && isempty(curT2nI)
                    curT2nIu=PHLinit{PHMap(iPHL,2),1}(1)+size(PHLines{iPHL,1},1)-j;
                else
                    curT2nIu=[];
                end
            end
            CIDsP2b=[];
            if ~isnan(PHMap(iPHL,3)) 
                curT2nDL=nan;
                curT2nIu=[];
                if size(PHLines{PHMap(iPHL,3),1},1)>(win-size(CIDsP1,1)+1)
                    CIDsP2b=PHLines{PHMap(iPHL,3),1}(1:win-size(CIDsP1,1),:);%was 1:win-size(CIDsP1,1)+1
                    if isnan(PHMap(PHMap(iPHL,3),2)) && isnan(PHMap(PHMap(iPHL,3),3))
                        curT2nDL=nan;
                    else
                        curT2nDL=size(PHLines{PHMap(iPHL,3),1},1)+size(PHLines{iPHL,1},1)-j;
                    end
                end
                if ~isempty(PHLinit{PHMap(iPHL,3),1}) && isempty(curT2nI)
                    curT2nIl=PHLinit{PHMap(iPHL,3),1}(1)+size(PHLines{iPHL,1},1)-j;
                else
                    curT2nIl=[];
                end
            end
            if isempty(curT2nI)
                if isempty(curT2nIu) && ~isempty(curT2nIl)
                    CIDs=[CIDsP1;CIDsP2b];
                    curT2nI=curT2nIl;
                    choose=curT2nDL;
                elseif ~isempty(curT2nIu) && isempty(curT2nIl)
                    CIDs=[CIDsP1;CIDsP2a];
                    curT2nI=curT2nIu;
                    choose=curT2nDU;
                elseif isempty(curT2nIu) && isempty(curT2nIl)
                    CIDs=[CIDsP1;CIDsP2b];
                    curT2nI=curT2nIl;
                elseif ~isempty(curT2nIu) && ~isempty(curT2nIl)
                    if curT2nIu<curT2nIl
                        CIDs=[CIDsP1;CIDsP2a];
                        curT2nI=curT2nIu;
                        choose=curT2nDU;
                    else
                        CIDs=[CIDsP1;CIDsP2b];
                        curT2nI=curT2nIl;
                        choose=curT2nDL;
                    end
                end
            end

        end
        curT2nD=size(PHLines{iPHL,1},1)-j;%was:-j-floor(win/2)+1
        if isempty(curT2nI); curT2nI=nan; end
        if curT2nD<0
            if ~isempty(choose)
                curT2nD=choose;
            else
                curT2nD=nan;
            end
        end
        curIDT=size(PHLines{iPHL,1},1);
        curAvgNCR=mean(cellfun(@(x) NCRmat(x(2), x(1)),num2cell(CIDs,2)),'all');

        uit=cat(1,uit,[curAvgNCR,curT2nI,curT2nD,curIDT]);
    end
end
end