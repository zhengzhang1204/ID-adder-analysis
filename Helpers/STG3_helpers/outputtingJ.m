%function [uit1,uit1Cell,uit3]=outputtingJ(oriLines, oriMap, PHLines,PHMap,SLines,SLMap)
function [uit1,uit1Cell,uit3,oriPerCell]=outputtingJ(oriLines, oriMap, PHLines,PHMap)
uit1=[];
uit1Cell={};
%1. retrieve oriFirePercentage in each division cycle.
%reshape the oriLines into ordinary order (map is fullized map.)
% oriLines0=oriLines;
% oriMap0=oriMap;

[oriLines, oriMap]=NormOriLinesAndMap(oriLines, oriMap);
divCycEnterOriIdx=cell(size(PHMap,1),1);
fireInfo=cell(length(oriLines),1);
oriPerCell=getOriPerCell(oriLines,PHLines);

for iPH=2:length(PHLines)%for each PH line, start from 2 intensionally.
%     if iPH==31
%         disp('debugging at Line 15')
%     end
    if isnan(PHMap(iPH,2)) && isnan(PHMap(iPH,3)) 
        continue;
    end
    stCIDs=PHLines{iPH}(1,:);
    enteredOriIdx=find(cellfun(@(x) ~isempty(find(x(:,1)==PHLines{iPH}(1,1) & x(:,2)==PHLines{iPH}(1,2),1)), oriLines));
    exitedOriIdx=find(cellfun(@(x) ~isempty(find(x(:,1)==PHLines{iPH}(end,1) & x(:,2)==PHLines{iPH}(end,2),1)), oriLines));
    if isempty(exitedOriIdx)%meaning there are analyzed PH cells without analyzed spot lines.
        continue;
    end
    OUT=zeros(length(enteredOriIdx),1);
    OUTcell=nan(length(enteredOriIdx),6);
    for j=1:length(enteredOriIdx)% for each entered origin
        collector=[];
        for k=1:length(exitedOriIdx)
            if enteredOriIdx(j)==exitedOriIdx(k)
                if stCIDs(1)==oriLines{enteredOriIdx(j)}(1,1)
                    collector=[collector;[1, exitedOriIdx(k)]];
                    fireInfo{enteredOriIdx(j)}=cat(1,fireInfo{enteredOriIdx(j)},[iPH,1]); 
                else
                    collector=[collector;[0, exitedOriIdx(k)]];
                    fireInfo{enteredOriIdx(j)}=cat(1,fireInfo{enteredOriIdx(j)},[iPH,0]);
                end
                continue;
            end
            uit=Tool_getMapLinkage(enteredOriIdx(j),exitedOriIdx(k),oriMap);
                        
            if uit==0%meaning the mom and daugter has no linkage.
                continue;
            else
                if isempty(uit)%meaning the mom and daughter are directly connected
                    if stCIDs(1)==oriLines{enteredOriIdx(j)}(1,1)
                        collector=[collector;[2,exitedOriIdx(k)];[2,enteredOriIdx(j)]];
                        fireInfo{enteredOriIdx(j)}=cat(1,fireInfo{enteredOriIdx(j)},[iPH,2]);
                    else
                        collector=[collector;[1,exitedOriIdx(k)]];
                        fireInfo{enteredOriIdx(j)}=cat(1,fireInfo{enteredOriIdx(j)},[iPH,1]);
                    end
                else
                    if stCIDs(1)==oriLines{enteredOriIdx(j)}(1,1)
                        collector=[collector;[length(uit)+2,exitedOriIdx(k)];[length(uit)+2,enteredOriIdx(j)]];
                        fireInfo{enteredOriIdx(j)}=cat(1,fireInfo{enteredOriIdx(j)},[iPH,length(uit)+2]);
                    else
                        collector=[collector;[length(uit)+1,exitedOriIdx(k)]];
                        fireInfo{enteredOriIdx(j)}=cat(1,fireInfo{enteredOriIdx(j)},[iPH,length(uit)+1]);
                    end
                end
            end
        end
        OUT(j)=max(collector(:,1),[],1);%error here: some adjacent spot lines has refilling problem.
        OUTcell(j,1)=OUT(j);
        
        if OUT(j)>=1
            sIdx=unique([oriMap(enteredOriIdx(j),2:3),collector(:,2)']);
            sIdx(isnan(sIdx))=[];
            initFrm2nd=unique(cellfun(@(x) x(1,1),oriLines(sIdx)));
            OUTcell(j,2:1+length(initFrm2nd))=initFrm2nd;
        end
        OUTcell(j,5:6)=[PHLines{iPH}(1,1),PHLines{iPH}(end,1)];

    end
    uit1=cat(1,uit1,[length(OUT),length(find(OUT==0))/length(OUT),length(find(OUT==1))/length(OUT),...
        length(find(OUT==2))/length(OUT),length(find(OUT==3))/length(OUT)]);
    divCycEnterOriIdx{iPH}=enteredOriIdx;
    uit1Cell=[uit1Cell;OUTcell];
end


%2. if current oriLine fired X times in a div cycle, how many times did its mother fire?
oriIdxPool=unique(cell2mat(divCycEnterOriIdx));
fireInfo=cellfun(@(x) ZZinitNumMax(x), fireInfo,'UniformOutput',false);
uit2=cell(length(oriLines),2);
for iOr=1:length(oriIdxPool)%<-----------------------------------------------------------------------change this back to 1
    if isnan(oriMap(oriIdxPool(iOr),1))
        continue;
    end
%     if oriIdxPool(iOr)==233
%         disp('D');
%     end
%if the entering oriLine starts from the first PH line, omit it.
    curOriLineStPHIdx= find(cellfun(@(x) ~isempty(find(x(:,1)==oriLines{oriIdxPool(iOr)}(1,1) & x(:,2)==oriLines{oriIdxPool(iOr)}(1,2),1)), PHLines));

    if curOriLineStPHIdx==1 && ~isequal(fireInfo{oriIdxPool(iOr)}(1,:),[2 0])
        continue;
    end
    if size(fireInfo{oriIdxPool(iOr)},1)==1
        if ~isequal(fireInfo{oriIdxPool(iOr)},[2 0])
            uit2{oriIdxPool(iOr),2}=fireInfo{oriIdxPool(iOr)};
            prevPHIdx=PHMap(fireInfo{oriIdxPool(iOr)}(1),1);%mom-divCyc of the divCyc that the current oriLine fires
            momOriLineIdxCand=divCycEnterOriIdx{prevPHIdx,1};%oriLines that enters this mom-divCyc
            for j=1:length(momOriLineIdxCand)
                uit=Tool_getMapLinkage(momOriLineIdxCand(j),oriIdxPool(iOr),oriMap);
                if uit==0%meaning the mom and daugter has no linkage.
                    continue;
                else
                    if isempty(uit)
                        if oriLines{momOriLineIdxCand(j)}(1,1)==PHLines{prevPHIdx}(1,1)% && oriLines{oriIdxPool(iOr),1}(1,1)~=PHLines{fireInfo{oriIdxPool(iOr)}(1),1}(1,1)
                            uit2{oriIdxPool(iOr),1}=[prevPHIdx,2];
                        else
                            uit2{oriIdxPool(iOr),1}=[prevPHIdx,1];
                        end
                    else
                        if oriLines{momOriLineIdxCand(j)}(1,1)==PHLines{prevPHIdx}(1,1)% && oriLines{oriIdxPool(iOr),1}(1,1)~=PHLines{fireInfo{oriIdxPool(iOr)}(1),1}(1,1)
                            uit2{oriIdxPool(iOr),1}=[prevPHIdx,length(uit)+2];
                        else%<---put the additional things to here.
                            uit2{oriIdxPool(iOr),1}=[prevPHIdx,length(uit)+1];
                        end
                    end
                    if oriLines{oriIdxPool(iOr),1}(1,1)==PHLines{fireInfo{oriIdxPool(iOr)}(1),1}(1,1)
                        uit2{oriIdxPool(iOr),1}=uit2{oriIdxPool(iOr),1}-[0,1];
                    end
                end
            end
            if isempty(momOriLineIdxCand)
                uit2{oriIdxPool(iOr),2}=[];
            end
        end
    else
        %decipher Row 1
        if ~isequal(fireInfo{oriIdxPool(iOr)}(1,:),[2 0])
            uit2{oriIdxPool(iOr),2}=fireInfo{oriIdxPool(iOr)}(1,:);
            prevPHIdx=PHMap(fireInfo{oriIdxPool(iOr)}(1,1),1);
            momOriLineIdxCand=divCycEnterOriIdx{prevPHIdx,1};
            for j=1:length(momOriLineIdxCand)
                uit=Tool_getMapLinkage(momOriLineIdxCand(j),oriIdxPool(iOr),oriMap);
                if uit==0%meaning the mom and daugter has no linkage.
                    continue;
                else
%                     mat=fireInfo{momOriLineIdxCand(j)};
%                     uit2{oriIdxPool(iOr),1}=mat(mat(:,1)==prevPHIdx,:);
                    if isempty(uit)
                        if oriLines{momOriLineIdxCand(j)}(1,1)==PHLines{prevPHIdx}(1,1)% && oriLines{oriIdxPool(iOr),1}(1,1)~=PHLines{fireInfo{oriIdxPool(iOr)}(1),1}(1,1)
                            uit2{oriIdxPool(iOr),1}=[prevPHIdx,2];
                        else
                            uit2{oriIdxPool(iOr),1}=[prevPHIdx,1];
                        end
                    else
                        if oriLines{momOriLineIdxCand(j)}(1,1)==PHLines{prevPHIdx}(1,1)% && oriLines{oriIdxPool(iOr),1}(1,1)~=PHLines{fireInfo{oriIdxPool(iOr)}(1),1}(1,1)
                            uit2{oriIdxPool(iOr),1}=[prevPHIdx,length(uit)+2];
                        else
                            uit2{oriIdxPool(iOr),1}=[prevPHIdx,length(uit)+1];
                        end
                    end
                    if oriLines{oriIdxPool(iOr),1}(1,1)==PHLines{fireInfo{oriIdxPool(iOr)}(1),1}(1,1)
                        uit2{oriIdxPool(iOr),1}=uit2{oriIdxPool(iOr),1}-[0,1];
                    end
                end
            end
            if isempty(momOriLineIdxCand) || prevPHIdx==1
                uit2{oriIdxPool(iOr),2}=[];
            end

           
        end
        %decipher Row 2
        addit=zeros(size(fireInfo{oriIdxPool(iOr)},1)-1,4);
        for iR=2:size(fireInfo{oriIdxPool(iOr)},1)
            addit(iR-1,3:4)=fireInfo{oriIdxPool(iOr)}(iR,:);
            prevPHIdx=PHMap(fireInfo{oriIdxPool(iOr)}(iR,1),1);
            iRowx=find(fireInfo{oriIdxPool(iOr)}(:,1)==prevPHIdx,1);
            addit(iR-1,1:2)=fireInfo{oriIdxPool(iOr)}(iRowx,:);
        end
        uit2{oriIdxPool(iOr),1}=cat(1,uit2{oriIdxPool(iOr),1},addit(:,1:2));
        uit2{oriIdxPool(iOr),2}=cat(1,uit2{oriIdxPool(iOr),2},addit(:,3:4));
    end
    %for debugging....remove after debugging<---------------------------------------------
    if ~isempty(uit2{oriIdxPool(iOr),1}) && ~isempty(uit2{oriIdxPool(iOr),2})
        if uit2{oriIdxPool(iOr),1}(1,2)==2 && uit2{oriIdxPool(iOr),2}(1,2)==2
            disp('C1');
        end
        if size(uit2{oriIdxPool(iOr),1},1)==2
            if uit2{oriIdxPool(iOr),1}(2,2)==2 && uit2{oriIdxPool(iOr),2}(2,2)==2
                disp('C2');
            end
        end
        if uit2{oriIdxPool(iOr),1}(1,2)==0 && uit2{oriIdxPool(iOr),2}(1,2)==0
            disp('D1');
        end
        if size(uit2{oriIdxPool(iOr),1},1)==2
            if uit2{oriIdxPool(iOr),1}(2,2)==0 && uit2{oriIdxPool(iOr),2}(2,2)==0
                disp('D2');
            end
        end
    end
end
uit3=cell2mat(uit2);%[prevDivCyc(PHindex), mother oriLineInitTime, current DivCyc(PHindex), current oriLineInitTime]
% Error here means there are singlets entries in uit2
% find those singles. The value in singlet uit2 is [PHpartIdx, NumberOfFire], the row number is OriLineIdx
% run oriMap([OriLineIdx], 1) to retrieve the momOriLineIdx,
% run oriLines{[momOriLineIdx]} to see if there are wrongly refilled spot line.

end

function [LinesOut,MapOut]=NormOriLinesAndMap(LinesIn,MapIn)
%internalMap=[(1:1:size(Lines,1))';ones(size(Lines,1),1)];[(1:1:size(Lines,1))';2*ones(size(Lines,1),1)];
sz=[size(LinesIn,1),2];
LinesOut=[LinesIn(:,1);LinesIn(:,2)];
MapOut=cellfun(@(x) ZZsub2ind(x,sz), MapIn);
[~,IX]=sort(MapOut(:,1));
MapOut=MapOut(IX,2:4);%fulloriMap
end

function [out]=ZZsub2ind(in,sz)
if isnan(in(1))
    out=nan;
else
    out=sub2ind(sz,in(1),in(2));
end
end

function [out]=ZZinitNumMax(in)
if isempty(in)
    out=[];
    return;
end
in=unique(in,'rows');
if size(in,1)==1
    out=in;
    return;
end

R=unique(in(:,1));
out=[];
for i=1:length(R)
    idx=find(in(:,1)==R(i));
    if idx==1
        out=[out;in(idx,:)];
    else
        Z=in(idx,:);
        X=sort(Z(:,2),'descend');
        out=[out;[R(i),X(1)]];
    end
end
end


function [out]=getOriPerCell(oriLines,PHLines)
SZ=[max(cellfun(@(x) max(x(:,1)), PHLines)),max(cellfun(@(x) max(x(:,2)), PHLines))];
out=zeros(SZ);
for i=1:length(oriLines)
    curIdx=sub2ind(SZ,oriLines{i}(:,1),oriLines{i}(:,2));
    out(curIdx)=out(curIdx)+1;
end
end

