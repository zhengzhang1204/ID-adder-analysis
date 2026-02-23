% clc;
% clear;
% load('E:\tempZ.mat');
% MTL_ACjudgeCore(A{1,1}(nFound(inF),:),curCat,C);
% 
function [uit,marker]=MTL_ACjudge(AInf,~,C)
%output format:
%uit: category, possibility, cellID in Cell C, row of spot in Cell C (spotID in cell)
marker=false(1,1);

if size(C,1)==2
% has two Cell Cs. First judge which cell the spot enters, and then the psb.

    uit=cell(1,1);%category, possibility, cellID, row of spot in Cell B (spotID in cell)
%1. which daughter did the spot enter? 
    [Prob1,funCell,curCat]=CatPos4div2(AInf(1,3));%Prob1 stores the psb that which cell the spot enters.(Format:[psb,cellID])
    uitX=[];
    for iC=1:2%for each cell in B, calculate all spots and psb.
        matr=C{Prob1(iC,2),1};
        selFun=funCell{Prob1(iC,2),1};
        uitXtemp=Prob1(iC,1)*selFun(matr(:,3),ones(size(matr,1),1)*AInf(1,3));%calc probability for this spot in Cell B
        uitX=cat(1,uitX,[uitXtemp,ones(size(matr,1),1)*Prob1(iC,2),[1:1:length(uitXtemp)]']);%[prob,cellID,spotID]
    end

%2. judge the possiblity that the spot is in each cell.
    [~,idxC]=sort(uitX(:,1),'descend');
    uitX=uitX(idxC,:);
    %judging:
    if ~isempty(uitX) && uitX(1,1)>0.83 %One spot maps to One spot.
        uit{1}=[curCat,uitX(1,1),uitX(1,2:3)];%category, possibility, cellID, row of spot in Cell B (spotID in cell)
        if size(uitX,1)>1 && uitX(2,1)>0.85 %One spot maps to Two spot.
            uit{1}=[curCat,uitX(1,1),uitX(1,2:3);curCat,uitX(2,1),uitX(2,2:3)];
        end
        marker=true;
    elseif ~isempty(uitX) && sum(find(uitX(:,1)>0.8))>1 %One spot splits into 2 or more.
        selRow=find(uitX(:,1)>0.8);
        uit{1}=[curCat*ones(length(selRow),1),uitX(selRow,1:3)];%category, possibility, cellID, row of spot in Cell B (spotID in cell)
        marker=true;
    else%lost
    % in case ter->ini happens.Then, One cell must be referred as the descendant cell. In this case, pack things and sent to 'B->C'
        if ~isempty(find(Prob1(:,1)>0.85, 1))% && ~isempty(C{find(Prob1(:,1)>0.85, 1),1})
            selRow=find(Prob1(:,1)>0.85, 1);%CellID that the spot likely enters
            CatPsbC=zeros(size(C{selRow,1},1),2);
            if size(C{selRow,1},1)==0
                uit=cell(1,1);
            else
                for iSR=1:size(C{selRow,1},1)
                    
                    CatPsbCTemp=CatPos4(C{selRow,1}(iSR,3));
                    if CatPsbCTemp(1,1)>0.7
                        CatPsbC(iSR,:)=CatPsbCTemp(1,1:2);
                    end
                end
                if isempty(find(CatPsbC(:,2)==2, 1)) && curCat(1,1)==2%confirmed a ter-ini event. Was curCat(1,1)==2 @220116
                    uit{1,1}=[2,1,-1,-1];
                    marker=true;
                end
            end
        else
            %<----------------add here a general block for 'simple-lost' case indicator
            uit=cell(1,1);
        end

    end


elseif size(C,1)==1 && ~isempty(C{1, 1})
% has only one Cell C, judge if the spots in this cell are descendants of the referred spot.
    uit=cell(1,1);
    if ~isempty(C{1, 1})% if one or both of CellA and CellB is 'Short'.
        isContainShort=~isempty(find([AInf(9),C{1,1}(1,9)]==1,1));
    else
        isContainShort=false;
    end
    curCat=CatPos4(AInf(1,3));% Category: from top to bottom: Cat1->3
    [SrchRng,funs]=PredictCat4v2(curCat,AInf(1,3));%Option Two
    
% calculate possibilities between the referred spot in Cell A and target spots in Cell C.
    psb=zeros(3,size(C{1,1},1));% consider psb of category
    psb2=zeros(1,size(C{1,1},1));% does not consider psb of category
    for iC=1:size(C{1,1},1)%for each spot in this Cell C.
        newLrt=C{1,1}(iC,3);%Length Ratio of the iC spot in Cell C.
        jud1=(newLrt*ones(3,1)>SrchRng(:,1) & newLrt*ones(3,1)<=SrchRng(:,2));
        if find(jud1)
            jud2=zeros(3,1);
            for iJud=1:3%intended to be 1:3
                if jud1(iJud)
                    selFun=funs{iJud};
                    %jud2(iJud)=selFun(newLrt);%Option One
                    jud2(iJud)=selFun(newLrt,AInf(1,3));%Option Two
                end
            end
            psb(:,iC)=curCat(:,1).*jud2;
            psb2(1,iC)=max(jud2);
        end
    end
    [psb,psb2]=modPsb12(curCat(:,1),psb,psb2);%Modify 'psb' for the case that the spot in A is 0.5-0.5 in two categories.
    psbU=sort(reshape(psb,[],1),'descend');
    psbU2=sort(reshape(psb2,[],1),'descend');
    if size(C{1,1},1)==1
        distAB=abs(AInf(1,8)-C{1,1}(1,8));
    else
        distAB=10000;
    end

% Judging by the possibility
    if psbU(1)>0.5 && psbU2(1)>0.75 && isContainShort %1->1 mapping. Was 0.28 for Option One
        [sRow,sCol]=find(psb==psbU(1));
        uit{1,1}=[curCat(sRow(1),2),psbU(1),1,sCol(1)];%category, possibility, cellID, row of spot in Cell B (spotID in cell)
        marker=true;
    elseif (psbU(1)>0.65 && psbU2(1)>0.92) || (psbU2(1)>0.56 && (isContainShort || distAB<=10))
        [sRow,sCol]=find(psb==psbU(1),1);
        uit{1,1}=[curCat(sRow,2),psbU(1),1,sCol]; 
        marker=true;
    elseif psbU(1)<psbU(2)+0.37 && psbU(1)>0.8 %1->2 mapping
        [sRow1,sCol1]=find(psb==psbU(1),1);
        [sRow2,sCol2]=find(psb==psbU(2),1);
        uitTmp=[curCat(sRow1,2),psbU(1),1,sCol1;curCat(sRow2,2),psbU(2),1,sCol2];
        if psbU(1)<psbU(3)+0.1 %1->3 mapping
            [sRow3,sCol3]=find(psb==psbU(3));
            uitTmp=cat(1,uitTmp,[curCat(sRow3,2),psbU(3),1,sCol3]);
        end
        %reset the order of the spots in cell:
        [~,idUT]=sort(uitTmp(:,4),'ascend');
        uitTmp=uitTmp(idUT,:);
        uit{1,1}=uitTmp;
        marker=true;
    elseif psbU(1)<=0.75 % no spot in cell B is likely to be the desendant. change:0.65->0.75
 %Two cases here: 1. really lost, 2. ter-ini
 %2.ter-in:
        if  curCat(1,1)>0.65 && curCat(1,2)==2 && ~isempty(C{1,1})%if catA is 2 and catC is not 2, likely to be ter->ini.
            CatPsbC=zeros(size(C{1,1},1),2);
            for iSC=1:size(C{1,1},1)%for each spot in Cell C, retrieve their category and psb
                CatPsbCTemp=CatPos4(C{1,1}(iSC,3));
                if CatPsbCTemp(1,1)>0.7
                    CatPsbC(iSC,:)=CatPsbCTemp(1,1:2);
                end
            end
            if isempty(find(CatPsbC(:,2)==2, 1)) && curCat(1,2)==2%confirmed a ter-ini event
                uit{1,1}=[2,1,-1,-1];
            end
            marker=true;
        else
            uit=cell(1,1);
%1. really lost<-------------------------------------------------------------need to rewrite!
        end
    elseif psbU(1)>1.5*psbU(2) && psbU(1)>0.5%means likely correct. <- need to re-write: according to category punctuality!!!
        %disp('Entered weired range!')
        marker=true;
        [sRow,sCol]=find(psb==psbU(1));
        uit{1,1}=[curCat(sRow,2),psbU(1),1,sCol];%category, possibility, row in B (spotID in cell)
    end
elseif isempty(C)
% No Cell C
    uit=cell(1,1);

elseif isempty(C{1,1})
% No Cell C
    uit=cell(1,1);
end
end
%% nested functions:
function [psb1u,psb2u]=modPsb12(cat1,psb1,psb2)
psb1u=psb1;
psb2u=psb2;
if cat1(2,1)>0.35%means almost equal possibility in two categories.
    
    for iPSB=1:size(psb1,2)%for each spots in Cell C
        if psb1(2,iPSB)==0; continue;end
        if psb1(2,iPSB)<psb1(1,iPSB)
            psb1u(:,iPSB)=psb2(1,iPSB).*[1;0;0];
        else
            psb1u(:,iPSB)=psb2(1,iPSB).*[0;1;0];
        end
        psb2u=psb1u;
    end
end
end