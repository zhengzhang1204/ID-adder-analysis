function [uit,is3rdCell]=getSpotLin_FAST(cellIn,lineNumber)
% locate all spots in Cell A, compare with Cell B
A=cellIn{1,1};
B=cellIn{1,2};
C=cellIn{1,3};
divAB=length(A)~=length(B);%if divided, divAB=1.

% possible division cases:
%A->B->C
%1->1->1
%1->1->2
%1->2->2
% possible replication cases:
%ter or lost
%ter -> ini
%1-1 continue
%2x2 continue(two very close spots, may be recognized as 'split;split');

%SList: [Frame, Idx, PercentageInCell, oriN, cellLength, SpotY(from image roof),SpotX(local tile), CoorInCell_Y), Short cell marker]

% A -> B
if ~divAB %if no division from A to B (i.e., length(A) =length(B) = 1)
    uit=cell(size(A{1,1},1),1);
    is3rdCell=false(size(A{1,1},1),1);
    res=cell(size(A{1,1},1),1);
    PassToNxtStp=cell(size(A{1,1},1),1);

    % initiation
    isFound=false(size(A{1,1},1),1);
    lostCnt=zeros(size(A{1,1},1),1);
    isSti=false(size(A{1,1},1),1);
    if ~isempty(B{1, 1})% if one or both of CellA and CellB is 'Short'.
        isContainShort=~isempty(find([A{1,1}(9),B{1,1}(9)]==1,1));
    else
        isContainShort=false;
    end

    %SList: [Frame, Idx, PercentageInCell, oriN, cellLength, SpotY(from image roof),SpotX(local tile), CoorInCell_Y), Short cell marker]
    %starting
    for iA=1:size(A{1,1},1)%for each spot in A
        %starting
        curCat=CatPos4(A{1,1}(iA,3));% Category: from top to bottom: Cat1->3
        %[SrchRng,funs]=PredictCat4(curCat,[1,0]);%Option One
        [SrchRng,funs]=PredictCat4v2(curCat,A{1,1}(iA,3));%Option Two

        if isempty(B{1,1})
            %1. really lost because no spots in Cell B
            isFound(iA,1)=false;
            lostCnt(iA,1)=lostCnt(iA,1)+1;

        else
            %2. B cell has one or more spots
            %retrieve possibility list for all spots in B cell. Store in psb.
            psb=zeros(3,size(B{1,1},1));
            psb2=zeros(3,size(B{1,1},1));%do not consider category possibility
            for iB=1:size(B{1,1},1)%found many spots in this cell, for each.
                if lostCnt(iA,1)==0
                    %[SrchRng,funs]=PredictCat4(curCat,[1,0]);%Option One
                    [SrchRng,funs]=PredictCat4v2(curCat,A{1,1}(iA,3));%Option Two
                end
                newLrt=B{1,1}(iB,3);
                jud1=(newLrt*ones(3,1)>SrchRng(:,1) & newLrt*ones(3,1)<=SrchRng(:,2));
                if find(jud1)
                    jud2=zeros(3,1);
                    for iJud=1:3%intended to be 1:3
                        if jud1(iJud)
                            selFun=funs{iJud};
                            %jud2(iJud)=selFun(newLrt);%Option One
                            jud2(iJud)=selFun(newLrt,A{1,1}(iA,3));%Option Two
                        end
                    end
                    psb(:,iB)=curCat(:,1).*jud2;
                    psb2(:,iB)=jud2;
                end
            end
            [psb,psb2]=modPsb12(curCat(:,1),psb,psb2);%Modify 'psb' for the case that the spot in A is 0.5-0.5 in two categories.
            psbU=sort(reshape(psb,[],1),'descend');
            psbU2=sort(reshape(psb2,[],1),'descend');
            if size(A{1,1},1)==1 && size(B{1,1},1)==1
                distAB=abs(A{1,1}(1,8)-B{1,1}(1,8));
            else
                distAB=1000;
            end

            if (psbU(1)<psbU(2)+0.37 && psbU(1)>0.8) || ((sum(psbU>0.5)>=2) && isContainShort) || is2LikelyInB(psb2)%1->2 mapping
                [sRow1,sCol1]=find(psb==psbU(1));
                [sRow2,sCol2]=find(psb==psbU(2));
                uitTmp=[curCat(sRow1(1),2),psbU(1),1,sCol1(1);curCat(sRow2(1),2),psbU(2),1,sCol2(1)];
                if psbU(1)<psbU(3)+0.1 %1->3 mapping
                    [sRow3,sCol3]=find(psb==psbU(3));
                    %added on 22.04.06 to avoid bugs
                    if length(sRow3)>1
                        sRow3=sRow3(1);
                        sCol3=sCol3(1);
                    end
                    %added ended
                    uitTmp=cat(1,uitTmp,[curCat(sRow3,2),psbU(3),1,sCol3]);
                end
                %reset the order of the spots in cell:
                [~,idUT]=sort(uitTmp(:,4),'ascend');
                uitTmp=uitTmp(idUT,:);
                if size(uitTmp,1)>1
                    [~,iss]=max(uitTmp(:,2));
                    if length(iss)~=1
                        iss=iss(1);
                    end
                    uitTmp=uitTmp(iss,:);
                end
                uit{iA,1}=uitTmp;
                isFound(iA,1)=true;
                
            elseif (psbU(1)>0.5 && psbU2(1)>0.75 && isContainShort) ||...
                    (psbU(1)>0.57 && psbU2(1)>0.82) || ...%0.88->0.82
                    (psbU2(1)>0.56 && (isContainShort || distAB<=10))%1->1 mapping  %0.58->0.56

                %(psbU(1)>10*psbU(2) && psbU(1)>0.65) || removed at 2022.02.14 09:16
                [sRow,sCol]=find(psb==psbU(1));
                if length(sRow)==2; sRow=sRow(1); sCol=sCol(1);end
                uit{iA,1}=[curCat(sRow,2),psbU(1),1,sCol];%category, possibility, cellID, row of spot in Cell B (spotID in cell)
                isFound(iA,1)=true;
%                 if (psbU(1)>10*psbU(2) && psbU(1)>0.65)
%                     disp('Via Here');
%                 end


            elseif psbU(1)<=0.65 && ~isContainShort && psbU2(1)<0.685% no spot in cell B is likely to be the desendant.
                %Two cases here: 1. really lost, 2. ter-ini
                %2. ter-in:
                if  curCat(1,1)>0.85 && curCat(1,2)==2 && ~isempty(C)%if catA is 2 and catB is not 2, likely to be ter->ini.
                    PassToNxtStp{iA,1}=[B;C];
                    isSti(iA,1)=true;
                end
                %1. really lost
                isFound(iA,1)=false;
                lostCnt(iA,1)=lostCnt(iA,1)+1;

            elseif psbU(1)>1.5*psbU(2) && (psbU(1)>0.5 || psbU2(1)>0.7)%means likely correct. <- need to re-write: according to category punctuality!!!
                if psbU2(1)>0.7
                    [sRow,sCol]=find(psb==psbU(1));
                    uit{iA,1}=[curCat(sRow,2),psbU(1),1,sCol];%category, possibility, row in B (spotID in cell)
                else
                    if psbU2(1)==psbU2(2)
                        [sRow,sCol]=find(psb==psbU(1));
                    else
                        [sRow,sCol]=find(psb2==psbU2(1));
                    end
                    uit{iA,1}=[curCat(sRow,2),psbU2(1),1,sCol];%category, possibility, row in B (spotID in cell) 
                end
                isFound(iA,1)=true;
            end
        end
    end
else%if A divides into two Bs (i.e., length(A) =1, length(B) = 2)
    uit=cell(size(A{1,1},1),1);
    is3rdCell=false(size(A{1,1},1),1);
    PassToNxtStp=cell(size(A{1,1},1),1);
    isFound=false(size(A{1,1},1),1);
    lostCnt=zeros(size(A{1,1},1));
    isSti=false(size(A{1,1},1),1);
    isContainShort=true;%because Cell Bs are new-born, they must be short cells
    %Since Cell Bs are new born, they are short, category may not accurate,
    %spots that are too close may be a pair.
    for iA=1:size(A{1,1},1)%for each spot in Cell A
        %initiation.
        isFound(iA,1)=false;
    
        %starting: which daughter did the spot enter? Prob1
        [Prob1,funCell,curCat]=CatPos4div2(A{1,1}(iA,3));
        uitX=[];

        for iB=1:2%for each cell in B
            matr=B{Prob1(iB,2),1};
            selFun=funCell{Prob1(iB,2),1};
            uitXtemp=Prob1(iB,1)*selFun(matr(:,3),ones(size(matr,1),1)*A{1,1}(iA,3));%calc probability for this spot in Cell B
            uitX=cat(1,uitX,[uitXtemp,ones(size(matr,1),1)*Prob1(iB,2),[1:1:length(uitXtemp)]']);%[prob,cellID,spotID]
        end
        %judge the possiblity that the spot is in each cell.
        [~,idxB]=sort(uitX(:,1),'descend');
        uitX=uitX(idxB,:);
        %judging:
        if ~isempty(uitX) && sum(find(uitX(:,1)>0.8))>1 %One spot splits into 2 or more.
            selRow=find(uitX(:,1)>0.8);
            if length(selRow)>1
                [~,iss]=max(uitX(:,1));
                selRow=iss;
            end
            uit{iA}=[curCat*ones(length(selRow),1),uitX(selRow,1:3)];%category, possibility, cellID, row of spot in Cell B (spotID in cell)
            isFound(iA,1)=true;
        elseif ~isempty(uitX) && (uitX(1,1)>0.85 || (uitX(1,1)>0.65 && isContainShort))%One spot maps to One spot.
            uit{iA}=[curCat,uitX(1,1),uitX(1,2:3)];%category, possibility, cellID, row of spot in Cell B (spotID in cell)
            isFound(iA,1)=true;

        else%in case ter->ini happens.Then, One cell must be referred as the descendant cell. In this case, pack things and sent to 'B->C'
            isFound(iA,1)=false;
            lostCnt(iA,1)=lostCnt(iA,1)+1;
            iProb1=find(Prob1(:,1)>0.85,1);
            selCellIDB=Prob1(iProb1,2);
            if ~isempty(iProb1) && ~isempty(B{selCellIDB,1})
                selRow=find(Prob1(:,1)>0.85, 1);
                PassToNxtStp{iA,1}=[B(Prob1(selRow,2),1);C(Prob1(selRow,2),1)];
                isSti(iA,1)=true;
            end
%<----------------add here a general block for 'simple-lost' case indicator
            


        end
    end
end

%% B -> C
% this works if isSti is true
nFound=find(~isFound);
if ~isempty(nFound) 
    for inF=1:length(nFound)%for each 'not-found' (i.e., really lost) spot
        if lostCnt(nFound(inF))>0 && isSti(nFound(inF))
            % look for spots on the 3rd frame. If the category of 2nd and 3rd frame
            % both lacks Cat2, then this is a ter->ini.
            for iTPNS=1:length(PassToNxtStp)
                if isempty(PassToNxtStp{iTPNS,1});continue;end
        
                matB=PassToNxtStp{iTPNS,1}{1,1};
                CatPsbB=zeros(size(matB,1),1);
                for iB=1:size(matB,1)
                    CatPsbBTmp=CatPos4(matB(iB,3));
                    if CatPsbBTmp(1,1)>0.8
                        CatPsbB(iB,1)=CatPsbBTmp(1,2);
                    end
                end
                
                if size(PassToNxtStp{iTPNS,1},1)==3
                    uit{iTPNS,1}=[2,1,-1,-1];
                else
                    matC=PassToNxtStp{iTPNS,1}{2,1};
                    CatPsbC=zeros(size(matC,1),1);
                    for iC=1:size(matC,1)
                        CatPsbCTmp=CatPos4(matC(iC,3));
                        if CatPsbCTmp(1,1)>=0.7
                            CatPsbC(iC,1)=CatPsbCTmp(1,2);
                        end
                    end
                    if isempty(find([CatPsbB;CatPsbC]==2, 1))
                        %this indicates a ter->ini case.
                        uit{iTPNS,1}=[2,1,-1,-1];%category, possibility, cellID (marked as -1), row of spot in Cell B (spotID in cell)(marked as -1)
                    end
                end
            end
        elseif lostCnt(nFound(inF))>0 && ~isSti(nFound(inF))
% Case 'Really lost': should check 'A -> C'<----------------working on here
            % look for spots on the 3rd frame. If the category of this spot
            % is the same as any one or two spot(s) on the 3rd frame, calculate the
            % possibility between these spots.
            % nFound(inF) : the index of not-found spot in Cell A.
            
            [uit(nFound(inF),1),is3rdCell(nFound(inF),1)]=MTL_ACjudge(A{1,1}(nFound(inF),:),curCat,C);%the info of this not-found spot in Cell A
        end
    end
end

%<----------------add here a general block for 'simple-lost' case resolver


%% write 'result' as text output:
sptIDtmp=[];%Case: 2->1
sptIDtmp2=cell(length(uit),1);% Case: 2=x=2 
isSptTemp2=false;
for iA=1:length(uit)
    curUit=uit{iA,1};

    if isempty(curUit)
        res{iA,1}='lost';
    end
    if size(curUit,1)>1
        res{iA,1}='splitted';
        if size(curUit,1)==2
            sptIDtmp2{iA,1}=curUit(:,end-1:end);
            isSptTemp2=true;
        end
    end
    if size(curUit,1)==1
        sptIDtmp=cat(1,sptIDtmp,[curUit(3:4),iA]);%<-original
       % sptIDtmp=cat(1,sptIDtmp,[curUit{1,1}(3:4),iA]);
        res{iA,1}='continue';
        if isequal(curUit(3:4),[-1,-1])
            res{iA,1}='TerIni';
        end
    end
end
if size(sptIDtmp,1)>1 % in case 2->1
    Dup=findDup(sptIDtmp);
    for iDup=1:length(Dup)
        sel=Dup{iDup,1}(3:end);
        for iSel=1:length(sel)
            res{sptIDtmp(sel(iSel),3),1}=strcat('kis',num2str(iSel));
        end
    end
end
if isSptTemp2
    Dup=findDupCell(sptIDtmp2);
    if ~isempty(Dup) && ~isempty(sptIDtmp2{Dup{1,1}(1)})
        for iDup=1:length(Dup)
            Atmp=Dup{iDup};%means they splited to the same set of spots in Cell B
            %now check if these spots in Cell A is close enough
            Xloc=A{1,1}(Atmp,3);
            if max(Xloc)-min(Xloc)<0.25
                %assign descendants by location
                %locA=[];
                locB=uit{Atmp(1),1}(:,end-1:end);
                LrtB=zeros(length(Atmp),3);
                for iAtmp=1:length(Atmp)
                    res{Atmp(iAtmp)}='continue';
                    if iAtmp>size(locB,1)
                        %added on 22.04.06 to avoid bug
                        if iAtmp>=size(LrtB,1)+1
                            LrtB(end,:)=[];
                        else
                            LrtB(iAtmp,:)=[];
                        end
                    else
                        if isempty(B{1,1})
                            LrtB(iAtmp,:)=cat(2,C{locB(iAtmp,1),1}(locB(iAtmp,2),3),locB(iAtmp,:));%lrt of selected spot in C.
                        else
                            LrtB(iAtmp,:)=cat(2,B{locB(iAtmp,1),1}(locB(iAtmp,2),3),locB(iAtmp,:));%lrt of selected spot in B.
                        end
                    end
%                     uitTmp1=uit{Atmp(iAtmp)};
%                     [~,idUT1]=min(uitTmp1(:,2));
%                     uitTmp1(idUT1,:)=[];
%                     uit{Atmp(iAtmp)}=uitTmp1;
                end
                uit=reShuffle(uit,Atmp,Xloc,LrtB);
            end
        end
    end
end

end
%% nested function
function [judge]=is2LikelyInB(in)
uit=false(size(in,2),1);
for i=1:size(in,2)
    if find(in(:,i)>0.8,1)
        uit(i)=true;
    end
end
if sum(uit)>=2
    judge=true;
else
    judge=false;
end
end
function [uit]=reShuffle(uit,Atmp,LrtA,LrtB)
if length(LrtA)==2
    LrtA=[LrtA,[1:1:size(LrtA,1)]'];
    LrtB=[LrtB,[1:1:size(LrtB,1)]'];
    Map=zeros(size(LrtA,1),1);
    [~,idA]=sort(LrtA(:,1));
    LrtAx=LrtA(idA,:);
    [~,idB]=sort(LrtB(:,1));
    LrtBx=LrtB(idB,:);
    Map(LrtAx(:,end),1)=LrtBx(:,end);
    for iAtp=1:length(Atmp)
        uit{Atmp(iAtp),1}=uit{Atmp(iAtp),1}(Map(iAtp),:);
    end
elseif length(LrtA)==3 && size(LrtB,1)==2
    idA=nchoosek(1:1:3,2);
    collec=zeros(size(idA,1),1);
    for ix=1:size(idA,1)
        collec(ix,1)=abs(LrtA(idA(ix,1))-LrtA(idA(ix,2)));
    end
    [~,selR]=min(collec);
    AtmpOri=Atmp;
    combine=idA(selR,:);%combined index in Atmp.
    Atmp(combine(2))=[];
    LrtA(combine(2))=[];

%repeat case length(LrtA)==2
    LrtA=[LrtA,[1:1:size(LrtA,1)]'];
    LrtB=[LrtB,[1:1:size(LrtB,1)]'];
    Map=zeros(size(LrtA,1),1);
    [~,idA]=sort(LrtA(:,1));
    LrtAx=LrtA(idA,:);
    [~,idB]=sort(LrtB(:,1));
    LrtBx=LrtB(idB,:);
    Map(LrtAx(:,end),1)=LrtBx(:,end);
    for iAtp=1:length(Atmp)
        uit{Atmp(iAtp),1}=uit{Atmp(iAtp),1}(Map(iAtp),:);
    end
    uit{AtmpOri(combine(2)),1}=uit{AtmpOri(combine(1)),1};
end
end

function [psb1u,psb2u]=modPsb12(cat1,psb1,psb2)
psb1u=psb1;
psb2u=psb2;
if cat1(2,1)>0.35%means almost equal possibility in two categories.
    
    for iPSB=1:size(psb1,2)
        if psb1(2,iPSB)==0; continue;end
        if psb1(2,iPSB)<psb1(1,iPSB)
            psb1u(:,iPSB)=psb2(:,iPSB).*[1;0;0];
            psb2u(:,iPSB)=psb2(:,iPSB).*[1;0;0];
        else
            psb1u(:,iPSB)=psb2(:,iPSB).*[0;1;0];
            psb2u(:,iPSB)=psb2(:,iPSB).*[0;1;0];
        end
        
    end
%     psb2u=psb1u;
%     [max_psb2,idx]=max(psb2);
%     if size(psb1,2)==1 && max_psb2>0.7
%         mult=zeros(3,1);
%         mult(idx)=1;
%         psb1u(:,1)=psb2(:,1).*mult;
%     end
end
end
function [res]=findDupCell(sptIDtmp2)
Mat1=zeros(length(sptIDtmp2),2);
Mat2=zeros(length(sptIDtmp2),2);
for iDup=1:length(sptIDtmp2)
    if isempty(sptIDtmp2{iDup,1})
        Mat1(iDup,:)=[rand,rand+1];
        Mat2(iDup,:)=[rand,rand-1];
        continue;
    end
    %fill up Mat1,2
    Mat1(iDup,:)=sptIDtmp2{iDup,1}(1,:);
    Mat2(iDup,:)=sptIDtmp2{iDup,1}(2,:);
end
%for Mat1
[~,~,ib1]= unique(Mat1, 'rows', 'stable');
ret1={};
for its1=1:max(ib1)
    sel=find(ib1==its1);
    if length(sel)>1
        ret1=[ret1;{sel'}];
    end
end

%for Mat2
[~,~,ib2]= unique(Mat2, 'rows', 'stable');
ret2={};
for its2=1:max(ib2)
    sel2=find(ib2==its2);
    if length(sel2)>1
        ret2=[ret2;{sel2'}];
    end
end
res={};
for k1=1:length(ret1)
    for k2=1:length(ret2)
        cs1=ret1{k1,1};
        cs2=ret2{k2,1};
        A=intersect(cs1,cs2);
        if length(A)>1
            res=[res;{A}];
        end
    end
end
end

function [ret]=findDup(sptid)
A=sptid(:,1:2);
[~,~,ib]= unique(A, 'rows', 'stable');
ret={};
for its=1:max(ib)
    sel=find(ib==its);
    if length(sel)>1
        ret=[ret;{[A(sel(1),:),sel']}];
    end
end
end

function [pred,funcs]=PredictCat4(curCat,ctrl)

    %ctrl: [1,0] not divided; [2,1] divided to mother end;
    %[2,2] divided to daughter end;
    pred=zeros(3,2);%intended to be (3,2):[rng1,rng2];
    funcs=cell(3,1);
    if ctrl(1,1)==1
        for iCat=1:3%intended to be 3
            if curCat(iCat,1)==0; continue; end
            switch curCat(iCat,2)
                case 3
                    pred(iCat,1:2)=[0.5, 1];
                    funcs{iCat,1}=@(x) sin(2*pi*x-pi);
                case 2
                    pred(iCat,1:2)=[0.25, 0.75];
                    funcs{iCat,1}=@(x) sin(2*pi*x-0.5*pi);
                case 1
                    pred(iCat,1:2)=[0, 0.5];
                    funcs{iCat,1}=@(x) sin(2*pi*x);
            end
        end
    elseif ctrl(1,1)==2 && ctrl(1,2)==2
        for iCat=1:3%intended to be 3
            if curCat(iCat,1)==0; continue; end
            switch curCat(iCat,2)
                case 3
                    pred(iCat,1:2)=[0, 0.15];
                case 2
                    pred(iCat,1:2)=[0, 0.5];
                case 1
                    pred(iCat,1:2)=[0, 1];
            end
        end
    elseif ctrl(1,1)==2 && ctrl(1,2)==1
        for iCat=1:3%intended to be 3
            if curCat(iCat,1)==0; continue; end
            switch curCat(iCat,2)
                case 3
                    pred(iCat,1:2)=[0, 1];
                case 2
                    pred(iCat,1:2)=[0.5, 1];
                case 1
                    pred(iCat,1:2)=[0.85,1];
            end
        end
    end
end
    