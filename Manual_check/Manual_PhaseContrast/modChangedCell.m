function [Msk,modFrames]=modChangedCell(prevMsk,Msk,delUnderRow)
%Version: v1.1 @ 2023.04.20 optimized for mesh fitting.
%Version: v1.02 @ 2023.02.10 fixed error when 3 cells are connected.

prevMsk=logical(prevMsk);
dif=xor(prevMsk,Msk);
if delUnderRow~=0
    dif(delUnderRow:end,:)=false;
end
L=regionprops(dif,"Centroid");
procFrm=[];         %for all modified frames, register processed ones.
modFrames=[];       
for iL=1:length(L)  %for each frame that contains dif region:
    frame=floor(L(iL).Centroid(1)/32)+1;                %frame number
    prevTile=prevMsk(:,(frame-1)*32+1:frame*32);        %prev. mask of this frame
    newTile=Msk(:,(frame-1)*32+1:frame*32);             %curr. Mask of this frame
    %core functions
    if ~ismember(frame,procFrm)
        newTile=refineModTile(prevTile,newTile,frame);%do the modification to the changed cellID.
        Msk(:,(frame-1)*32+1:frame*32)=newTile;
        procFrm=cat(1,procFrm,frame);%register which frame is processed. Don't redo it.
    end
    modFrames=cat(1,modFrames,frame);
end
modFrames=unique(modFrames);

function [prevX]=refineModTile(prev,new,iFrm)
prevX=prev;
dif1=xor(prev,new);
Lprev=regionprops(prev,'PixelIdxList','Centroid');
Lprev=sortProps(Lprev);
Lnew=regionprops(new,'PixelIdxList','Centroid');
Lnew=sortProps(Lnew);
Ldif=regionprops(dif1,'PixelIdxList');
modX=[];%format of modX [cellID, modification method idx].

for i=1:length(Ldif)%for each modified region.
%1. judge if the modification is adding or erasing.
    %for this modified region, which regionprop in Lprev and Lnew may contain it?
    isInPrev=false;
    isInNew=false;
    ratInPrev=0;
    ratInNew=0;
    selILP=0;
    selILN=0;
    KP=cellfun(@(x) sum(ismember(Ldif(i).PixelIdxList,x))/length(Ldif(i).PixelIdxList),Lprev);
    if find(KP,1)
        selILP=find(KP,1);%by definition, find(KP) or find(KN) can only be a single value.
        isInPrev=true;
        ratInPrev=sum(ismember(Ldif(i).PixelIdxList,Lprev{selILP,1}))/length(Lprev{selILP,1});
    end
    KN=cellfun(@(x) sum(ismember(Ldif(i).PixelIdxList,x))/length(Ldif(i).PixelIdxList),Lnew);
    if find(KN,1)
        selILN=find(KN,1);
        isInNew=true;
        ratInNew=sum(ismember(Ldif(i).PixelIdxList,Lnew{selILN,1}))/length(Lnew{selILN,1});
    end
    clear('iLp','iLn','isInPrevX','isInNewX');

% 2. see which case does this modification belong.     
    if ~isInPrev && isInNew && ratInNew<0.98 % Case 'join' or 'partially expand'
        %Overlapped cells In Prev with the given cell in new
        oip=cellfun(@(x) sum(ismember(Lnew{selILN},x))/length(x),Lprev);
        selPrevIdx=find(oip);
    % 2.1 Case 'join'
        if length(selPrevIdx)==2
        %look for original two masks in PrevMsk:
            prevX(Lprev{selPrevIdx(1),1})=false;%clear candidate regions in prevmask.
            prevX(Lprev{selPrevIdx(2),1})=false;%clear candidate regions in prevmask.
            modX=cat(1,modX,[selILN,1]);
            disp(['On Frame ',num2str(iFrm),', Cell ',num2str(selPrevIdx(1)),' and ',num2str(selPrevIdx(2)),' is joined into Cell ',num2str(selILN),'.']);
    % 2.2 Case 'partially expand'
        elseif length(selPrevIdx)==1
            prevX(Lprev{selPrevIdx,1})=false;
            disp(['On Frame ',num2str(iFrm),', Cell ',num2str(selILN),' is expanded.']);
            modX=cat(1,modX,[selILN,1]);
        else
            for iDel=1:length(selPrevIdx)
                prevX(Lprev{selPrevIdx(iDel),1})=false;
            end
            modX=cat(1,modX,[selILN,1]);
            %error(['The Modification ', num2str(i),' on Frame ',num2str(iFrm),'covers on >3 cells.']);
        end

    elseif isInPrev && ~isInNew && ratInPrev<0.98 % Case 'split' or 'partially erase'
        %Overlapped cells In Prev with the given cell in new
        oin=cellfun(@(x) sum(ismember(x,Lprev{selILP,1}))/length(x),Lnew);
        selNewIdx=find(oin);
    % 2.3 Case 'split'
        if length(selNewIdx)==2
            prevX(Lprev{selILP,1})=false;%clear candidate regions in prevmask.
            modX=cat(1,modX,[selNewIdx(1),2;selNewIdx(2),2]);
            disp(['On Frame ',num2str(iFrm),', Cell ',num2str(selILP),' is splitted into Cell ',num2str(selNewIdx(1)),' and ',num2str(selNewIdx(2)),'.']);
    % 2.4 Case 'partially erase'
        elseif length(selNewIdx)==1
            prevX(Lprev{selILP,1})=false;
            modX=cat(1,modX,[selNewIdx,2]);
            disp(['On Frame ',num2str(iFrm),', Cell ',num2str(selNewIdx),' is partially erased.']);
        else
            %error(['The Modification ', num2str(i),' on Frame ',num2str(iFrm),'covers on >3 cells.']);
        end
    % 2.5 Case 'delete'
    elseif isInPrev && ~isInNew && ratInPrev>=0.98 
        prevX(Lprev{selILP,1})=false;
        %disp(['On Frame ',num2str(iFrm),', Cell ',num2str(selILP),' is deleted.']);
    % 2.6 Case 'add'
    elseif ~isInPrev && isInNew && ratInNew>=0.98 
        modX=cat(1,modX,[selILN,2]);
        disp(['On Frame ',num2str(iFrm),', Cell ',num2str(selILN),' is added.']);
    else % Case 'unknown'
        disp('the modification belongs to NO category!')
    end
end

% 3. perform modification.
if ~isempty(modX)
    %remove duplicated register.
    modX=remDupReg(modX);
    for iRow=1:size(modX,1)
        newX=false(256,32);
        newX(Lnew{modX(iRow,1)})=true;
        if modX(iRow,2)==1%for Case join and partially expand
            newX=imopen(newX,strel("disk",2));
        elseif modX(iRow,2)==2%for Case split and patially shrink
            newX=imopen(newX,strel("disk",4));
        end
        prevX(newX)=true;
    end
    prevX=bwmorph(prevX,"hbreak");%added on 2024.03.04
end
end
function [uit]=remDupReg(in)
% find out duplicated first column.
[in1,~,ib1]= unique(in(:,1), 'rows', 'stable');
in2=zeros(size(in1));
ret1={};
for its1=1:max(ib1)
    sel=find(ib1==its1);
    if length(sel)>1
        ret1=[ret1;{sel'}];
    end
end
% for each duplicated cellID, designate a type indicator.
if ~isempty(ret1)
    for ix=1:length(ret1)
        modX=in(ret1{ix,1},2);
        if ~isempty(find(modX==2, 1))
            num=in(ret1{ix,1}(1),1);
            in2(in1==num,1)=2;
        else
            num=in(ret1{ix,1}(1),1);
            in2(in1==num,1)=1;
        end

    end
end
% for the rest of cellIDs, use its own type indicator.
for i=1:length(in2)
    if in2(i)==0
        num=in1(i);
        in2(i)=in(in(:,1)==num,2);
    end
end
uit=[in1,in2];
end
function [uit]=sortProps(L)
Yraw=cellfun(@(x) x(2),{L.Centroid});
[~,idx]=sort(Yraw,'descend');
K={L.PixelIdxList}';
%uit=struct('PixelIdxList',K(idx));
uit=K(idx);
end
end