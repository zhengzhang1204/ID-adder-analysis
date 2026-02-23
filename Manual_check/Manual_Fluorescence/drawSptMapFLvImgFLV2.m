%Version 2.0 @23.02.20  optimized for time cost.
function [out,canvLineRawImg,CanvIdxSpot,CanvIdxLine,CanvIdxLineStEnd]=drawSptMapFLvImgFLV2(locList,Map,RJ,PHMask)
% Preparation and initation:
Contour=uint8(imdilate(PHMask,strel('diamond',2))-PHMask)*255*0.4;
canvas=uint8(zeros(size(PHMask)));
sz=[size(canvas,1),size(canvas,2)];

% remove bad entries for spots.
del=find(isnan(locList(:,1)));
sptNum=(1:1:size(locList,1))';
if ~isempty(del)
    locList(del,:)=[];
    Map(del,:)=[];
    sptNum(del)=[];
end

% draw marker and spot row ID:
canvas = insertMarker(canvas,locList,'x-mark','color','green');
canvas = insertText(canvas,locList, sptNum,'BoxOpacity',0,'FontSize',9,'TextColor','green');

% create yellow lines
StartEndLineInd=cell2mat(cellfun(@(x,y) createStartEndLineMat(x,y),num2cell(Map,2),num2cell([1:1:length(Map)]'),'UniformOutput',false));
StartEndLineMat=cell2mat(arrayfun(@(x,y) [locList(x,:),locList(y,:)],StartEndLineInd(:,1),StartEndLineInd(:,2),'UniformOutput',false));
canvas=insertShape(canvas,'Line',StartEndLineMat,'LineWidth',1,'Color','red');

% cover the yellow lines near spots so that yellow lines do not overlap with green x.
sptCent=num2cell(locList,2);
CanvIdxSpot=cellfun(@(x) coverSpotRec(x,sz,4),sptCent,'UniformOutput',false);
tape=canvas(:,:,1);
tape(cell2mat(CanvIdxSpot))=0;
canvas(:,:,1)=tape;

% create thick lines and get their pixel index
CanvIdxLine=getThickLinePixelIds(StartEndLineMat,sz);

% draw frame numbers:
FrmNum=size(RJ,2)/32;
TxtPos=[((0:1:FrmNum-1)*32+6)',235*ones(FrmNum,1)];
txt=num2cell((1:1:FrmNum)');
txt=cellfun(@(x) num2str(x),txt,'UniformOutput',false);
canvas = insertText(canvas,TxtPos,txt,'BoxOpacity',0,'FontSize',16,'TextColor','white');
out=cat(3,RJ+Contour,RJ+2*canvas(:,:,2),RJ);
canvLineRawImg=canvas(:,:,1);%[uint8] black canvas with white lines only

% prepare output:
CanvIdxLineStEnd=[StartEndLineInd,StartEndLineMat];
end

function [ind]=coverSpotRec(coor,sz,value)
[X,Y]=meshgrid(coor(1)-value:coor(1)+value,coor(2)-value:coor(2)+value);

X(X<1 | X> sz(2))=0;
Y(X<1 | X> sz(2))=0;
Y(Y<1 | Y> sz(1))=0;
X(Y<1 | Y> sz(1))=0;
Xnew=reshape(X,[],1);
Ynew=reshape(Y,[],1);

delX=find(Xnew==0);
delY=find(Ynew==0);
del=unique([delX;delY]);
Xnew(del)=[];
Ynew(del)=[];

ind=sub2ind(sz,Ynew,Xnew);
end

function [uit]=getThickLinePixelIds(mat,sz)
uit=cell(size(mat,1),1);
for i=1:size(mat,1)%for each line
    xq=mat(i,1)+4:1:mat(i,3)-4;
    if isempty(xq);continue;end
    yq=round(interp1([mat(i,1),mat(i,3)],[mat(i,2),mat(i,4)],xq));
    yqExt=repmat(yq,7,1)+(3:-1:-3)'*ones(1,length(yq));
    xqExt=repmat(xq,7,1);
    yqExt(yqExt>256|yqExt<1)=nan;
    xqExt(yqExt>256|yqExt<1)=nan;
    yqExt=reshape(yqExt,[],1);
    xqExt=reshape(xqExt,[],1);
    del=isnan(yqExt);
    yqExt(del)=[];
    xqExt(del)=[];
    uit{i}=sub2ind(sz,yqExt,xqExt);
end
end

function [uit]=createStartEndLineMat(map,num)
    map(isnan(map))=[];
    if isempty(map)
        uit=[];
        return;
    end
    if map(1)==-1
        uit=[];
        return;
    end
    uit=[num*ones(length(map),1),map'];
end






