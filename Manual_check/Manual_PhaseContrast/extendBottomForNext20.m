function [Msk]=extendBottomForNext20(prevMsk,Msk,delUnderRow)
prevMsk=logical(prevMsk);
dif=xor(prevMsk,Msk);
if delUnderRow~=0
    dif(delUnderRow:end,:)=false;
end
L=regionprops(dif,"Centroid");

%1. identify the last modified frame:
lastModFrm=zeros(length(L),1);
for iL=1:length(L)  %for each frame that contains dif region:
    lastModFrm(iL)=floor(L(iL).Centroid(1)/32)+1;                %frame number
end
lastModFrm=max(lastModFrm);

%2. spread the bottom location to the next 20 frames
refTile=Msk(:,(lastModFrm-1)*32+1:lastModFrm*32);             %curr. Mask of this frame
newBotLoc=find(sum(refTile,2)~=0,1,'last');
nTiles=size(Msk,2)/32;
frmToBeProcessed=[lastModFrm+1:1:min([lastModFrm+20,nTiles])];
for iL=1:length(frmToBeProcessed)  %for each frame that contains dif region:
    oldTile=Msk(:,(frmToBeProcessed(iL)-1)*32+1:frmToBeProcessed(iL)*32);             %raw Mask of this frame
    L=regionprops(oldTile,"Centroid","PixelIdxList");
    Cent=zeros(length(L),2);
    for k=1:length(L)
        Cent(k,:)=L(k).Centroid;
    end
    [~,idx]=max(Cent(:,2));
    oldTile(L(idx).PixelIdxList)=0;
    newTile=false(size(oldTile));
    newTile(L(idx).PixelIdxList)=1;
    K=bwboundaries(newTile);K=K{1};
    k=(newBotLoc-min(K(:,1)))/(max(K(:,1))-min(K(:,1)));
    b=newBotLoc-k*max(K(:,1));
    K(:,1)=k*K(:,1)+b;
    V=poly2mask(K(:,2),K(:,1),256,32);
    V=imdilate(V, strel('line', 3, 0)) | oldTile;
    Msk(:,(frmToBeProcessed(iL)-1)*32+1:frmToBeProcessed(iL)*32)=V;
end
end