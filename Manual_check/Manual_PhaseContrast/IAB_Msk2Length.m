function [outL,outMsk,outIdxList]=IAB_Msk2Length(in,minAL)
%extract length from Mask
%fill contours that is smaller than 180
nTile=size(in,2)/32;
outMsk=in;
outL=cell(nTile,1);
outIdxList=cell(nTile,1);
for i=1:nTile
    curTile=in(:,(i-1)*32+1:i*32);
    bwL=bwlabel(curTile,4);
    L=regionprops(bwL,'Area','MajorAxisLength','PixelIdxList','PixelList');
    RegL=[];
    RegMinIdx=[];
    RegIdx={};
    for k=1:length(L)
        MinReq=min(L(k).PixelList(:,2),[],'all');
        if L(k).Area<minAL(1) || MinReq<=minAL(3)%<---arbitary parameter detected....
            curTile(L(k).PixelIdxList)=false;
            continue;
        else
            RegL=cat(1,RegL,L(k).MajorAxisLength);
            RegMinIdx=cat(1,RegMinIdx,min(L(k).PixelList(:,2),[],'all'));
            RegIdx=[RegIdx;{L(k).PixelIdxList}];
        end
    end
    [~,idx]=sort(RegMinIdx,'descend');
    %RegL=RegL(idx,:);
    outMsk(:,(i-1)*32+1:i*32)=curTile;
    outL{i,1}=RegL(idx,:);%This is for lineage linking
    outIdxList{i,1}=RegIdx(idx,:);
end
end