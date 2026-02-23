% load('T:\test.mat');
% [uList,ICv2,Umsk]=IAB_testFittingErr2(Msk,PHparts);
function [uList,ICv2,Umsk]=IAB_testFittingErrv2(Msk,PHparts)
uList=[];
validCID=cell2mat(PHparts);%frame, cid

% preparing labelled image L for each frame:
nFrm=size(Msk,2)/32;
L=uint8(zeros(256,32,nFrm));
maxCIDonEachFrm=cell2mat(cellfun(@(x) max(validCID(validCID(:,1)==x,2)), num2cell((1:1:nFrm)'),'UniformOutput',false));
for iFrm=1:length(maxCIDonEachFrm)%for each Frame
    curMSK=Msk(:,(iFrm-1)*32+1:iFrm*32);
    curL=uint8(zeros(256,32));
    CC = bwconncomp(curMSK,4);
    K=regionprops(CC,'PixelIdxList','Centroid');
    Yraw=cellfun(@(x) x(2),{K.Centroid});
    [~,idx]=sort(Yraw,'descend');
    Z={K.PixelIdxList}';
    Z=Z(idx);
    for ix=1:maxCIDonEachFrm(iFrm)
        curL(Z{ix})=ix;
    end
    L(:,:,iFrm)=curL;
end

%test fitting problem by each frame
cellMesh=cell(size(validCID,1),1);
WL=cell(size(validCID,1),1);
clc;
disp('Creating meshes for all valid cells...')
parfor iSel=1:size(validCID,1)%PARFOR each valid CID entry
    inse=L(:,:,validCID(iSel,1))==validCID(iSel,2);
    [cellMesh{iSel},~,WL{iSel},~]=msk2mesh(inse);%<--------------core
end
disp('Creating meshes is done.')

emptyEntry=find(cellfun(@(x) isempty(x), cellMesh));%an empty cell mesh means fitting has problem.
if ~isempty(emptyEntry)
    uList=validCID(emptyEntry,:);
    Umsk=false(256,nFrm*32);
    for iUlist=1:size(uList,1)
        tempUmsk=Umsk(:,(uList(iUlist,1)-1)*32+1:uList(iUlist,1)*32);
        xxx=L(:,:,uList(iUlist,1))==uList(iUlist,2);
        tempUmsk(xxx)=true;
        Umsk(:,(uList(iUlist,1)-1)*32+1:uList(iUlist,1)*32)=tempUmsk;
        WL{emptyEntry(iUlist)}=[0,0];
    end
else
    Umsk=false(256,nFrm*32);
    uList=[];
end

%prepare Infcellv2 for the first time:
ICv2=cell(nFrm,1);
for iFrm=1:length(maxCIDonEachFrm)%for each frame
    curICv2=cell(1,3);
    % storing labelled image.
    curICv2{1}=L(:,:,iFrm);
    newMat=zeros(maxCIDonEachFrm(iFrm),8);
    meshTemp=cell(maxCIDonEachFrm(iFrm),1);
    meshLength=zeros(maxCIDonEachFrm(iFrm),1);
    for iCID=1:maxCIDonEachFrm(iFrm)%number of cells on this frame
        entIdx=find(validCID(:,1)==iFrm & validCID(:,2)==iCID);
        indx=L(:,:,iFrm)==iCID;
        [Area,Volume]=getAreaVolume(cellMesh{entIdx});

        % storing newMat [binarySize, backboneLength, width, OuftiV, OuftiS, ModelS, ModelV, MajorAxisLength]
        newMat(iCID,1)=sum(indx,'all');% binary size
        newMat(iCID,2)=WL{entIdx}(2);% backbone length (pole-to-pole)
        newMat(iCID,3)=WL{entIdx}(1);% width
        newMat(iCID,4)=Volume;% Oufti Volume from mesh
        newMat(iCID,5)=Area;  % Oufti Area from mesh
        newMat(iCID,6)=WL{entIdx}(1)*(WL{entIdx}(2)-WL{entIdx}(1))+pi*(WL{entIdx}(1)^2)/4; %Model area
        newMat(iCID,7)=(WL{entIdx}(2)-WL{entIdx}(1))*pi*(WL{entIdx}(1)^2)/4 + 4/3*pi*(WL{entIdx}(1)^3)/8;% Model volume
        MAL=regionprops(indx,'MajorAxisLength');
        newMat(iCID,8)=MAL.MajorAxisLength;
        meshTemp{iCID}=cellMesh{entIdx};
        meshLength(iCID)=size(cellMesh{entIdx},1);
    end
    curICv2{2}=newMat;

    % storing newMesh
    newMesh=single(zeros(max(meshLength),4,maxCIDonEachFrm(iFrm)));
    for iCID=1:maxCIDonEachFrm(iFrm)
        newMesh(:,:,iCID)=[meshTemp{iCID};zeros(max(meshLength)-meshLength(iCID),4)];
    end
    curICv2{3}=newMesh;
    ICv2{iFrm}=curICv2;
end
end

%% neseted functions
function [Area,Volume]=getAreaVolume(mesh)
    if isempty(mesh)
        Area=0;Volume=0;
        return;
    end

   steplength = edist(mesh(2:end,1)+mesh(2:end,3),mesh(2:end,2)+mesh(2:end,4),...
                              mesh(1:end-1,1)+mesh(1:end-1,3),mesh(1:end-1,2)+mesh(1:end-1,4))/2;
   % area
   steparea = zeros(size(mesh,1)-1,1);
   for i=1:size(mesh,1)-1
       steparea(i,1) = polyarea([mesh(i:i+1,1);mesh(i+1:-1:i,3)],[mesh(i:i+1,2);mesh(i+1:-1:i,4)]); 
   end
   Area = sum(steparea);

   % volume
   d = edist(mesh(:,1),mesh(:,2),mesh(:,3),mesh(:,4));
   stepvolume = (d(1:end-1).*d(2:end) + (d(1:end-1)-d(2:end)).^2/3).*steplength*pi/4;
   Volume = sum(stepvolume);
end

function d=edist(x1,y1,x2,y2)
    % complementary for "getextradata", computes the length between 2 points
    d=sqrt((x2-x1).^2+(y2-y1).^2);
end