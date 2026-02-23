

%This is to create mesh fitting for selected frames. The input is entire set of tile with selected frame number.
% Version 1.0   @23.04.25   The purpose is to save time during manual correction of segmentations. 
%                           Previously, all cells will do mesh fitting once pressing 'T'. Now, only the modified frames need to do mesh fitting.

function [uList,prevICv2,Umsk]=IAB_testFittingErrSelFrms(Msk,PHparts,selFrms,prevICv2,prevUlist2,prevUmsk2)
if isempty(selFrms)
    uList=[];
    Umsk=[];
    return;
end

uList=[];
validCID=cell2mat(PHparts);%frame, cid

%1. Valid CIDs on selected Frames and the max ID number on those frames.
selVCID=[];
maxCIDonEachSelFrm=nan(length(selFrms),1);
for i=1:length(selFrms)
    curVCID=validCID(validCID(:,1)==selFrms(i),:);
    maxCIDonEachSelFrm(i)=max(curVCID(:,2));
    selVCID=cat(1,selVCID,curVCID);
end

% 2. preparing labelled image L for each selected frame:
L=uint8(zeros(256,32,length(selFrms)));
for iFrm=1:length(selFrms)%for each Frame
    curMSK=Msk(:,(selFrms(iFrm)-1)*32+1:selFrms(iFrm)*32);
    curL=uint8(zeros(256,32));
    CC = bwconncomp(curMSK,4);
    K=regionprops(CC,'PixelIdxList','Centroid');
    Yraw=cellfun(@(x) x(2),{K.Centroid});
    [~,idx]=sort(Yraw,'descend');
    Z={K.PixelIdxList}';
    Z=Z(idx);
    for ix=1:maxCIDonEachSelFrm(iFrm)
        curL(Z{ix})=ix;
    end
    L(:,:,iFrm)=curL;
end

%3. par for each L.
cellMesh=cell(size(selVCID,1),1);
WL=cell(size(selVCID,1),1);
sVCID2=selVCID(:,2);sVCID1=selVCID(:,1);
for iSel=1:size(selVCID,1)%PARFOR each valid CID entry
    inse=L(:,:,(selFrms==sVCID1(iSel)))==sVCID2(iSel);
    [cellMesh{iSel},~,WL{iSel},~]=msk2mesh(inse);%<--------------core
end

%4. on the selected frames, are meshes all created correctly?
emptyEntry=find(cellfun(@(x) isempty(x), cellMesh));
%4.1 remove entries on selected frames from the previous ulist2:


if ~isempty(prevUlist2)
    del=ismember(prevUlist2(:,1),selFrms);
    prevUlist2(del,:)=[];
end

if ~isempty(emptyEntry)% fitting errors occured for one or more cells.
%update Ulist2:
    newUlist2=selVCID(emptyEntry,:);
    newUlist2=[newUlist2;prevUlist2];
    newUlist2=sortrows(newUlist2,[1,2]);

% find out if there is any frame in selFrms that all cells are correct.
    correctSelFrms=selFrms(~ismember(selFrms,unique(sVCID1(emptyEntry))));%absolute frame number.

%patch previous Umsk2
    for iSelFrm=1:length(selFrms)
        curTile=false(256,32);
        atFrm=find(newUlist2(:,1)==selFrms(iSelFrm));
        if ~isempty(atFrm)
            for k=1:length(atFrm)
                curTile(L(:,:,(selFrms==newUlist2(atFrm(k),1)))==newUlist2(atFrm(k),2))=true;
            end
        end
        prevUmsk2(:,(selFrms(iSelFrm)-1)*32+1:selFrms(iSelFrm)*32)=curTile;
    end
else% current tested fittings are all good.
    newUlist2=prevUlist2;
    for iSelFrm=1:length(selFrms)

        prevUmsk2(:,(selFrms(iSelFrm)-1)*32+1:selFrms(iSelFrm)*32)=false(256,32);
    end
    correctSelFrms=selFrms;
end


% update ICv2 for those correct frames.
if ~isempty(correctSelFrms)
    for i=1:length(correctSelFrms)
        iselFrm=find(selFrms==correctSelFrms(i),1);
        curICv2=cell(1,3);

        % storing labelled image.
        curICv2{1}=L(:,:,iSelFrm);

        % storing results.
        newMat=zeros(maxCIDonEachSelFrm(iselFrm),8);
        meshTemp=cell(maxCIDonEachSelFrm(iselFrm),1);        
        meshLength=zeros(maxCIDonEachSelFrm(iselFrm),1);
        for iCID=1:maxCIDonEachSelFrm(iselFrm)%number of cells on this frame
            entIdx=find(selVCID(:,1)==selFrms(iselFrm) & selVCID(:,2)==iCID);
            indx=L(:,:,iselFrm)==iCID;
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
        newMesh=single(zeros(max(meshLength),4,maxCIDonEachSelFrm(iselFrm)));
        for iCID=1:maxCIDonEachSelFrm(iselFrm)
            newMesh(:,:,iCID)=[meshTemp{iCID};zeros(max(meshLength)-meshLength(iCID),4)];
        end
        curICv2{3}=newMesh;

        prevICv2{correctSelFrms(i)}=curICv2;
    end
end




%4.2 see if newUlist2 is empty:
if ~isempty(newUlist2)%newUlist2 is not empty. Still need adjust fittings manually.
    Umsk=prevUmsk2;
    uList=newUlist2;
else
%     for iselFrm=1:length(selFrms)%for each selected frame
%         curICv2=cell(1,3);
%         % storing labelled image.
%         curICv2{1}=L(:,:,iselFrm);
% 
%         % storing results.
%         newMat=zeros(maxCIDonEachSelFrm(iselFrm),8);
%         meshTemp=cell(maxCIDonEachSelFrm(iselFrm),1);
%         meshLength=zeros(maxCIDonEachFrm(iselFrm),1);
% 
%         for iCID=1:maxCIDonEachSelFrm(iselFrm)%number of cells on this frame
%             entIdx=find(selVCID(:,1)==selFrms(iselFrm) & selVCID(:,2)==iCID);
%             indx=L(:,:,iselFrm)==iCID;
%             [Area,Volume]=getAreaVolume(cellMesh{entIdx});
% 
%             % storing newMat [binarySize, backboneLength, width, OuftiV, OuftiS, ModelS, ModelV, MajorAxisLength]
%             newMat(iCID,1)=sum(indx,'all');% binary size
%             newMat(iCID,2)=WL{entIdx}(2);% backbone length (pole-to-pole)
%             newMat(iCID,3)=WL{entIdx}(1);% width
%             newMat(iCID,4)=Volume;% Oufti Volume from mesh
%             newMat(iCID,5)=Area;  % Oufti Area from mesh
%             newMat(iCID,6)=WL{entIdx}(1)*(WL{entIdx}(2)-WL{entIdx}(1))+pi*(WL{entIdx}(1)^2)/4; %Model area
%             newMat(iCID,7)=(WL{entIdx}(2)-WL{entIdx}(1))*pi*(WL{entIdx}(1)^2)/4 + 4/3*pi*(WL{entIdx}(1)^3)/8;% Model volume
%             MAL=regionprops(indx,'MajorAxisLength');
%             newMat(iCID,8)=MAL.MajorAxisLength;
%             meshTemp{iCID}=cellMesh{entIdx};
%             meshLength(iCID)=size(cellMesh{entIdx},1);
%         end
%         curICv2{2}=newMat;
% 
%         % storing newMesh
%         newMesh=single(zeros(max(meshLength),4,maxCIDonEachSelFrm(iselFrm)));
%         for iCID=1:maxCIDonEachSelFrm(iselFrm)
%             newMesh(:,:,iCID)=[meshTemp{iCID};zeros(max(meshLength)-meshLength(iCID),4)];
%         end
%         curICv2{3}=newMesh;
% 
%         prevICv2{selFrm(iselFrm)}=curICv2;

        Umsk=[];
        uList=[];
%     end
end
end

%% neseted functions
function [Area,Volume]=getAreaVolume(mesh)
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