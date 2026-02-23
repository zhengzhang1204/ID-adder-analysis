function []=Neo_TakeInYmetWarp(Target,numPos,rot,uitCIy,uitCIx)
% init for general.
prm=rot.prm;RotDeg=rot.deg;isInv=rot.isMouthDown;
load([Target,'\temp\p',num2str(numPos,'%0.2i'),'_FileNameList.mat']);
FileNameList=FileNameList(:,:,2);
gsKern = ndGauss([3 3],[18 18]);

% init for warp
load([Target,'\temp\align.mat']);
if shift(1)<=0 && shift(2)>0
    func=@(x,y) fc1(x,y);
elseif shift(1)>0 && shift(2)>0
    func=@(x,y) fc2(x,y);
elseif shift(1)<=0 && shift(2)<=0
    func=@(x,y) fc3(x,y);
elseif shift(1)>0 && shift(2)<=0
    func=@(x,y) fc4(x,y);
end

%% 1. Read in Y and get ready for LLR
uitY=uint16(zeros(256,prm(4)-prm(3)+1,size(FileNameList,1)));
TSstack=uint16(zeros(256,size(FileNameList,1)*32,size(uitCIx,2)));
for iFrm=1:size(FileNameList,1)%for each frame
    % 1. Read in raw Y, warp
    J=imread(FileNameList(iFrm,1));
    C=func(imwarp(J,tform,'OutputView',SAI),shift);
    
    % 2. rotate warpped Y
    if isInv
        C=flipud(C);
    end
    C=imrotate(C,RotDeg,'bilinear');
    C=C(prm(1):prm(2),prm(3):prm(4),:);
    
    % 3. rotate warpped Y. uitY, TSstackA are ready for use.
    %TSstackA: tilized Y, the 3rd dimension is frame
    %uitY: 256 * x uint16 images, in 3D.
    imwrite(C(uitCIy(iFrm,3):uitCIy(iFrm,4),:),FileNameList(iFrm,2));%for storage;
    uitY(:,:,iFrm)=C(uitCIy(iFrm,1):uitCIy(iFrm,2),:);%for use
    for iSchN=1:size(uitCIx,2)%for each side channel
        if isnan(uitCIx(1,iSchN,iFrm))
            continue;
        end
        TSstack(:,(iFrm-1)*32+1:iFrm*32,iSchN)=C(uitCIy(iFrm,1):uitCIy(iFrm,2),uitCIx(1,iSchN,iFrm):uitCIx(2,iSchN,iFrm));
    end
end
clear('C');

%% 2. get Y_res.mat
curLLRY=cell(1,size(TSstack,3));
Est=cell(1,size(TSstack,3));
meTSstack=zeros(size(TSstack,3),2);
for j4=1:size(uitCIx,2)%parfor each side channel
    [curLLRY{1,j4},Est{j4},meTSstack(j4,:)]=DL_getLLR(double(TSstack(:,:,j4)),gsKern,1.6*ones(size(TSstack(:,:,j4))));
end

save([Target,'\Y_res\p',num2str(numPos,'%0.2i'),'_resY.mat'],'Est','TSstack','meTSstack','-v7.3');%Est:A1,LLR; FLTSstack; meTSstack: median and mean intensity.
clear('Est','TSstack','meTSstack');

%% 3. get threshold for LLR
Thr=getThr4Schn(curLLRY);

%% 4. Saving results:
save([Target,'\Y_res\p',num2str(numPos,'%0.2i'),'_resY.mat'],'Thr','-append');
end





%% nested functions
function [uit]=getThr4Schn(curLLR)
%get Thr for each side channel
nTile=length(curLLR);
uit=zeros(nTile,1);
for iSchn=1:nTile
    B=curLLR{1,iSchn};
    C=ones(32,1)*max(B,[],1);
    uit(iSchn,1)=multithresh([B;C]);
end
end
function [x]=fc1(x,shift)
if shift(1)~=0
    V=x(:,1:-shift(1));
    x(:,1:-shift(1))=[];
    x=[x,V];
end
if shift(2)~=0
    W=x(end-shift(2):end,:);
    x(end-shift(2):end,:)=[];
    x=[W;x];
end
end
function [x]=fc2(x,shift)
if shift(1)~=0
    V=x(:,end-shift(1):end);
    x(:,end-shift(1):end)=[];
    x=[V,x];
end
if shift(2)~=0
    W=x(end-shift(2):end,:);
    x(end-shift(2):end,:)=[];
    x=[W;x];
end
end
function [x]=fc3(x,shift)
if shift(1)~=0
    V=x(:,1:-shift(1));
    x(:,1:-shift(1))=[];
    x=[x,V];
end
if shift(2)~=0
    W=x(1:-shift(2),:);
    x(1:-shift(2),:)=[];
    x=[x;W];
end
end
function [x]=fc4(x,shift)
if shift(1)~=0
    V=x(:,end-shift(1):end);
    x(:,end-shift(1):end)=[];
    x=[V,x];
end
if shift(2)~=0
    W=x(1:-shift(2),:);
    x(1:-shift(2),:)=[];
    x=[x;W];
end
end