cllc;

Source='G:\20250214_AC2_M13_4min';
Target='E:\20250214_AC2_M13_4minProc';
Channels={'M','Y'};             %'M': brightfield or PH, 'Y': spots, 'R': hupA-mChr
interval=1;int=4;   %
ValidFrame=[1, 5];            %Valid frames. Put empty([]) if all frames are valid.
ValidPos=[4];isPosStartFrom0=false;        %Valid positions. If from Fov_00, ValidPos=10 means processing Fov_10, generating P11.

narrow.isNarrowChn=true; narrow.width=44; narrow.restoreOffset=-2; narrow.optimMethod=1;   %By Stg0Neo_test_DL_width.m
isRotate=true; rotateAngle=-0.2; %By imageJ, add '-' to the rotation degree in ImageJ if mouths are up. Add nothing to the degree if mouths are down.
isFacingDown=false;                 %set to 0 if the side channels are facing upwards. Need to automatically shift up side down if true.
isUseGPU=false;
isWarp=true;     %is channel shifted?
XcrsShift=103;                  %By Stg0Neo_test_DL_width.m
FromStage=0;%debugging only

% %-----JinMicroscope only-----(seqAmChr used this)
JinNameFolderPrefix='field';JinNameFolderPrefixDecimal=4;
JinNameChn1FolderPrefix='BF1';JinNameChn1ImgPrefix='imageBF1';%suffix before number.
JinNameChn2FolderPrefix='Venus';JinNameChn2ImgPrefix='imageVenus';%suffix before number.
JinNameImgIndexDecimal=5;

%-----JinMicroscope fov_format-----
% JinNameFolderPrefix='fov_';JinNameFolderPrefixDecimal=0;
% JinNameChn1FolderPrefix='phase';JinNameChn1ImgPrefix='t';%suffix before number.
% JinNameChn2FolderPrefix='YFP';JinNameChn2ImgPrefix='t';%suffix before number.
% JinNameImgIndexDecimal=0;


%do not touch from here:
dlmodelLocation='D:\Work-Program\Neo2\stills\deltaNet.mat';
PFSloc='D:\Work-Program\Neo2\stills\DL_PSF.mat';
Xpatternloc='D:\Work-Program\Neo2\stills\CrossPatternFLN.mat';
TwoDNXCloc='D:\Work-Program\Neo2\stills\CoreNew.tif';
JinName={{JinNameFolderPrefix,JinNameFolderPrefixDecimal},{JinNameChn1FolderPrefix,JinNameChn1ImgPrefix},{JinNameChn2FolderPrefix,JinNameChn2ImgPrefix},JinNameImgIndexDecimal};
gsKern = ndGauss([3 3],[18 18]);%was [3][18] works for 80 min samples (1-2 spots).
load(PFSloc);
minArea=140;
majorL=14;


%% 0. generate file name list and initiate folders
% The FileNameList.mat stores the input and output file paths. It is a 'nFrm x 2 x nChn' string matrix. 2 means input and output paths.
if FromStage<=0%<----change this back to 1
    mkdir(strcat(Target,'\temp'));
    type=DL_validateFileFormatv2(Source,JinName);
    if type==1; DL_InitDirAndFileList_NIS(Source,Target,Channels,ValidFrame,ValidPos,isRotate,rotateAngle,isFacingDown);%-tested OK
    elseif type==2 || type==3; DL_InitDirAndFileList_CPatHuangV2(Source,Target,Channels,ValidFrame,ValidPos,JinName,isRotate,rotateAngle,isFacingDown,interval,isPosStartFrom0);
%     elseif type==3
%         JinName={'fov',{'phase','t'},{'YFP','t'},0};
%         DL_InitDirAndFileList_CPatHuang(Source,Target,Channels,ValidFrame,ValidPos,JinName,isRotate,rotateAngle,isFacingDown,interval,isPosStartFrom0);%-tested OK
    else; error('Unrecognized file format');
    end
end
if ~exist([Target,'\Y_res'],'dir')
    mkdir([Target,'\Y_res']);
end
if ~exist([Target,'\M_res'],'dir')
    mkdir([Target,'\M_res']);
end
clear('JinName*');

%% 1. Read and process M (par process)
% Read in images, rotate, crop Y coordinates, crop X coordinates, PH tiles.
% In: path of input image, rotation and cropping parameters, 
% Out: 
% ChnYRng, 
% M2Store (cropped M images for storaging), 
% storePH (tilized PH image with ), 
% cropPararam(X range for side channel.), 
% cropping (Y range for side channels)


if isRotate
    load([Target,'\temp\rot.mat']);
    prm=rot.prm;RotDeg=rot.deg;isInv=rot.isMouthDown;
end
Core=imread(TwoDNXCloc);
if size(Core,2)>(prm(4)-prm(3)+1)
    Core=Core(:,1:(prm(4)-prm(3)+1));
end
if isempty(ValidPos)
    F1x=dir([Target,'\temp\p*_FileNameList.mat']);
    ValidPos=cellfun(@(x) str2double(x(2:3)),{F1x.name},'UniformOutput',true);
    clear("F1x");
end
if isPosStartFrom0
    ValidPos=ValidPos+1;
end

for iPos=1:length(ValidPos)%for each valid position
    load([Target,'\temp\p',num2str(ValidPos(iPos),'%0.2i'),'_FileNameList.mat']);
%% 1. Read in PH images:
    disp(['1. Reading p',num2str(ValidPos(iPos),'%0.2i'),'...']);
    [storePH,cropParam,uitCIy,tK]=Neo_TakeInM(Target,ValidPos(iPos),rot,Core,XcrsShift);
    disp(['1. p',num2str(ValidPos(iPos),'%0.2i'),' is read in ',num2str(tK(1),'%0.2f'),' seconds, wrote in ',num2str(tK(2),'%0.2f'),' seconds']);

%% 2. Align and extract side channels
    [uitPH4seg,uitPH,uitCIx,narrowParams,nFrm,nSchN]=Neo_AlignAndExtractSchN(cropParam,storePH,narrow,psf_c);
    clear('storePH','cropParam');
    disp(['2. p',num2str(ValidPos(iPos),'%0.2i'),' is tilized.']);
% 
%% 3. do DL
    %optimized for gpu usage. since 22.08.18 v2.05.
    %input: uitPH4seg [256 x 32 x 1 x [nFrm*nSchN]];
    %output: Masks  [256 x 32 x 1 x [nFrm*nSchN]];

    if ~exist('net','var'); load(dlmodelLocation); end
    tic
    if ~isUseGPU
        Masks=predict(net,uitPH4seg,'MiniBatchSize',64);
    else
        Masks=predict(net,uitPH4seg);
    end
    
    clear("uitPH4seg");
    if narrow.optimMethod==1
        parfor iFrm=1:nFrm
            pick=(0:1:nSchN-1)*nFrm+iFrm;
            curFrm=squeeze(mat2cell(squeeze(Masks(:,:,1,pick)),256,32,ones(1,nSchN)));
            PosMask{1,iFrm}=cellfun(@(x) DL_ModMasksMet1(x,minArea,majorL,narrow.isNarrowChn,narrowParams,narrow.restoreOffset), curFrm,'UniformOutput',false);
        end
    elseif optimMethod==2
        parfor iFrm=1:nFrm
            pick=(0:1:nSchN-1)*nFrm+iFrm;
            curFrm=squeeze(mat2cell(squeeze(Masks(:,:,1,pick)),256,32,ones(1,nSchN)));
            PosMask{1,iFrm}=cellfun(@(x) DL_ModMasksMet2(x,minArea,majorL,narrow.isNarrowChn,narrowParams,narrow.restoreOffset), curFrm,'UniformOutput',false);
        end
    end
    clear('Masks');
    tK2=toc;
    disp(['3. p',num2str(ValidPos(iPos),'%0.2i'),' is masked in ',num2str(tK2,'%0.2f'),' seconds']);
    

%% 4. saving Mres.mat
    % write M.res including:
    % uitCIx (X coordinates for side channels on each frame)
    % uitCIy (Y coordinates for side channels on each frame, 1:2 for 256 cropping, 3:4 for storage cropping)
    % uitPH (tilized PH images, to replace images in M_crop);
    % PosMask (masks)
    save([Target,'\M_res\p',num2str(ValidPos(iPos),'%0.2i'),'_resM.mat'],'uitCIx','uitCIy','uitPH','PosMask','narrow');
    clear('PosMask','uitPH');
    disp(['4. p',num2str(ValidPos(iPos),'%0.2i'),' is saved to M_res.']);

%% 5. Process YFP images, warp if necessary.
    tic
    if isWarp
        Neo_TakeInYmetWarp(Target,ValidPos(iPos),rot,uitCIy,uitCIx);
    else
        Neo_TakeInY(Target,ValidPos(iPos),rot,uitCIy,uitCIx);%<-----------------not done.
    end
    tK3=toc;
    disp(['5. p',num2str(ValidPos(iPos),'%0.2i'),' is saved to Y_res in ',num2str(tK3,'%0.2f'),' seconds']);
end




%% Nested functions.
function [SchNRng]=foundSchNRng(J,Core)
CoreSize2=round(size(Core)/2);
RawSize2=size(J);
C = normxcorr2(Core,J);
C=C(CoreSize2(1):CoreSize2(1)+RawSize2(1)-1,...
    CoreSize2(2):CoreSize2(2)+RawSize2(2)-1);
%[p,l]=findpeaks(sum(C,2),'SortStr','descend');
%p=p(1:min([3,length(p)]));
[~,l]=findpeaks(sum(C,2),'SortStr','descend');
l(l>(size(C,1)-100) | l<100)=[];
l=l(1:min([3,length(l)]));
cnt=1;
found=[];
while cnt<=length(l) && isempty(found)
    if l(cnt)<5 || l(cnt)>RawSize2(1)-5
        cnt=cnt+1;
        continue
    else
        found=l(cnt);
        cnt=cnt+1;
    end
end
if ~isempty(found)
    SchNRng=[found-CoreSize2(1), found+CoreSize2(1)];
else
    SchNRng=[];
end
SchNRng(SchNRng<1)=1;
SchNRng(SchNRng>RawSize2(1))=RawSize2(1);
end

function [storePH,cropParam]=getSchNinPH(J)
xEdge=size(J,2);
[~,b,~,d]=findpeaks(movmean(sum(J,1),20),'MinPeakWidth',10);
minProm=0.3*quantile(d,0.9);
b=b(d>minProm);%channel positions from findpeaks:

% correct channel positions
chk1=find(abs(diff(diff(b)))<18);
yy=[0,find(diff(chk1)>1),length(diff(chk1))];
[~,idx]=max(diff(yy));
sel=[chk1(yy(idx)+1):1:chk1(yy(idx+1))]+2;

f=polyfit(sel,b(sel),1);
X=-5:1:length(b)+5;
Y_pred=f(1)*X+f(2);

idxX=Y_pred<20 |Y_pred>xEdge-20;
Y_pred(idxX)=[];
b_out=Y_pred;
for j=1:length(Y_pred)
    xxxx=find(abs(Y_pred(j)-b)<15);
    if length(xxxx)==1
        b_out(j)=b(xxxx);
    elseif length(xxxx)==2
        b_out(j)=b(xxxx(1));
    elseif isempty(xxxx)
        b_out(j)=round(Y_pred(j));
    end
end

% crop out channels:
cropParam=nan(length(b_out),2);
% storePH=cell(length(b_out),1);
storePH=uint16(zeros(256,length(b_out)*32));

for i=1:length(b_out)
    if b_out(i)+15>xEdge
        offset = b_out(i)+15-xEdge;
        b_out(i)=b_out(i)-offset;
    end
    if b_out(i)-16<1
        offset = 1-b_out(i)+16;
        b_out(i)=b_out(i)+offset;
    
    end
    if b_out(i)+15<=size(J,2)
        Jtemp=J(:,b_out(i)-16:b_out(i)+15);
        %imwrite(Jtemp,'E:\test.tif','WriteMode','append');%delete this<----------
    else
        if i==length(b_out)
            %Jtemp=J(:,b_out(i)-16:size(J,2));
            return;
        else
            error('bad cropping in the middle.');
        end
    end
%     storePH{i,1}=Jtemp;
    storePH(:,(i-1)*32+1:i*32)=Jtemp;
    if b_out(i)+15<=size(J,2)
        cropParam(i,:)=[b_out(i)-16,b_out(i)+15];
    else
        cropParam(i,:)=[b_out(i)-16,size(J,2)];
    end
end
end
