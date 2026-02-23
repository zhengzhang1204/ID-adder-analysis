function []=DL_InitDirAndFileList_CPatHuangV3(sourceEXP,targetEXP,ValidTh,ValidFov,JinName,isRot,RotDeg,isMouthDown,interval,isPosFromZero)
%Version: v3.00 @24.7.23

F=dir([sourceEXP,'\',JinName{1}{1},'*']);
nFrm=nan(size(JinName{3},1),1);
for i=1:size(JinName{3},1)
    str1=JinName{3}(i,3);
    mkdir( [targetEXP,'\',convertStringsToChars(str1)]);
    Z=dir([F(1).folder,'\',F(1).name,'\',convertStringsToChars(JinName{3}(i,1)),'\*.tif*']);
    nFrm(i)=length(Z);%total file numbers for each channel in each fov
    [~,~,ext] = fileparts(fullfile(Z(1).folder,Z(1).name));%.tiff or .tif
end
if sum(nFrm)==size(JinName{3},1)*nFrm(1)%all channels has the same number of frames.
    nFrm=nFrm(1);
else
    error('Frame numbers in each channel is not the same!');
end

% select valid fov
if isempty(ValidFov)
    ValidFov=(0:1:length(F)-1)';
    if ~isPosFromZero
        ValidFov=ValidFov+1;
    end
end
% select valid range of time-lapse
if ValidTh(2)<nFrm%process all frames
    frmEnd=ValidTh(2);
else
    frmEnd=nFrm-1;
end
frmSt=ValidTh(1);
posStrSt=strlength(JinName{1})+1;
rawPosNum=cellfun(@(x) str2double(x(posStrSt:end)),{F.name}');%raw fov number, in the order as in F struct.
validInd=find(ismember(rawPosNum,ValidFov));%index of valid pos with the order as in F struct

% create IO filename list:
for i=1:length(validInd)%for each valid fov 
    FileNameList=[];
    if isPosFromZero
        decPosName=['p',num2str(str2double(F(validInd(i)).name(posStrSt:end))+1,'%0.2i')];%'pxx' in output path
    else
        decPosName=['p',num2str(str2double(F(validInd(i)).name(posStrSt:end)),'%0.2i')];%'pxx' in output path
    end

    for j=1:size(JinName{3},1)%for each imaging channel:
        mkdir([targetEXP,'\',convertStringsToChars(JinName{3}{j,3}),'\',decPosName]);
        if JinName{4}==0
            iFrm=convertCharsToStrings(arrayfun(@(x) num2str(x),(frmSt:interval:frmEnd)','UniformOutput',false));
        else
            iFrm=convertCharsToStrings(arrayfun(@(x) num2str(x,['%0.',num2str(JinName{4}),'i']),(frmSt:interval:frmEnd)','UniformOutput',false));
        end
        readName=strcat(F(validInd(i)).folder,'\',F(validInd(i)).name,'\',JinName{3}{j,1},'\',JinName{3}{j,2},iFrm,ext);
        oFrm=convertCharsToStrings(arrayfun(@(x) num2str(x,'%0.4i'), (1:1:length(iFrm)),'UniformOutput',false));
        writeName=strcat(targetEXP,'\',JinName{3}{j,3},'\', decPosName,'\',decPosName,'_t',oFrm,'.tif');
        FileNameList=cat(3,FileNameList,[readName,writeName']);
    end
    save(strcat(targetEXP,'\temp\', decPosName,'_FileNameList.mat'),'FileNameList')
end

% create rotation parameters
K=dir([sourceEXP,'\',F(1).name,'\',convertStringsToChars(JinName{3}(1,1)),'\*',ext]);
if isRot
    rot.prm=getRotCropParam(imread([K(1).folder,'\',K(1).name]),RotDeg);
    rot.deg=RotDeg;
    rot.isMouthDown=isMouthDown;
else
    X=imfinfo([K(1).folder,'\',K(1).name]);
    rot.isMouthDown=isMouthDown;
    rot.prm=[1,X.Height,1,X.Width];
end
save(strcat(targetEXP,'\temp\rot.mat'),"rot");

% initialize output folders:
for i=1:size(JinName{3},1)
    if ~exist([targetEXP,'\',convertStringsToChars(JinName{3}(i,3)),'_res'],'dir')
        mkdir([targetEXP,'\',convertStringsToChars(JinName{3}(i,3)),'_res']);
    end
end
end

function [param]=getRotCropParam(J,deg)
   % deg=abs(deg);
   Jy=imrotate(J,deg,'bilinear');
   nW=size(Jy,2);
   nH=size(Jy,1);
   W=size(J,2);
   H=size(J,1);
   tg=abs(tand(deg));
   a1=(H*tg-W)/(tg^2-1);
   a2=((tg^2)*W-H*tg)/(tg^2-1);
   x=a1/cosd(deg);
   y=a2/abs(sind(deg));
%   DD=Jy(ceil((nH-y)/2)+1:floor((nH-y)/2+y),ceil((nW-x)/2)+1:floor((nW-x)/2+x));
   param=[ceil((nH-y)/2)+1,floor((nH-y)/2+y),ceil((nW-x)/2)+1,floor((nW-x)/2+x)];
end

