function []=DL_InitDirAndFileList_CPatHuangV2(sourceEXP,targetEXP,chNx,ValidTh,ValidFov,JinName,isRot,RotDeg,isMouthDown,interval,isPosFromZero)
%Version: v1.00 @23.6.2
% preparing directories.
F=dir([sourceEXP,'\',JinName{1}{1},'*']);
decChnName=cell(length(chNx),1);
nFrm=nan(length(chNx),1);
for i=1:length(chNx)
    str1=chNx{i};
    mkdir( [targetEXP,'\',str1]);
    switch str1
        case 'M'
            decChnName{i}=JinName{i+1}{1};
        case 'Y'
            decChnName{i}=JinName{i+1}{1};
        case 'R'
            decChnName{i}=JinName{i+1}{1};
    end
    nFrm(i)=length(dir([F(1).folder,'\',F(1).name,'\',decChnName{i},'\*.tif*']));
end
if sum(nFrm)==length(chNx)*nFrm(1)%all channels has the same number of frames.
    nFrm=nFrm(1);
else
    error('Frame numbers in each channel is not the same!');
end

% initialize output folders
if JinName{4}~=0
    if exist([F(1).folder,'\',F(1).name,'\',JinName{2}{1},'\',JinName{2}{2},num2str(1,['%0.',num2str(JinName{4}),'i']),'.tif'],'file')
        ext='tif';
    end
    if exist([F(1).folder,'\',F(1).name,'\',JinName{2}{1},'\',JinName{2}{2},num2str(1,['%0.',num2str(JinName{4}),'i']),'.tiff'],'file')
        ext='tiff';
    end
else
    if exist([F(1).folder,'\',F(1).name,'\',JinName{2}{1},'\',JinName{2}{2},'1.tif'],'file')
        ext='tif';
    end
    if exist([F(1).folder,'\',F(1).name,'\',JinName{2}{1},'\',JinName{2}{2},'1.tiff'],'file')
        ext='tiff';
    end
end



if isempty(ValidFov)
    ValidFov=(0:1:length(F)-1)';
    if ~isPosFromZero
        ValidFov=ValidFov+1;
    end
end

if ValidTh(2)<nFrm%process all frames
    frmEnd=ValidTh(2);
else
    frmEnd=nFrm-1;
end
frmSt=ValidTh(1);
% shift=ValidTh(1)-1;

posStrSt=length(JinName{1}{1})+1;
rawPosNum=cellfun(@(x) str2double(x(posStrSt:end)),{F.name}');%raw fov number, in the order as in F struct.
validInd=find(ismember(rawPosNum,ValidFov));%index of valid pos with the order as in F struct

% validInd=find(cellfun(@(x) ismember(str2double(x(posStrSt:end)),ValidFov),{F.name}')); %index in F struct
% decPosName=cell(length(validInd),1);
for i=1:length(validInd)
    FileNameList=[];
    if isPosFromZero
        decPosName=['p',num2str(str2double(F(validInd(i)).name(posStrSt:end))+1,'%0.2i')];%'pxx' in output path
    else
        decPosName=['p',num2str(str2double(F(validInd(i)).name(posStrSt:end)),'%0.2i')];%'pxx' in output path
    end

    for j=1:length(chNx)
        mkdir([targetEXP,'\',chNx{j},'\',decPosName]);
        if JinName{4}==0
            iFrm=convertCharsToStrings(arrayfun(@(x) num2str(x),(frmSt:interval:frmEnd)','UniformOutput',false));
        else
            iFrm=convertCharsToStrings(arrayfun(@(x) num2str(x,['%0.',num2str(JinName{4}),'i']),(frmSt:interval:frmEnd)','UniformOutput',false));
        end
        readName=strcat(F(validInd(i)).folder,'\',F(validInd(i)).name,'\',decChnName{j},'\',JinName{j+1}{2},iFrm,'.',ext);
        oFrm=convertCharsToStrings(arrayfun(@(x) num2str(x,'%0.4i'), (1:1:length(iFrm)),'UniformOutput',false));
        writeName=strcat(targetEXP,'\',chNx{j},'\', decPosName,'\',decPosName,'_t',oFrm,'.tif');
        FileNameList=cat(3,FileNameList,[readName,writeName']);
    end

    save(strcat(targetEXP,'\temp\', decPosName,'_FileNameList.mat'),'FileNameList')

end
 


if isRot
    K=dir([sourceEXP,'\',F(1).name,'\',JinName{2}{1},'\*.',ext]);
    rot.prm=getRotCropParam(imread([K(1).folder,'\',K(1).name]),RotDeg);
    rot.deg=RotDeg;
    rot.isMouthDown=isMouthDown;
end
save(strcat(targetEXP,'\temp\rot.mat'),"rot");
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

