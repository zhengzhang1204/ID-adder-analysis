function []=DL_InitDirAndFileList_Jin(sourceEXP,targetEXP,chNx,ValidTh,ValidFov,JinName,isRot,RotDeg,isMouthDown,interval)
%Version: v1.0 @23.5.31
isFrom0=false;
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
F=dir([sourceEXP,'\',JinName{1}{1},'*']);

%test if fov starts from 00?
if str2double(F(1).name(end-JinName{5}+1:end))==0
    isFrom0=true;
end

decPosName=cell(length(F),1);
if isempty(ValidFov)
    for i=1:length(F)
        ValidFov=[ValidFov,str2double(F(i).name(end-JinName{5}+1:end))];
    end
end
for i=1:length(F)
    if ~ismember(str2double(F(i).name(end-JinName{5}+1:end)),ValidFov)
        continue;
    end
    for j=1:length(chNx)
        if isFrom0
            mkdir([targetEXP,'\',chNx{j},'\p',num2str(str2double(F(i).name(end-JinName{5}+1:end))+1,'%0.2i')]);
        else
            mkdir([targetEXP,'\',chNx{j},'\p',num2str(str2double(F(i).name(end-JinName{5}+1:end)),'%0.2i')]);
        end
    end
    if isFrom0
        decPosName{i}=['p',num2str(str2double(F(i).name(end-JinName{5}+1:end))+1,'%0.2i')];
    else
        decPosName{i}=['p',num2str(str2double(F(i).name(end-JinName{5}+1:end)),'%0.2i')];
    end
end

% prepare transfer:
if ValidTh(2)<nFrm%process all frames
    frmEnd=ValidTh(2);
else
    frmEnd=nFrm-1;
end
frmSt=ValidTh(1);

% do transfer:
for pos=1:length(F)%for each position
    FileNameList=[];
    if ~ismember(str2double(F(pos).name(end-JinName{5}+1:end)),ValidFov)
        continue;
    end
    for iChn=1:length(chNx)%for each channel
        curFileNameList=getReadWritePaths(strcat(sourceEXP,'\',F(pos).name,'\',decChnName{iChn},'\',JinName{iChn+1}{2},'*.tif*'),[frmSt,frmEnd],...
            strcat(targetEXP,'\',chNx{iChn},'\', decPosName{pos},'\',decPosName{pos}),interval);
        FileNameList=cat(3,FileNameList,curFileNameList);
    end
    save(strcat(targetEXP,'\temp\', decPosName{pos},'_FileNameList.mat'),'FileNameList')
end

if isRot
    K=dir([sourceEXP,'\',F(1).name,'\',JinName{2}{1},'\*.tiff']);
    rot.prm=getRotCropParam(imread([K(1).folder,'\',K(1).name]),RotDeg);
    rot.deg=RotDeg;
    rot.isMouthDown=isMouthDown;
    save(strcat(targetEXP,'\temp\rot.mat'),"rot");
end
end

function [outPaths]=getReadWritePaths(in, inRng, out,interv)
F=dir(in);
inPaths=cellfun(@(x,y) strcat(x,'\',y),{F.folder}',{F.name}','UniformOutput',false);
inPaths=inPaths(inRng(1):interv:inRng(2)-1);

outPaths=arrayfun(@(x) [out,'_t',num2str(x,'%0.4i'),'.tif'],(1:1:length(inPaths))','UniformOutput',false);
outPaths=[inPaths,outPaths];


outPaths=convertCharsToStrings(outPaths);
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

