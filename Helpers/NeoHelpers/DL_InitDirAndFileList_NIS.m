function []=DL_InitDirAndFileList_NIS(sourceEXP,targetEXP,chNx,ValidTh,ValidXY,isRot,RotDeg,isMouthDown)
%Version: v1.01 @23.5.31
% Core:
F=dir([sourceEXP,'\*.tif']);
ch=zeros(10,1);
shifts=1-ValidTh(1);%normally 0;
% test how many channels are involved: (optional)
for tst=1:10
    tsn=F(tst).name;
    usl3=strfind(tsn,'c');
    ch(tst,1)=str2double(tsn(usl3+1));
end
uls4=chNx;
% for each channel, make a new folder to hold corresponding images.
for color=1:length(uls4)
    mkdir( [targetEXP,'\',chNx{color}]);
end

% decipher filename, extract idx of XY, T, C
Nlist={'c';'xy';'t';'.'};
idx=zeros(1,3);
idx(1)=strfind(tsn,Nlist{1});
idx(2)=strfind(tsn,Nlist{2});
idxT=strfind(tsn,Nlist{3});
idx(3)=idxT(1);
idx(4)=strfind(tsn,Nlist{4});
[idx,Sid]=sort(idx,'ascend');
Nxlist=Nlist(Sid);
for i=1:3
    Xid=find(strcmp(Nxlist,Nlist{i})==1);
    if i==1
        Cidx=idx(Xid)+1;
    elseif i==2
        XYidx=idx(Xid)+2:1:idx(Xid+1)-1;
    elseif i==3
        Tidx=idx(Xid)+1:1:idx(Xid+1)-1;
    end
end
suffix=convertCharsToStrings(tsn(1));

%determine frame range and selected pos
if isempty(ValidTh)
    stFrm=1;endFrm=max(cellfun(@(x) str2double(x(Tidx)), {F.name}),[],"all");
else
    stFrm=ValidTh(1);endFrm=min(ValidTh(2),max(cellfun(@(x) str2double(x(Tidx)), {F.name}),[],"all"));
end

if isempty(ValidXY)
    ValidXY=unique(cellfun(@(x) str2double(x(XYidx)), {F.name}));
end


test=["Cstr","XYstr","Tstr"];
[~,order]=sort([Cidx;XYidx(1);Tidx(1)]);

test=test(order);

%generate In_imageName list:
CoreTstr=convertCharsToStrings(arrayfun(@(x) ['t',num2str(x,['%0.',num2str(Tidx(2)-Tidx(1)+1),'i'])], (stFrm:1:endFrm)','UniformOutput',false));
CoreCstr=convertCharsToStrings(arrayfun(@(x) ['c',num2str(x,'%0.1i')], (1:1:length(chNx))','UniformOutput',false));
CoreXYstr=convertCharsToStrings(arrayfun(@(x) ['xy',num2str(x,'%0.2i')], ValidXY','UniformOutput',false));
CoreXYnum=ValidXY';
CoreCnum=(1:1:length(chNx))';

Tstr=repmat(CoreTstr,length(CoreCstr)*length(CoreXYstr),1);
XYstr=reshape(repmat(CoreXYstr',length(CoreTstr),1),[],1);
XYnum=reshape(repmat(CoreXYnum',length(CoreTstr),1),[],1);
XYstr=repmat(XYstr,length(CoreCstr),1);
XYnum=repmat(XYnum,length(CoreCstr),1);
Cstr=repmat(CoreCstr',length(CoreXYstr)*length(CoreTstr),1);
Cnum=repmat(CoreCnum',length(CoreXYstr)*length(CoreTstr),1);
Cstr=reshape(Cstr,[],1);
Cnum=reshape(Cnum,[],1);

M=struct();
M.Tstr=Tstr;
M.Cstr=Cstr;
M.XYstr=XYstr;

A=M.(test(1));
B=M.(test(2));
C=M.(test(3));
Q122=strcat(sourceEXP,"\",suffix,A,B,C,".tif");

%generate Out_imageName list:
CoreTstr_out=convertCharsToStrings(arrayfun(@(x) ['_t',num2str(x,'%0.4i'),'.tif'], (1:1:endFrm-stFrm+1)','UniformOutput',false));
CorePos_out=convertCharsToStrings(arrayfun(@(x) ['\p',num2str(x,'%0.2i')], ValidXY','UniformOutput',false));
CoreCstr_out=convertCharsToStrings(cellfun(@(x) ['\',x],chNx','UniformOutput',false));

Tstr_out=repmat(CoreTstr_out,length(CoreCstr)*length(CoreXYstr),1);
XYstr_out=reshape(repmat(CorePos_out',length(CoreTstr),1),[],1);
XYstr_out=repmat(XYstr_out,length(CoreCstr),1);
Cstr_out=repmat(CoreCstr_out',length(CoreXYstr)*length(CoreTstr),1);
Cstr_out=reshape(Cstr_out,[],1);
P122=strcat(targetEXP,Cstr_out,XYstr_out,XYstr_out,Tstr_out);

for Pos=ValidXY%for each position:
    FNLtemp=[Q122(XYnum==Pos),P122(XYnum==Pos)];
    Ctemp=Cnum(XYnum==Pos);
    FileNameList=[];
    for iC=1:length(chNx)
        FileNameList=cat(3,FileNameList,FNLtemp(Ctemp==iC,:));
    end

    for color=1:length(chNx)
        mkdir( [targetEXP,'\',chNx{color},'\p',num2str(Pos,'%0.2i')]);
    end
    save([targetEXP,'\temp\p',num2str(Pos,'%0.2i'),'_FileNameList.mat'],'FileNameList');
end

if isRot
    rot.prm=getRotCropParam(imread(strcat(sourceEXP,'\',F(1).name)),RotDeg);
    rot.deg=RotDeg;
    rot.isMouthDown=isMouthDown;
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
