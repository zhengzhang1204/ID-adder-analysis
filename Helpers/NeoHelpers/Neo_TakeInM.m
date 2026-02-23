function [storePH,cropParam,uitCIy,tK]=Neo_TakeInM(Target,numPos,rot,Core,XcrossShift)%XcrossShift is negative.
if isfield(rot,'deg')
    prm=rot.prm;RotDeg=rot.deg;isInv=rot.isMouthDown;tK=zeros(1,2);
    isRot=true;
else
    prm=rot.prm;isInv=rot.isMouthDown;tK=zeros(1,2);RotDeg=0;
    isRot=false;
end
load([Target,'\temp\p',num2str(numPos,'%0.2i'),'_FileNameList.mat']);
% XcrossShift=-23;

% Initialization:
inPath=FileNameList(:,1,1);
uitCIy=zeros(length(inPath),2);
uitCIy4storage=zeros(length(inPath),2);
M4Store=cell(length(inPath),1);
cropParam=cell(length(inPath),1);
storePH=cell(length(inPath),1);

% Core loop:
% PP=parpool('Processes');
%tic
for i=1:length(inPath)%PARFOR, 'PROCESSES', for each image:
% read and reshape
    J=imread(inPath(i));
    if isInv
        J=flipud(J);
    end
    if isRot
        J=imrotate(J,RotDeg,'bilinear');
        J=J(prm(1):prm(2),prm(3):prm(4));
    end
% look for cropping range in Y
    ChnYRng=foundSchNRng(J,Core);
    uitCIy(i,:)=[ChnYRng(2)-(255-XcrossShift),ChnYRng(2)+XcrossShift];%for 4080 wafer, use these two arbitary values. The sum of the two numbers should be 255!

% identify each side channels and extract PH side channel pictures.
    [storePH{i},cropParam{i}]=getSchNinPH(J(ChnYRng(2)-(255-XcrossShift):ChnYRng(2)+XcrossShift,:));% the sum of the two numbers should be 255!

% Crop image for storaging.
    ChnYRng(1)=max([1,ChnYRng(1)-10]);%was 50 before 240723
    ChnYRng(2)=min([size(J,1),ChnYRng(2)+XcrossShift+20]);%was 50 before 240723. [ChnYRng] is the Y range for storaging pictures.
    uitCIy4storage(i,:)=ChnYRng;
    M4Store{i}=J(ChnYRng(1):ChnYRng(2),:);%this is for storaging.
end
%tK(1)=toc;
uitCIy=[uitCIy,uitCIy4storage];
% save M for storage.
tic
for iFrm=1:length(inPath)
    imwrite(M4Store{iFrm}, [Target,'\M\p',num2str(numPos,'%0.2i'),'\p',num2str(numPos,'%0.2i'),'_t',num2str(iFrm,'%0.4i'),'.tif']);
end
tK(2)=toc;
end

%% nested functions:
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
b=b(d>0.3*quantile(d,0.9));%channel positions from findpeaks:

% do prediction
X=cumsum(round(diff(b)./75));
Y=b(2:end);
f=polyfit(X,Y,1);
X_pred=-15:1:length(b)+15;
Y_pred=f(1)*X_pred+f(2);
Y_pred(Y_pred<20 |Y_pred>xEdge-20)=[];

% correct channel positions
b_out=Y_pred;
for j=1:length(Y_pred)
    cand=find(abs(Y_pred(j)-b)<15);
    if isscalar(cand)
        b_out(j)=b(cand);
    elseif length(cand)==2
        b_out(j)=b(cand(1));
    elseif isempty(cand)
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
