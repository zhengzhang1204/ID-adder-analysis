%output simp:
clear;
clc;
close('all');

%% User specified region
inDir='D:\mChr_RescaleL260224';

%% Stills:
nSD=3;
ThrPercent=0.005;%remove the bin if sample size in the bin is lower than this ratio.
edgesX=0.775:0.05:1.275;% bin edges, may center at 1;%: 0.875:0.05:1.25
scaleX=150;
scaleY=150;
isRemOutlier=true;
isExample1K=true;

%% codes start here:
F=dir([inDir,'\*.mat']);
[~,idx]=sort(cellfun(@(x) str2double(extractBetween(x,'M','_data.mat')),{F.name}),'ascend');
F=F(idx);
mediaID=cellfun(@(x) extractBefore(x,'_data.mat'),{F.name},'UniformOutput',false);
outMat=zeros(length(F),269);
for i=1:length(F)
    if ~exist([inDir,'\',mediaID{i}],'dir')
        mkdir([inDir,'\',mediaID{i}]);
    end

% Prepare adder and lineage data:
    load(fullfile(F(i).folder,F(i).name));
    [ID,linID]=ZZ_deMap(tab_ID{:,:},int);
    [DD,linDD]=ZZ_deMap(tab_DD{:,:},int);
    [II,linII]=ZZ_deMap(tab_II{:,:},int);

%Row 1: li_Did, li_Tid, li_Ld
    li_id=ID(:,1);
    did=ID(:,2);
    T_id=ID(:,3);
    Ld_id=ID(:,4);
%Row 2: li_Dii, li_Tii, li(n)_li(n+1)
    li_ii=II(:,1);
    dii=II(:,2);
    T_ii=II(:,3);
    lri=II(:,4);
%Row 3: Lb_Ddd, Lb_Tdd, Lb_Ddd
    lb=DD(:,1);
    ddd=DD(:,2);
    Ld_dd=DD(:,4);
    T_dd=DD(:,3);

% pass data
%Row 1: li_Did, li_Tid, li_Ld 
    [outMat(i,1:22)]=drawFig250727b(edgesX,li_id,did,nSD,[inDir,'\',mediaID{i},'\1_li_Did.tif'],ThrPercent,isExample1K);
    [outMat(i,23:44)]=drawFig250727b(edgesX,li_id,T_id,nSD,[inDir,'\',mediaID{i},'\2_li_Tid.tif'],ThrPercent,isExample1K);
    [outMat(i,45:66)]=drawFig250727b(edgesX,li_id,Ld_id,nSD,[inDir,'\',mediaID{i},'\3_li_Ld.tif'],ThrPercent,isExample1K);
    
%Row 2: li_Dii, li_Tii, li(n)_li(n+1) 
    [outMat(i,67:88)]=drawFig250727b(edgesX,li_ii,dii,nSD,[inDir,'\',mediaID{i},'\4_li_Dii.tif'],ThrPercent,isExample1K);
    [outMat(i,89:110)]=drawFig250727b(edgesX,li_ii,T_ii,nSD,[inDir,'\',mediaID{i},'\5_li_Tii.tif'],ThrPercent,isExample1K);
    [outMat(i,111:132)]=drawFig250727b(edgesX,li_ii,lri,nSD,[inDir,'\',mediaID{i},'\6_li(n)_li(n+1).tif'],ThrPercent,isExample1K);
    
%Row 3: Lb_Ddd, Lb_Tdd, Lb_Ld 
    [outMat(i,133:154)]=drawFig250727b(edgesX,lb,ddd,nSD,[inDir,'\',mediaID{i},'\7_Lb_Ddd.tif'],ThrPercent,isExample1K);
    [outMat(i,155:176)]=drawFig250727b(edgesX,lb,T_dd,nSD,[inDir,'\',mediaID{i},'\8_Lb_Tdd.tif'],ThrPercent,isExample1K);
    [outMat(i,177:198)]=drawFig250727b(edgesX,lb,Ld_dd,nSD,[inDir,'\',mediaID{i},'\9_Lb_Ld.tif'],ThrPercent,isExample1K);

%Row 4: Dii_Dii, Did_Did, Ddd_Ddd
    [outMat(i,199:220)]=drawFig250727b(edgesX,linII(:,1),linII(:,2),nSD,[inDir,'\',mediaID{i},'\10_Dii_Dii.tif'],ThrPercent,isExample1K);
    [outMat(i,221:242)]=drawFig250727b(edgesX,linID(:,1),linID(:,2),nSD,[inDir,'\',mediaID{i},'\11_Did_Did.tif'],ThrPercent,isExample1K);
    [outMat(i,243:264)]=drawFig250727b(edgesX,linDD(:,1),linDD(:,2),nSD,[inDir,'\',mediaID{i},'\12_Ddd_Ddd.tif'],ThrPercent,isExample1K);

%Row 5: avgDdd, avgDii, avgDid, CVDdd, CVDii, CVDid, GR, avgTid
    outMat(i,265)=mean(ddd,'all','omitnan')*0.065;%mean(ddd)
    outMat(i,266)=mean(dii,'all','omitnan')*0.065;%mean(dii)
    outMat(i,267)=mean(did,'all','omitnan')*0.065;%mean(did)

    outMat(i,268)=median(DD_GR,"all","omitnan");%GR
    outMat(i,269)=median(T_id,"all","omitnan");%Tid
    disp(i);
end

%% prepare table:
N1=["li_Did_","li_Tid_","li_Ld_","li_Dii_","li_Tii_","li(n)_li(n+1)_","Lb_Ddd_","Lb_Tdd_","Lb_Ld_","Dii_Dii_","Did_Did_","Ddd_Ddd_"];
N2=["ano", "ano_low", "ano_up", "awo", "awo_low", "awo_up",...
    "bavgno", "bavgno_low", "bavgno_up", "bmedno", "bmedno_low", "bmedno_up",...
    "bavgwo", "bavgwo_low", "bavgwo_up", "bmedwo", "bmedwo_low", "bmedwo_up",...
    "covno", "covwo", "nno", "nwo"];
Nms=paringTitle(N1,N2);
Nms=[Nms,"avgDdd(\mum)","avgDii(\mum)","avgDid(\mum)","GR(h^{-1})","avgTid(min)"];

T1 = array2table(outMat,'VariableNames',Nms);
T0 = array2table(mediaID','VariableNames',{'Media'});
T=[T0,T1];
writetable(T,[inDir,'\results.xlsx']);

function [out]=paringTitle(N1,N2)
out=strings(1,length(N1)*length(N2));
cnt=1;
for i=1:length(N1)
    for j=1:length(N2)
        out(cnt)=strcat(N1{i},N2{j});
        cnt=cnt+1;
    end
end
end

function [XY,LinXY]=ZZ_deMap(in,int)
% adder:
if size(in,2)==5
    in1=in(:,4);in2=in(:,5);
    in1(cellfun(@isempty, in1))=[];
    in2(cellfun(@isempty, in2))=[];
    XY=[cell2mat(cellfun(@(x) [x(1),x(end)-x(1),(length(x)-1)*int,x(end)],in1,'UniformOutput',false));...
        cell2mat(cellfun(@(x) [x(1),x(end)-x(1),(length(x)-1)*int,x(end)],in2,'UniformOutput',false))];
else
    XY=cell2mat(cellfun(@(x) [x(1),x(end)-x(1),(length(x)-1)*int,x(end)],in(:,4),'UniformOutput',false));
end

% map:
map = string(in(:,1:3));
newMap=nan(size(map,1),2);
for i=1:size(map,1)
    [row2, ~] = find(strcmp(map(:,1),map(i,2)));
    if ~isempty(row2)
        newMap(i,1)=row2;
    end
    [row3, ~] = find(strcmp(map(:,1),map(i,3)));
    if ~isempty(row3)
        newMap(i,2)=row3;
    end
end

% lineage:
if size(in,2)==5
    LinXY=cell(size(newMap,1),1);
    for i=1:size(newMap,1)
           Zu1=[];
           Zl1=[];
        if ~isnan(newMap(i,1))
           Z1=in{i,4}(end)-in{i,4}(1);

           if ~isempty(in{newMap(i,1),4})
            Zu1=[Z1,in{newMap(i,1),4}(end)-in{newMap(i,1),4}(1)];
           end

           if ~isnan(in{newMap(i,1),5})
            Zl1=[Z1,in{newMap(i,1),5}(end)-in{newMap(i,1),5}(1)];
           end
        end
           Zu2=[];
           Zl2=[];
        if ~isnan(newMap(i,2))
           Z2=in{i,5}(end)-in{i,5}(1);

           if ~isempty(in{newMap(i,2),4})
            Zu2=[Z2,in{newMap(i,2),4}(end)-in{newMap(i,2),4}(1)];
           end

           if ~isempty(in{newMap(i,2),5})
            Zl2=[Z2,in{newMap(i,2),5}(end)-in{newMap(i,2),5}(1)];
           end
        end
        LinXY{i}=[Zu1;Zl1;Zu2;Zl2];
    end
    LinXY=cell2mat(LinXY);
else
    [Zrow,Zcol]=find(~isnan(newMap));
    LinXY=cell2mat(arrayfun(@(x,y) [XY(x,2),XY(newMap(x,y),2)],Zrow,Zcol,'UniformOutput',false));
end
end