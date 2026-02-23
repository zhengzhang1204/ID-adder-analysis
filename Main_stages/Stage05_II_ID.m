cllc;

%% change here:
load('G:\ATCpath2.mat');
% sel=[1:1:3,6,8,10,12:15];
% Rpath=Rpath(sel,:);
Rpath=ATCpath;
path=cellfun(@(x,y) [x,y], repmat({'g:\'},size(Rpath,1),1),Rpath(:,1),'UniformOutput',false);
exportDir='D:\mChr_RescaleL260131\';%'D:\mChr_Rescale\';'D:\yPet_Rescale\'
nSD=3;
ThrPercent=0.005;%remove the bin if sample size in the bin is lower than this ratio.
edgesX=0.775:0.05:1.275;% bin edges, may center at 1;%: 0.875:0.05:1.25
% edgesX2=0.65:0.1:1.35;
scaleX=150;
scaleY=150;
isRemOutlier=true;
isExample1K=true;

%% codes start here:
outMat=zeros(length(path),269);

for i=1:length(path)%<-----CHANGE THIS BACK TO 1
%in path:
    M=load([path{i},'\M_res\Tab_GRIDTADDERRv3.mat']);
    Y=load([path{i},'\resV2.mat'],'outA','outGHI','int');
    ADLc=load([path{i},'\AdderLinV3.mat']);

    if ~exist([exportDir,Rpath{i,2}],'dir')
        mkdir([exportDir,Rpath{i,2}]);
    end
%Row 1: li_Did, li_Tid, li_Ld
    li_id=Y.outGHI(:,10);
    did=Y.outGHI(:,6)-Y.outGHI(:,10);
    T_id=Y.outGHI(:,1)*Y.int;
    Ld_id=Y.outGHI(:,6);
%Row 2: li_Dii, li_Tii, li(n)_li(n+1)
    li_ii=[Y.outA(:,4);Y.outA(:,4)];
    dii=[Y.outA(:,5);Y.outA(:,6)];
    T_ii=[Y.outA(:,10);Y.outA(:,11)]*Y.int;
    lri=li_ii+dii;%reinitiation length. Equivalent to Ld.
%Row 3: Lb_Ddd, Lb_Tdd, Lb_Ddd
    lb=M.tab.BirthLength;
    ddd=M.tab.DivAdderByCombDaug;
    Ld_dd=lb+ddd;
    T_dd=M.tab.IDT;
%Row 4: 
%Dii_Dii, Did_Did, Ddd_Ddd
%Row 5:
    li_ii_half=Y.outA(:,4);
    Ld_closest=Y.outA(:,18);
    

%% labels:
%Row 1: li_Did, li_Tid, li_Ld   -   outMat(1~12)
    lab11.xlab='Rescaled li';
    lab11.ylab='Rescaled \delta_{id}';

    lab12.xlab='Rescaled li';
    lab12.ylab='Rescaled \tau_{id}';

    lab13.xlab='Rescaled li';
    lab13.ylab='Rescaled L_d';

%Row 2: li_Dii, li_Tii, li(n)_li(n+1)   -   outMat(13~24)
    lab21.xlab='Rescaled li';
    lab21.ylab='Rescaled \delta_{ii}';

    lab22.xlab='Rescaled li';
    lab22.ylab='Rescaled \tau_{ii}';

    lab23.xlab='Rescaled li';
    lab23.ylab='Rescaled lreinit';

%Row 3: Lb_Ddd, Lb_Tdd, Lb_Ld   -   outMat(25~36)
    lab31.xlab='Rescaled Lb';
    lab31.ylab='Rescaled \delta_{dd}';

    lab32.xlab='Rescaled Lb';
    lab32.ylab='Rescaled \tau_{id}';    % IDT

    lab33.xlab='Rescaled Lb';
    lab33.ylab='Rescaled Ld';


%Row 4: Dii_Dii, Did_Did, Ddd_Ddd   -   outMat(37~48)
    lab41.xlab='Rescaled current \delta_{ii}';
    lab41.ylab='Rescaled next \delta_{ii}';

    lab42.xlab='Rescaled current \delta_{id}';
    lab42.ylab='Rescaled next \delta_{id}';

    lab43.xlab='Rescaled current \delta_{dd}';
    lab43.ylab='Rescaled next \delta_{dd}';

% %Row 4a: Li_Dit, Li_Tit, Li_Lt   -   outMat(37~48)-debugging 20250612
%     lab41.xlab='Rescaled Li';
%     lab41.ylab='Rescaled \delta_{it}';
%     Dit=Y.outA(:,20).*Y.outA(:,21)-Y.outA(:,4);
%     Li_ii=Y.outA(:,4);
% 
%     lab42.xlab='Rescaled Li';
%     lab42.ylab='Rescaled next \tau_{it}'; 
%     Tit=Y.outA(:,9);
% 
%     lab43.xlab='Rescaled Li';
%     lab43.ylab='Rescaled next Lt';
%     Lt=Y.outA(:,20).*Y.outA(:,21);


%Row 5: misc [liiiN_liiiN+1]
    lab51.xlab='Rescaled current l_{i}';
    lab51.ylab='Rescaled next l_{i}';
    lab52.xlab='Rescaled l_{i}';
    lab52.ylab='Rescaled nearest \Delta_{id}';
    lab53.xlab='Rescaled \lambdaii';
    lab53.ylab='Rescaled \Deltaii';

%Row 6: CVDdd, CVdDii, CVDid
    lab61.xlab='\Deltadd';
    lab61.ylab='Count';
    lab62.xlab='\Deltaii';
    lab62.ylab='Count';
    lab63.xlab='\Deltaid';
    lab63.ylab='Count';
%% send:
%Row 1: li_Did, li_Tid, li_Ld 
    [outMat(i,1:22)]=drawFig250727b(edgesX,li_id,did,nSD,[exportDir,Rpath{i,2},'\1_li_Did.tif'],ThrPercent,isExample1K);
    [outMat(i,23:44)]=drawFig250727b(edgesX,li_id,T_id,nSD,[exportDir,Rpath{i,2},'\2_li_Tid.tif'],ThrPercent,isExample1K);
    [outMat(i,45:66)]=drawFig250727b(edgesX,li_id,Ld_id,nSD,[exportDir,Rpath{i,2},'\3_li_Ld.tif'],ThrPercent,isExample1K);
    
%Row 2: li_Dii, li_Tii, li(n)_li(n+1) 
    [outMat(i,67:88)]=drawFig250727b(edgesX,li_ii,dii,nSD,[exportDir,Rpath{i,2},'\4_li_Dii.tif'],ThrPercent,isExample1K);
    [outMat(i,89:110)]=drawFig250727b(edgesX,li_ii,T_ii,nSD,[exportDir,Rpath{i,2},'\5_li_Tii.tif'],ThrPercent,isExample1K);
    [outMat(i,111:132)]=drawFig250727b(edgesX,li_ii,lri,nSD,[exportDir,Rpath{i,2},'\6_li(n)_li(n+1).tif'],ThrPercent,isExample1K);
    
%Row 3: Lb_Ddd, Lb_Tdd, Lb_Ld 
    [outMat(i,133:154)]=drawFig250727b(edgesX,lb,ddd,nSD,[exportDir,Rpath{i,2},'\7_Lb_Ddd.tif'],ThrPercent,isExample1K);
    [outMat(i,155:176)]=drawFig250727b(edgesX,lb,T_dd,nSD,[exportDir,Rpath{i,2},'\8_Lb_Tdd.tif'],ThrPercent,isExample1K);
    [outMat(i,177:198)]=drawFig250727b(edgesX,lb,Ld_dd,nSD,[exportDir,Rpath{i,2},'\9_Lb_Ld.tif'],ThrPercent,isExample1K);

%Row 4: Dii_Dii, Did_Did, Ddd_Ddd
    [outMat(i,199:220)]=drawFig250727b(edgesX,ADLc.IIAdderLin(:,1),ADLc.IIAdderLin(:,2),nSD,[exportDir,Rpath{i,2},'\10_Dii_Dii.tif'],ThrPercent,isExample1K);
    [outMat(i,221:242)]=drawFig250727b(edgesX,ADLc.IDAdderLin(:,1),ADLc.IDAdderLin(:,2),nSD,[exportDir,Rpath{i,2},'\11_Did_Did.tif'],ThrPercent,isExample1K);
    [outMat(i,243:264)]=drawFig250727b(edgesX,ADLc.DDAdderLin(:,1),ADLc.DDAdderLin(:,2),nSD,[exportDir,Rpath{i,2},'\12_Ddd_Ddd.tif'],ThrPercent,isExample1K);
%Row 5: misc (li_ii(n) vs. li_ii(n+1); li_ii vs closest ld)
    % [outMat(i,97:104)]=drawFig250727(edgesX,li_ii,(li_ii+dii)/2,nSD,[exportDir,Rpath{i,2},'\13_liiiN_liiiNPlusOne.tif'],ThrPercent,scaleX,scaleY,isRemOutlier);

%Row 6: avgDdd, avgDii, avgDid, CVDdd, CVDii, CVDid, GR, avgTid
    outMat(i,265)=mean(ddd,'all','omitnan')*0.065;%mean(ddd)
    outMat(i,266)=mean(dii,'all','omitnan')*0.065;%mean(dii)
    outMat(i,267)=mean(did,'all','omitnan')*0.065;%mean(did)
    % outMat(i,108)=drawFigHistCV(lb,ddd,nSD,[exportDir,Rpath{i,2},'\15_Hist_Ddd.tif'],Rpath{i,2},lab61);
    % outMat(i,109)=drawFigHistCV(li_ii,dii,nSD,[exportDir,Rpath{i,2},'\16_Hist_Dii.tif'],Rpath{i,2},lab62);
    % outMat(i,110)=drawFigHistCV(li_id,did,nSD,[exportDir,Rpath{i,2},'\17_Hist_Did.tif'],Rpath{i,2},lab63);
    outMat(i,268)=median(M.tab.GR,"all","omitnan");%GR
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
T0 = array2table(Rpath(:,2),'VariableNames',{'Media'});
T=[T0,T1];
writetable(T,[exportDir,'\results.xlsx']);

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