cllc;close('all');

%% change here:
% load('H:\mChrRpath.mat');
% sel=[1:3,6,8,10,12:15];
% Rpath=Rpath(sel,:);scl2=scl2(sel,:);
load('G:\ATCpath2.mat');

ATCpath=[{'20260412_M18_RdnaA-DnaNyPet-atc_0.8'},{'atc0.8'},{3}];
scl2(1,1)=300;scl2(2,1)=800;
Rpath=ATCpath;
path=cellfun(@(x,y) [x,y], repmat({'G:\'},size(Rpath,1),1),Rpath(:,1),'UniformOutput',false);
exportDir='D:\mChr_RescaleL260409b';%'D:\mChr_Rescale\';'D:\yPet_Rescale\'

for i=1:length(path)%<-----CHANGE THIS BACK TO 1
    M=load([path{i},'\M_res\Tab_GRIDTADDERRv3simp.mat']);
    Y=load([path{i},'\resV2simp.mat'],'outGHI');
    int=Rpath{i,3};
    li_ID=Y.outGHI(:,10)*0.065;             % Li
    time_li_ID=Y.outGHI(:,16)*int;          % Time Li

    lb_all=M.tab.BirthLength*0.065;         % Lb
    ddd_all=M.tab.DivAdderByCombDaug*0.065; % added dd
    GR=M.tab.GR;                            % GR
    time_dd=M.tab.CycInitTime;              % time for Lb, dd or GR (unit:min)

    drawFig250402(time_dd,GR,scl2(i,1:4),[exportDir,'\880\',Rpath{i,2},'_ss_GR.tif']);
    drawFig250402(time_dd,lb_all,scl2(i,[1,2,5,6]),[exportDir,'\880\',Rpath{i,2},'_ss_Lb.tif'])
    drawFig250402(time_li_ID,li_ID,[scl2(i,[1,2]),scl2(i,[5,6])*2],[exportDir,'\880\',Rpath{i,2},'_ss_li_ID.tif'])
    disp(i);
end


%% nested functions:
function []=drawFig250402(t,y,scl,svPath)
del=(t<scl(1) | t>scl(2));
t(del)=[];y(del)=[];
t=t-scl(1);
% [t,y];
% MX=max(t,[],'all');
edges=0:50:diff(scl(1:2));
[ctg, be]= discretize(t,edges);
cent=mean([(be(2:end))',(be(1:end-1))'],2);
[a,~,c]=unique(ctg);
a(isnan(a))=[];
cent=cent(a);
avg=zeros(length(a),1);
for j=1:length(a)
    sel=c==a(j);
    avg(j)=mean(y(sel),'all');
end

% F=figure('Position',[100,100,247.6,247.6],'Visible','off');%can not change this
F=figure('Position',[100,100,247,237],'Visible','off');%can not change this

ax=gca;
hold on
scatter(t,y,4,[0.5 0.5 0.5],'filled','MarkerFaceAlpha',0.4);
plot(cent,smoothdata(avg,'gaussian',5),'-','LineWidth',1.5,'Color','k');
% ax.XTickLabel=[];
% ax.YTickLabel=[];
Bx=[cent,smoothdata(avg,'gaussian',5)];
Ax=[t,y];


ax.TickDir='out';
ax.LineWidth=1;
xlim([0 diff(scl(1:2))]);
ylim(scl(3:4));
yL=ylim;xL=xlim;plot(xL(1),yL(1),'.k');
% box on
% set(gca, 'color', 'none');
% set(gca,'XColor', 'none','YColor','none')
hold off
%saving
if ~isempty(svPath)
    saveas(F,svPath);
end
close('all');
end







