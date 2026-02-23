% Version 2: Screening samples with nSD. Remove bins that containing too
% less points. @230803

function [slope,slope2,N,rho]=drawFigDensV2(edgesX,X0,Y0,nSD,svPath,lgdLabel,med,lab,thrPct)
% remove outliers according to rescaling
xM=mean(X0,'omitnan');
yM=mean(Y0,'omitnan');
xLow=xM-nSD*std(X0);
xUp=xM+nSD*std(X0);
yLow=yM-nSD*std(Y0);
yUp=yM+nSD*std(Y0);

del=find(X0<xLow | X0>xUp | Y0<yLow | Y0>yUp | isnan(X0) | isnan(Y0));

refMeanX0=mean(X0,'omitnan');
X0(del)=[];
Y0(del)=[];
Nsample=length(X0);%<------------------%

% Binning
edgesX=edgesX.*refMeanX0;
values=(edgesX(1:end-1)+edgesX(2:end))/2;
bin1=discretize(X0,edgesX,values);
scatY=zeros(length(values),1);
SEM=zeros(length(values),1);
delBin=[];
Ntemp=0;%<---------------------------------
for i=1:length(values)
    pool=Y0(bin1==values(i));
    scatY(i)=mean(pool);%<--------------------mean
    SEM(i)=std(pool)./length(pool);
    if length(pool)<Nsample*thrPct
        delBin=[delBin;i];
    else%<---------------------------------
        Ntemp=Ntemp+length(pool);%<---------------------------------
    end
end
values(delBin)=[];
scatY(delBin)=[];
SEM(delBin)=[];

% Drawing
F=figure('Position',[100,100,350,410],'Visible','off');
ax=gca;

hold on

% -density scatter:
DensScat(X0,Y0,'TargetAxes',ax);%,'ColorMap'
colormap(flipud(gray));
caxis([-0.01 0.3]);
cb=colorbar;
set(cb,'Limits',[0 0.15]);
set(cb,'YTick',[])



% -regression_bin:
X = [ones(size(values')) values'];
[b,~] = regress(scatY,X) ;
t=[0.75:0.5:1.3]*refMeanX0;
yq=b(2)*t+b(1);
plot(t,yq,'-k','LineWidth',2.5);


% -regression_all
X_all=[ones(size(X0)) X0];
[b_all,~] = regress(Y0,X_all) ;


% -bin values:
% scatter(values,scatY,15,'filled','MarkerFaceAlpha',0.9,'MarkerFaceColor',[0.8500, 0.3250, 0.0980]);
errorbar(values,scatY,SEM,'Color',[0.8500, 0.3250, 0.0980],'Marker','o','MarkerFaceColor',[0.8500, 0.3250, 0.0980],'MarkerSize',5,'CapSize',0,'LineStyle','none');
% axis([0 2 0 2])
ax=gca;
ax.FontSize=12;
% xticks([0 1 2]);
% yticks([0 1 2]);
hold off
% correlation index:
[rho,pval] = corr(X0,Y0,'type','Pearson');


% annotating
slope=b(2);slope2=b_all(2);
N=length(X0);

set(ax, 'Position',[0.15, 0.02, 0.65, 0.8]);

dim_vis = [0.15 0.74 0.75 0.2];
dim_left = [0.29 0.74 0.33 0.2];
dim_right = [0.62 0.74 0.28 0.2];
dim_med = [0.15 0.74 0.14 0.2];

% str = {[lgdLabel,'Bin = ',num2str(slope,'%0.2f')],...
%     [lgdLabel,'All = ',num2str(slope2,'%0.2f')],...
%     ['N = ',num2str(N)],...
%     ['\rho = ',num2str(rho,'%0.2f')]};

str_left = {[lgdLabel,'Bin = ',num2str(slope,'%0.2f')],...
    [lgdLabel,'All = ',num2str(slope2,'%0.2f')]};
str_right = {['N = ',num2str(N)],...
    ['\rho = ',num2str(rho,'%0.2f')]};

annotation('textbox',dim_vis);
annotation('textbox',dim_left,'String',str_left,'Interpreter','tex','FontSize',12,'EdgeColor','none');
annotation('textbox',dim_right,'String',str_right,'Interpreter','tex','FontSize',12,'EdgeColor','none');
annotation('textbox',dim_med,'String',med,'Interpreter','tex','FontSize',18,'FontName','Arial','FontWeight','bold','EdgeColor','none');
xlabel(lab.xlab);
ylabel(lab.ylab);
%saving
saveas(F,svPath);
close('all');
end