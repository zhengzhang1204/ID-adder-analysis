% Version 2.1: added Nmax so that scatter plots only randomly selected spots.
% Version 2: Screening samples with nSD. Remove bins that containing too less points. @230803


function [k_bin_All,k_all_All,nAll,rho_all]=drawFigDensRescaleV2(edgesX,X0,Y0,nSD,svPath,lgdLabel,med,lab,thrPct,Nmax)
% Rescaling
xM=mean(X0,'omitnan');
yM=mean(Y0,'omitnan');
xLow=xM-nSD*std(X0);
xUp=xM+nSD*std(X0);
yLow=yM-nSD*std(Y0);
yUp=yM+nSD*std(Y0);

del=find(X0<xLow | X0>xUp | Y0<yLow | Y0>yUp | isnan(X0) | isnan(Y0));
xR=X0./xM;
yR=Y0./yM;
xR(del)=[];
yR(del)=[];

xAll=xR;
yAll=yR;

% choose 1050 samples:
if length(xR)>Nmax%debugging
    selIdx=randi(length(xR),[Nmax,1]);
    xR=xR(selIdx);%selected to plot and 
    yR=yR(selIdx);
end
nAll=length(xAll);%<------------------% have to think over again.

% Binning for plotting (selected dots)
values=(edgesX(1:end-1)+edgesX(2:end))/2;
binPlot=discretize(xR,edgesX,values);
scatYPlot=zeros(length(values),1);
SEM=zeros(length(values),1);
delBin=[];
Ntemp=0;%<---------------------------------
for i=1:length(values)
    pool=yR(binPlot==values(i));
    scatYPlot(i)=mean(pool);%<--------------------mean
    SEM(i)=std(pool)./length(pool);
    if length(pool)<nAll*thrPct
        delBin=[delBin;i];
    else%<---------------------------------
        Ntemp=Ntemp+length(pool);%<---------------------------------
    end
end
values(delBin)=[];
scatYPlot(delBin)=[];
SEM(delBin)=[];

% Binning for calculating (all dots)
valuesAll=(edgesX(1:end-1)+edgesX(2:end))/2;
binAll=discretize(xAll,edgesX,valuesAll);
scatYAll=zeros(length(valuesAll),1);
% SEM=zeros(length(values),1);
delBinAll=[];
Ntemp=0;%<---------------------------------
for i=1:length(valuesAll)
    pool=yAll(binAll==valuesAll(i));
    scatYAll(i)=mean(pool);%<--------------------mean
    %SEM(i)=std(pool)./length(pool);
    if length(pool)<nAll*thrPct
        delBinAll=[delBinAll;i];
    else%<---------------------------------
        Ntemp=Ntemp+length(pool);%<---------------------------------
    end
end
valuesAll(delBinAll)=[];
scatYAll(delBinAll)=[];


% Drawing
F=figure('Position',[100,100,350,410],'Visible','off');
ax=gca;
xlabel(lab.xlab);
ylabel(lab.ylab);
hold on

% -density scatter:
DensScat(xR,yR,'TargetAxes',ax);%,'ColorMap'
colormap(flipud(gray));
caxis([-0.01 0.3]);
cb=colorbar;
set(cb,'Limits',[0 0.15]);
set(cb,'YTick',[])

% -regression_bin:
X = [ones(size(values')) values'];
[b,~] = regress(scatYPlot,X) ;
t=0.5:0.5:1.5;
yq=b(2)*t+b(1);
plot(t,yq,'-k','LineWidth',2.5);

% -regression_all
X_all=[ones(size(xR)) xR];
[b_all,~] = regress(yR,X_all) ;


% Calculation for all dots:
X_bin_All = [ones(size(valuesAll')) valuesAll'];
[b_bin_All,~] = regress(scatYAll,X_bin_All);k_bin_All=b_bin_All(2);%all samples, k_bin

X_all_All=[ones(size(xAll)) xAll];
[b_all_All,~] = regress(yAll,X_all_All);k_all_All=b_all_All(2);%all samples, k_all

[rho_all,~] = corr(xAll,yAll,'type','Pearson');%all samples, rho

% -bin values:
errorbar(values,scatYPlot,SEM,'Color',[0.8500, 0.3250, 0.0980],'Marker','o','MarkerFaceColor',[0.8500, 0.3250, 0.0980],'MarkerSize',5,'CapSize',0,'LineStyle','none');
axis([0 2 0 2])
ax=gca;
ax.FontSize=12;
xticks([0 1 2]);
yticks([0 1 2]);
hold off
% correlation index:
% [rho,pval] = corr(X(:,2),scatY,'type','Pearson');
[rho,~] = corr(xR,yR,'type','Pearson');

% annotating
slope=b(2);slope2=b_all(2);
N=length(xR);
set(ax, 'Position',[0.15, 0.02, 0.65, 0.8]);

dim_vis = [0.15 0.74 0.75 0.2];
dim_left = [0.29 0.74 0.33 0.2];
dim_right = [0.62 0.74 0.28 0.2];
dim_med = [0.15 0.74 0.14 0.2];

str_left = {[lgdLabel,'Bin = ',num2str(slope,'%0.2f')],...
    [lgdLabel,'All = ',num2str(slope2,'%0.2f')]};
str_right = {['N = ',num2str(N)],...
    ['\rho = ',num2str(rho,'%0.2f')]};

annotation('textbox',dim_vis);
annotation('textbox',dim_left,'String',str_left,'Interpreter','tex','FontSize',12,'EdgeColor','none');
annotation('textbox',dim_right,'String',str_right,'Interpreter','tex','FontSize',12,'EdgeColor','none');
annotation('textbox',dim_med,'String',med,'Interpreter','tex','FontSize',18,'FontName','Arial','FontWeight','bold','EdgeColor','none');

%saving
if ~isempty(svPath)
    saveas(F,svPath);
end
close('all');
end