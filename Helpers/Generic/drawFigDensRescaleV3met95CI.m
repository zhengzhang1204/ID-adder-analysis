% Version 2.1: added Nmax so that scatter plots only randomly selected spots.
% Version 2: Screening samples with nSD. Remove bins that containing too less points. @230803
% Version 3: calculate CI95 for each k. remove sem error bar


function [UIT]=drawFigDensRescaleV3met95CI(edgesX,X0,Y0,nSD,svPath,lgdLabel,med,lab,thrPct,Nmax)
%% preparing data:
% Rescaling
xM=mean(X0,'omitnan');
yM=mean(Y0,'omitnan');
xLow=xM-nSD*std(X0);
xUp=xM+nSD*std(X0);
yLow=yM-nSD*std(Y0);
yUp=yM+nSD*std(Y0);

% Removing outliers (nSD>3)
del=find(X0<xLow | X0>xUp | Y0<yLow | Y0>yUp | isnan(X0) | isnan(Y0));
xR=X0./xM;
yR=Y0./yM;
xR(del)=[];
yR(del)=[];
xAll=xR;
yAll=yR;
nAll=length(xAll);

%%%%%%%%%%%%%%%%%%%%%%%%
% DELETE THE FOLLOWING %
%%%%%%%%%%%%%%%%%%%%%%%%
%This is to make sure the randomly selected samples yield similar results as the entire dataset.
CHECK=1;cnt=1;thresh=0.09;
while CHECK>thresh

%%%%%%%%%%%%%%%%%%%%
% DELETE TILL HERE %
%%%%%%%%%%%%%%%%%%%%

if nAll>Nmax
    selIdx=randi(nAll,[Nmax,1]);
    xR=xAll(selIdx);%selected to plot 
    yR=yAll(selIdx);
end

%xR and yR are ramdomly selected samples for display.
%xAll and yAll are the entire dataset for calculation.

%% binning data:
% Binning for plotting (dots for display)
values=(edgesX(1:end-1)+edgesX(2:end))/2;
binPlot=discretize(xR,edgesX,values);
scatYPlot=zeros(length(values),1);
SEM=zeros(length(values),1);
delBin=[];
Ntemp=0;
for i=1:length(values)
    pool=yR(binPlot==values(i));
    scatYPlot(i)=mean(pool);%mean Y value for each bin
    SEM(i)=std(pool)./length(pool);
    if length(pool)<nAll*thrPct%if the bin has too few samples, remove the bin
        delBin=[delBin;i];
    else%
        Ntemp=Ntemp+length(pool);
    end
end
values(delBin)=[];
scatYPlot(delBin)=[];
SEM(delBin)=[];

% Binning for calculating (all dots, for calculation)
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

%% Calculation (consider only dots for display)
F=figure('Position',[100,100,350,410],'Visible','off');
ax=gca;
xlabel(lab.xlab);
ylabel(lab.ylab);
hold on

% -density scatter:
DensScat(xR,yR,'TargetAxes',ax);%display for selected dots (dots for display).
colormap(flipud(gray));
caxis([-0.01 0.3]);
cb=colorbar;
set(cb,'Limits',[0 0.15]);
set(cb,'YTick',[])

% -regression_bin (dots for display)
X = [ones(size(values')) values'];
[b,~] = regress(scatYPlot,X) ;


% -regression_all
X_all=[ones(size(xR)) xR];
[b_all,bint] = regress(yR,X_all) ;

t=0.5:0.5:1.5;
yq=b_all(2)*t+b_all(1);
% ylow = bint(1,1)+bint(2,1)*t;%<---original CI95 patching 
% yupp = bint(1,2)+bint(2,2)*t;%<---original CI95 patching 
passPoint=[1,b_all(2)+b_all(1)];
ylow = bint(2,1)*t+passPoint(2)-bint(2,1);%<---SJ-type CI95 patching 
yupp = bint(2,2)*t+passPoint(2)-bint(2,2);%<---SJ-type CI95 patching 

plot(t,yq,'-k','LineWidth',1);
patch([t fliplr(t)], [ylow fliplr(yupp)], 'k','EdgeColor','white')
alpha(0.3);


% Pearson correlation index:
[rho,~] = corr(xR,yR,'type','Pearson');

%% Calculation (consider all dots for calculation)
X_bin_All = [ones(size(valuesAll')) valuesAll'];%regression for binned values
[b_bin_All,bint_bin_All] = regress(scatYAll,X_bin_All);k_bin_All=b_bin_All(2);%all samples, k_bin
CI95_binL=bint_bin_All(2,1);CI95_binU=bint_bin_All(2,2);

X_all_All=[ones(size(xAll)) xAll];%regression for all dots
[b_all_All,bint_all_All] = regress(yAll,X_all_All);k_all_All=b_all_All(2);%all samples, k_all
CI95_allL=bint_all_All(2,1);CI95_allU=bint_all_All(2,2);

[rho_all,~] = corr(xAll,yAll,'type','Pearson');%all samples, rho

%%%%%%%%%%%%%%%%%%%%%%%%
% DELETE THE FOLLOWING %
%%%%%%%%%%%%%%%%%%%%%%%%
%This is to make sure the randomly selected samples yield similar results as the entire dataset.
CHECK=abs(k_all_All-b_all(2));
cnt=cnt+1;
if cnt>100
    thresh=thresh+0.01;
    cnt=1;
    disp(['New thresh: ',num2str(thresh)]);
end

if nAll<=Nmax
    CHECK=0;
end
end
%%%%%%%%%%%%%%%%%%%%
% DELETE TILL HERE %
%%%%%%%%%%%%%%%%%%%%

%% Draw (consider only dots for display)
errorbar(values,scatYPlot,SEM,'Color',[0.8500, 0.3250, 0.0980],'Marker','o','MarkerFaceColor',[0.8500, 0.3250, 0.0980],'MarkerSize',5,'CapSize',0,'LineStyle','none');
axis([0 2 0 2])
ax=gca;
ax.FontSize=12;
xticks([0 1 2]);
yticks([0 1 2]);
hold off

%% annotating (consider only dots for display)
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

%% OUTPUT SIMPLIFY:
UIT=[k_all_All,CI95_allU,CI95_allL,k_bin_All,CI95_binU,CI95_binL,nAll,rho_all];
end