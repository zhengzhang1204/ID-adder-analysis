% Version 2.1: added Nmax so that scatter plots only randomly selected spots.
% Version 2: Screening samples with nSD. Remove bins that containing too less points. @230803
% Version 3: calculate CI95 for each k. remove sem error bar


function [UIT]=drawFig250401(edgesX,X0,Y0,nSD,svPath,lgdLabel,med,lab,thrPct,scaleX,scaleY)
%% preparing data:
% Rescaling
xM=mean(X0,'omitnan');
yM=mean(Y0,'omitnan');
xLow=xM-nSD*std(X0,'omitnan');
xUp=xM+nSD*std(X0,'omitnan');
yLow=yM-nSD*std(Y0,'omitnan');
yUp=yM+nSD*std(Y0,'omitnan');

% Removing outliers (nSD>3)
del=find(X0<xLow | X0>xUp | Y0<yLow | Y0>yUp | isnan(X0) | isnan(Y0));
xR=X0./xM;
yR=Y0./yM;
xR(del)=[];
yR(del)=[];
xAll=xR;
yAll=yR;
nAll=length(xAll);

%% binning data:
% Binning for calculating (all dots, no outliers, for calculation)
valuesAll=(edgesX(1:end-1)+edgesX(2:end))/2;
binAll=discretize(xAll,edgesX,valuesAll);
scatYAll=zeros(length(valuesAll),1);
delBinAll=[];
Ntemp=0;%<---------------------------------
for i=1:length(valuesAll)
    pool=yAll(binAll==valuesAll(i));
    scatYAll(i)=mean(pool);
    if length(pool)<nAll*thrPct
        delBinAll=[delBinAll;i];
    else
        Ntemp=Ntemp+length(pool);
    end
end
valuesAll(delBinAll)=[];    % bin center
scatYAll(delBinAll)=[];     % bin values

%% Calculation
F=figure('Position',[100,100,247.6,247.6],'Visible','on');%can not change this
ax=gca;
hold on

%%%Option Three: contour
% ZZdensityContourF(xAll,yAll,colormap(flipud(gray)),ax,scaleX,scaleY);
scatter(xR.*scaleX/2,yR.*scaleY/2,4,'filled','MarkerFaceAlpha',0.2);


% -regression_bin (slope for binned data)
X = [ones(size(valuesAll')) valuesAll'];
[b,~] = regress(scatYAll,X) ;

% -regression_all
X_all=[ones(size(xAll)) xAll*scaleX/2];
[b_all,bint] = regress(yAll*scaleY/2,X_all) ;

t=(0.5:0.5:1.5)*scaleX/2;
yq=b_all(2)*t+b_all(1);
% ylow = bint(1,1)+bint(2,1)*t;%<---original CI95 patching 
% yupp = bint(1,2)+bint(2,2)*t;%<---original CI95 patching 
passPoint=[scaleX/2,b_all(2)*scaleX/2+b_all(1)];
ylow = bint(2,1)*t+passPoint(2)-bint(2,1)*passPoint(1);%<---SJ-type CI95 patching 
yupp = bint(2,2)*t+passPoint(2)-bint(2,2)*passPoint(1);%<---SJ-type CI95 patching 

plot(t,yq,'-k','LineWidth',1);
acp=patch([t fliplr(t)], [ylow fliplr(yupp)], 'k','EdgeColor','white');
% alpha(0.3);
set(acp,'EdgeColor','none','FaceAlpha',0.3)

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



%% Draw
%errorbar(values,scatYPlot,SEM,'Color',[0.8500, 0.3250, 0.0980],'Marker','o','MarkerFaceColor',[0.8500, 0.3250, 0.0980],'MarkerSize',5,'CapSize',0,'LineStyle','none');
%scatter(ax,valuesAll*scaleX/2,scatYAll*scaleY/2,15,[0.8500, 0.3250, 0.0980],'filled','Marker','o','MarkerFaceColor',[0.8500, 0.3250, 0.0980]);
axis([0 scaleX 0 scaleY])
ax=gca;
ax.FontSize=12;
xticks([0 scaleX/2 scaleX]);
yticks([0 scaleY/2 scaleY]);
ax.XTickLabel=[];
ax.YTickLabel=[];
ax.TickDir='out';
ax.LineWidth=1;
box on
hold off

%saving
if ~isempty(svPath)
    saveas(F,svPath);
end
close('all');

%% OUTPUT SIMPLIFY:
UIT=[k_all_All,CI95_allU,CI95_allL,k_bin_All,CI95_binU,CI95_binL,nAll,rho_all];
end