% Version 2.1: added Nmax so that scatter plots only randomly selected spots.
% Version 2: Screening samples with nSD. Remove bins that containing too less points. @230803
% Version 3: calculate CI95 for each k. remove sem error bar
% Version 4 (250727): No outliers, gray scatters, red bin mean, fit bin mean black dash line

function [UIT]=drawFig250727(edgesX,X0,Y0,nSD,svPath,thrPct,~,~,isRemOutlier)

%% preparing data:
%1. keep or remove outliers
xM=mean(X0,'omitnan');
yM=mean(Y0,'omitnan');
if isRemOutlier% Removing outliers (nSD>3)
    xLow=xM-nSD*std(X0,'omitnan');
    xUp=xM+nSD*std(X0,'omitnan');
    yLow=yM-nSD*std(Y0,'omitnan');
    yUp=yM+nSD*std(Y0,'omitnan');
    del=find(X0<xLow | X0>xUp | Y0<yLow | Y0>yUp | isnan(X0) | isnan(Y0));
    xR=X0./xM;
    yR=Y0./yM;
    xR(del)=[];
    yR(del)=[];
    xAll=xR;
    yAll=yR;
    nAll=length(xAll);
else
    del=find(isnan(X0) | isnan(Y0));
    xR=X0./xM;
    yR=Y0./yM;
    xR(del)=[];
    yR(del)=[];
    xAll=xR;
    yAll=yR;
    nAll=length(xAll);
end

%2. Binning
% Binning for calculating (all dots, no outliers, for calculation)
binCentAll=(edgesX(1:end-1)+edgesX(2:end))/2;
binNumAll=discretize(xAll,edgesX,binCentAll);
binYmean=zeros(length(binCentAll),1);
% binYmedian=zeros(length(binCentAll),1);
binSEM=zeros(length(binCentAll),1);
delBinAll=[];

for i=1:length(binCentAll)
    pool=yAll(binNumAll==binCentAll(i));
    binYmean(i)=mean(pool);
    % binYmedian(i)=median(pool);
    binSEM(i)=std(pool)/sqrt(length(pool));   
    if length(pool)<nAll*thrPct
        delBinAll=[delBinAll;i];
    end
end
binCentAll(delBinAll)=[];    % bin center (x-axis)
binYmean(delBinAll)=[];     % bin mean values (y-axis)
% binYmedian(delBinAll)=[];     % bin median values (y-axis)
binSEM(delBinAll)=[];   % bin SEM values (y-error-bar)

mdl_bin_avg = fitlm(binCentAll,binYmean);
k_bin_avg=mdl_bin_avg.Coefficients.Estimate(2);

% mdl_bin_med = fitlm(binCentAll,binYmedian);
% k_bin_med=mdl_bin_med.Coefficients.Estimate(2);

%% Calculation
F=figure('Position',[100,100,247.6,247.6],'Visible','on');%can not change this
ax=gca;
hold on

%%%Option Three: contour
% scatter(xR.*scaleX/2,yR.*scaleY/2,4,'k','filled','MarkerFaceAlpha',0.2);
scatter(xR,yR,4,'k','filled','MarkerFaceAlpha',0.2);
% plot(binCentAll,binYmean,'.r');
% errorbar(binCentAll,binYmean,binSEM)
tx=[binCentAll(1)-0.1,binCentAll(end)+0.1];
ty=feval(mdl_bin_avg,tx);
plot(tx,ty,'-k','LineWidth',1.5);
errorbar(binCentAll,binYmean,binSEM,"o","MarkerSize",5,"MarkerEdgeColor","none","MarkerFaceColor","red",'CapSize',3,"Color","r","LineStyle","none");
% plot([1.5 1.5],[0.5 1.5],'-k','LineWidth',1);
% plot([0.5 1.5],[1.5 1.5],'-k','LineWidth',1);
xlim([0.5 1.5]);ylim([0.5 1.5]);
xline(ax,ax.XLim(2),'linewidth',1)
yline(ax,ax.YLim(2),'linewidth',1)
hold off

set(ax,'linewidth',1)
xticks([0.5 1 1.5]);
yticks([0.5 1 1.5]);
ax.TickDir='out';
ax.LineWidth=1;
ax.XTickLabel=[];
ax.YTickLabel=[];
%saving
if ~isempty(svPath)
    saveas(F,svPath);
end
close('all');

%% OUTPUT SIMPLIFY:
UIT=[k_all_All,CI95_allU,CI95_allL,k_bin_All,CI95_binU,CI95_binL,nAll,rho_all];
end