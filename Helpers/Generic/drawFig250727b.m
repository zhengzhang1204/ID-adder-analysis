% Version 2.1: added Nmax so that scatter plots only randomly selected spots.
% Version 2: Screening samples with nSD. Remove bins that containing too less points. @230803
% Version 3: calculate CI95 for each k. remove sem error bar
% Version 4 (250727b): No outliers, gray scatters, red bin mean, fit bin mean black dash line

function [UIT]=drawFig250727b(edgesX,X0,Y0,nSD,svPath,thrPct,isEg1K)
%% preparing data:
%1. keep or remove outliers
xM=mean(X0,'omitnan');
yM=mean(Y0,'omitnan');
xR=X0./xM;
yR=Y0./yM;

xLow=xM-nSD*std(X0,'omitnan');
xUp=xM+nSD*std(X0,'omitnan');
yLow=yM-nSD*std(Y0,'omitnan');
yUp=yM+nSD*std(Y0,'omitnan');
del_OL=X0<xLow | X0>xUp | Y0<yLow | Y0>yUp | isnan(X0) | isnan(Y0);
x_noOL=xR(~del_OL);
y_noOL=yR(~del_OL);
n_noOL=length(y_noOL);

del_Basic=isnan(X0) | isnan(Y0);
x_withOL=xR(~del_Basic);
y_withOL=yR(~del_Basic);
n_withOL=length(y_withOL);

[slps_bwo,~,~]=dealBinning(edgesX,x_withOL,y_withOL,n_withOL,thrPct);
[slps_bno,plt_bno1,plt_bno2]=dealBinning(edgesX,x_noOL,y_noOL,n_noOL,thrPct);

[slps_awo,~]=dealAll(x_withOL,y_withOL);
[slps_ano,plt_ano]=dealAll(x_noOL,y_noOL);


cov_wo=dealCov(x_withOL,y_withOL);
cov_no=dealCov(x_noOL,y_noOL);
%uit (1x22 <double>): 
% Col.1-6: [ano, ano_low, ano_up, awo, awo_low, awo_up]
% Col.7-12: [bavgno, bavgno_low, bavgno_up, bmedno, bmedno_low, bmedno_up];
% Col.13-18: [bavgwo, bavgwo_low, bavgwo_up, bmedwo, bmedwo_low, bmedwo_up];
% Col.19-20: [covno, covwo];
% Col.21,22: [nno, nwo];
UIT=[slps_ano,slps_awo, slps_bno, slps_bwo, cov_no, cov_wo, n_noOL, n_withOL];

%% Draw:
if length(x_noOL)>=1000 && isEg1K
    rngIdx=randperm(length(x_noOL),1000);
    x_noOL=x_noOL(rngIdx);y_noOL=y_noOL(rngIdx);
end
F=figure('Position',[100,100,247,237],'Visible','off');%can not change this
ax1 = axes;
hold on
scatter(x_noOL,y_noOL,4,[0.5 0.5 0.5],'filled','MarkerFaceAlpha',0.4);
plot(plt_ano(:,1),plt_ano(:,2),'-k','LineWidth',1.5);
errorbar(plt_bno2(:,1),plt_bno2(:,2),plt_bno2(:,3),"o","MarkerSize",3,"MarkerEdgeColor","none","MarkerFaceColor","red",'CapSize',3,"Color","r","LineStyle","none");

xlim([0 2]);ylim([0 2]);
yL=ylim;xL=xlim;
% plot(xL(1),yL(1),'.k');

% AD's frame:
hold off
box(ax1, 'off');
ax1.XAxisLocation = 'bottom';
ax1.YAxisLocation = 'left';
ax1.TickDir = 'out';
xticks([0 1 2]);
yticks([0 1 2]);
ax1.TickDir='out';
ax1.XTickLabel=[];
ax1.YTickLabel=[];
ax2 = axes('Position', ax1.Position, ...
           'Color', 'none', ...
           'XAxisLocation', 'top', ...
           'YAxisLocation', 'right', ...
           'XTick', [], ...
           'YTick', [], ...
           'Box', 'on');
ax2.XTickLabel=[];
ax2.YTickLabel=[];
linkaxes([ax1 ax2]);

set(ax1, 'color', 'none');
set(ax1,'XColor', 'k','YColor','k')
set(ax2, 'color', 'none');
set(ax2,'XColor', 'k','YColor','k')

if ~isempty(svPath)
    saveas(F,svPath);
end
close('all');

end

%% nested functions:
function [uit1,uit2,uit3]=dealBinning(edgesX,xAll,yAll,nAll,thrPct)
%2. Binning
% Binning for calculating (all dots, no outliers, for calculation)
binCentAll=(edgesX(1:end-1)+edgesX(2:end))/2;
binNumAll=discretize(xAll,edgesX,binCentAll);
binYmean=zeros(length(binCentAll),1);
binYmedian=zeros(length(binCentAll),1);
binSEM=zeros(length(binCentAll),1);
delBinAll=[];

for i=1:length(binCentAll)
    pool=yAll(binNumAll==binCentAll(i));
    binYmean(i)=mean(pool);
    binYmedian(i)=median(pool);
    binSEM(i)=std(pool)/sqrt(length(pool));   
    if length(pool)<nAll*thrPct
        delBinAll=[delBinAll;i];
    end
end
binCentAll(delBinAll)=[];    % bin center (x-axis)
binYmean(delBinAll)=[];     % bin mean values (y-axis)
binYmedian(delBinAll)=[];     % bin median values (y-axis)
binSEM(delBinAll)=[];   % bin SEM values (y-error-bar)

mdl_bin_avg = fitlm(binCentAll,binYmean);
k_bin_avg=mdl_bin_avg.Coefficients.Estimate(2);

CI = coefCI(mdl_bin_avg); % 95% confidence intervals
CI_bin_avg = CI(2,:); % CI for slope (second coefficient)

mdl_bin_med = fitlm(binCentAll,binYmedian);
k_bin_med=mdl_bin_med.Coefficients.Estimate(2);

CI = coefCI(mdl_bin_med); % 95% confidence intervals
CI_bin_med = CI(2,:); % CI for slope (second coefficient)

tx=[0.5,1.5];
ty=feval(mdl_bin_avg,tx);
uit1=[k_bin_avg, CI_bin_avg, k_bin_med, CI_bin_med];
uit2=[tx',ty'];
uit3=[binCentAll',binYmean,binSEM];
end

function [uit,uit2]=dealAll(X,Y)
mdl_all = fitlm(X,Y);
k_all=mdl_all.Coefficients.Estimate(2);
CI = coefCI(mdl_all); % 95% confidence intervals
CI_all = CI(2,:); % CI for slope (second coefficient)
uit=[k_all,CI_all];
tx=[0.5,1.5];
ty=feval(mdl_all,tx);
uit2=[tx',ty'];
end

function [k]=dealCov(X,Y)
C = cov(X,Y);
k=C(1,2);
end