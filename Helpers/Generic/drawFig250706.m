% Version 2.1: added Nmax so that scatter plots only randomly selected spots.
% Version 2: Screening samples with nSD. Remove bins that containing too less points. @230803
% Version 3: calculate CI95 for each k. remove sem error bar


function [UIT]=drawFig250706(edgesX,X0,Y0,nSD,svPath,lgdLabel,med,lab,thrPct,isRemOutlier)

%% preparing data:
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

%% binning data:
% Binning for calculating (all dots, for calculation)
binCentAll=(edgesX(1:end-1)+edgesX(2:end))/2;
binNumAll=discretize(xAll,edgesX,binCentAll);
binYmean=zeros(length(binCentAll),1);
binYmedian=zeros(length(binCentAll),1);
delBinAll=[];
Ntemp=0;%<---------------------------------
for i=1:length(binCentAll)
    pool=yAll(binNumAll==binCentAll(i));
    binYmean(i)=mean(pool);
    binYmedian(i)=median(pool);
    if length(pool)<nAll*thrPct
        delBinAll=[delBinAll;i];
    else%<---------------------------------
        Ntemp=Ntemp+length(pool);%<---------------------------------
    end
end
binCentAll(delBinAll)=[];
binYmean(delBinAll)=[];

%% Calculation (consider only dots for display)
mdl_all = fitlm(xAll,yAll);
k_all=mdl_all.Coefficients.Estimate(2);

mdl_bin_avg = fitlm(binCentAll,binYmean);
k_bin_avg=mdl_bin_avg.Coefficients.Estimate(2);

mdl_bin_med = fitlm(binCentAll,binYmedian);
k_bin_med=mdl_bin_med.Coefficients.Estimate(2);

C = cov(xAll,yAll);
k_cov=C(1,2);


UIT=[k_all,k_bin_avg,k_bin_med,k_cov];
end