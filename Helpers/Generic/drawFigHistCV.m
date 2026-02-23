% Version 2.1: added Nmax so that scatter plots only randomly selected spots.
% Version 2: Screening samples with nSD. Remove bins that containing too less points. @230803


function [CV]=drawFigHistCV(X0,Y0,nSD,svPath,med,lab)
% Rescaling
xM=mean(X0,'omitnan');
yM=mean(Y0,'omitnan');
xLow=xM-nSD*std(X0);
xUp=xM+nSD*std(X0);
yLow=yM-nSD*std(Y0);
yUp=yM+nSD*std(Y0);

del=find(X0<xLow | X0>xUp | Y0<yLow | Y0>yUp | isnan(X0) | isnan(Y0));
% xR=X0./xM;
% yR=Y0./yM;
X0(del)=[];
Y0(del)=[];


% Drawing
F=figure('Position',[100,100,350,410],'Visible','off');
xlabel(lab.xlab);
ylabel(lab.ylab);
hold on

% -density scatter:
histogram(Y0)
CV=std(Y0)./mean(Y0,'all','omitnan');


% Calculation for all dots:

xL=xlim;yL=ylim;
text((xL(2)-xL(1))*0.7+xL(1),(yL(2)-yL(1))*0.9+yL(1),['CV = ',num2str(CV,'%0.2f')],'FontSize',12);


%saving
saveas(F,svPath);
close('all');
end