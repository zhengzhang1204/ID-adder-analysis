cllc;

%% batch mode:
% load('G:\mChrRpath.mat');
% for iExp=1:size(Rpath,1)
%     Exp=['G:\',Rpath{iExp,1}];
%     int=Rpath{iExp,3};
%     coreQC(Exp,int,iExp);
% end

%% single mode:
coreQC('G:\RdnaA-seqAR-atc 3to10-M30 medium 20260130',3,1);

function []=coreQC(Exp,int,iexp)
F=dir([Exp,'\M_res\p*_lkgM.mat']);
GRIDTclt=cell(length(F),1);%GR,IDT,collector. Format: Pos, SChn, {Idx-Frame lineage}, {[idt, gr, gof, CycInitTime]}
for i=1:length(F)%PARFOR each position, read its Mres.mat
	M=load([F(i).folder,'\',F(i).name],'PHparts', 'Infcellv2', 'PHMap');
    % disp(['Processing P',F(i).name(2:3),'(',num2str(i,'%0.2i'),')...']);
    data=[];
    for iSchN=1:length(M.PHparts)%for each side channel
        if iSchN>size(M.Infcellv2,2);continue;end
        if isempty(M.Infcellv2{1,iSchN});continue;end
        curPHparts=M.PHparts{iSchN};
        if isempty(curPHparts);continue;end
        curPHMap=M.PHMap{iSchN};
        curInfCell=M.Infcellv2(:,iSchN);
        
        %get indices of complete cycles:
        idxEntry=find(~isnan(curPHMap(:,2)) & ~isnan(curPHMap(:,3)));
        idxEntry=idxEntry(2:end);
        
        %'IDT','GR','CycInitTime','Gof','BirthLength','DivAdder_by_dividingCell','Pos','SchN','PHpartsID', 'DivAdder_by_twoDaughters'
        curSchNData=zeros(length(idxEntry),10);
        try
            for j=1:length(idxEntry)%for each PH parts
                daugIDs=curPHMap(idxEntry(j),2:3);
                daugIDs=[M.PHparts{iSchN}{daugIDs(1)}(1,:);M.PHparts{iSchN}{daugIDs(2)}(1,:)];
                Ld_2daug=curInfCell{daugIDs(1,1)}{2}(daugIDs(1,2),2)+curInfCell{daugIDs(2,1)}{2}(daugIDs(2,2),2);
                curSchNData(j,1)=size(curPHparts{idxEntry(j)},1)*int;%'IDT'
                curSchNData(j,3)=curPHparts{idxEntry(j)}(1,1)*int;%'CycInitTime'
                Ls=cellfun(@(x) curInfCell{x(1)}{2}(x(2),2), num2cell(curPHparts{idxEntry(j)},2));
                x=(1:1:length(Ls))'*int/60;
                [f,gof] = fit(x,Ls,'exp1');
                curSchNData(j,2)=f.b;%'GR'
                curSchNData(j,4)=gof.rsquare;%'Gof'
                curSchNData(j,5)=Ls(1);%'BirthLength'
                curSchNData(j,6)=Ld_2daug-Ls(1);%adder where Ld is Lb1 + Lb2
                curSchNData(j,7)=i;
                curSchNData(j,8)=iSchN;
                curSchNData(j,9)=idxEntry(j); 
                curSchNData(j,10)=Ls(end)-Ls(1);%adder where Ld is dividing cell
            end
            del=curSchNData(:,6)<=0;
            curSchNData(del,:)=[];
            data=cat(1,data,curSchNData);
        catch
            disp(['E',num2str(iexp,'%0.2i'),'P',num2str(i,'%0.2i'),'S',num2str(iSchN,'%0.2i'),'F',num2str(daugIDs(1),'%0.3i')]);
        end

    end
    GRIDTclt{i}=data;
end

%% select complete parts only: 
GRIDTclt=cell2mat(GRIDTclt);
tab=array2table(GRIDTclt,'VariableNames',{'IDT','GR','CycInitTime','Gof','BirthLength','DivAdderByCombDaug','Pos','SchN','PHpartsID','DivAdderBySd-Sb'});
save([Exp,'\M_res\Tab_GRIDTADDERRv3.mat'],'tab');
end

%% nested functions
function [spuit,abuit]=getSpeAbsER(L,int)
spuit=zeros(length(L)-1,1);
abuit=zeros(length(L)-1,1);
dt=int/60;
for i=1:length(L)-1
    spuit(i)=(log(L(i+1))-log(L(i)))/dt;
    abuit(i)=(L(i+1)-L(i))/dt;
end
spuit={spuit};
abuit={abuit};
end
function []=plot3in1QCfig(tab,edgesX,int,ti)
F=figure('Position',[100,100,360,800]);
subplot(3,1,1)
X=(tab.BirthLength)./mean(tab.BirthLength);
Y=(tab.DivAdder)./mean(tab.DivAdder);
values=(edgesX(1:end-1)+edgesX(2:end))/2;
bin1=discretize(X,edgesX,values);
%F1=figure;
hold on
scatter(X,Y,10,'filled','MarkerFaceAlpha',0.1,'MarkerFaceColor','b');
b1=boxchart(bin1,Y,'MarkerStyle','none','BoxWidth',0.07);
hold off
xlim([0 2]);
ylim([0 2]);
xt=xlim;
xt=xt(1)+(xt(2)-xt(1))*0.75;
yt=ylim;
yt=yt(1)+(yt(2)-yt(1))*0.05;
xlabel('Rescaled birth length','FontSize',15);
ylabel('Rescaled added length','FontSize',15);
text(xt,yt,['n=',num2str(length(X))],'FontSize', 15);
X1 = [ones(size(X)) X];
[b,bint] = regress(Y,X1) ;
text(0.01,1.8,['k_{dd}: ',num2str(b(2),'%0.2f'),' [',num2str(bint(2,1),'%0.2f'),',',num2str(bint(2,2),'%0.2f'),']']);
title(ti,'FontSize',15);

%F1.Position=[100,100,300,300];

%plot GR vs Time:
subplot(3,1,2)
plot(tab.CycInitTime,tab.GR,'.');
text(0.05,0.9,['Elongation rate: ',num2str(mean(tab.GR),'%0.2f'),char(177),...
    num2str(std(tab.GR),'%0.2f'),' h^{-1}'],'Units','normalized');
xlabel('Cyc. Init. Time (min)','FontSize',15);
ylabel('Growth rate (h^{-1})','FontSize',15);
grid on;
grid minor;
% ylim([0 2]);
%plot Lb vs Time
subplot(3,1,3)
% plot(tab.CycInitTime*int,tab.BirthLength,'.');
% xlabel('Cyc. Init. Time (min)','FontSize',15);
% ylabel('Birth Length (px)','FontSize',15);
plot(tab.CycInitTime,60*log(2)./tab.IDT,'.');
text(0.05,0.9,['Doubling rate: ',num2str(mean(60*log(2)./tab.IDT),'%0.2f'),char(177),...
    num2str(std(60*log(2)./tab.IDT),'%0.2f'),' h^{-1}'],'Units','normalized');
xlabel('Cyc. Init. Time (min)','FontSize',15);
ylabel('Doubling rate (h^{-1})','FontSize',15);
grid on;
grid minor;
%ylim([0 2]);
end
function []=plotDivAdderOnly(tab,edgesX,int)
X=(tab.BirthLength)./mean(tab.BirthLength);
Y=(tab.DivAdder)./mean(tab.DivAdder);
values=(edgesX(1:end-1)+edgesX(2:end))/2;
bin1=discretize(X,edgesX,values);
F1=figure;
hold on
scatter(X,Y,10,'filled','MarkerFaceAlpha',0.1,'MarkerFaceColor','b');
b1=boxchart(bin1,Y,'MarkerStyle','none','BoxWidth',0.07);
hold off
xlim([0 2]);
ylim([0 2]);
xt=xlim;
xt=xt(1)+(xt(2)-xt(1))*0.75;
yt=ylim;
yt=yt(1)+(yt(2)-yt(1))*0.05;
xlabel('Rescaled birth length','FontSize',15);
ylabel('Rescaled added length','FontSize',15);
text(xt,yt,['n=',num2str(length(X))],'FontSize', 15);
F1.Position=[100,100,300,300];
end
function [uit]=selCompleteParts(in)
uit=cell(size(in,1),4);
uit(:,1:2)=in(:,1:2);
for i=1:size(in,1)
    col3=in{i,3};
    col4=in{i,4};
    del=[];
    for j=1:size(in{i,5},1)
        %select effective entries:
        if isnan(in{i,5}(j,2)) && isnan(in{i,5}(j,3))
            del=[del;j];
        end

        if isnan(in{i,5}(j,1))
            del=[del;j];
        end
    end
    col3(del)=[];
    col4(del)=[];
    uit{i,3}=col3;
    uit{i,4}=col4;
end
del2=cellfun(@(x) isempty(x),uit(:,3));
uit(del2,:)=[];
end
