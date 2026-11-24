%output simp: 
clear;
clc;
close('all');

%% User specified region
% Option 1: process a group of experiments
load('G:\mChrRpath10.mat');
exportDir='D:\WT10Exp\';

% Option 2: process a specified experiment
% Rpath=[{'20260622_M30_atc1to5'},{'M30'},{3}];
% Rpath=[{'20260301_RdnaA-DnaNyPet-atc_1to5'},{'M30'},{3}];
% exportDir='D:\DnaAupshift260301\';%'D:\mChr_Rescale\';'D:\yPet_Rescale\'

%% codes start here:
path=cellfun(@(x,y) [x,y], repmat({'G:\'},size(Rpath,1),1),Rpath(:,1),'UniformOutput',false);
for i=1:length(path)%<-----CHANGE THIS BACK TO 1
    M=load([path{i},'\M_res\Tab_GRIDTADDERRv3simp.mat']);
    Y=load([path{i},'\resV2simp.mat']);
    int=Y.int;
    pth1=[exportDir,'\',Rpath{i,2},'_data.xlsx'];
    pth2=[exportDir,'\',Rpath{i,2},'_data.mat'];
    
    StEd_DD=[M.tab.CycInitTime,M.tab.IDT+M.tab.CycInitTime];
    [tab_DD]=ZZ_writeXLSX(M.GRIDTData,pth1,'DD',StEd_DD);
    StEd_II=Y.outA(:,22:24)*int;
    StEd_ID=Y.outGHI(:,16:17)*int;
    [tab_II]=ZZ_writeXLSX(Y.outAData,pth1,'II',StEd_II);
    [tab_ID]=ZZ_writeXLSX(Y.outGHIData,pth1,'ID',StEd_ID);
    DD_GR=M.tab.GR;
    save(pth2,'tab_ID','tab_II','tab_DD','int','DD_GR');
end

function [tab]=ZZ_writeXLSX(in,pth,sheetName,StEd)
in1=in(:,1:3);
if size(in,2)==4 && size(StEd,2)==2     %DD and ID
    in2=in(:,4);
    in=[in(:,1:3),num2cell(StEd),in(:,4)];
    tt={'current_ID','daughter1_ID','daughter2_ID','Starting time (min)','ending time (min)','cur_data'};
elseif size(in,2)==5                    %II
    in2=cellfun(@(x,y) [x;nan;y],in(:,4),in(:,5),'UniformOutput',false);
    in=[in(:,1:3),num2cell(StEd),in(:,4:5)];
    tt={'current_ID','daughter1_ID','daughter2_ID','starting time_1 (min)', 'ending time_1 (min)',...
        'ending time_2 (min)','cur-to-d1_data','cur-to-d2_data'};
end
maxH=max(cellfun(@length, in2),[],'all');
uit=[in1,num2cell(StEd),num2cell(cell2mat(cellfun(@(x) [x',nan(1,maxH-length(x)+1)], in2,'UniformOutput',false)))];
writecell(uit,pth,'Sheet', sheetName);
tab=cell2table(in,"VariableNames",tt);
end