%output simp: 
clear;
clc;
close('all');

%% User specified region
load('G:\mChrRpath.mat');
sel=[1:1:3,6,8,10,12:15];
Rpath=Rpath(sel,:);
path=cellfun(@(x,y) [x,y], repmat({'G:\'},size(Rpath,1),1),Rpath(:,1),'UniformOutput',false);
exportDir='D:\mChr_RescaleL260224\';%'D:\mChr_Rescale\';'D:\yPet_Rescale\'

%% codes start here:
for i=1:length(path)%<-----CHANGE THIS BACK TO 1
    M=load([path{i},'\M_res\Tab_GRIDTADDERRv3simp.mat']);
    Y=load([path{i},'\resV2simp.mat']);
    int=Y.int;
    pth1=[exportDir,'\',Rpath{i,2},'_data.xlsx'];
    pth2=[exportDir,'\',Rpath{i,2},'_data.mat'];
    [tab_DD]=ZZ_writeXLSX(M.GRIDTData,pth1,'DD');
    [tab_II]=ZZ_writeXLSX(Y.outAData,pth1,'II');
    [tab_ID]=ZZ_writeXLSX(Y.outGHIData,pth1,'ID');
    DD_GR=M.tab.GR;
    save(pth2,'tab_ID','tab_II','tab_DD','int','DD_GR');
end


function [tab]=ZZ_writeXLSX(in,pth,sheetName)
in1=in(:,1:3);
if size(in,2)==4
    in2=in(:,4);
    tt={'current_ID','daughter1_ID','daughter2_ID','cur_data'};
else
    in2=cellfun(@(x,y) [x;nan;y],in(:,4),in(:,5),'UniformOutput',false);
    tt={'current_ID','daughter1_ID','daughter2_ID','cur-to-d1_data','cur-to-d2_data'};
end
maxH=max(cellfun(@length, in2),[],'all');
uit=[in1,num2cell(cell2mat(cellfun(@(x) [x',nan(1,maxH-length(x)+1)], in2,'UniformOutput',false)))];
writecell(uit,pth,'Sheet', sheetName);
tab=cell2table(in,"VariableNames",tt);
end