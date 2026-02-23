%Version v2.0 @ 23.03.03    removed Stage 4.
clear;
clc;

%% User-specified region:
% specify folders and channels for image sorting:
Target='G:\20260124 RdnaA_M30_atc1_dnaNyPet2';    %Folder storing all processed stuffs. Output folder.
COE=7.5;minArea4L92=35;%<-----Use MoMa_DL_kil_Stage15_test_LLR_params.m to determine these two parameters.
isPicOut=false;      %Always false. True to export AsPicOut images for debugging.

F=dir([Target,'\M_res\p*_resM.mat']);
for i=14%length(F)%for each position, read its Mres.mat <------------------change this back to 1
    pos=str2double(F(i).name(2:3));
    YcropFolder=[Target,'\Y_crop\p',num2str(pos,'%0.2i'),'\p',num2str(pos,'%0.2i')];
    load([Target,'\M_res\p',num2str(pos,'%0.2i'),'_resM.mat'],'PosMask');
    load([Target,'\Y_res\p',num2str(pos,'%0.2i'),'_resY.mat'],'Est','Thr','TSstack','Cand');
    if ~exist("Cand",'var')
        Cand=[];
    end
    [Cand,Est,PHMaskTS]=DL_stage16v2(Thr,Est,PosMask,TSstack,minArea4L92,COE,isPicOut,Cand);%,CropInfo,[256,ImageSize],YcropFolder);%Fred's codes;
    save([Target,'\Y_res\p',num2str(pos,'%0.2i'),'_resY.mat'],'Est','Cand','PHMaskTS','-append');
    clear('Est','Thr','Cand','CropInfo','PosMaskN','PosMask');
    disp(['Pos ',num2str(pos,'%0.2i'), '(',num2str(i),') is done!']);
end