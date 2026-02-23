function [type]=DL_validateFileFormatv2(EXP,~)
F=dir(EXP);
i=3;isFov=false;
while i<13 && ~isFov
    if strcmp(F(i).name(1:3),'fov') || strcmp(F(i).name(1:3),'fie')
        isFov=true;
    end
    i=i+1;
end
i=i-1;
if isFov%means in Jin format
    G=dir([EXP,'\',F(i).name]);
    H=dir([EXP,'\',F(i).name,'\',G(3).name]);
    if strcmp(H(3).name,'t0.tiff')
        type = 3;% CP format on Huangscope
    else
        type = 2;% Jin format on Jinscope
    end
else%means in NIS format
    type = 1; return;
end

