function [uit,newCanvLineIdx]=ImgFLN_createNewLine(img,pos)
uit=[];
J=cat(3,img,img,img);
stPoint=pos(1:2);
endPoint=pos(3:4);
canvas=insertShape(J,'Line',[stPoint,endPoint],'LineWidth',1,'Color','red');

canvLine=uint8(zeros(size(img)));
canvLine=insertShape(canvLine,'Line',[stPoint,endPoint],'LineWidth',5,'Color','red');

canvas(stPoint(1,2)-4:stPoint(1,2)+4,stPoint(1,1)-4:stPoint(1,1)+4,1)=0;
canvas(endPoint(1,2)-4:endPoint(1,2)+4,endPoint(1,1)-4:endPoint(1,1)+4,1)=0;
canvLine(stPoint(1,2)-4:stPoint(1,2)+4,stPoint(1,1)-4:stPoint(1,1)+4,1)=0;
canvLine(endPoint(1,2)-4:endPoint(1,2)+4,endPoint(1,1)-4:endPoint(1,1)+4,1)=0;

L=regionprops(logical(canvLine(:,:,1)),'PixelIdxList');
if length(L)>1
    error('debugging-more than one region detected during line.');
else
    newCanvLineIdx={L(1).PixelIdxList};
end
uit=canvas(:,:,1);
end