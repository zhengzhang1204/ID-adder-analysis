function [uit,newCanvLineIdx]=ImgFLN_createNewLineBatch(img,route,sptCent)
uit=[];
canvas=cat(3,img,img,img);
newCanvLineIdx={};
for i=1:size(route,1)
    stPoint=sptCent(route(i,1),:);
    endPoint=sptCent(route(i,2),:);
    canvas=insertShape(canvas,'Line',[stPoint,endPoint],'LineWidth',1,'Color','red');
    canvLine=uint8(zeros(size(img)));
    canvLine=insertShape(canvLine,'Line',[stPoint,endPoint],'LineWidth',5,'Color','red');
    canvas(stPoint(1,2)-4:stPoint(1,2)+4,stPoint(1,1)-4:stPoint(1,1)+4,1)=0;
    canvas(endPoint(1,2)-4:endPoint(1,2)+4,endPoint(1,1)-4:endPoint(1,1)+4,1)=0;
    canvLine(stPoint(1,2)-4:stPoint(1,2)+4,stPoint(1,1)-4:stPoint(1,1)+4,1)=0;
    canvLine(endPoint(1,2)-4:endPoint(1,2)+4,endPoint(1,1)-4:endPoint(1,1)+4,1)=0;
    L=regionprops(logical(canvLine(:,:,1)),'PixelIdxList');
    newCanvLineIdx=[newCanvLineIdx;{L(1).PixelIdxList}];
end
uit=canvas(:,:,1);
end