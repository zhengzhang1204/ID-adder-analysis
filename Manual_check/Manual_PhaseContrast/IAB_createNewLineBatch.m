function [canvasISO,canvasLNG]=IAB_createNewLineBatch(sz,route,isISO)
canvasISO=zeros(cat(2,sz,3));
canvasLNG=zeros(cat(2,sz,3));
for i=1:size(route,1)%for each PHpart
    curRoute=route{i,1};
    if size(curRoute,1)==1; continue; end
    pos=nan(size(curRoute,1)-1,4);
    for j=1:size(curRoute,1)-1
        pos(j,1:2)=fliplr(curRoute(j,:));
        pos(j,3:4)=fliplr(curRoute(j+1,:));
    end
    if isISO(i)
        canvasISO=insertShape(canvasISO,'Line',pos,'LineWidth',1,'Color','red');
    else
        canvasLNG=insertShape(canvasLNG,'Line',pos,'LineWidth',1,'Color','red');
    end
end
canvasISO=logical(canvasISO(:,:,1));
canvasLNG=logical(canvasLNG(:,:,1));
end