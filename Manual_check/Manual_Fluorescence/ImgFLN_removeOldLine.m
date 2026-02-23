function [uit]=ImgFLN_removeOldLine(img,info)
J=cat(3,img,img,img);

%draw a dark line:
stPoint=info(3:4);
endPoint=info(5:6);
canvas=insertShape(J,'Line',[stPoint,endPoint],'LineWidth',1,'Color','black');
uit=canvas(:,:,1);
end