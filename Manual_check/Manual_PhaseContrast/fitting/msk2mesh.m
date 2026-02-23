% cllc;
% J=imread('E:\testExtCurve.tif');
% BW = imbinarize(J);
% % BW = imbinarize(imread('E:\testExtCurve.tif'));

function [cellMesh,skel,WL,cellContour]=msk2mesh(BW)

%1. initial contour from FFTsmoothed raw mask:
[ii,jj]=find(bwperim(BW),1,'first');
RRx=bwtraceboundary(BW,[ii,jj],'n',4,inf,'counterclockwise');
cellContour =fftSmooth(RRx,13);%<---param can be 10, was 20

%2. calculate MajorAxisLength:
stats = regionprops(BW,'MajorAxisLength');

%adapted contour2mesh from Morphometrics
[boneX,boneY]=contour2meshZZ(cellContour(:,2),cellContour(:,1),stats.MajorAxisLength);
if isempty(boneX)
    cellMesh=[];skel=[];WL=[];cellContour=[];
    return;
end

%debugging-----
% contour2mesh(cellContour(:,2),cellContour(:,1),'plot');
% figure;hold on;plot(cellContour(:,2),cellContour(:,1),'-');plot(boneX,boneY,'xr');hold off;axis equal

%polyfit:
timepoint=floor(min(cellContour(:,1),[],'all'))-1:1:ceil(max(cellContour(:,1),[],'all'))+1;
[p,S,mu]=polyfit(boneY,boneX,2);

if (S.normr>14 && stats.MajorAxisLength<100 ) || (S.normr>17 && stats.MajorAxisLength>100 )
%debugging-----
%     disp([num2str(S.normr),' - ',num2str(stats.MajorAxisLength)]);
%     baseline = polyval(p,timepoint,[],mu);%baseline=>x, timepoint=>y
%     F=figure("Visible","on");hold on;plot(cellContour(:,2),cellContour(:,1),'-');plot(baseline',timepoint','-xr');
%     plot(boneX,boneY,'xb');
%     hold off;axis equal
    cellMesh=[];skel=[];WL=[];cellContour=[];
    return;
end
baseline = polyval(p,timepoint,[],mu);%baseline=>x, timepoint=>y

%debugging-----
% F=figure("Visible","on");hold on;plot(cellContour(:,2),cellContour(:,1),'-');plot(baseline',timepoint','-xr');
% plot(boneX,boneY,'xb');
% hold off;axis equal
% name2=[num2str(name(1)),'-',num2str(name(2)),'-',num2str(name(3)),'-',num2str(name(4))];
% saveas(F,['E:\test\',num2str(name2,'%0.3i'),'.tif']);
% close(F);
%-----debugging
cellContour=flipud(cellContour);

%get mesh
[cellMesh,skel,WL] = model2meshZZSimp([baseline',timepoint'],cellContour,0.5,0.01,14);%coordinatePoints,stepSize,tolerance,meshWidth)
%cellMesh = model2mesh(cellContour,1,0.01,14);%coordinatePoints,stepSize,tolerance,meshWidth)
% 
% toPlot=[cellMesh(:,1:2);flipud(cellMesh(:,3:4))];
% F=figure;hold on
% zx=imshow(J);
% AX=get(zx,'Parent');
% plot(AX,cellContour(:,2),cellContour(:,1),'-r');
% plot(AX,skel(:,1),skel(:,2),'.b');
end


