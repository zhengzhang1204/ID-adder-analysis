function [cellMesh,skel,WL] = model2meshZZSimp(skel,coordinatePoints,stepSize,~,meshWidth)
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%function res = model2mesh(coordinatePoints,stp,tolerance,lng)
%oufti.v0.3.0
%@author:  oleksii sliusarenko
%@copyright 2012-2014 Yale University
%==========================================================================
%**********output********:
%cellMesh:  created mesh from coordinate points.
%**********Input********:
%coordinatePoints:  coordinate vector for a cell contour.
%stepSize:  steps to void between each segment in a mesh.
%tolerance: ?
%meshWidth: width of the mesh to be created.
%=========================================================================
% PURPOSE:
% This function performs a medial axis transform to a non-branching axis
% Takes the outline coordinates, step size on the final centerline,
% tolerance to non-centrality, and the length of the ribs. The last
% parameter should be longer than the ribs, but should not be too long to
% intersect the countour agan: though most cases of >1 intersection will
% be resolved by the function, some will not. Outputs the coordinates of
% the centerline points.
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

%% ZZ start
% interpolate skeleton to equal step, extend outside of the cell and smooth
% tolerance=0.001;
coordinatePoints = [coordinatePoints;coordinatePoints(1,:)];
%figure;hold on;plot(coordinatePoints(:,2),coordinatePoints(:,1),'-');plot(skel(:,1),skel(:,2),'xr');hold off;
[intersectX,intersectY,ii]=polyxpoly(skel(:,1),skel(:,2),coordinatePoints(:,2),coordinatePoints(:,1));
if length(intersectY)~=2
    cellMesh=[];skel=[];WL=[];return;
end
[~,idx]=min(intersectY);
interXst=intersectX(idx);
interYst=intersectY(idx);
ii=ii(1,idx);
skel=[[interXst,interYst];skel(ii+1:end,:)];

d=diff(skel,1,1);
l=cumsum([0;sqrt((d.*d)*[1 ;1])]);
if size(skel,1)<3
    cellMesh=[];skel=[];WL=[];return;
end
skel = [interp1(l,skel(:,1)',0:stepSize:l(end))' interp1(l,skel(:,2)',0:stepSize:l(end))'];

[intersectX,intersectY,ii]=polyxpoly(skel(2:end,1),skel(2:end,2),coordinatePoints(:,2),coordinatePoints(:,1));
if length(intersectY)~=1
    cellMesh=[];skel=[];WL=[];return;
end

skel(ii(1)+2:end,:)=[];
skel=[skel;[intersectX,intersectY]];
if ((skel(end-1,1)-skel(end,1))^2+(skel(end-1,2)-skel(end,2))^2)<stepSize/2
    skel(end-1,:)=[];
end
coordinatePoints = [coordinatePoints;coordinatePoints(1,:)];

% recenter and smooth again the skeleton
[cellMesh,WL]=skel2meshZZSimp(fliplr(skel),coordinatePoints,meshWidth,stepSize);
if isempty(cellMesh)
    cellMesh=[];skel=[];WL=[];return;
end
% figure;hold on
% plot(coordinatePoints(:,2),coordinatePoints(:,1),'-');
% % plot(skel(:,1),skel(:,2),'-xr');
% plot(skel(:,1),skel(:,2),'-xb');
% tobeplot=[cellMesh(:,1:2);flipud(cellMesh(:,3:4))];
% plot(tobeplot(:,1),tobeplot(:,2),'-r');
% hold off
% axis equal

% if numel(cellMesh)==1 || length(cellMesh)<=4, cellMesh=0; disp('Unable to create mesh'); end
% if length(cellMesh)>1 && (cellMesh(1,1)~=cellMesh(1,3) || cellMesh(end,1)~=cellMesh(end,3)), cellMesh=0; disp('Mesh creation error! Cell rejected'); end

end 
function [meshUit,WL]=skel2meshZZSimp(sk,coordinatePoints,meshWidth,stepSize)
% This function finds intersections of ribs with the contour
if isempty(sk), pintx=[]; pinty=[]; q=false; error('what is wrong at Line 97'); end

pintx = zeros(size(sk,1),2);
pinty = zeros(size(sk,1),2);

d=diff(sk,1,1);
len=cumsum([0;sqrt((d.*d)*[1 ;1])]);%note: the last value of len is the length
plinesx1 = repmat(sk(1:end-1,1),1,2)+meshWidth/stepSize*d(:,2)*[0 1];% left wing Y
plinesy1 = repmat(sk(1:end-1,2),1,2)-meshWidth/stepSize*d(:,1)*[0 1];% left wing X
plinesx2 = repmat(sk(1:end-1,1),1,2)+meshWidth/stepSize*d(:,2)*[0 -1];% right wing Y
plinesy2 = repmat(sk(1:end-1,2),1,2)-meshWidth/stepSize*d(:,1)*[0 -1];% right wing X

% Define the first pair of intersections as the prevpoint
pintx(1,:) = [sk(1,2) sk(1,2)];
pinty(1,:) = [sk(1,1) sk(1,1)];
pintx(end,:) = [sk(end,2) sk(end,2)];
pinty(end,:) = [sk(end,1) sk(end,1)];
    
% figure;hold on
% plot(coordinatePoints(:,2),coordinatePoints(:,1),'-');
% plot(sk(:,2),sk(:,1),'-xr')
% 
% for i=1:size(plinesx1,1)
%     plot([plinesy1(i,1),plinesy1(i,2)],[plinesx1(i,1),plinesx1(i,2)],'-b');
%     plot([plinesy2(i,1),plinesy2(i,2)],[plinesx2(i,1),plinesx2(i,2)],'-g');
% end
% hold off
% axis equal

ax=quantile(1:1:(size(plinesx1,1)+1),[0.25,0.35,0.65,0.75]);
sel=[floor(ax(1)):1:ceil(ax(2)), floor(ax(3)):1:ceil(ax(4))];
widthPool=nan(length(sel),1);iWDP=1;
rescaledLeft=nan(size(plinesx1,1)-1,2);
rescaledRight=nan(size(plinesx1,1)-1,2);
for i=2:size(plinesx1,1)%i=2~46, all=1-47
    [res1x,~]=polyxpoly([plinesy1(i,1),plinesy1(i,2)],[plinesx1(i,1),plinesx1(i,2)],coordinatePoints(:,2),coordinatePoints(:,1));
    [res2x,~]=polyxpoly([plinesy2(i,1),plinesy2(i,2)],[plinesx2(i,1),plinesx2(i,2)],coordinatePoints(:,2),coordinatePoints(:,1));%right
    if isempty(res1x) ||isempty(res2x) || length(res1x)~=1 || length(res2x)~=1
        meshUit=[];WL=[];return
    end
    [pintx(i,1),pinty(i,1)]=polyxpoly([plinesy1(i,1),plinesy1(i,2)],[plinesx1(i,1),plinesx1(i,2)],coordinatePoints(:,2),coordinatePoints(:,1));
    [pintx(i,2),pinty(i,2)]=polyxpoly([plinesy2(i,1),plinesy2(i,2)],[plinesx2(i,1),plinesx2(i,2)],coordinatePoints(:,2),coordinatePoints(:,1));%right
    useWidth=sqrt((pintx(i,1)-pintx(i,2))^2 + (pinty(i,1)-pinty(i,2))^2)/2;
    vModL=sqrt((plinesy1(i,2)-plinesy1(i,1))^2+(plinesx1(i,2)-plinesx1(i,1))^2);
    rescaledLeft(i-1,1)=useWidth/vModL*(plinesy1(i,2)-plinesy1(i,1))+plinesy1(i,1);
    rescaledLeft(i-1,2)=useWidth/vModL*(plinesx1(i,2)-plinesx1(i,1))+plinesx1(i,1);

%     vModR=sqrt((plinesy2(i,2)-plinesy2(i,1))^2+(plinesx2(i,2)-plinesx2(i,1))^2);
    rescaledRight(i-1,1)=useWidth/vModL*(plinesy2(i,2)-plinesy2(i,1))+plinesy2(i,1);
    rescaledRight(i-1,2)=useWidth/vModL*(plinesx2(i,2)-plinesx2(i,1))+plinesx2(i,1);
    if ismember(i,sel)
        widthPool(iWDP)=sqrt((pintx(i,1)-pintx(i,2))^2 + (pinty(i,1)-pinty(i,2))^2);
        iWDP=iWDP+1;
    end
end

% figure;hold on
% plot(coordinatePoints(:,2),coordinatePoints(:,1),'-');
% plot(sk(:,2),sk(:,1),'-r')
% for i=1:size(plinesx1,1)
%     plot([plinesy1(i,1),plinesy1(i,2)],[plinesx1(i,1),plinesx1(i,2)],'-b');
%     plot([plinesy2(i,1),plinesy2(i,2)],[plinesx2(i,1),plinesx2(i,2)],'-g');
% end
% plot(rescaledLeft(:,1),rescaledLeft(:,2),'xb');
% plot(rescaledRight(:,1),rescaledRight(:,2),'xb');
% hold off

pintx(2:end-1,:)=[rescaledLeft(:,1),rescaledRight(:,1)];
pinty(2:end-1,:)=[rescaledLeft(:,2),rescaledRight(:,2)];
meshUit=[pintx(:,1),pinty(:,1),pintx(:,2),pinty(:,2)];

% plot(pintx(:,1),pinty(:,1),'ob');
% plot(pintx(:,2),pinty(:,2),'ob');
% hold off
WL=[mean(widthPool),len(end)];

% [~,idx1]=min(abs(len-(len(end)-wid)));
% [~,idx2]=min(abs(len-wid));


end % function skel2mesh

