function [drift,tForm,sameAsInput]=DL_NeoManualAlignChannels(FLpath,PHpath,Xpath)
%This program registers dual camera with a rotation matrix as output.
%% Usage:
%1. From horizontalized images, choose one and open its PH and FL images into ImageJ
%2. For each FL channel, register a few points over the FL image and the corresponding PH image.
%3. Retrieve the point locations and save in an Excel sheet with row format: [xFL, yFL, xPH, yPH].
%4. Type mapSpots=[];
%5. Copy all point locations (in the format as mentioned above) into mapSpots and Run this codes.
%6. Obtain drift values.

%% Main
%1. read rotated images, including PH and a fluorecent channel:
R=imread(FLpath);
M=imread(PHpath);M=imadjust(M);
Rmin=floor(double(min(R,[],'all'))/25)*25;
Rmax=ceil(double(max(R,[],'all'))/25)*25;
load(Xpath);
%2. In GUI, use 'a' 's' 'd' 'w' to drag the fluorecent image, left click to select aligned spots.
tool = NeoImageAlign(M,R,pat,[Rmin,Rmax]);
selSpots=tool.tabData;

%3. calculate the Mxx by selected spots.
[tForm,Mxx]=getTForm(selSpots);
p=getParamsFromTMatrix(median(Mxx,3));
disp(['dX: ',num2str(p.dx)]);
disp(['dY: ',num2str(p.dy)]);
disp(['theta: ',num2str(p.theta)]);
disp(['sX: ',num2str(p.sx)]);
disp(['sX: ',num2str(p.sx)]);


%4. generate a fluorescent image using Mxx.
sameAsInput = affineOutputView(size(R),tForm,'BoundsStyle','SameAsInput');
Rr = imwarp(imadjust(R),tForm,'OutputView',sameAsInput);

%5. align this affined image with the original PH image. Register the 2D drift values.
tool = ImageAlign(M,Rr,pat);
drift=tool.driftX;
end
%% nested functions.
function [tform,UIT]=getTForm(mapSpots)
%% Main
%1. Obtain transformation matrix by selected points.
v = 1:1:size(mapSpots,1);
Cs = nchoosek(v,3);
UIT=zeros(3,3,size(Cs,1));
for Q=1:size(Cs,1)
    R=Cs(Q,:);
    xa0=mapSpots(R(1),3);
    ya0=mapSpots(R(1),4);
    za0=1;
    xa1=mapSpots(R(1),1);
    ya1=mapSpots(R(1),2);
    za1=1;
    
    xb0=mapSpots(R(2),3);
    yb0=mapSpots(R(2),4);
    zb0=1;
    xb1=mapSpots(R(2),1);
    yb1=mapSpots(R(2),2);
    zb1=1;   
    
    xc0=mapSpots(R(3),3);
    yc0=mapSpots(R(3),4);
    zc0=1;
    xc1=mapSpots(R(3),1);
    yc1=mapSpots(R(3),2);
    zc1=1;    
    
    syms aS bS cS dS eS fS gS hS kS
    eqns=[
        aS*xa0+bS*ya0+cS*za0==xa1;
        dS*xa0+eS*ya0+fS*za0==ya1;
        gS*xa0+hS*ya0+kS*za0==za1;
        
        aS*xb0+bS*yb0+cS*zb0==xb1;
        dS*xb0+eS*yb0+fS*zb0==yb1;
        gS*xb0+hS*yb0+kS*zb0==zb1;     
        
        aS*xc0+bS*yc0+cS*zc0==xc1;
        dS*xc0+eS*yc0+fS*zc0==yc1;
        gS*xc0+hS*yc0+kS*zc0==zc1;          

    ];
    S=vpasolve(eqns,[aS,bS,cS,dS,eS,fS,gS,hS,kS]);
    A=double(S.aS);
    B=double(S.bS);
    C=double(S.cS);
    D=double(S.dS);
    E=double(S.eS);
    F=double(S.fS);
    G=double(S.gS);
    H=double(S.hS);
    I=double(S.kS);
    UIT(:,:,Q)=[A,D,G;B,E,H;C,F,I];
end
tform = affine2d(median(UIT,3));%NOTE: transformation matrix is median(UIT,3).

end

function [p]=getParamsFromTMatrix(Mxx)
% The following can calculate drift vector (x,y), rotation angle (theta, in degree), zooming factor(Sx, Sy)
% from Mxx. This process is optional.
syms x y theta sx sy
eqns=[...
sx*cosd(theta)==Mxx(1,1);...
sy*sind(theta)==Mxx(1,2);...
-sx*sind(theta)==Mxx(2,1);...
%sy*cosd(theta)==Mxx(2,2);...
sx*(cosd(theta)*x-sind(theta)*y)==Mxx(3,1);
sy*(sind(theta)*x+cosd(theta)*y)==Mxx(3,2)];
S = vpasolve(eqns,[x y theta sx sy]);
p.dx=double(S.x);
p.dy=double(S.y);
p.theta=double(S.theta);
p.sx=double(S.sx);
p.sy=double(S.sy);
end
