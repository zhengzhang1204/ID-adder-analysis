function [uitX,uitY]=contour2meshZZ(X,Y,refL)

% input ready
[~,I,~]=unique([X,Y],'rows');
if length(I)~=length(X)
    error('There can not be any repeat values in (X,Y).')
end
X=X';Y=Y';

%****** BRANCHING *****
%get voronoi vertices
[Vx,Vy]=voronoi(X,Y);

%sig figs
sig=4;
Vx=round(Vx*10^sig)/10^sig;
Vy=round(Vy*10^sig)/10^sig;

%get voronoi cells
% [V,~]=voronoin([X',Y']);
% V_inside=find(inpolygon(V(:,1),V(:,2),X,Y))';
% 
% V=round(V*10^sig)/10^sig;

%determine which points lie inside the polygon
pts_inside=inpolygon(Vx,Vy,X,Y);

%find the segements that lie entirely inside the polygon
segs_inside=find(sum(pts_inside,1)==2);


%************************************
%   PARSE NEUTRAL AXES
%************************************

%get internal segments
brch_X=Vx(:,segs_inside);
brch_Y=Vy(:,segs_inside);

%calculate number of distinct internal points and point connection matrix
points0=unique([reshape(brch_X,[],1),reshape(brch_Y,[],1)],'rows');
npts=size(points0,1);
point_mat=zeros(npts,npts);
for i=1:npts
    %get current point
    curr_pt=points0(i,:);
    
    %generate list of segments in voronoi space that share this point
    listx=brch_X==curr_pt(1);
    listy=brch_Y==curr_pt(2);
    list_both=listx.*listy;
    
    %indices of shared segments
    seg_list=find(sum(list_both,1)>0);
    
    %list of connected points in X Y space
    for j=1:length(seg_list)
        %ignore zero-length segments created by round-off error
        if ~and(brch_X(1,seg_list(j))==brch_X(2,seg_list(j)),brch_Y(1,seg_list(j))==brch_Y(2,seg_list(j)))
            %get distinct row
            row1=find(list_both(:,seg_list(j))==0);
            
            %record connected points
            xconn=brch_X(row1,seg_list(j));
            yconn=brch_Y(row1,seg_list(j));
            
            %find this point in the condensed list
            listx2=points0(:,1)==xconn;
            listy2=points0(:,2)==yconn;
            listboth2=find((listx2.*listy2)');
            
            %get list of points that are connected to current point
            point_mat(i,listboth2)=1;
        end
    end
end

%get point connection degrees
%degrees > 2 indicate a branch point
degrees=sum(point_mat,1);

if max(degrees(:))>3
    warning('Analysis found a branch point with more than three connections.')
end

%first start by handling end points of the network
ends=find(degrees==1)';
nodes=find(degrees>2);
node_num=length(nodes);

%create first grouping matrix (remove nodes from the regular diagram)
group_mat=point_mat+eye(npts);
group_mat(nodes',:)=0;
group_mat(:,nodes)=0;

%find matrix representation of nodes
node_mat=zeros(npts,npts);
node_mat(nodes',:)=1;
node_mat(:,nodes)=1;
node_mat=point_mat.*node_mat;

%find ends connected directly to nodes
lones=find(degrees==1); %these are the ends

q=0;
for i=1:length(lones)
    curr_vec=zeros(npts,1);
    curr_vec(lones(i))=1;
    
    end1=find(node_mat*curr_vec);
    
    if ~isempty(end1)
        q=q+1;
        branches(q).degree=[degrees(lones(i)),degrees(end1)];
        branches(q).Xpos=[points0(lones(i),1),points0(end1,1)];
        branches(q).Ypos=[points0(lones(i),2),points0(end1,2)];
    end
end

%find other groups
[groups0,~]=graph_analysis_morphometrics(group_mat);

%find connecting nodes for these groups
for i=1:length(groups0)
    q=q+1;
    
    curr_vec=zeros(npts,1);
    curr_vec(groups0(i).elements)=1;
    
    %find elements with nodes in this branch group
    elements=find((node_mat*curr_vec+curr_vec)>0);
    nodes1=find((node_mat*curr_vec)>0);
    ends1=intersect(ends,groups0(i).elements);
    
    %get points in this branch
    X0=points0(elements,1);
    Y0=points0(elements,2);
    
    %choose starting point
    if ~isempty(ends1)
        P=find(elements==ends1(1));
    else
        P=find(elements==nodes1(1));
    end
    
    %get ordering of points into a coherent branch
    if length(X0)<=2
        uitX=[];uitY=[];return;
    end
    [Xout,Yout,I]=points2contour(X0,Y0,P,'cw');
    
    branches(q).degree=degrees(I);
    branches(q).Xpos=Xout;
    branches(q).Ypos=Yout;
end

%find nodes connected directly to other nodes
if node_num>1
    for i=1:node_num
        %get connections to first node
        list1=[nodes(i),find(point_mat(nodes(i),:))];
        for j=1:node_num
            if j>i
                %get connections to second node
                list2=[nodes(j),find(point_mat(nodes(j),:))];
                
                %find overlapping nodes
                overlaps=intersect(list2,list1);
                
                if length(overlaps)>2
                    disp(['Something is up with node-node branch determination at node: ' num2str(i)])
                end
                
                %record overlaps (should be two)
                if ~isempty(overlaps)
                    q=q+1;
                    
                    %handle node-node connection
                    if length(overlaps)==2
                        branches(q).degree=degrees(overlaps);
                        branches(q).Xpos=points0(overlaps,1);
                        branches(q).Ypos=points0(overlaps,2);
                    else %handle node-point-node connection
                        p1=[nodes(i),overlaps(1),nodes(j)];
                        
                        branches(q).degree=degrees(p1);
                        branches(q).Xpos=points0(p1,1);
                        branches(q).Ypos=points0(p1,2);
                    end
                end
            end
        end
    end
end

% %figure out which branches are neighbors
% branchN=length(branches);
% if branchN>1
%     branchL=zeros(branchN,branchN);
%     branchR=zeros(branchN,branchN);
%     for i=1:branchN
%         %compare LHS of branch to other branches
%         xLi=branches(i).Xpos(1);
%         yLi=branches(i).Ypos(1);
%         
%         %compare RHS of branch to other branches
%         xRi=branches(i).Xpos(end);
%         yRi=branches(i).Ypos(end);
%         
%         for j=1:branchN
%             if i~=j
%                 %compare LHS of branch to other branches
%                 xj=[branches(j).Xpos(1),branches(j).Xpos(end)];
%                 yj=[branches(j).Ypos(1),branches(j).Ypos(end)];
%                 
%                 %branch touching conditions
%                 cond1=and(xLi==xj(1),yLi==yj(1));
%                 cond2=and(xLi==xj(2),yLi==yj(2));
%                 cond3=and(xRi==xj(1),yRi==yj(1));
%                 cond4=and(xRi==xj(2),yRi==yj(2));
%                 
%                 if or(cond1,cond2)
%                     branchL(i,j)=1;
%                 end
%                 
%                 if or(cond3,cond4)
%                     branchR(i,j)=1;
%                 end
%             end
%         end
%         
%         %get list of neighbors
%         branches(i).neighbors_left=find(branchL(i,:));
%         branches(i).neighbors_right=find(branchR(i,:));
%     end
% end

%calculate branch lengths
BrchL=zeros(length(branches),1);
for i=1:length(branches)
    Xcurr=branches(i).Xpos;
    Ycurr=branches(i).Ypos;
    D=sqrt((Xcurr(1:end-1)-Xcurr(2:end)).^2+(Ycurr(1:end-1)-Ycurr(2:end)).^2);
    BrchL(i)=sum(D);
end

[~,idx]=max(BrchL,[],'all');
[test,idSort]=sort(BrchL,'descend');
if (length(BrchL)>=3 && (test(1)-test(2))/test(1) < 0.5) || (BrchL(idx) < refL*0.65 && refL>31 && length(BrchL)>=3)
    k1=[branches(idSort(1)).Ypos(end)-branches(idSort(1)).Ypos(1),branches(idSort(1)).Xpos(end)-branches(idSort(1)).Xpos(1)];
    k2=[branches(idSort(2)).Ypos(end)-branches(idSort(2)).Ypos(1),branches(idSort(1)).Xpos(end)-branches(idSort(1)).Xpos(1)];
    sigma = acosd(dot(k1,k2)/(norm(k1)*norm(k2)));
    if (sigma<50 || sigma>130 )&& length(branches(idSort(2)).Xpos)>=5
        uitX=[branches(idSort(1)).Xpos,branches(idSort(2)).Xpos];
        uitY=[branches(idSort(1)).Ypos,branches(idSort(2)).Ypos];
%         skL=BrchL(idSort(1))+BrchL(idSort(2));
    else
        uitX=branches(idx).Xpos;
        uitY=branches(idx).Ypos;
%         skL=BrchL(idx);
    end
else
    uitX=branches(idx).Xpos;
    uitY=branches(idx).Ypos;
%     skL=BrchL(idx);
end
% if skL/refL<0.65
%     uitX=[];uitY=[];
% %     disp(num2str(skL/refL));
% end
end

















