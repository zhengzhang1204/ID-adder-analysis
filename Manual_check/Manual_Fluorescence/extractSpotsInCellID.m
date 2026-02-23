function [SpotsInfoInCell,CellIdNxt3,SelUit]=extractSpotsInCellID(curID,List,Lxg,Int)
%Note: curID is 'current frame + cur CellID (from bottom)'
%  Format: 
%1. cell frame
%2. cellID (start from channel bottom)
%3. spot relative location in cell (LRT)
%4. oric# (optional)
%5. cellLength
%6. SpotY (length from image roof)
%7. SpotX (local tile)
%8. CoorInCell_Y
%9. Marker. 0: big cell, 1: small cell (<30% quantile), 2: small cell with
%squeezed spot LRT
%10. mean intensity
%11. median intensity
%12. 98% quantile
%13. cell length
%14. spot mean location -row
%15. spot mean location -column

%1. return if cellID is at the end of movie.
if curID(1)>size(Lxg,2)-2 || curID(2)==0%added on 23.2.13
    SpotsInfoInCell=[];
    CellIdNxt3{1,2}=[];
    SelUit=[];
    return;
end

%2. get next 3 cell from current cellID. Format: cellID, frame
CellIdNxt3=cell(1,3);
CellIdNxt3{1,1}=fliplr(curID);

NxtID=Lxg{curID(2),curID(1)};
if length(NxtID)==2
    CellIdNxt3{1,2}=[NxtID;[curID(1)+1,curID(1)+1]]';
elseif length(NxtID)==1
    CellIdNxt3{1,2}=[NxtID,curID(1)+1];
elseif isempty(NxtID)
    %do not fill anything
end

if size(CellIdNxt3{1,2},1)==2
    uitX=zeros(2,2);
    for iCIN3=1:2%intended to 2
        if CellIdNxt3{1,2}(iCIN3,1)==-1%case divLost
            NxtID3=nan;
            uitX(iCIN3,:)=[NxtID3,curID(1)+2];
            continue;
        end
        NxtID3=Lxg{CellIdNxt3{1,2}(iCIN3,1),CellIdNxt3{1,2}(iCIN3,2)};
        if length(NxtID3)==2
            error('Continuous division detected!');
        end
        uitX(iCIN3,:)=[NxtID3,curID(1)+2];
    end
    CellIdNxt3{1,3}=uitX;
elseif size(CellIdNxt3{1,2},1)==1
    NxtID3=Lxg{CellIdNxt3{1,2}(1),CellIdNxt3{1,2}(2)};
    if length(NxtID3)==2
        CellIdNxt3{1,3}=[NxtID3;[curID(1)+2,curID(1)+2]]';
    elseif length(NxtID3)==1
        CellIdNxt3{1,3}=[NxtID3,curID(1)+2];
    end
elseif isempty(CellIdNxt3{1,2})
    %do not fill anything
end
clear 'uitX';

%3. collect all spots (row) in these cellID
SpotsRowInCell=cell(1,3);
SpotsInfoInCell=cell(1,3);
SelUit=cell(1,3);
for iFrm=1:3%intended to be 1:3
    if size(CellIdNxt3{1,iFrm},1)==2
        uitX=cell(2,1);
        uitY=cell(2,1);
        uitZ=cell(2,1);

        uitX{1,1}=find(List(:,1)==CellIdNxt3{1,iFrm}(1,2) & List(:,2)==CellIdNxt3{1,iFrm}(1,1));
        sel1=find(List(:,1)==CellIdNxt3{1,iFrm}(1,2) & List(:,2)==CellIdNxt3{1,iFrm}(1,1));
        uitY{1,1}=cat(2,List(sel1,:),Int(sel1,:));
        uitZ{1,1}=sel1;

        uitX{2,1}=find(List(:,1)==CellIdNxt3{1,iFrm}(2,2) & List(:,2)==CellIdNxt3{1,iFrm}(2,1));
        sel2=find(List(:,1)==CellIdNxt3{1,iFrm}(2,2) & List(:,2)==CellIdNxt3{1,iFrm}(2,1));
        uitY{2,1}=cat(2,List(sel2,:),Int(sel2,:));
        uitZ{2,1}=sel2;
        
        SpotsRowInCell{1,iFrm}=uitX;
        SpotsInfoInCell{1,iFrm}=uitY;
        SelUit{1,iFrm}=uitZ;
    elseif size(CellIdNxt3{1,iFrm},1)==1
        SpotsRowInCell{1,iFrm}={find(List(:,1)==CellIdNxt3{1,iFrm}(2) & List(:,2)==CellIdNxt3{1,iFrm}(1))};
        sel=find(List(:,1)==CellIdNxt3{1,iFrm}(2) & List(:,2)==CellIdNxt3{1,iFrm}(1));
        SpotsInfoInCell{1,iFrm}={cat(2,List(sel,:),Int(sel,:))};
        SelUit{1,iFrm}{1,1}=sel;


    elseif isempty(CellIdNxt3{1,iFrm})
        %do not fill anything
    end
end
end