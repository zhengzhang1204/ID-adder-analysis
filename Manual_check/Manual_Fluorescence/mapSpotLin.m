function [Map,CellReg]=mapSpotLin(uit,Nxt3CellID,isSkip2nd,Map,SelNxt3,CellReg)
%1. cell frame
%2. cellID (from channel bottom)
%3. spot relative location in cell
%4. oric# (optional)
%5. ?
%6. ?
%7. spot centroid -row
%8. spot centroid -column
%9. mean intensity
%10. median intensity
%11. 98% quantile
%12. cell length
%13. spot mean location -row
%14. spot mean location -column

% Format of uit:
% category of spot in Cell A, PSB of spot in Cell B, CellID, SpotID
A=SelNxt3{1,1}{1,1};

for iSN3=1:length(A)
    for iMap=1:size(uit{iSN3},1)
        CellID=uit{iSN3}(iMap,end-1);
        SpotID=uit{iSN3}(iMap,end);
        if CellID>0
            if isSkip2nd(iSN3)
                Map(A(iSN3,1),iMap)=SelNxt3{1,3}{CellID,1}(SpotID,1);
            else
                Map(A(iSN3,1),iMap)=SelNxt3{1,2}{CellID,1}(SpotID,1);
            end
        else
            Map(A(iSN3,1),iMap)=-1;
        end
    end
end
CellReg=cat(1,CellReg,Nxt3CellID{1,1});

end