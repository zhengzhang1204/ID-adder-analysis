function [Map,ini]=imgFL_Spt2Lxgv3_simp(ini,lxg,intInf)
%inputs are: (SList,lxg,IntInfo,TSCenCoor)
% prepare mapping:
Map=nan(size(ini,1),3);%[cat, next1, next2]
Mark=true(size(ini,1),1);
CellIdReg=[];

% mapping
for i=1:size(ini,1)%for each row in ini<-----------------change back to 1
%     if i==4893
%         disp('C');
%     end
    if ~Mark(i)
        continue;
    end
%3.2 collect all spots in this cellID.
    if ~RemDupCellID(ini(i,1:2),CellIdReg)
        continue;
    end
%3.3 collect the current Cell (A) and the next 2 cells (B and C)
    [SpotInfoInNxt3,Nxt3CellID,SelNxt3]=extractSpotsInCellID(ini(i,1:2),ini,lxg,intInf);%get all spots (row numbers) in the following 3 frames.
    if isempty(Nxt3CellID{1, 2})
        continue;
    end
%3.4 Main. For each spot, judge its connenction in the 3-consecutive cells. 
    [uit,isSkip2ndCell]=getSpotLin_FAST(SpotInfoInNxt3,i);
    %format: [CellID, frame]
%3.5 Main. Write the Mapping table.
    [Map,CellIdReg]=mapSpotLin(uit,Nxt3CellID,isSkip2ndCell,Map,SelNxt3,CellIdReg);
end

Map=Map(:,1:3);
Map=RemSinglets(Map);

end

%% Nested functions:

function [uit]=RemDupCellID(cur,CellIdReg)
    %if cur is on the reg-list, false.
    %if cur is not on the reg-list, true.
    % NOTE: cur: [frame, cellID]
    %       CellIdReg: [cellID, frame]
    uit=true;
    if isempty(CellIdReg)
        uit=true;
        return;
    end
    
    for iRow=1:size(CellIdReg,1)
        if isequal(fliplr(cur),CellIdReg(iRow,:))
            uit=false;
            return;
        end
    end
end


