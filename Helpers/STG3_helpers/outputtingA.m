%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%function [outA]=outputtingA(SptMaps,SptLines,PHMaps,PHLines,divNumByGen,oriLines,oriNummat,InfCell,Spt2PHMap,int)
%Stg3Neo_V2
%Author:  Zheng Zhang
%==========================================================================
%**********Output********
%outA:          Results about spot lines
%**********Input*********
%SptLines:      CIDs of all spot lines
%SptMaps:       Map of spot lines
%PHMaps:        Map of PH lines
%PHLines:       CIDs of all PH lines
%divNumByGen:   Gen number for each valid cell
%oriLines:      CIDs of all origin lines
%oriNummat:     ori number for each valid cell
%InfCell:       2D info of each valid cell
%Spt2PHMap:     Map of each spot line and its corresponding PH line
%int:           Acquisition interval in minutes.
%=========================================================================
% PURPOSE:
%generate information for II adder and others
%=========================================================================
% OTHER INFO:
% Format of outputA: 
% Col. 1~6: [Li/ori_m, Δii_U(by ori), Δii_L(by ori), Li/gen_m, Δii_U(by gen), Δii_L(by gen)] 
% Col. 7~14: [GRII_DaugU, GRII__DaugL, C-period (by spot line length, unit:Frm), Tii_U (unit:Frm), Tii_L (unit:Frm), ΔTii (min, initUp-initLow), Lbup,Lblow] 
% Col.15-17: [TInitdaughterUp,TinitdaughterLow,Tdiv](all are min)
% Col.18-19: [closest Ld/gen_m, gen_m] (Upon Chenli request)
% Col.20-21: [Lt, gen_t] real cell length at termination, generation number for Lt.@23.12.07
% Col.22-24: [t_start,t_end_U,t_end_L] Frame number for each spot line
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

function [outA]=outputtingA(SptMaps,SptLines,PHMaps,PHLines,divNumByGen,oriLines,oriNummat,InfCell,Spt2PHMap,int)
outA=[];
for i=1:size(SptMaps,1)%for each spot line:
    col=nan(1,24);
    if ~isnan(SptMaps(i,1)) || ~isnan(SptMaps(i,2))%if the current spot line has daughter(s).
        CIDm=SptLines{i,1}(1,2:3);%current spot line's initial CID (frame+CID).
        if CIDm(1)==1
            continue;
        end
        col(1)=InfCell{CIDm(1),1}{1,2}(CIDm(2),2)/oriNummat(CIDm(2),CIDm(1));   %Li/ori_m.
        col(4)=InfCell{CIDm(1),1}{1,2}(CIDm(2),2)/divNumByGen(i);               %Li/gen_m.

        if ~isnan(SptMaps(i,1))%if the upper daughter spot line exists
            CIDdup=SptLines{SptMaps(i,1),1}(1,2:3);
            curInitCycCIDup=oriLines{i,1};
            if CIDdup(2)<=size(InfCell{CIDdup(1),1}{1,2},1)
                col(2)=(InfCell{CIDdup(1),1}{1,2}(CIDdup(2),2)/oriNummat(CIDdup(2),CIDdup(1)))*2-col(1);    %Δii_U(by ori)
                col(5)=(InfCell{CIDdup(1),1}{1,2}(CIDdup(2),2)/divNumByGen(SptMaps(i,1)))*2-col(4);       %Δii_U(by gen)
            end
            [col(7),~]=getGRthrSpotLine(curInitCycCIDup,InfCell,PHLines,int);    %GR of upper ori
            col(10)=SptLines{SptMaps(i,1),1}(1,2)-SptLines{i,1}(1,2);     %Tii_U (unit:Frm)
            col(22)=CIDm(1);           %starting frame of current spot line
            col(23)=CIDdup(1)-1;       %ending frame of upper oriLine
        end
        if ~isnan(SptMaps(i,2))%if the upper daughter spot line exists
            CIDdlow=SptLines{SptMaps(i,2),1}(1,2:3);
            curInitCycCIDlow=oriLines{i,1};
            if CIDdlow(2)<=size(InfCell{CIDdlow(1),1}{1,2},1)
                col(3)=(InfCell{CIDdlow(1),1}{1,2}(CIDdlow(2),2)/oriNummat(CIDdlow(2),CIDdlow(1)))*2-col(1); %Δii_L(by ori)
                col(6)=(InfCell{CIDdlow(1),1}{1,2}(CIDdlow(2),2)/divNumByGen(SptMaps(i,2)))*2-col(4);      %Δii_L(by gen)
            end
            [col(8),~]=getGRthrSpotLine(curInitCycCIDlow,InfCell,PHLines,int);   %GR of upper ori
            col(11)=SptLines{SptMaps(i,2),1}(1,2)-SptLines{i,1}(1,2);     %Tii_U (unit:Frm)
            col(22)=CIDm(1);           %starting frame of current spot line
            col(24)=CIDdlow(1)-1;      %ending frame of lower oriLine
        end
        if ~isnan(SptMaps(i,1)) && ~isnan(SptMaps(i,2))%the current spot line has two descendants
            col(12)=int*(SptLines{SptMaps(i,1),1}(1,2)-SptLines{SptMaps(i,2),1}(1,2));%Δt in min. (T_initUp-T_initLow)
            col(15)=int*SptLines{SptMaps(i,1),1}(1,2);
            col(16)=int*SptLines{SptMaps(i,2),1}(1,2);
        end
        if ~isnan(SptMaps(i,1)) || ~isnan(SptMaps(i,2))%the current spot line has at least one descendant
            col(20)=InfCell{SptLines{i}(end,2),1}{1,2}(SptLines{i}(end,3),2);
            %any division during the current spot line?
            [~,PHroute]=extractCIDsBtwAB(SptLines{i}(1,2:3),SptLines{i}(end,2:3),PHLines,PHMaps);
            col(21)=(2^(length(PHroute)-1))*divNumByGen(i);%
        end

        col(9)=SptLines{i}(end,2)-SptLines{i}(1,2);% C period raw
        if ~isnan(Spt2PHMap(i,4))
            col(13)=InfCell{Spt2PHMap(i,5),1}{1,2}(Spt2PHMap(i,6),2);%upper daughter Lb
            col(17)=int*Spt2PHMap(i,5);%upper daughter initTime
        end
        if ~isnan(Spt2PHMap(i,7))
            col(14)=InfCell{Spt2PHMap(i,8),1}{1,2}(Spt2PHMap(i,9),2);%lower daughter Lb
            col(17)=int*Spt2PHMap(i,8);%lower daughter initTime
        end
        if ~isequal(SptLines{i}(1,2:3),[1,1])
            %look for the division of PHpart that contains this CID:
            indCIDinParts=outs_findCIDpartContainingCID(SptLines{i}(1,2:3),PHLines);
            selCID_A=PHLines{indCIDinParts}(end,1:2);
            col(18)=InfCell{selCID_A(1),1}{1,2}(selCID_A(2),2)/divNumByGen(i);
            col(19)=divNumByGen(i);%the closest Ld to this CID.
        end
    end
    if sum(isnan(col))~=11
        outA=cat(1,outA,col);
    end
end
end