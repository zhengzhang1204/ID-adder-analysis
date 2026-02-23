cllc;
clc;

%% User specified region
% Option 1: process a specified experiment
EXP='G:\RdnaA-seqAR-atc 1to3-M30 medium 20260130';int=3;WinK=3;WinL=10;%
coreF(EXP,int);

% % Option 2: process a list of experiments
% load('F:\mChrRpath.mat');%'G:\mChrRpath.mat';'G:\yPetRpath.mat'
% path=cellfun(@(x,y) [x,y], repmat({'F:\'},size(Rpath,1),1),Rpath(:,1),'UniformOutput',false);
% int=cell2mat(Rpath(:,3));
% for i=13:length(path)%for each experiment
%     coreF(path{i},int(i));
% end

%% core
function []=coreF(EXP,int)
Cperiod=41;%unit min C=41min for M6.
stEntry=1;%start from which entry in 'data'?
mappingMethod=4;
if exist([EXP,'\selChannels.mat'],'file')
    load([EXP,'\selChannels.mat']);
else
    K=dir([EXP,'\Y_res\p*_lkgY.mat']);
    data=[];
    for iK=1:length(K)
        load([K(iK).folder,'\',K(iK).name]);
        pos=str2double(K(iK).name(2:3));
        for iX=1:length(SLines)
            if ~isempty(SLinesMap{1,iX})
                data=cat(1,data,[pos,iX]);
            end
        end
    end
end
Entry=1:1:size(data,1);
if stEntry~=1
    data=[data(stEntry:end,:);data(1:stEntry-1,:)];
    Entry=[Entry(stEntry:end),Entry(1:stEntry-1)]';
end
% outAsy=[];
% outSyn=[];
% outOriN=[];
% 
% CVOutCDadder=zeros(4,2);
% CVOutDivSize=zeros(4,2);
% CVOutLi=zeros(4,2);
IDAdderLin=cell(size(data,1),1);
IIAdderLin=cell(size(data,1),1);
DDAdderLin=cell(size(data,1),1);
nAll=num2str(size(data,1));%str format
parfor i=1:size(data,1) % parfor each entry
    MPath=[EXP,'\M_res\p',num2str(data(i,1),'%0.2i'),'_lkgM.mat'];
    FPath=[EXP,'\Y_res\p',num2str(data(i,1),'%0.2i'),'_lkgY.mat'];
    disp(['Processing ',num2str(Entry(i)),'/',nAll]);
    [IDAdderLin{i},IIAdderLin{i},DDAdderLin{i}]=genResult0611(MPath,FPath,data(i,2),mappingMethod,round(Cperiod/int),int);
end
IDAdderLin=cat(1,IDAdderLin{:});
IIAdderLin=cat(1,IIAdderLin{:});
DDAdderLin=cat(1,DDAdderLin{:});

%save the above 3 Lin variables.
save([EXP,'\AdderLinV3.mat'],'IDAdderLin','IIAdderLin','DDAdderLin','EXP');
end



%% nested functions:
function [deltaIDLin,deltaIILin,deltaDDLin]=genResult0611(MPath,FPath,SchN,mapMet,Cperiod,int)
load(FPath,'SLinesparts','SLinesMap','SLines');
load(MPath,'Infcellv2','lxg','PHparts','PHLines','PHMap');%<-----------

% Initialization
InfCell=Infcellv2(:,SchN(1));%<-----------
if exist('Infcellv2','var')
    for iM=1:length(InfCell)
        if isempty(InfCell{iM})
            continue;
        end
        InfCell{iM}=InfCell{iM}(2:3);
    end
end
lxg=lxg{:,SchN(1)}{1,2};
SPTparts=SLinesparts{1,SchN(1)};
selSptMap=SLinesMap{1,SchN(1)};
CIDinParts=PHparts{1,SchN(1)};
PHparts=PHLines{1,SchN(1)};
HorLineMap=PHMap{1,SchN(1)};
selSpotLines=SLines{1,SchN(1)};
% SpotEndFrm=max(cellfun(@(x) x(end,2), selSpotLines),[],1);
%error here means: [tree plot is drawn wrongly]

% Generation of Spt2PHMap. Format of Spt2PHMap:
% For each spot line, Column 1 is the idx of the PH line that this spot line attaches to.
% Column 2 and 3 is the PH line's ending CID.
% Column 4-9 are for daughter's PH index and PH initial CID. 
% [ref cell's rowIdx+PHendingCID, upper cell's rowIdx+PHinitCID, lower cell's rowIdx+PHinitCID];
Spt2PHMap=nan(size(SPTparts,1),9);%main output.
for i=1:size(SPTparts,1)%for each spot line
    sptY=SPTparts(i,2);
    if isnan(sptY) && i==size(SPTparts,1)% if the last row in SPTparts is nan, omit it.
        Spt2PHMap=Spt2PHMap(1:i-1,:);
        SPTparts=SPTparts(1:i-1,:);
        selSptMap=selSptMap(1:i-1,:);
        selSpotLines=selSpotLines(1:i-1,:);
        selCIDCollect_beg=selCIDCollect_beg(1:i-1,:);
        selCIDCollect_end=selCIDCollect_end(1:i-1,:);
        continue;
    end
    idx=find(PHparts(:,2)==sptY);
    if ~isempty(idx)
        Spt2PHMap(i,1:3)=[idx,CIDinParts{idx,1}(end,:)];
        if ~isnan(HorLineMap(idx,2))
            Spt2PHMap(i,4:6)=[HorLineMap(idx,2),CIDinParts{HorLineMap(idx,2),1}(1,:)];
        end
        if ~isnan(HorLineMap(idx,3))
            Spt2PHMap(i,7:9)=[HorLineMap(idx,3),CIDinParts{HorLineMap(idx,3),1}(1,:)];
        end
    end
end

% Creating splitAndDivInfo:
% format of splitAndDivInfo:[UpSpotReInitFrm, LowSpotReInitFrm, PHdivFrm, PHmomEndLength, UpCellLength, LowCellLength]
splitAndDivInfo=nan(size(SPTparts,1),6);%per spot line
full_selSptMap=fullizeMap(selSptMap);
virtualSelSptLines=cell(size(selSpotLines));

%format of numOfSpotInEachCID: [CID matrix]
% virtualDNAmat=nan(max(cellfun(@(x) x(end,3), selSpotLines)),max(cellfun(@(x) x(end,2), selSpotLines)));REPLACED @220616
virtualDNAmat=nan(max(cellfun(@(x) x(end,2), CIDinParts)),max(cellfun(@(x) x(end,1), CIDinParts)));
% actualDNAmat=virtualDNAmat;
oriNummat=virtualDNAmat;
for i=1:size(SPTparts,1)%for each spot line.
% calculation of splitAndDivInfo:
    % UpSpotReInitFrm
    if mapMet==4
        if ~isnan(selSptMap(i,1))%THis is fast GR, count the 2nd init of upper spot.
            if ~isnan(selSptMap(selSptMap(i,1),1)) 
                splitAndDivInfo(i,1)=SPTparts(selSptMap(selSptMap(i,1),1),1);
            end
            if ~isnan(selSptMap(selSptMap(i,1),2))
                splitAndDivInfo(i,1)=SPTparts(selSptMap(selSptMap(i,1),2),1);
            end
        end
    
        % LowSpotReInitFrm
        if ~isnan(selSptMap(i,2))
            if ~isnan(selSptMap(selSptMap(i,2),1))
                splitAndDivInfo(i,2)=SPTparts(selSptMap(selSptMap(i,2),1),1);
            end
            if ~isnan(selSptMap(selSptMap(i,2),2))
                splitAndDivInfo(i,2)=SPTparts(selSptMap(selSptMap(i,2),2),1);
            end  
        end
    elseif mapMet==5
        if ~isnan(selSptMap(i,1))%THis is slow GR, count the 2nd init of upper spot.
            splitAndDivInfo(i,1)=selSpotLines{selSptMap(i,1),1}(1,2);
        end
        if ~isnan(selSptMap(i,2))
            splitAndDivInfo(i,2)=selSpotLines{selSptMap(i,2),1}(1,2);
        end
    end

    % PHdivFrm
    % PHmomEndLength
    if ~isnan(Spt2PHMap(i,1))
        if ~isnan(HorLineMap(Spt2PHMap(i,1),2)) || ~isnan(HorLineMap(Spt2PHMap(i,1),3))
            splitAndDivInfo(i,3)=PHparts(Spt2PHMap(i,1),3);
            if Spt2PHMap(i,3)>size(InfCell{Spt2PHMap(i,2),1}{1,1},1)
                disp(['Line',num2str(i),'Mom cell exceeding error']);
            else
                splitAndDivInfo(i,4)=InfCell{Spt2PHMap(i,2),1}{1,1}(Spt2PHMap(i,3),2);
            end
        else
            splitAndDivInfo(i,3)=nan;
            splitAndDivInfo(i,4)=nan;
        end
    end
    % UpCellLength - upper new born cell length
    if ~isnan(Spt2PHMap(i,4))
        if Spt2PHMap(i,6)>size(InfCell{Spt2PHMap(i,5),1}{1,1},1)
            disp(['Line',num2str(i),'Upper cell exceeding error']);
        else
            splitAndDivInfo(i,5)=InfCell{Spt2PHMap(i,5),1}{1,1}(Spt2PHMap(i,6),2);
        end
    end
    % LowCellLength - lower new born cell length
    if ~isnan(Spt2PHMap(i,7))
        if Spt2PHMap(i,9)>size(InfCell{Spt2PHMap(i,8),1}{1,1},1)
            disp(['Line',num2str(i),'Lower cell exceeding error']);
        else
            splitAndDivInfo(i,6)=InfCell{Spt2PHMap(i,8),1}{1,1}(Spt2PHMap(i,9),2);
        end
    end

% calculation of actualDNAmat:
% a). create virtualSelSptLines persuming constant C and accurate initiation event determination.
    if isnan(Spt2PHMap(i,1))
        curSptLineMomIdx=full_selSptMap(i,1);
        if ~isnan(curSptLineMomIdx)%this spot line has no mother and it lasts till end
            CIDs=getFullCIDs(Spt2PHMap(curSptLineMomIdx,1:3),CIDinParts,HorLineMap,selSpotLines{i,1}(1,2));
        else
            CIDs=[];
        end
        if isempty(CIDs)
            CIDs=findCIDpartContainingCID_SPE(selSpotLines{i,1},CIDinParts,HorLineMap);
        end
    else
        CIDs=getFullCIDs(Spt2PHMap(i,1:3),CIDinParts,HorLineMap,selSpotLines{i,1}(1,2));
    end
    if ~isnan(full_selSptMap(i,1)) && ~isempty(CIDs)
        virtualSelSptLines{i,1}=completeSptLines(selSpotLines{i,1},CIDs,Cperiod,i);
    else
        virtualSelSptLines{i,1}=selSpotLines{i,1};
    end
    
end

% c) For each SLine, assign a division number when this line is initiated.
divNumByGen=getDivNum4SLinesByGen(selSptMap,selSpotLines,CIDinParts,HorLineMap,Spt2PHMap);

%% genRes220611 - Summary
% A)init Adder by ori# method.
% B)DNAc, NCR before and after initiation.
% C)for each div, DNA content vs divSite: Δ(DNAcont) vs Δ(divSite)
% D)NCratio (2D) and 1D 
% E)1. DNA content at birth vs GR of divCyc;
%   2. NCratio vs GR of divCyc
%   3. DNA content vs SpeGR
%   4. NCratio vs SpeGR
% F)1. NCratio at birth vs inter-init-time (IIT)
% 	2. NCratio at birth vs C(original Cperiod)
% 	3. NCratio at birth vs GR
% ---------LCL assignments x3-----------
% GHI) C to D params, incl. C+D. [each init-div pair]
% J) ori-firing time in each div cyc. [each div cycle]
% K) Moving Window NCR
% L) NCRs before division and initiation.

% trace all PH lines to summarize oriPerCell
% subs=cell2mat(CIDinParts);
% inds=sub2ind(size(oriNummat),subs(:,2),subs(:,1));
% outOriN=uint8(oriNummat(inds));
% ZZ=find(outOriN==6);
% if ~isempty(ZZ)
%     [Zrow,Zcol]=ind2sub(size(oriNummat),inds(ZZ));
% end
% outOriN=[Zrow,Zcol];


% ----- output A init Adder by ori#, gen, ter methods. 
% format of outputA: Col. 1~6: [Li/ori_m, Li/ori_DaugU, Li/ori_DaugL, Li/gen_m, Li/gen_DaugU, Li/gen_DaugL] 
% Col. 7~12: [lambdaII_DaugU, lambdaII__DaugL, C-period (by spot line length, unit:Frm), IIT_U, IIT_L, Δt (min, initUp-initLow)]

outA=nan(size(selSptMap,1),12);
for i=1:size(selSptMap,1)%for each spot line:
    col=nan(1,12);
    if ~isnan(selSptMap(i,1)) || ~isnan(selSptMap(i,2))
        CIDm=selSpotLines{i,1}(1,2:3);%mother's initial CID (frame+CID).
        if CIDm(1)==1
            continue;
        end
       % col(1)=InfCell{CIDm(1),1}{1,1}(CIDm(2),2)/oriNummat(CIDm(2),CIDm(1));%get Li/ori_m.
        col(4)=InfCell{CIDm(1),1}{1,1}(CIDm(2),2)/divNumByGen(i);%get Li/gen_m.
        if ~isnan(selSptMap(i,1))
            CIDdup=selSpotLines{selSptMap(i,1),1}(1,2:3);
%            curInitCycCIDup=oriLines{i,1};
            if CIDdup(2)<=size(InfCell{CIDdup(1),1}{1,1},1)
              %  col(2)=(InfCell{CIDdup(1),1}{1,1}(CIDdup(2),2)/oriNummat(CIDdup(2),CIDdup(1)))*2-col(1);%get Li/ori_DU.
                col(5)=(InfCell{CIDdup(1),1}{1,1}(CIDdup(2),2)/divNumByGen(selSptMap(i,1)))*2-col(4);%get Li/gen_DU.
            end
 %           [col(7),~]=getGRthrSpotLine(curInitCycCIDup,InfCell,CIDinParts,int);
            col(10)=selSpotLines{selSptMap(i,1),1}(1,2)-selSpotLines{i,1}(1,2);
        end
        if ~isnan(selSptMap(i,2))
            CIDdlow=selSpotLines{selSptMap(i,2),1}(1,2:3);
%            curInitCycCIDlow=oriLines{i,1};
            if CIDdlow(2)<=size(InfCell{CIDdlow(1),1}{1,1},1)
              %  col(3)=(InfCell{CIDdlow(1),1}{1,1}(CIDdlow(2),2)/oriNummat(CIDdlow(2),CIDdlow(1)))*2-col(1);%get Li/ori_DL.
                col(6)=(InfCell{CIDdlow(1),1}{1,1}(CIDdlow(2),2)/divNumByGen(selSptMap(i,2)))*2-col(4);%get Li/ori_DL.
            end
%            [col(8),~]=getGRthrSpotLine(curInitCycCIDlow,InfCell,CIDinParts,int);
            col(11)=selSpotLines{selSptMap(i,2),1}(1,2)-selSpotLines{i,1}(1,2);
        end
        if ~isnan(selSptMap(i,1)) && ~isnan(selSptMap(i,2))
            col(12)=selSpotLines{selSptMap(i,1),1}(1,2)-selSpotLines{selSptMap(i,2),1}(1,2);%Δt in Frm. (T_initUp-T_initLow)
            col(12)=col(12)*int;%Δt in Frm. (T_initUp-T_initLow)
        end
        col(9)=selSpotLines{i}(end,2)-selSpotLines{i}(1,2);
    end
    if sum(isnan(col))~=11
        outA(i,:)=col;
    end
end

% ----- output GHI C+D related things
outGHI=nan(size(Spt2PHMap,1),14); %format:[C+D(frm), rawC, Li, Li/ori#, Li/genN, Ld, ori#@division,lambda*(C+D),lambda,LiRat].

%to plot C+D adder, use C10 and C6-C10
stCIDs=cellfun(@(x) x(1,2:3),selSpotLines,'UniformOutput',false);
SLstInPHCIDidx=nan(length(stCIDs),1);
for i=1:length(stCIDs)
    SLstInPHCIDidx(i)=find(cellfun(@(x) ~isempty(find(x(:,1)==stCIDs{i}(1) & x(:,2)==stCIDs{i}(2), 1)),CIDinParts,'UniformOutput',true));
end
clear('stCIDs');

for i=1:size(Spt2PHMap,1)%for each spot line
    %omit if this cell does not have a PH line to attach.
    if isnan(Spt2PHMap(i,1))
        continue;
    end
    %omit if this spot line is an initial thread.
    if selSpotLines{i,1}(1,2)==1
        continue;
    end
    %check if the corresponding PH line has division.
    if ~isnan(HorLineMap(Spt2PHMap(i,1),2)) && ~isnan(HorLineMap(Spt2PHMap(i,1),3))
        lastCID=CIDinParts{Spt2PHMap(i,1),1}(end,:);
        daugCIDs=[Spt2PHMap(i,5:6);Spt2PHMap(i,8:9)];
        startCID=selSpotLines{i}(1,2:3);
        [~,PHroute]=extractCIDsBtwAB(startCID,lastCID,CIDinParts,HorLineMap);
        CIDidDivEvent=getDivCIDs(PHroute,CIDinParts,HorLineMap);%the CIDidDivEvent is an array of cells, each cell is {[CID]}
        CIDrat=cellfun(@(x) InfCell{x(1),1}{1,1}(x(2),2),CIDidDivEvent);
        divRats=CIDrat(1:2:(length(CIDrat)-1))./(CIDrat(1:2:(length(CIDrat)-1))+CIDrat(2:2:length(CIDrat)));
        Li=InfCell{selSpotLines{i,1}(1,2),1}{1,1}(selSpotLines{i,1}(1,3),2);
        CD=Spt2PHMap(i,2)-selSpotLines{i,1}(1,2);%time of C+D
        [CIDcd,~]=extractCIDsBtwAB(selSpotLines{i,1}(1,2:3),lastCID,CIDinParts,HorLineMap);
        [lambda,~]=getGRthrSpotLine([selSpotLines{i,1}(1,2:3);CIDcd;lastCID],InfCell,CIDinParts,int);
        stIdxInPHpart=find(CIDinParts{SLstInPHCIDidx(i)}==selSpotLines{i}(1,2));
        Time2=(size(CIDinParts{SLstInPHCIDidx(i)},1)-stIdxInPHpart)*int;
        Ld=InfCell{daugCIDs(1,1),1}{1,1}(daugCIDs(1,2),2)+InfCell{daugCIDs(2,1),1}{1,1}(daugCIDs(2,2),2);
        outGHI(i,:)=[CD,...%C+D
                size(selSpotLines{i,1},1),...%rawC
                Li,...%Li
                Li./oriNummat(selSpotLines{i,1}(1,3),selSpotLines{i,1}(1,2)),...%Li/ori#
                Li./divNumByGen(i)...%Li/genN
                Ld,...%Ld was InfCell{lastCID(1),1}{1,1}(lastCID(2),2)
                oriNummat(lastCID(2),lastCID(1)),...%ori#@division
                lambda*CD,...%lambda*(C+D)
                lambda,...%lambda
                Li*prod(divRats),...%LiRats
                stIdxInPHpart*int,...
                Time2,...
                stIdxInPHpart./size(CIDinParts{SLstInPHCIDidx(i)},1),...
                (Time2/int)./size(CIDinParts{SLstInPHCIDidx(i)},1)];
    end
end

%deltaIDLin:
deltaIDLin=[];
for i=1:size(Spt2PHMap,1)
    if ~isnan(selSptMap(i,1))
        deltaIDLin=[deltaIDLin;[outGHI(i,6)-outGHI(i,10),outGHI(selSptMap(i,1),6)-outGHI(selSptMap(i,1),10)]];
    end
    if ~isnan(selSptMap(i,2))
        deltaIDLin=[deltaIDLin;[outGHI(i,6)-outGHI(i,10),outGHI(selSptMap(i,2),6)-outGHI(selSptMap(i,2),10)]];
    end
end
del=[];
for i=1:size(deltaIDLin,1)
    if isnan(deltaIDLin(i,1)) || isnan(deltaIDLin(i,2))
        del=[del;i];
    end
end
deltaIDLin(del,:)=[];

%deltaIILin:
deltaIILin=[];
for i=1:size(selSptMap,1)
    if isnan(outA(i,4))
        continue;
    end
    if ~isnan(outA(i,5))
        Add_mu=outA(i,5);%Δu
        idx_U=selSptMap(i,1);
        if ~isnan(outA(idx_U,5))
            Add_du=outA(idx_U,5);
            deltaIILin=cat(1,deltaIILin,[Add_mu,Add_du]);%Δu and Δuu
        end
        if ~isnan(outA(idx_U,6))
            Add_dl=outA(idx_U,6);
            deltaIILin=cat(1,deltaIILin,[Add_mu,Add_dl]);%Δu and Δul
        end
    end
    if ~isnan(outA(i,6))
        Add_ml=outA(i,6);%Δl
        idx_L=selSptMap(i,2);
        if ~isnan(outA(idx_L,5))
            Add_du=outA(idx_L,5);%Δlu
            deltaIILin=cat(1,deltaIILin,[Add_ml,Add_du]);%Δl and Δlu
        end
        if ~isnan(outA(idx_L,6))
            Add_dl=outA(idx_L,6);%Δll
            deltaIILin=cat(1,deltaIILin,[Add_ml,Add_dl]);%Δl and Δll
        end
    end
end

%deltaDDLin:
deltaDDLin=[];
for i=2:size(HorLineMap,1)
    if isnan(HorLineMap(i,2)) || isnan(HorLineMap(i,3)); continue; end
    daugCIDs=[CIDinParts{HorLineMap(i,2)}(end,:);CIDinParts{HorLineMap(i,3)}(end,:)];
    cur_Ld=InfCell{daugCIDs(1,1),1}{1,1}(daugCIDs(1,2),2)+InfCell{daugCIDs(2,1),1}{1,1}(daugCIDs(2,2),2);
    Add_cur=cur_Ld...%InfCell{CIDinParts{i}(end,1),1}{1,1}(CIDinParts{i}(end,2),2)...
        -InfCell{CIDinParts{i}(1,1),1}{1,1}(CIDinParts{i}(1,2),2);
    if ~isnan(HorLineMap(i,2)) && (~isnan(HorLineMap(HorLineMap(i,2),2)) && ~isnan(HorLineMap(HorLineMap(i,2),3)))
        daugCIDsU=[CIDinParts{HorLineMap(HorLineMap(i,2),2)}(end,:);CIDinParts{HorLineMap(HorLineMap(i,2),3)}(end,:)];
        curLdU=InfCell{daugCIDsU(1,1),1}{1,1}(daugCIDsU(1,2),2)+InfCell{daugCIDsU(2,1),1}{1,1}(daugCIDsU(2,2),2);
        stCID=CIDinParts{HorLineMap(i,2)}(1,1:2);
        % endCID=CIDinParts{HorLineMap(i,2)}(end,1:2);
        Add_u=curLdU...%InfCell{endCID(1),1}{1,1}(endCID(2),2)...
            -InfCell{stCID(1),1}{1,1}(stCID(2),2);
        deltaDDLin=cat(1,deltaDDLin,[Add_cur,Add_u]);
    end
    if ~isnan(HorLineMap(i,3)) && (~isnan(HorLineMap(HorLineMap(i,3),2)) && ~isnan(HorLineMap(HorLineMap(i,3),3)))
        daugCIDsD=[CIDinParts{HorLineMap(HorLineMap(i,3),2)}(end,:);CIDinParts{HorLineMap(HorLineMap(i,3),3)}(end,:)];
        curLdD=InfCell{daugCIDsD(1,1),1}{1,1}(daugCIDsD(1,2),2)+InfCell{daugCIDsD(2,1),1}{1,1}(daugCIDsD(2,2),2);
        stCID=CIDinParts{HorLineMap(i,3)}(1,1:2);
        % endCID=CIDinParts{HorLineMap(i,3)}(end,1:2);
        Add_d=curLdD...%InfCell{endCID(1),1}{1,1}(endCID(2),2)...
            -InfCell{stCID(1),1}{1,1}(stCID(2),2);
        deltaDDLin=cat(1,deltaDDLin,[Add_cur,Add_d]);
    end
end
end

function [out]=doCutoffXY(in,cutoff)
xResc=in(:,1)./mean(in(:,1));
yResc=in(:,2)./mean(in(:,2));
del=xResc<cutoff(1) | xResc>cutoff(2) | yResc<cutoff(1) | yResc>cutoff(2);
disp([num2str(find(length(del))),' out of ',num2str(size(in,1)),' were removed during cutoff.']);
in(del,:)=[];
out=in;

end
function [DivNum]=getDivNum4SLinesByGen(SptMap,SLines,PHLines,PHMap,S2Pmap)
SptMap=fullizeMap(SptMap);
initThread=find(isnan(SptMap(:,1)));
%find if each initial thread has a corresponding PHLINE.
DivNum=zeros(size(SptMap,1),1);
for iInt=1:length(initThread)
    if isnan(S2Pmap(initThread(iInt),1))
        continue;
    end
    if S2Pmap(initThread(iInt),1)==1
        DivNum(iInt)=1;
    else
        lnk=getMapLinkage(1,S2Pmap(initThread(iInt),1),PHMap);%corresponding PH line.
        DivNum(iInt)=2^(length(lnk)+1);
    end

end

%mark divNum for other spot lines.
for iSL=length(initThread)+1:size(SLines,1)
    idxSLm=SptMap(iSL,1);%index of mother spot line of the current spot line
    idxPHm=find(cellfun(@(x) ~isempty(find(SLines{idxSLm,1}(1,2)==x(:,1) & SLines{idxSLm,1}(1,3)==x(:,2), 1)), PHLines));
    idxPHd=find(cellfun(@(x) ~isempty(find(SLines{iSL,1}(1,2)==x(:,1) & SLines{iSL,1}(1,3)==x(:,2), 1)), PHLines));
    if idxPHm==idxPHd
        DivNum(iSL)=DivNum(idxSLm)*2;
    else
        lnk=getMapLinkage(idxPHm,idxPHd,PHMap);
        DivNum(iSL)=DivNum(idxSLm)/(2^(length(lnk)));
    end
end

end
function [uit]=getMapLinkage(idxM,idxD,PHmap)
uit=[];
while idxD~=idxM
    target=PHmap(idxD,1);
    if target~=0
        uit=cat(2,target,uit);
        idxD=target;
    else
        uit=[];
        return;
    end
end
if ~isempty(uit)
    uit(1)=[];
end
end
function [uit,cnt]=getGRthrSpotLine(CIDin,InfCell,CIDinParts,int)
uit=[];
Ls=nan(size(CIDin,1),1);
idx=find(cellfun(@(x) ~isempty(find(x(:,1)==CIDin(1,1) & x(:,2)==CIDin(1,2),1)), CIDinParts));
remLength=size(CIDinParts{idx,1},1)-find(CIDinParts{idx,1}(:,1)==CIDin(1,1),1)+1;
cnt=0;
i=1;
curPos=1;
while i<=size(CIDin,1)
    if curPos>remLength
        idx=find(cellfun(@(x) ~isempty(find(x(:,1)==CIDin(i,1) & x(:,2)==CIDin(i,2),1)), CIDinParts));
        remLength=size(CIDinParts{idx,1},1)-find(CIDinParts{idx,1}(:,1)==CIDin(i,1),1)+1;
        cnt=cnt+1;
        curPos=1;
    end
    Ls(i)=2^(cnt)*InfCell{CIDin(i,1),1}{1,1}(CIDin(i,2),2);
    i=i+1;
    curPos=curPos+1;
end

x=(1:1:length(Ls))'*int/60;
if length(Ls)>=5
    [f,~] = fit(x,Ls,'exp1');
    uit=f.b;
else
    uit=nan;
end

end
function [out,out2]=extractCIDsBtwAB(A,B,PHParts,PHMap)
%extract CIDs Btw CID_A and CID_B BEG
idx1=find(cellfun(@(x) ~isempty(find(x(:,1)==A(1) & x(:,2)==A(2),1)),PHParts));%which CID in CIDParts contains CID_A?
idx2=find(cellfun(@(x) ~isempty(find(x(:,1)==B(1) & x(:,2)==B(2),1)),PHParts));%which CID in CIDParts contains CID_B?
if idx1==idx2
    cur=PHParts{idx1,1};
    out=cur(find(A(1)==cur(:,1),1)+1:find(B(1)==cur(:,1),1),:);
    out2=idx1;%linkIDs_full
else
    target=PHMap(idx2,1);
    linkIdx=[];
    cnt=1;
    while target~=idx1
        linkIdx=[target,linkIdx];
        target=PHMap(target,1);
        cnt=cnt+1;
    end
    out=[PHParts{idx1,1}(find(A(1)==PHParts{idx1,1}(:,1),1)+1:end,:);...
        cell2mat(PHParts(linkIdx,:));...
        PHParts{idx2,1}(1:find(B(1)==PHParts{idx2,1}(:,1),1),:)];
    out2=[idx1,linkIdx,idx2];%linkIDs_full
end
out=out(1:end-1,:);%CIDs in between.

%extract CIDs Btw CID_A and CID_B END
end
function [CIDs]=findCIDpartContainingCID_SPE(curSptLine,PHlines,PHMap)
%which PHline is the end of this spot line in
uit1=[];
cidEnd=curSptLine(end,2:3);
rngEng=find(cellfun(@(x) x(1,1)<=cidEnd(1) & x(end,1)>=cidEnd(1),PHlines));
for i=1:length(rngEng)
    idxEnd=find(PHlines{rngEng(i),1}(:,1)==cidEnd(1) & PHlines{rngEng(i),1}(:,2)==cidEnd(2), 1);
    if ~isempty(idxEnd)
        uit1=cat(2,uit1,rngEng(i));
    end
end
if length(uit1)~=1
    error('something wrong');
end

%which PHline is the begining of this spot line in
uit2=[];
cidBeg=curSptLine(1,2:3);
rngBeg=find(cellfun(@(x) x(1,1)<=cidBeg(1) & x(end,1)>=cidBeg(1),PHlines));

for i=1:length(rngBeg)
    idxBeg=find(PHlines{rngBeg(i),1}(:,1)==cidBeg(1) & PHlines{rngBeg(i),1}(:,2)==cidBeg(2), 1);
    if ~isempty(idxBeg)
        uit2=cat(2,uit2,rngBeg(i));
        selidxBeg=idxBeg;
    end
end
if length(uit2)~=1
    error('something wrong');
end
CIDs=PHlines{uit2,1}(selidxBeg:end,:);

% up the start phline and end phline
link=uit1;
if uit1~=uit2
    target=PHMap(uit1,1);
    cnt=0;
    while target~=uit2
        link=[target;link];
        target=PHMap(target,1);
        cnt=cnt+1;
        if target==0;break;end
        if cnt>50; error('linking failed!');end
    end
    link=[uit2;link];
end
if length(link)>1
    for iLk=2:length(link)%start from 2 intentionally
        CIDs=cat(1,CIDs,PHlines{link(iLk),1});
    end
end
end
function [uit]=getFullCIDs(endCIDInf,PHParts,PHMap,startFrm)
uit=[];

% get PH row index chain for the given index (i.e., endCIDInf(1,1)).
idxChain=endCIDInf(1,1);
if isnan(idxChain)
    return;
end
target=PHMap(endCIDInf(1,1),1);
while target~=0
    idxChain=[target,idxChain];
    target=PHMap(target,1);
end
% get CIDs from the beginning of the initiation till the end of its corr.
% div.
entireCID=cell2mat(PHParts(idxChain,1));
idx=find(entireCID(:,1)==startFrm,1);
uit=entireCID(idx:end,:);
end
function [uit]=fullizeMap(in)
uit=nan(size(in,1),1);
for i=1:size(in,1)
    [row,~]=find(in==i);
    if ~isempty(row)
        uit(i)=row;
    end
end
uit=[uit,in];
end
function [newSptLine]=completeSptLines(curSptLine,CIDs,Cperiod,i4disp)

idx1=find(CIDs(:,1)==curSptLine(1,2) & CIDs(:,2)==curSptLine(1,3));%must exist an index
idx2=find(CIDs(:,1)==curSptLine(end,2) & CIDs(:,2)==curSptLine(end,3));
if size(curSptLine,1)<=Cperiod
    if size(CIDs,1)-idx1+1<=Cperiod
        newPart=CIDs(idx2+1:size(CIDs,1),:);
        newSptLine=[curSptLine;[nan(size(newPart,1),1),newPart,nan(size(newPart,1),1)]];
    else
        newPart=CIDs(idx2+1:idx2+Cperiod-size(curSptLine,1),:);
        newSptLine=[curSptLine;[nan(size(newPart,1),1),newPart,nan(size(newPart,1),1)]];
    end
else
    newSptLine=curSptLine(1:Cperiod,:);
end

end
function [uit]=getDivCIDs(route,PHLines,map)
uit=zeros(2*(length(route)-1),2);
%for each pair of daughter CIDs, the interested one always goes to the first row. [interested1; disgarded1; interested2; disgarded2...]
for i=2:length(route)%from 2 intentionally
    uit(2*(i-2)+1,:)=PHLines{route(i)}(1,1:2);
    interestIdx=find(map(route(i-1),:)==route(i));
    disgardx=map(route(i-1),5-interestIdx);
    if ~isnan(disgardx)
        uit(2*(i-1),:)=PHLines{disgardx}(1,1:2);
    else
        uit(2*(i-1),:)=PHLines{route(i)}(1,1:2);
    end
end
uit=mat2cell(uit,ones(size(uit,1),1));
end

function []=drawFig(edgesX,X0,Y0,cutLow,cutUp,svPath,lgdLabel,mediaN)
% CVX=std(X0,[],"all","omitnan")./mean(X0,'all','omitnan');
% CVY=std(Y0,[],"all","omitnan")./mean(Y0,'all','omitnan');
x1=X0./mean(X0,'omitnan');
y=Y0./mean(Y0,'omitnan');
del=find(x1<cutLow | x1>cutUp | isnan(x1) | isnan(y));
x1(del)=[];
y(del)=[];
values=(edgesX(1:end-1)+edgesX(2:end))/2;
bin1=discretize(x1,edgesX,values);
scatY=zeros(length(values),1);
for i=1:length(values)
    scatY(i)=median(y(bin1==values(i)));
end
F=figure('Position',[100,100,200,200],'Visible','off');
hold on
scatter(x1,y,10,'filled','MarkerFaceAlpha',0.1,'MarkerFaceColor',[0, 0.4470, 0.7410]);
scatter(values,scatY,15,'filled','MarkerFaceAlpha',0.9,'MarkerFaceColor',[0.8500, 0.3250, 0.0980]);
axis([0 2 0 2])
ax=gca;
ax.FontSize=12;
xticks([0 1 2]);
yticks([0 1 2]);
X = [ones(size(values')) values'];
[b,~] = regress(scatY,X) ;

t=0.5:0.5:1.5;
yq=b(2)*t+b(1);
plot(t,yq,'-r','LineWidth',1.5);
hold off

text(0.1,0.2,[lgdLabel,'=',num2str(b(2),'%0.2f')],'FontSize',10);
text(0.1,0.37,['Media: ',mediaN],'FontSize',10);
text(0.1,0.54,['N=',num2str(length(x1))],'FontSize',10);

saveas(F,svPath);
close(F);
end
