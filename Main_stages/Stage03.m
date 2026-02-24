clear;
clc;
%v4.3 @25.12.01 Ld is now defined as sum of two daughters.
%v4.2 @23.08.16 Added closest Ld/gen corresponding to each li in outA.
%v4.1 @22.11.21 Calculates Lb of daughters. This is to illustrate how Δt affects division site. Not very useful.
%v4.0 @22.11.21 Calculates deltaT in outA
%Output: resAv4-L2.mat

%% User specified region
% Option 1: process a specified experiment
% EXP='F:\RdnaA-seqAR-atc 1to3-M30 medium 20260130';int=3;WinK=3;WinL=10;%
% Stg3Neo_outs_with_closest_LdCore(EXP,int,WinK,WinL);

% Option 2: process a list of experiments
WinK=4;WinL=3;
load('G:\mChrRpath.mat');
path=cellfun(@(x,y) [x,y], repmat({'G:\'},size(Rpath,1),1),Rpath(:,1),'UniformOutput',false);
int=cell2mat(Rpath(:,3));
for i=1:length(path)%for each experiment
    Stg3Neo_outs_with_closest_LdCore(path{i},int(i),WinK,WinL);
end

%% Core function starts here.
function []=Stg3Neo_outs_with_closest_LdCore(EXP,int,WinK,WinL)

%WinK=2;%Unit is frame. Consider fast growth[5 min], slow growth [10 min]
%WinL=15;%Unit is frame. Consider fast growth[20-30 min]; slow is [75 min]
Cperiod=48;%unit frame C=48min for M18
% % int=5;%1 for fast growth, 5 for slow growth.
minIDT=20;  %unit: min. 10 for fast growth RDM+GLY
            %  20-30 for slow growth MOPS+GLU
mappingMethod=5;    %[4] fast growth
                    %[5] slow growth<-was always here
stEntry=1;

%% core
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

outA=cell(size(data,1),1);
outAData=cell(size(data,1),1);
outGHI=cell(size(data,1),1);
outGHIData=cell(size(data,1),1);
nAll=num2str(size(data,1));%str format
parfor i=1:size(data,1)%parfor each entry
    MPath=[EXP,'\M_res\p',num2str(data(i,1),'%0.2i'),'_lkgM.mat'];
    FPath=[EXP,'\Y_res\p',num2str(data(i,1),'%0.2i'),'_lkgY.mat'];
    disp(['Processing ',num2str(Entry(i)),'/',nAll]);
    try
        [outA{i},outAData{i},outGHI{i},outGHIData{i}]=genResult0611(MPath,FPath,data(i,:),mappingMethod,round(Cperiod/int),int,minIDT,WinK,WinL);
    catch
        disp(['P',num2str(data(i,1),'%0.2i'),'S',num2str(data(i,2),'%0.2i')]);
    end
end
outA=cat(1,outA{:});
outAData=cat(1,outAData{:});
outGHI=cat(1,outGHI{:});
outGHIData=cat(1,outGHIData{:});
% clear('outa','outb','outc','outd1','outd2','oute1','oute2','outf','outgh','outghi','outj1','outj1C','outj2','outk','outl1','outl2','i','MPath','FPath');
clear('mappingMethod','minIDT','WinK','WinL','Cperiod','Entry','stEntry','nAll','SLinesMap','SLines');
save([EXP,'\resV2simp.mat'],"int","outA","outAData","outGHI","outGHIData");
end

function [outA,outAData,outGHI,outGHIData]=genResult0611(MPath,FPath,SchN,mapMet,Cperiod,int,minIDT,WindowK,WindowL)
%% core - ini
load(FPath,'SLinesparts','SLinesMap','SLines');
load(MPath,'Infcellv2','lxg','PHparts','PHLines','PHMap');

% Initialization
InfCell=Infcellv2(:,SchN(2));
lxg=lxg{:,SchN(2)}{1,2};
SPTparts=SLinesparts{1,SchN(2)};
selSptMap=SLinesMap{1,SchN(2)};
CIDinParts=PHparts{1,SchN(2)};
PHparts=PHLines{1,SchN(2)};
HorLineMap=PHMap{1,SchN(2)};
selSpotLines=SLines{1,SchN(2)};
SpotEndFrm=max(cellfun(@(x) x(end,2), selSpotLines),[],1);
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
        disp(['last yellow line is nan: i=',num2str(i)]);
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
actualDNAmat=virtualDNAmat;
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
            if Spt2PHMap(i,3)>size(InfCell{Spt2PHMap(i,2),1}{1,2},1)
                disp(['Line',num2str(i),'Mom cell exceeding error']);
            else
                splitAndDivInfo(i,4)=InfCell{Spt2PHMap(i,2),1}{1,2}(Spt2PHMap(i,3),2);
            end
        else
            splitAndDivInfo(i,3)=nan;
            splitAndDivInfo(i,4)=nan;
        end
    end
    % UpCellLength - upper new born cell length
    if ~isnan(Spt2PHMap(i,4))
        if Spt2PHMap(i,6)>size(InfCell{Spt2PHMap(i,5),1}{1,2},1)
            disp(['Line',num2str(i),'Upper cell exceeding error']);
        else
            splitAndDivInfo(i,5)=InfCell{Spt2PHMap(i,5),1}{1,2}(Spt2PHMap(i,6),2);
        end
    end
    % LowCellLength - lower new born cell length
    if ~isnan(Spt2PHMap(i,7))
        if Spt2PHMap(i,9)>size(InfCell{Spt2PHMap(i,8),1}{1,2},1)
            disp(['Line',num2str(i),'Lower cell exceeding error']);
        else
            splitAndDivInfo(i,6)=InfCell{Spt2PHMap(i,8),1}{1,2}(Spt2PHMap(i,9),2);
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

% b) Calculate DNA content in each cell - create an additional table
% registering DNA content where mother and daughter spot lines are not
% connected.
    %identify unconnected links:

% obtain the ori lines:
[oriLines,oriMap]=getOriInfo(CIDinParts,HorLineMap,selSpotLines,full_selSptMap,SpotEndFrm);
[oriNummat]=getOriLines(oriLines,oriNummat);
[virtualDNAmat,~]=identifyUnconnectedLink(virtualSelSptLines,full_selSptMap,CIDinParts,HorLineMap,virtualDNAmat,oriNummat);
for iPH=1:length(CIDinParts)%for each PH part. <----change this back to 1.
    for iCID=1:size(CIDinParts{iPH,1},1)%for each cell in this PH part. <----change this back to 1.
        CID=CIDinParts{iPH,1}(iCID,:);%format of CID: [frame, cid];
        if CID(1)>size(virtualDNAmat,2)
            continue;
        end
        [ou,~]=ImgFLN_findDNABunchInCID(CID,virtualSelSptLines,full_selSptMap,Cperiod);
        if isnan(virtualDNAmat(CID(2),CID(1)))
            actualDNAmat(CID(2),CID(1))=ou;
        else
            actualDNAmat(CID(2),CID(1))=ou+virtualDNAmat(CID(2),CID(1));
        end
%         if isnan(oriNummat(CID(2),CID(1)))
%             oriNummat(CID(2),CID(1))=nori;
%         else
%             oriNummat(CID(2),CID(1))=nori+oriNummat(CID(2),CID(1));
%         end
    end
end

% c) For each SLine, assign a division number when this line is initiated.
divNumByGen=getDivNum4SLinesByGen(selSptMap,selSpotLines,CIDinParts,HorLineMap,Spt2PHMap);

%% genRes220611 - Summary(see docx for more info):
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


% ----- output A init Adder by ori#, gen, ter methods. 
% Format of outputA: 
% Col. 1~6: [Li/ori_m, Δii_U(by ori), Δii_L(by ori), Li/gen_m, Δii_U(by gen), Δii_L(by gen)] 
% Col. 7~14: [GRII_DaugU, GRII__DaugL, C-period (by spot line length, unit:Frm), Tii_U (unit:Frm), Tii_L (unit:Frm), ΔTii (min, initUp-initLow), Lbup,Lblow] 
% Col.15-17: [TInitdaughterUp,TinitdaughterLow,Tdiv](all are min)
% Col.18-19: [closest Ld/gen_m, gen_m] (Upon Chenli request)
% Col.20-21: [Lt, gen_t] real cell length at termination, generation number for Lt.@23.12.07
% Col.22-24: [t_start,t_end_U,t_end_L] Frame number for each spot line
outA=nan(size(selSptMap,1),24);
outAData=cell(size(selSptMap,1),2);
outAMap=nan(size(selSptMap,1),3);
for i=1:size(selSptMap,1)%for each spot line:
    col=nan(1,24);
    colData=cell(1,2);
    colMap=cell(1,2);
    if ~isnan(selSptMap(i,1)) || ~isnan(selSptMap(i,2))
        CIDm=selSpotLines{i,1}(1,2:3);%mother's initial CID (frame+CID).
        if CIDm(1)==1
            continue;
        end
        col(1)=InfCell{CIDm(1),1}{1,2}(CIDm(2),2)/oriNummat(CIDm(2),CIDm(1));   %Li/ori_m.
        col(4)=InfCell{CIDm(1),1}{1,2}(CIDm(2),2)/divNumByGen(i);               %Li/gen_m.
        
        if ~isnan(selSptMap(i,1))
            CIDdup=selSpotLines{selSptMap(i,1),1}(1,2:3);
            curInitCycCIDup=oriLines{i,1};
            if CIDdup(2)<=size(InfCell{CIDdup(1),1}{1,2},1)
                col(2)=(InfCell{CIDdup(1),1}{1,2}(CIDdup(2),2)/oriNummat(CIDdup(2),CIDdup(1)))*2-col(1);    %Δii_U(by ori)
                col(5)=(InfCell{CIDdup(1),1}{1,2}(CIDdup(2),2)/divNumByGen(selSptMap(i,1)))*2-col(4);       %Δii_U(by gen)
            end
            [col(7),LsRaw1]=getGRthrSpotLine(curInitCycCIDup,InfCell,CIDinParts,int);    %GR of upper ori
            col(10)=selSpotLines{selSptMap(i,1),1}(1,2)-selSpotLines{i,1}(1,2);     %Tii_U (unit:Frm)
            col(22)=CIDm(1);           %starting frame of current spot line
            col(23)=CIDdup(1)-1;       %ending frame of upper oriLine
            colData{1}=[LsRaw1/divNumByGen(i);col(5)+col(4)];
        end
        if ~isnan(selSptMap(i,2))
            CIDdlow=selSpotLines{selSptMap(i,2),1}(1,2:3);
            curInitCycCIDlow=oriLines{i,1};
            if CIDdlow(2)<=size(InfCell{CIDdlow(1),1}{1,2},1)
                col(3)=(InfCell{CIDdlow(1),1}{1,2}(CIDdlow(2),2)/oriNummat(CIDdlow(2),CIDdlow(1)))*2-col(1); %Δii_L(by ori)
                col(6)=(InfCell{CIDdlow(1),1}{1,2}(CIDdlow(2),2)/divNumByGen(selSptMap(i,2)))*2-col(4);      %Δii_L(by gen)
            end
            [col(8),LsRaw2]=getGRthrSpotLine(curInitCycCIDlow,InfCell,CIDinParts,int);   %GR of upper ori
            col(11)=selSpotLines{selSptMap(i,2),1}(1,2)-selSpotLines{i,1}(1,2);     %Tii_U (unit:Frm)
            col(22)=CIDm(1);           %starting frame of current spot line
            col(24)=CIDdlow(1)-1;      %ending frame of lower oriLine
            colData{2}=[LsRaw2/divNumByGen(i);col(6)+col(4)];
        end
        if ~isnan(selSptMap(i,1)) && ~isnan(selSptMap(i,2))%the current spot line has two descendants
            col(12)=int*(selSpotLines{selSptMap(i,1),1}(1,2)-selSpotLines{selSptMap(i,2),1}(1,2));%Δt in min. (T_initUp-T_initLow)
            col(15)=int*selSpotLines{selSptMap(i,1),1}(1,2);
            col(16)=int*selSpotLines{selSptMap(i,2),1}(1,2);
        end
        if ~isnan(selSptMap(i,1)) || ~isnan(selSptMap(i,2))%the current spot line has at least one descendant
            col(20)=InfCell{selSpotLines{i}(end,2),1}{1,2}(selSpotLines{i}(end,3),2);
            %any division during the current spot line?
            [~,PHroute]=extractCIDsBtwAB(selSpotLines{i}(1,2:3),selSpotLines{i}(end,2:3),CIDinParts,HorLineMap);
            col(21)=(2^(length(PHroute)-1))*divNumByGen(i);%
        end

        col(9)=selSpotLines{i}(end,2)-selSpotLines{i}(1,2);% C period raw
        if ~isnan(Spt2PHMap(i,4))
            col(13)=InfCell{Spt2PHMap(i,5),1}{1,2}(Spt2PHMap(i,6),2);%upper daughter Lb
            col(17)=int*Spt2PHMap(i,5);%upper daughter cell birth time in min
        end
        if ~isnan(Spt2PHMap(i,7))
            col(14)=InfCell{Spt2PHMap(i,8),1}{1,2}(Spt2PHMap(i,9),2);%lower daughter Lb
            col(17)=int*Spt2PHMap(i,8);%lower daughter cell birth time in min
        end
        if ~isequal(selSpotLines{i}(1,2:3),[1,1])
            %look for the division of PHpart that contains this CID:
            indCIDinParts=outs_findCIDpartContainingCID(selSpotLines{i}(1,2:3),CIDinParts);
            selCID_A=CIDinParts{indCIDinParts}(end,1:2);
            col(18)=InfCell{selCID_A(1),1}{1,2}(selCID_A(2),2)/divNumByGen(i);
            col(19)=divNumByGen(i);%the closest Ld to this CID.
        end
    end
    if sum(isnan(col))~=11
        outA(i,:)=col;
        outAData(i,:)=colData;
        outAMap(i,:)=[i,selSptMap(i,:)];
        % outA=cat(1,outA,col);
        % outAData=cat(1,outAData,colData);
        % outAMap=cat(1,outAMap,colMap);
    end
end

% ----- output B: %format:[idx, DNAc, L, NCR, DNAc/ori, L/ori];
% ----- output C: for each div, DNA content vs divSite:  Δ(divSite) vs Δ(DNAcont)
% ----- output D: NCratio (2D) and 1D 
% ----- output E -OK
% ----- output F -OK
% ----- output GHI C+D related things
outGHI=nan(size(Spt2PHMap,1),17); 
outGHIData=cell(size(Spt2PHMap,1),1);
% outGHIMap=nan(size(Spt2PHMap,1),3);
% Format of outputGHI: 
% Col. 1~6: [C+D(frm), rawC, Li, B-period, Li/genN, Ld]. (Note: Col.4 is now B(frm))
% Col. 7-11:[ori number at division,lambda*(C+D),lambda,Li*Rat,Ti].
% Col. 12-15:[spot line index in Spt2PHMap, Lb of the corresponding division, number of division during C+D, previous Ld].
% Col. 16-17:[t_id_start, t_id_end] Starting and ending frame of this ID.
%to plot C+D adder, use C10 and C6-C10
del=[];
for i=1:size(Spt2PHMap,1)%for each spot line
    %omit if this cell does not have a PH line to attach.
    if isnan(Spt2PHMap(i,1))
        del=[del;i];
        continue;
    end
    %omit if this spot line is an initial thread.
    if selSpotLines{i,1}(1,2)==1
        del=[del;i];
        continue;
    end
    %check if the corresponding PH line has division.
    if ~isnan(HorLineMap(Spt2PHMap(i,1),2)) && ~isnan(HorLineMap(Spt2PHMap(i,1),3))% changed from || to &&
        lastCID=CIDinParts{Spt2PHMap(i,1),1}(end,:);
        daugCIDs=[Spt2PHMap(i,5:6);Spt2PHMap(i,8:9)];
        startCID=selSpotLines{i}(1,2:3);
        [~,PHroute]=extractCIDsBtwAB(startCID,lastCID,CIDinParts,HorLineMap);
        CIDidDivEvent=getDivCIDs(PHroute,CIDinParts,HorLineMap);%the CIDidDivEvent is an array of cells, each cell is {[CID]}
        CIDrat=cellfun(@(x) InfCell{x(1),1}{1,2}(x(2),2),CIDidDivEvent);
        divRats=CIDrat(1:2:(length(CIDrat)-1))./(CIDrat(1:2:(length(CIDrat)-1))+CIDrat(2:2:length(CIDrat)));
        Li=InfCell{selSpotLines{i,1}(1,2),1}{1,2}(selSpotLines{i,1}(1,3),2);
        CD=Spt2PHMap(i,2)-selSpotLines{i,1}(1,2)+1;%time of C+D UPDATED
        startCIDofCorrPHpart=CIDinParts{Spt2PHMap(i,1),1}(1,:);%starting CID of the corresponding PH part.
        if isequal(startCIDofCorrPHpart,[1,1]);startCIDofCorrPHpart=nan(1,2);end
        [CIDcd,~]=extractCIDsBtwAB(selSpotLines{i,1}(1,2:3),lastCID,CIDinParts,HorLineMap);
        [lambda,LsGHI]=getGRthrSpotLine([selSpotLines{i,1}(1,2:3);CIDcd;lastCID],InfCell,CIDinParts,int);
        [idxLine,~]=find(HorLineMap(:,2:3)==Spt2PHMap(i,1));
        if HorLineMap(idxLine,1)~=0
            prevCID=CIDinParts{idxLine}(end,:);
            prevLd=InfCell{prevCID(1)}{1,2}(prevCID(2),2);
        else
            prevLd=nan;
        end
        Ld=InfCell{daugCIDs(1,1),1}{1,2}(daugCIDs(1,2),2)+InfCell{daugCIDs(2,1),1}{1,2}(daugCIDs(2,2),2);% upgraded
        uix=[CD,...%1. C+D
            size(selSpotLines{i,1},1),...%2. rawC(frm)
            Li,...%3. Li
...Li./oriNummat(selSpotLines{i,1}(1,3),selSpotLines{i,1}(1,2)),...%Li/ori#%removed on 24.03.14
            startCID(1)-startCIDofCorrPHpart(1),...%4. B-period (frm). Added on 24.03.14
            Li./divNumByGen(i)...%5. Li/genN
            Ld,...%6. Ld. was InfCell{lastCID(1),1}{1,2}(lastCID(2),2)
            oriNummat(lastCID(2),lastCID(1)),...%7. ori#@division
            lambda*CD,lambda,...%8-9. lambda*(C+D), lambda(h-1)
            Li*prod(divRats),...%10. Li*Rats
            selSpotLines{i,1}(1,2),i,...%11-12. Ti, index of spotline in Spt2PHmap
            InfCell{CIDinParts{Spt2PHMap(i,1)}(1,1)}{1,2}(CIDinParts{Spt2PHMap(i,1)}(1,2),2),...%13. Lb of the corresponding division
            length(divRats),...%14. number of division during C+D,
            prevLd,...%15. previous Ld
            startCID(1),lastCID(1)];%16-17. Starting and ending frame of this ID.
        outGHI(i,:)=uix;
        outGHIData{i}=[LsGHI*prod(divRats);Ld];
        % outGHIMap(i,:)=[i,Spt2PHMap(i,:)];
    else
        del=[del;i];
    end
end

delA=isnan(outA(:,1));
delR=isnan(outAMap);
outAMap=arrayfun(@(x) ['P',num2str(SchN(1),'%0.2i'),'S',num2str(SchN(2),'%0.2i'),'N',num2str(x,'%0.2i')], outAMap,'UniformOutput',false);
outAMap(delR)={''};
outGHIMap=outAMap;
outAMap(delA,:)=[];
outAData(delA,:)=[];
outA(delA,:)=[];
outAData=[outAMap,outAData];

outGHI(del,:)=[];
outGHIData(del,:)=[];
outGHIMap(del,:)=[];
outGHIData=[outGHIMap,outGHIData];
% ----- output J ori fire in div cyc.
% outJ1=[]; %format:see below.
% outJ2={}; %format:see below.
%[outJ1, outJ1Cell, outJ2]=outputtingJ(oriLines, oriMap, CIDinParts,HorLineMap,selSpotLines,full_selSptMap);
%format of outJ1: [enterOriC number, init 0 number percentage, init once number percentage, ...2%, ...3% ]
%format of outJ2: [PHidxOfMomPHpart, NumOfFiringInMomPHpart, PHidxOfCurPHpart, NumOfFiringInCurPHpart]
% ----- output K Moving Window NCR: [NCR, T2nI, T2nD]
% ----- output L NCR before div and init: NCRDiv and NCRInit (N x Window double)

end
%% nested function
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
function [uit,Ls]=getGRthrSpotLine(CIDin,InfCell,CIDinParts,int)
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
    Ls(i)=2^(cnt)*InfCell{CIDin(i,1),1}{1,2}(CIDin(i,2),2);
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
function [out]=getInfoByCID(CIDs,actualDNAmat,InfCell,oriNummat)
out=zeros(size(CIDs,1),5);
for i=1:size(CIDs,1)
   out(i,1)=actualDNAmat(CIDs(i,2),CIDs(i,1));
   out(i,2)=InfCell{CIDs(i,1),1}{1,2}(CIDs(i,2),2);
   out(i,3)=out(i,1)./out(i,2);
   out(i,4)=out(i,1)./oriNummat(CIDs(i,2),CIDs(i,1));
   out(i,5)=out(i,2)./oriNummat(CIDs(i,2),CIDs(i,1));
   %format:[DNAc, L, NCR, DNAc/ori, L/ori];
end
end
function [meanNCR]=getMeanNCRatiothrSpotLine(CIDs,NCRmat)
NCR=nan(size(CIDs,1),1);
for i=1:size(CIDs,1)
    NCR(i)=NCRmat(CIDs(i,2),CIDs(i,1));
end
meanNCR=mean(NCR);
end
function [uit1,uit2]=getPHPartData(PHparts,PHMap,InfCell,NCratmat,actDNAmat,int,minFrm)
%output1: [DNA content at birth, NCratio at birth, meanNCR, medianNCR, maxNCR, GR, IDT(frm)]
%output2: {[DNA content, NCratio, SpeGR]}
uit1=zeros(length(PHparts),7);
uit2=cell(length(PHparts),1);
for i=2:length(PHparts)%For each PH part, omit the 1st one.
    curCIDs=PHparts{i,1};
    if size(curCIDs,1)<minFrm% omit too short PH parts
        continue;
    end
    if isnan(PHMap(i,2)) && isnan(PHMap(i,3))
        continue;
    end
    Ls=zeros(size(curCIDs,1),1);
    DnaCont=zeros(size(curCIDs,1),1);
    NCrat=zeros(size(curCIDs,1),1);

    for j=1:size(curCIDs,1)
        Ls(j)=InfCell{curCIDs(j,1),1}{1,2}(curCIDs(j,2),2);
        if curCIDs(j,1)<=size(actDNAmat,2)
            DnaCont(j)=actDNAmat(curCIDs(j,2),curCIDs(j,1));
            NCrat(j)=NCratmat(curCIDs(j,2),curCIDs(j,1));
        else
            DnaCont(j)=nan;
            NCrat(j)=nan;
        end

    end
    x=(1:1:length(Ls))'*int/60;
    [f,~] = fit(x,Ls,'exp1');
    uit1(i,:)=[actDNAmat(curCIDs(1,2),curCIDs(1,1)),NCratmat(curCIDs(1,2),curCIDs(1,1)),mean(NCrat,'all','omitnan'),median(NCrat,'all','omitnan'),max(NCrat,[],'all','omitnan'),f.b,length(Ls)];
    [SpeGR,~]=getSpeAbsER(Ls,int);
    uit2{i,1}=[DnaCont,NCrat,[nan;SpeGR{1,1}]];
end
end
function [spuit,abuit]=getSpeAbsER(L,int)
spuit=zeros(length(L)-1,1);
abuit=zeros(length(L)-1,1);

dt=int/60;

for i=1:length(L)-1
    spuit(i)=(log(L(i+1))-log(L(i)))/dt;
    abuit(i)=(L(i+1)-L(i))/dt;
    
end
spuit={spuit};
abuit={abuit};

end
function [vDNAmat,vOriNummat]=identifyUnconnectedLink(lines,map,PHParts,PHMap,vDNAmat,vOriNummat)
%NOTE: THIS WORKS FOR THE FAST GR ONLY!!!
%1. find out those unconnected mother-daughter pairs.
uit=[];%row number. format: [iRow in spot map, iCol in spot map].
uitC=[];
for iLine=1:length(lines)%for each spot line
    curEnd=lines{iLine,1}(end,2);
    if ~isnan(map(iLine,2))
        curDaugBeg1=lines{map(iLine,2),1}(1,2);
        if curDaugBeg1>curEnd+1
            uit=[uit;[iLine,2]];% current spot line and its upper daughter are not connected
        else
            uitC=[uitC;[iLine,2]];% current spot line and its upper daughter are connected
        end
    end 
    if ~isnan(map(iLine,3))
        curDaugBeg2=lines{map(iLine,3),1}(1,2);
        if curDaugBeg2>curEnd+1
            uit=[uit;[iLine,3]];% current spot line and its lower daughter are not connected
        else
            uitC=[uitC;[iLine,3]];% current spot line and its lower daughter are connected
        end
    end
    if isnan(map(iLine,2)) && isnan(map(iLine,3))% if current has no daughter:
        sptEndCID=lines{iLine,1}(end,2:3);%then, after its ending CID, all cells has to assign a DNAc
        idx=cellfun(@(x) ~isempty(find(x(:,1)==sptEndCID(1) & x(:,2)==sptEndCID(2),1)), PHParts);%PHpart that contains this CID.
        if sptEndCID(1)~=PHParts{idx,1}(end,1)% if this CID is not the last cell of this PH part:
            idx1=find(PHParts{idx,1}(:,1)==sptEndCID(1) & PHParts{idx,1}(:,2)==sptEndCID(2),1);%which cell is this CID in this PH part?
            curOut=PHParts{idx,1}(idx1+1:end,:);%CIDs in this PH part that to be assigned DNAc.
            for icO=1:size(curOut,1)% for each of these CIDs:
                if size(vDNAmat,2)>=curOut(icO,1)
                    if isnan(vDNAmat(curOut(icO,2),curOut(icO,1)))%add 2 to this CID.
                        vDNAmat(curOut(icO,2),curOut(icO,1))=2;
                    else
                        vDNAmat(curOut(icO,2),curOut(icO,1))=vDNAmat(curOut(icO,2),curOut(icO,1))+2;
                    end
                    if isnan(vOriNummat(curOut(icO,2),curOut(icO,1)))
                        vOriNummat(curOut(icO,2),curOut(icO,1))=2;
                    else
                        vOriNummat(curOut(icO,2),curOut(icO,1))=2+vOriNummat(curOut(icO,2),curOut(icO,1));
                    end
                end
            end
        end
        % check if this cell has decendants. If yes, fill it with half DNAc
        desIdx=PHMap(idx,2:3);desIdx(isnan(desIdx))=[];
        if ~isempty(desIdx)
            curOut1=cell2mat(PHParts(desIdx));%CIDs in this PH part that to be assigned DNAc.
            DNAcInDes=vDNAmat(PHParts{idx,1}(end,2),PHParts{idx,1}(end,1))/2;%DNAc to be filled in CIDs of descendants'.
            for icO=1:size(curOut1,1)
                if size(vDNAmat,2)>=curOut1(icO,1)
                    if isnan(vDNAmat(curOut1(icO,2),curOut1(icO,1)))%add 2 to this CID.
                        vDNAmat(curOut1(icO,2),curOut1(icO,1))=DNAcInDes;
                    else
                        vDNAmat(curOut1(icO,2),curOut1(icO,1))=vDNAmat(curOut1(icO,2),curOut1(icO,1))+DNAcInDes;
                    end
                end
            end
        end

    end
end

%2. for those unconnected pairs:
Reg=[];
Reg3=[];
% uit([2,6,11],:)=[];%<-------------------------testting!!!!!!!!!!!!!!!!!!!!!!
for iUc=1:size(uit,1)
%2.1 find linkage between the mother ending CID to daughter starting CID.
    %lines{uit(iUc,1)}(end,2:3);%mother ending cid
    %lines{map(uit(iUc,1),uit(iUc,2))}(1,2:3);%daughter starting cid
    CIDME=lines{uit(iUc,1)}(end,2:3);
    CIDDS=lines{map(uit(iUc,1),uit(iUc,2))}(1,2:3);

%2.2 extract CIDs in between and complete the virtual DNA contant matrix;
    [selCIDs,~]=extractCIDsBtwAB(CIDME,CIDDS,PHParts,PHMap);
    for iCID=1:size(selCIDs,1)
        if isnan(vDNAmat(selCIDs(iCID,2),selCIDs(iCID,1)))
            vDNAmat(selCIDs(iCID,2),selCIDs(iCID,1))=1;
        else
            vDNAmat(selCIDs(iCID,2),selCIDs(iCID,1))=1+vDNAmat(selCIDs(iCID,2),selCIDs(iCID,1));
        end
        if isnan(vOriNummat(selCIDs(iCID,2),selCIDs(iCID,1)))
            vOriNummat(selCIDs(iCID,2),selCIDs(iCID,1))=1;
        else
            vOriNummat(selCIDs(iCID,2),selCIDs(iCID,1))=1+vOriNummat(selCIDs(iCID,2),selCIDs(iCID,1));
        end
    end

%2.3 register for singlet:
    if isempty(Reg)
        Reg=uit(iUc,:);
    else
        if uit(iUc,1)~=Reg(1,1)
            if ~isempty(uitC)
                if isempty(find(uitC(:,1)==Reg(1) & uitC(:,2)==(5-Reg(2)),1))
                    Reg3=[Reg3;Reg];
                end
            end
            Reg=uit(iUc,:);
        else
            Reg=[Reg;uit(iUc,:)];
            if isequal(Reg(:,2),[2;3]) || isequal(Reg(:,2),[3;2])
                Reg=[];
            else
                error('very wrong');
            end
        end
    end
end
idxN=find((~isnan(map(:,2)) & isnan(map(:,3))) | (isnan(map(:,2)) & ~isnan(map(:,3))));
if ~isempty(Reg3)
    addIdx=idxN(~ismember(idxN,Reg3(:,1)));
    for iAi=1:length(addIdx)
        Reg3=[Reg3;[addIdx(iAi),find(~isnan(map(addIdx(iAi),2:3)))+1]];
    end
else
    for iAi=1:length(idxN)
        Reg3=[Reg3;idxN(iAi),find(~isnan(map(idxN(iAi),2:3)))+1];
    end
end

%3 for those singlets in the uit:
for iReg=1:size(Reg3,1)
    CIDME=lines{Reg3(iReg,1)}(end,2:3);
    CIDDS=lines{map(Reg3(iReg,1),Reg3(iReg,2))}(1,2:3);
    idx1=find(cellfun(@(x) ~isempty(find(x(:,1)==CIDME(1) & x(:,2)==CIDME(2),1)),PHParts));%which CIDidx in CIDParts contains CID_mom_end?
    idx2=find(cellfun(@(x) ~isempty(find(x(:,1)==CIDDS(1) & x(:,2)==CIDDS(2),1)),PHParts));%which CIDidx in CIDParts contains CID_daug_beg?
    if idx1==idx2
        CIDDE=lines{map(Reg3(iReg,1),Reg3(iReg,2))}(end,2:3);
        idx3=find(cellfun(@(x) ~isempty(find(x(:,1)==CIDDE(1) & x(:,2)==CIDDE(2),1)),PHParts));%which CIDidx in CIDParts contains CID_daug_end?
        if idx3==idx2
            cur=PHParts{idx2,1};
            out=cur(find(CIDME(1)==cur(:,1),1)+1:end,:);
        else
            target=PHMap(idx3,1);
            linkIdx=idx3;
            cnt=1;
            while target~=idx1
                linkIdx=[target,linkIdx];
                target=PHMap(target,1);
                cnt=cnt+1;
            end
            linkIdx=[idx1,linkIdx];
            pX=PHMap(linkIdx(1),2:3);
            pX(pX==linkIdx(2))=[];
            cur=PHParts{idx1,1};
            if ~isnan(pX)%<-------------changed on 22.06.23
                out=[cur(find(CIDME(1)==cur(:,1),1)+1:end,:);PHParts{pX}];%<------original line (others were added on 22.06.23)
            else
                out=cur(find(CIDME(1)==cur(:,1),1)+1:end,:);
            end
        end
    else
        target=PHMap(idx2,1);
        linkIdx=idx2;
        cnt=1;
        while target~=idx1 && target~=0
            linkIdx=[target,linkIdx];
            target=PHMap(target,1);
            cnt=cnt+1;
        end
        linkIdx=[idx1,linkIdx];
        pX=PHMap(linkIdx(1),2:3);
        pX(pX==linkIdx(2))=[];
        cur=PHParts{idx1,1};
        if ~isnan(pX)
            out=[cur(find(CIDME(1)==cur(:,1),1)+1:end,:);PHParts{pX}];
        else
            out=[cur(find(CIDME(1)==cur(:,1),1)+1:end,:)];
        end
    end
%refill the virtual DNA content map
    for iCID=1:size(out,1)
        if isnan(vDNAmat(out(iCID,2),out(iCID,1)))
            vDNAmat(out(iCID,2),out(iCID,1))=1;
        else
            vDNAmat(out(iCID,2),out(iCID,1))=1+vDNAmat(out(iCID,2),out(iCID,1));
        end
        if isnan(vOriNummat(out(iCID,2),out(iCID,1)))
            vOriNummat(out(iCID,2),out(iCID,1))=1;
        else
            vOriNummat(out(iCID,2),out(iCID,1))=1+vOriNummat(out(iCID,2),out(iCID,1));
        end
    end
end

%4 at the initiation cell, exclude newly formed origin while calculating:
for iLine=1:length(lines)%for each spot line
    if lines{iLine,1}(1,2)==1
        continue;
    end

    if ~isnan(vOriNummat(lines{iLine,1}(1,3),lines{iLine,1}(1,2)))
        vOriNummat(lines{iLine,1}(1,3),lines{iLine,1}(1,2))=vOriNummat(lines{iLine,1}(1,3),lines{iLine,1}(1,2))-1;
    else
        vOriNummat(lines{iLine,1}(1,3),lines{iLine,1}(1,2))=-1;
    end
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
        target=PHMap(target,1);%if error here, check PHpart between these two CIDs [A, B]
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
function [NR]=getOriLines(oLines,oriNummat)%the output 'NR' is for per-ori calculations
oriNummat=zeros(size(oriNummat));
debOriNummat=oriNummat;
for i=1:size(oLines,1)
    for j=1:2%intentionally 1:2
        curOL=oLines{i,j};
        for k=1:size(curOL,1)%from 2 intentionally.
            oriNummat(curOL(k,2),curOL(k,1))=oriNummat(curOL(k,2),curOL(k,1))+1;
        end
    end
    debOriNummat(oLines{i,1}(1,2),oLines{i,1}(1,1))=debOriNummat(oLines{i,1}(1,2),oLines{i,1}(1,1))+1;
end
debOriNummat(1,1)=0;
NR=oriNummat-debOriNummat;
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
