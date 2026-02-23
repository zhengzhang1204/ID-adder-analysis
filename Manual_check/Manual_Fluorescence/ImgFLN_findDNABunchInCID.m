function [uit,uitOri]=ImgFLN_findDNABunchInCID(CID,SptLines,SptMap,Cperiod)

initThreads=find(isnan(SptMap(:,1)));
%1 retrieve IDXs of spot lines that exist in this CID
IDXs=find(cellfun(@(x) ~isempty(find(x(:,2)==CID(1) & x(:,3)==CID(2),1)),SptLines));
stFrm=cellfun(@(x) find(x(:,2)==CID(1) & x(:,3)==CID(2),1),SptLines,'UniformOutput',false);


%2 check if there are mother-daughter pairs in these idxs.
if length(IDXs)==1
    C=IDXs;
else
    C=nchoosek(IDXs,2);
end
MDcollector={};%mother-daughter collector
bunchOut=[];
mono=IDXs;
for iC=1:size(C,1)
    if isMomDaug(C(iC,:),SptMap)
        MDcollector=[MDcollector;{C(iC,:)}];
        bunchOut=[bunchOut,C(iC,:)];
    end
end
bunchOut=unique(bunchOut);
for iB=1:length(bunchOut)
    mono(mono==bunchOut(iB))=[];
end

%3 calculate the DNA content in this new born cell:
uit=0;
uitOri=0;
%mono spot lines:
for iM=1:length(mono)
    if ismember(mono(iM),initThreads)
        endFrm=SptLines{mono(iM),1}(end,2);
        if endFrm>Cperiod
            uit=uit+2-(endFrm-stFrm{mono(iM)})/endFrm;
        else
            uit=uit+2-(endFrm-stFrm{mono(iM)})/Cperiod;
        end
        uitOri=uitOri+2;
    else
        if size(SptLines{mono(iM),1},1)>=Cperiod%added this if-end here on 22.06.23
                                                %if the given C period is too short, use Cperiod instead.
            uit=uit+stFrm{mono(iM)}/size(SptLines{mono(iM),1},1)+1;%original
            uitOri=uitOri+2;%original
        else
            uit=uit+stFrm{mono(iM)}/Cperiod+1;
            uitOri=uitOri+2;
        end
    end
end
if isempty(MDcollector)
    return;
end

%DNA bunch:
BunchMap=[];% collect all mother-daugter pair:
for iB=1:length(MDcollector)
    curB=MDcollector{iB,1};
    if ~isempty(BunchMap)
        idx=find(BunchMap(:,1)==curB(1),1);
        if ~isempty(idx)
            BunchMap(idx,3)=curB(2);
        else
            BunchMap=[BunchMap;[curB,nan]];
        end
    else
        BunchMap=[curB,nan];
    end
end

%see how many bunches exist in this map?
[~,idx1]=sort(BunchMap(:,1),'ascend');
BunchMap=BunchMap(idx1,:);
occupy=zeros(size(BunchMap,1),1);
Bunch={};
cnt=1;
for iBM=1:size(BunchMap,1)
    idx2=find(BunchMap(:,2)==BunchMap(iBM,1) | BunchMap(:,3)==BunchMap(iBM,1));
    if ~isempty(idx2)
        cont=Bunch{occupy(idx2)};
        cont=[cont;BunchMap(iBM,:)];
        Bunch{occupy(idx2)}=cont;
        occupy(iBM,1)=occupy(idx2);
    else
        occupy(iBM,1)=cnt;
        cnt=cnt+1;
        Bunch=[Bunch;{BunchMap(iBM,:)}];
    end
end

for iBunch=1:size(Bunch,1)
    curBLength=size(Bunch{iBunch},1);
    curBInfo=Bunch{iBunch};
    if curBLength==1%3Y-shaped
        A1=ratioOfCompletion(CID,SptLines{curBInfo(1)});%mother spot line
        A2=ratioOfCompletion(CID,SptLines{curBInfo(2)});%daughter
        AO1=0;
        AO2=2;
        if isnan(curBInfo(3))%2Y-shaped
            A3=0;
            AO3=1;
        else
            A3=ratioOfCompletion(CID,SptLines{curBInfo(3)});%daughter
            AO3=2;
        end
        uit=uit+1+A1+A2+A3;
        uitOri=uitOri+AO1+AO2+AO3;
    elseif curBLength==2%5Y-shaped
        %reshape the bunch so that the mother line is the first row.
        idx3=find(curBInfo(2,2:3)==curBInfo(1,1),1);
        if ~isempty(idx3)
            curBInfo=flipud(curBInfo);
        end
        idx3=find(curBInfo(1,:)==curBInfo(2,1),1);
        if idx3==3
            curBInfo(1,2:3)=fliplr(curBInfo(1,2:3));
        end
        A1=ratioOfCompletion(CID,SptLines{curBInfo(1,1)});%mother
        A2=ratioOfCompletion(CID,SptLines{curBInfo(1,2)});%daughter which has daughters
        AO1=0;
        AO2=0;
        if isnan(curBInfo(1,3))
            A3=0;
            AO3=1;
        else
            A3=ratioOfCompletion(CID,SptLines{curBInfo(1,3)});%daughter which does not have a daughter
            AO3=2;
        end
        A4=ratioOfCompletion(CID,SptLines{curBInfo(2,2)});%granddaughter
        AO4=2;
        if isnan(curBInfo(2,3))
            A5=0;
            AO5=1;
        else
            A5=ratioOfCompletion(CID,SptLines{curBInfo(2,3)});%granddaughter
            AO5=2;
        end
        uit=uit+1+A1+A2+A3+A4+A5;
        uitOri=uitOri+AO1+AO2+AO3+AO4+AO5;
    elseif curBLength==3%7Y-shaped
        idx3=find(curBInfo(:,2:3)==curBInfo(1,1),1);
        idx4=find(curBInfo(:,2:3)==curBInfo(2,1),1);
        if ~isempty(idx3) && ~isempty(idx4)
            curBInfo=curBInfo([3;1;2],:);
        elseif isempty(idx4)
            curBInfo=curBInfo([2;1;3],:);
        end
        idx3=find(curBInfo(1,:)==curBInfo(2,1),1);
        if idx3==3
            curBInfo(2:3,:)=flipud(curBInfo(2:3,:));
        end
        A1=ratioOfCompletion(CID,SptLines{curBInfo(1,1)});%mother
        A2=ratioOfCompletion(CID,SptLines{curBInfo(1,2)});%daughter A
        A3=ratioOfCompletion(CID,SptLines{curBInfo(1,3)});%daughter B
        A4=ratioOfCompletion(CID,SptLines{curBInfo(2,2)});%granddaughter A1
        A5=ratioOfCompletion(CID,SptLines{curBInfo(2,3)});%granddaughter A2
        %error here means: [broken SL],run: SptLines{curBInfo(1,1)} to see which SL errors
        A6=ratioOfCompletion(CID,SptLines{curBInfo(3,2)});%granddaughter B1
        A7=ratioOfCompletion(CID,SptLines{curBInfo(3,3)});%granddaughter B2
        uit=uit+1+A1+A2+A3+A4+A5+A6+A7;
        uitOri=uitOri+8;
    end
end
end
function [uit]=ratioOfCompletion(CID,curSpotLine)
idx=find(curSpotLine(:,2)==CID(1),1);
uit=(idx-1)/size(curSpotLine,1);
if isempty(idx)
    error('wrong')
end
end
function [uit]=isMomDaug(in,Map)
uit=false;
if length(in)==1
    return;
end
%1. if in(1) is mother:
rr1=find(Map(in(1),:)==in(2),1);
if ~isempty(rr1)
    uit=true;
end

%2. if in(1) is daughter:
rr2=find(Map(in(2),:)==in(1),1);
if ~isempty(rr2)
    uit=true;
end
end