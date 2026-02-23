function [parts_CID,Map]=IAB_simpDecipherLineagev2(lxg)

% Part 1 (formerlly imgPH_simpDecipherLineage)
% 1. Initialization
if ~isempty(find(~cellfun(@(x) isempty(x),lxg(:,end)), 1))
    lxg=[lxg ,cell(size(lxg,1),1)];
end

%remove entries that are not originated from Entry {1,1}.
A=lxg;%{1,nSCh}{1,2};%change here for testing. Take 5 for instance, lxg{1,5}{1,2}
Paridx={};%OUTPUT

init=A{1,1};
inTime=1;
stkDiv=[];%dynamic stack
stkCid=[];%dynamic stack
stkDrCo=[];%dynamic stack
stkDivT=[];%dynamic stack

curdrawCoor=[0,0];%dynamic stack [x,y]
curdivTime=0;%dynamic stack

isDivCyc=false;
isFirstRun=true;

LP=[];% Y of previous ending point.
LC=[];% [x0,y0,x1,y1].
cellID=[];
curLP=[];

while ~isempty(stkDiv) || isFirstRun
    %delete the first row of tempStk (that was used in previous round). Begin
    if ~isFirstRun
        stkDiv=stkDiv(2:end,:);
        stkDivT=stkDivT(2:end,:);
        stkDrCo=stkDrCo(2:end,:);
    end
    %delete the first row of tempStk. End
    for t=inTime:size(A,2)%for each time point
        if length(init)==2
            disp(['SChn',num2str(nSCh),' divided at 1st frame!']);
            error('out!')
        end

        if length(A{init,t})==1 && ~isDivCyc
            %<----stopping here indicates that the cell divided at the 1st frame. Use ManCheck to correct.
            init=A{init,t};
            LP=cat(1,LP,curdrawCoor(1,2));
            LC=cat(1,LC,[curdrawCoor(1,1),curdrawCoor(1,2),t,curdrawCoor(1,2)]);     %<---------=-----------=-----
            cellID=cat(1,cellID,[t,init]);
            curdrawCoor=[t,curdrawCoor(1,2)];
            continue;
        elseif length(A{init,t})==2 && ~isDivCyc
            isDivCyc=true;
            stkDiv=[stkDiv;[init,t]];
            stkCid=[stkCid;[init,t+1]];
            init=A{init,t}(1);
            LC=cat(1,LC,[curdrawCoor(1,1),curdrawCoor(1,2),t,curdrawCoor(1,2)]);     %<---------=-----------=-----
            cellID=cat(1,cellID,[t,init]);
            LP=cat(1,LP,curdrawCoor(1,2));
            stkDrCo=[stkDrCo;[t,curdrawCoor(1,2)]];
            stkDivT=[stkDivT;curdivTime];
            curLP=curdrawCoor(1,2);
            curdivTime=curdivTime+1;
            curdrawCoor=[t,curdrawCoor(1,2)+1/(2^curdivTime)];
        elseif length(A{init,t})==1 && isDivCyc
            init_old=init;
            init=A{init,t}(1);
            stkCid=[stkCid;[init,t+1]];
            LP=cat(1,LP,curLP);
            LC=cat(1,LC,[curdrawCoor(1,1),curdrawCoor(1,2),t,curdrawCoor(1,2)]);     %<---------=-----------=-----
            cellID=cat(1,cellID,[t,init_old]);
            curdrawCoor(1,:)=[t,curdrawCoor(1,2)];
            curLP=curdrawCoor(1,2);
            
        elseif length(A{init,t})==2 && isDivCyc
            stkDiv=[stkDiv;[init,t]];
            cellID=cat(1,cellID,[t,init]);
            init=A{init,t}(1);
            Paridx=[Paridx;{stkCid}];
            stkCid=[init,t+1];
            LP=cat(1,LP,curLP);
            LC=cat(1,LC,[curdrawCoor(1,1),curdrawCoor(1,2),t,curdrawCoor(1,2)]);     %<---------=-----------=-----
            
            stkDrCo=[stkDrCo;[t,curdrawCoor(1,2)]];
            stkDivT=[stkDivT;curdivTime];

            curdivTime=curdivTime+1;
            curdrawCoor(1,:)=[t,curdrawCoor(1,2)+1/(2^curdivTime)];
        elseif isempty(A{init,t})

            LP=cat(1,LP,curLP);
            LC=cat(1,LC,[curdrawCoor(1,1),curdrawCoor(1,2),t,curdrawCoor(1,2)]);     %<---------=-----------=-----
            cellID=cat(1,cellID,[t,init]);
            if isempty(stkDiv)
                break;
            end

            init=A{stkDiv(1,1),stkDiv(1,2)}(2);
            if init==-1% divLost case
                while init==-1 && ~isempty(stkDiv)
                    stkDiv(1,:)=[];
                    stkDrCo(1,:)=[];
                    stkDivT(1,:)=[];
                    if ~isempty(stkDiv)
                        init=A{stkDiv(1,1),stkDiv(1,2)}(2);
                        %add the follow to test BEG:
                        inTime=stkDiv(1,2)+1;
                        isDivCyc=true;
                        curLP=stkDrCo(1,2);
                        stkCid=[init,inTime];
                        curdrawCoor=[inTime-1,stkDrCo(1,2)-1/(2^(stkDivT(1)+1))];
                        curdivTime=stkDivT(1)+1;
                        %add the follow to test END.
                    end
                end
                break;
            end
            inTime=stkDiv(1,2)+1;
            isDivCyc=true;
            curLP=stkDrCo(1,2);
            stkCid=[init,inTime];
            curdrawCoor=[inTime-1,stkDrCo(1,2)-1/(2^(stkDivT(1)+1))];
            curdivTime=stkDivT(1)+1;
            break;
        end
    end
    isFirstRun=false;
end

% Part 2 (formerly drawTreeFromLL)

marker=true(size(LC,1),1);
marker(1)=false;
parts=[];
StY=[];
parts_CID={};

%1. retrieve horizontal lines [x0,y0,x1,y1];
%acquire: parts
curPart=LC(1,:);
curStY=LP(1);
curCID=cellID(1,:);

while ~isempty(curPart)
    for iIn=1:size(LC,1)
        if isequal(LC(iIn,1:2),curPart(1,3:4))
            curPart(1,3:4)=LC(iIn,3:4);
            marker(iIn)=false;
            curCID=cat(1,curCID,cellID(iIn,:));
            %continue;
        end
    end
    parts=cat(1,parts,curPart);
    StY=cat(1,StY,curStY);
    parts_CID=[parts_CID;{curCID}];
    newStRow=find(marker,1,"first");

    if ~isempty(newStRow)
        curPart=LC(newStRow,:);
        marker(newStRow)=false;
        curStY=LP(newStRow);
        curCID=cellID(newStRow,:);
    else
        curPart=[];
        curStY=[];
        curCID=[];
    end
end

%2. retrive vertical lines.
Map=zeros(size(parts,1),1);
for iH=2:size(parts,1)
    selRow=find(parts(:,4)==StY(iH));
    if ~isempty(selRow) && length(selRow)==1
        Map(iH)=selRow;
    else
        error('not unique in Y.')
    end
end

%3. reshuffle Y
[~,I]=sort(parts(:,2));
newParts=parts;
for iD=1:size(parts,1)
    newParts(I(iD),2)=iD;
    newParts(I(iD),4)=iD;
end

%4. draw lineage tree for PH cell contours
F=figure('visible','off');
set(F,'position',[10,10,900,1200]);
grid on
hold on
% invert the tree plot upside down.
Inv=max(newParts(:,2));
newParts(:,2)=Inv-newParts(:,2);
newParts(:,4)=newParts(:,2);

% Part 3 (formerly fullHLM)
Map=fullHLM(Map,newParts);

function [uit]=fullHLM(in,PHLines)
uit2=nan(size(in,1),2);
for iSL=1:size(in,1)
    xx=find(in==iSL);
    if isempty(xx)
        continue;
    end
    mom_Y=PHLines(iSL,2);
    if length(xx)==2
        if PHLines(xx(1),2)<PHLines(xx(2),2)
            xx=flipud(xx);
        end
        uit2(iSL,:)=xx';
    elseif length(xx)==1
        if PHLines(xx(1),2)>mom_Y
            uit2(iSL,1)=xx;
        else
            uit2(iSL,2)=xx;
        end

    end
end
uit=[in,uit2];
end
end