function [lkg,lng,Ulist,Umsk,OMsk,Pole]=IAB_Length2LNG(A,Umsk,idxList,OMsk)
%version v1.0 @23.02.10     Remove the top divLost case when number of cells in this side channel is >2.
R={};%cell(length(A),1);
for k=1:length(A)
    if isempty(A{k,1})
        break;
    else
        R=[R,A(k,1)];
    end
end
%linking begin, build up lkg,lng
endFrame=length(R);
maxRow=max(cellfun(@(x) length(x),R));
lkg=cell(maxRow,endFrame);
lng=cell(maxRow,endFrame);
Pole=zeros(maxRow,endFrame);
Pole=getPoleMarker(Pole,'0',R{1,1},1);
Ulist=[];

for curFrame=1:endFrame-1%for each frame
    ML=R{1,curFrame};
    DL=R{1,curFrame+1};
    [tdt,tnt,~]=DL_ZZlinkMoMaMain(ML,DL);
    lkg(1:length(tdt),curFrame)=tdt;%output linkage
    lng(1:length(tnt),curFrame)=tnt;%output lineage
    Pole(:,curFrame:curFrame+1)=getPoleMarker(Pole(:,curFrame),tnt,'0',2);
    curMsk=false(256,32);
    for iLKG=1:length(tdt)%for each regions in this side channel at this frame
       if strcmp(tdt{iLKG}(3:end),'tooFastGR') || (strcmp(tdt{iLKG}(3:end),'divLost') && iLKG~=length(tdt))
           Ulist=cat(1,Ulist,[curFrame,iLKG]);
           idxL=idxList{curFrame,1}{iLKG,1};
           curMsk(idxL)=true;
       elseif (strcmp(tdt{iLKG}(3:end),'divLost') && iLKG==length(tdt)) && length(tdt)>2%added this block on 23.02.10. 
           tempCurOMsk=OMsk(:,curFrame*32+1:curFrame*32+32);
           tempCurOMsk(idxList{curFrame+1,1}{lng{iLKG,curFrame}(1)})=false;
           OMsk(:,curFrame*32+1:curFrame*32+32)=tempCurOMsk;
           Pole(lng{iLKG,curFrame}(1),curFrame+1)=nan;
           
           A{curFrame+1,1}(lng{iLKG,curFrame}(1))=[];
           R{1,curFrame+1}(lng{iLKG,curFrame}(1))=[];
           idxList{curFrame+1,1}(lng{iLKG,curFrame}(1))=[];
           originalTxt=lkg{iLKG,curFrame};
           lkg{iLKG,curFrame}=[originalTxt(1:2),'lost'];
           lng{iLKG,curFrame}=[];
       end
    end
    Umsk(:,(curFrame-1)*32+1:(curFrame)*32)=curMsk;
end
%linking end

%Mark the first frame if division happens
idDiv=find(cellfun(@(x) length(x), lng(:,1))==2);
if ~isempty(idDiv)
    curMsk=false(256,32);
    for idv=1:length(idDiv)
        Ulist=cat(1,Ulist,[1,idDiv(idv)]);
        idxL=idxList{1,1}{idDiv(idv),1};
        curMsk(idxL)=true;
    end
    Umsk(:,1:32)=curMsk;
end




%nested fucntions：
function [tmplkg,tmplng,uit]=DL_ZZlinkMoMaMain(mothLen,daugLen)%MAIN
%init
global daughterOne count
tmplkg={};
tmplng={};
uit=[];
daughterOne=0;
daugID=1:1:length(daugLen);
for j=1:length(mothLen)%for each cell on mother frame
    
    if isempty(daugLen)
        tmplkg=[tmplkg;{strcat(num2str(j,'%0.2d'),'lost')}];
        tmplng=[tmplng;{}];
        break;%break means to start comparing next mother (j+1)
    end
	count=1;
	tmp=[];
    while count<=length(daugID)
        [next,txt,msg,chs]=ZZMoMaComp(mothLen(j,:),daugLen(count,:));%CORE.
        count=count+next;
        if strcmp(txt,'divided')
            if daughterOne~=0 && count>length(daugID)
                txt='divLost';
            else
                tmp=chs;
            end
        elseif strcmp(txt,'DivVerified')
            chs=[tmp,chs];
            break
        elseif strcmp(txt,'tooFastGR')
            break
        elseif strcmp(txt,'lost')%Case Lost
            if ~isempty(tmp)%Case Divide + Lost
                daughterOne=0;
                txt='divLost';%
                chs=tmp;%
                break;
            end
            count=count+1;
        elseif strcmp(txt,'elongated')
            break;
        end
    end
    %record data:
    tmplkg=[tmplkg;{strcat(num2str(j,'%0.2d'),txt)}];
    if strcmp(txt,'divLost')
        tmplng=[tmplng;{[daugID(chs),-1]}];%case divLost
    else
        tmplng=[tmplng;{daugID(chs)}];
    end
    if ~strcmp(txt,'tooFastGR')
        daugLen(chs)=[];
        daugID(chs)=[];
    end    
    
    
end
end
function [next,UIT,unexpected,chosen]=ZZMoMaComp(m,d)
%init
global daughterOne count
elo=[0.77, 1.32];
ediv=0.7;
unexpected=[];
chosen=[];

if m>elo(1)*d && m<elo(2)*d
    UIT='elongated';
elseif d<m*ediv
    UIT='divided';
elseif d>m*elo(2)
    UIT='tooFastGR';
else
    UIT='lost';    
end

%copy part 2
if daughterOne~=0 && strcmp(UIT, 'divided')
    if (daughterOne+d)<elo(2)*m && (daughterOne+d)>elo(1)*m
        %division is verified.
        UIT='DivVerified';            
    end
elseif daughterOne~=0 && ~strcmp(UIT, 'divided')
    %divided but the 2nd daughter lost
    %UIT='lost';
    UIT='lost';
end

%copy part 3
switch UIT
    case 'divided'
        next=1;
        if daughterOne==0
            daughterOne=d;
        end
        chosen=count;
    case 'lost'
        next=0;
    case 'divLost'
        next=0;
    case 'tooFastGR'
        next=0;
        chosen=count;
    case 'DivVerified'
        next=1;
        daughterOne=0;
        chosen=count;
    case 'elongated'
        chosen=count;
        next=1;
end
end
function [uit]=getPoleMarker(in,tnt,ML,sw)
%for the first frame:
if sw==1
    for i=1:length(ML)
        in(i,1)=rem(i,2);%facing up (old -> new) is 1.
    end
    uit=in;
    return;
end

newPS=nan(length(in),1);
for i=1:length(tnt)
    As=tnt{i,1};
    if length(As)==1 && As~=0
        newPS(As,1)=in(i,1);
    elseif length(As)==2 && As(2)~=-1
        newPS(As(1))=1;
        newPS(As(2))=0;
    elseif length(As)==2 && As(2)==-1
        newPS(As(1))=1;
    end
end
uit=[in,newPS];
end
end
