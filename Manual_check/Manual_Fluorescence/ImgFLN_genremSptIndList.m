function [uit]=ImgFLN_genremSptIndList(curM,map)
uit=[];
curCand=map(curM,:);
curCand(isnan(curCand))=[];

Reg=[];
oldReg=[];
usedReg=[];
while ~isempty(curCand) || ~isempty(Reg)
    if length(curCand)==1
        temp=[curM,curCand];
        curM=curCand;
        curCand=map(curM,:);
        curCand(isnan(curCand))=[];
    elseif length(curCand)>=2
        temp=[curM*ones(length(curCand),1),curCand'];
        idx=find(ismember(curCand,usedReg));
        if ~isempty(idx)
            temp(idx,:)=[];curCand(idx)=[];
        end

        curM=curCand(1);
        usedReg=curCand(1);
        if length(curCand)>1
            curCand=curCand(2:end);
        else
            curCand=[];
        end
        Reg=[Reg;curCand'];Reg=unique(Reg);
        oldReg=[oldReg;curCand'];oldReg=unique(oldReg);
        curCand=map(curM,:);
        curCand(isnan(curCand))=[];
    elseif ~isempty(Reg)
        curM=Reg(1);
        usedReg=[usedReg;Reg(1)];
        Reg=Reg(2:end);
        curCand=map(curM,:);
        curCand(isnan(curCand))=[];
        continue;
    end
    uit=[uit;temp];
end
uit = unique(uit,'rows');
end