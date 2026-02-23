function [uit,cnt]=getGRthrSpotLine(CIDin,InfCell,CIDinParts,int)

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