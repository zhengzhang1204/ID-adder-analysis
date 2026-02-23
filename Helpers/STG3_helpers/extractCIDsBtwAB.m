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