function [uit]=RemSinglets(Map)

%This program is to remove 'singlet arrows' from Map.
%There can be 'multi->1' arrows. If any original spots (i.e., spot where
%arrow is from) is orphan (i.e., it does not have original spot), remove
%this arrow.

A=unique(Map);
A(isnan(A))=[];
%1. look for 'multi -> 1' cases, record spotID, and rows that point to this ID
selSptID=[];
selLoc={};
for i = 1:length(A)
    if A(i)==-1;continue;end
    [r1,~]=find(Map(:,1)==A(i));
    [r2,~]=find(Map(:,2)==A(i));
    [r3,~]=find(Map(:,3)==A(i));

    if (length(r1)+length(r2)+length(r3))>1
        selSptID=cat(1,selSptID,A(i));
        selLoc=[selLoc;{[r1;r2;r3]}];
    end
end

%2. for all rows that points to One spotID, who is continuous and who is
%not?
for i=1:length(selLoc)%for each record in selLoc (spots that is 'multi->1')
    rows=selLoc{i,1};
    judge=false(length(rows),1);
    for iRow=1:length(rows)
        if ~isempty(find(Map==rows(iRow), 1))
            judge(iRow,1)=true;
        end
    end

%3. delete arrows.
    if sum(judge)==1%means there is one continuous arrow and one singlet.
        singlet=rows(~judge);
        singlet=unique(singlet);
        %remove this singlet from Map:
        for ix=1:length(singlet)
            Map(singlet(ix),Map(singlet(ix),:)==selSptID(i))=nan;
        end
%         idx=Map(singlet,:)==selSptID(i);
%         QQ=Map(singlet,idx);
%         QQ(idx)=nan;
%         Map(singlet,idx)=nan;
%         rep=Map(singlet,:);
%         rep(isnan(rep))=[];
%         Map(singlet(1),1:3)=[rep,nan(1,3-length(rep))];
    end
end

uit=Map;

end