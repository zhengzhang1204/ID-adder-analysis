function [uit]=CatPos4(pos)
    % as for an 0-, 1- or 2-spot cell, give the possibility of categories that the spot belongs to by its relative
    % location in the cell. Cat1: 1/4, Cat2: 1/2, Cat3: 3/4. 'uit(:,1)' is the
    % possibility of these categories (uit(:,2)) in descend order.
    if pos>0.75 || isnan(pos)%Cat3: 3/4
        uit=[0;0;1];
    elseif pos<=0.75 && pos>0.5%Cat2: 1/2, Cat3: 3/4
        uit=[0;(3-4*pos);(4*pos-2)];
    elseif pos<=0.5 && pos>0.25%Cat1: 1/4, Cat2: 1/2
        uit=[(2-4*pos);(4*pos-1);0];
    elseif pos>=0 && pos<=0.25%Cat1: 1/4
        uit=[1;0;0];
    end
    uit=[uit,[1;2;3]];
    [~,idxE]=sort(uit(:,1),'descend');
    uit=uit(idxE,:);
end