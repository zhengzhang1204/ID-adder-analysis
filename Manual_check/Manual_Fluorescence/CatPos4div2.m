function [uit,funcs,cate]=CatPos4div2(pos)
    % Suitable case: 0-, 1- or 2-spot cell, which is dividing.
    % For a given spot, this program provides the probability that this
    % spot enters each sister, and the corresponding prob-functions for
    % spots in each sister.
    % Output format:
    % 1. uit: 2x2 matrix. 1st column is probability that the spot enters
    % daughter, 2nd column is the daughter index. NOTE, UPPER DAUGHTER IS
    % 2, LOWER DAUGHTER IS 1.
    % 2. funcs: 2x1 cell for calculating probabilities for spots in each
    % cell.
    
    funcs=cell(2,1);
    funcs{2,1}=@(x,pos) cos(pi*(0.5*x+0.5-pos));
    funcs{1,1}=@(x,pos) cos(pi*(0.5*x-pos));

    if pos>0.75 || isnan(pos)%Cat3: 3/4
        cate=3;
        uit=[0;1];
    elseif pos<=0.75 && pos>0.25%Cat2: 1/2, Cat3: 3/4
        cate=2;
        uit=[(1.5-2*pos);(2*pos-0.5)];
    elseif pos>=0 && pos<=0.25%Cat1: 1/4
        cate=1;
        uit=[1;0];
    end
    uit=[uit,[1;2]];
end