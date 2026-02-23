function [pred,funcs]=PredictCat4v2(curCat,pos)
    pred=zeros(3,2);%intended to be (3,2):[rng1,rng2];
    funcs=cell(3,1);
    for iCat=1:3%intended to be 3
        if curCat(iCat,1)==0; continue; end
        switch curCat(iCat,2)
            case 3
                pred(iCat,1:2)=[0.5, 1];
                funcs{iCat,1}=@(x,pos) cos(1.25*pi*(x-pos));
            case 2
                pred(iCat,1:2)=[0.25, 0.75];
                funcs{iCat,1}=@(x,pos) cos(1.25*pi*(x-pos));
            case 1
                pred(iCat,1:2)=[0, 0.5];
                funcs{iCat,1}=@(x,pos) cos(1.25*pi*(x-pos));
        end
    end
end