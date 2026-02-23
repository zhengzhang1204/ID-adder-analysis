function [J]=DL_ModMasksMet2(J,minA,minL,isNarrow,p,Offset)
%0. if input is for narrow side channels, shrink it.
%1. Course find bottom and erase pixels under it.
%2. Smooth
%3. Erase small things

%narrowParams=[a,b,width];can be squeeze case: [-3,-3,26]
if isNarrow
    if p(1)>0%expanded for segmentation, now should be shrinked
        J=[single(zeros(256,p(1)-Offset)),J,single(zeros(256,p(2)-Offset))];
        J=imresize(J,[256,32]);
    elseif p(1)<0
        J=J(:,1-p(1)+Offset:32+p(2)-Offset);
        J=imresize(J,[256,32]);
    else
        if Offset>0
            J=J(Offset:32-Offset);
        elseif Offset<0
            J=[zeros(256,-Offset),J,zeros(256,-Offset)];
        end
        J=imresize(J,[256,32]);
    end
end
J=(J>0.1);%converted to logical.
% J(end-8:end,:)=zeros(9,32);
J=bwmorph(J,'hbreak');
J=bwmorph(J,'spur');
J=imerode(J,strel('disk',1));
J=imopen(J,strel('disk',3));
J=imopen(J,strel('diamond',2));
L=regionprops(J,'Area','MajorAxisLength','PixelIdxList');
for j=1:length(L)
    if L(j).Area<minA || L(j).MajorAxisLength<minL%test numbers
        J(L(j).PixelIdxList)=false;
    end
end
end