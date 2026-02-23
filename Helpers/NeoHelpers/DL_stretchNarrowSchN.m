function [uit]=DL_stretchNarrowSchN(in,p,psf_c)
%in: 256*32 uint16
%uit: 256*32 uint16
%Params: [a,b,width];
J=imresize(in,[256 p(3)]);
if p(1)<=0
    J=[repmat(J(:,1),1,-p(1)),J,repmat(J(:,end),1,-p(2))];
else
    J=J(:,p(1):p(1)+31);
end

J=deconvlucy(J,0.01*psf_c,4,uint16(100));
uit = adapthisteq(J,'NumTiles',[10 10]);
end