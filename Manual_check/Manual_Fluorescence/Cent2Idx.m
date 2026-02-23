function [idx]=Cent2Idx(in,sz)
X=in(1)-3:1:in(1)+3;
Y=in(2)-3:1:in(2)+3;
X=repmat(X,7,1);
Y=repmat(Y',1,7);
Xi=X(:);
Yi=Y(:);
del=Xi<1|Yi<1|(sz(2)-Xi)<0|(sz(1)-Yi)<0;
Xi(del)=[];
Yi(del)=[];
idx=sub2ind(sz,Yi,Xi);
end