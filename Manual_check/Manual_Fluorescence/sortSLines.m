function [uit1,uit2]=sortSLines(in1,in2)
out=cell2mat(cellfun(@(x) x(1,2:3), in1,'UniformOutput',false));
[~,idx]=sort(out(:,1),'ascend');
uit1=in1(idx,:);
uit2=in2(idx,:);
end