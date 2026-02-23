function [uit]=Neo_AdjustM(in,raw,tgt)
p=polyfit(raw,tgt,2);
F=@(x) p(1)*x.^2+p(2)*x+p(3);
uit=cellfun(@(x) uint16(floor(F(double(x)))), in,'UniformOutput',false);
% imwrite(uit{3},'E:\xxli.tif');
end