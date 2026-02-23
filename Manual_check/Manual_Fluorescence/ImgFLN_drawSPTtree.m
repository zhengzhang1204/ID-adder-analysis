function []=ImgFLN_drawSPTtree(SPTLine,SPTMap,sft,visStr,savePath,PathIn)
c1=[242, 141, 0]./255;
c2=[212, 212, 0]./255;
HLine_add=[];
VLine=[];
uit=SPTLine;
%1. collect data for drawing
for i=1:size(SPTMap,1) %for each spot line.
    [iRow,~]=find(SPTMap==i);%find the mother spot line of the current spot line
    if ~isempty(iRow)
        if SPTLine(i,1)<=uit(iRow,3)
            uit(iRow,3)=SPTLine(i,1);
        else
            if length(iRow)~=1
                error(['Spot line #',num2str(i),' is wrong!']);
            else
                HLine_add=cat(1,HLine_add,[uit(iRow,3:4),SPTLine(i,1),uit(iRow,4)]);%this may indicate: 1. wrong spot line mapping.
            end
        end
        VLine=cat(1,VLine,[SPTLine(i,1:2),SPTLine(i,1),uit(iRow,4)]);
    end
end

%2. draw spot line tree
%decompose input path
EXP=PathIn{2,1}(1:strfind(PathIn{2,1},'Y_res')-1);%EXP, with \
posStr=PathIn{2,1}(end-11:end-9);%e.g., p10
SchN=PathIn{3,1};%SchN
f=openfig([EXP,'Figs\',posStr,'S',num2str(SchN,'%0.2i'),'.fig'],visStr);%<changed on 23.02.20

% f=openfig('E:\DLstg32DivTree.fig',visStr);
ax = gca;
grid minor
for iA=1:size(uit,1)
    plot(ax, uit(iA,[1,3]),uit(iA,[2,4])+sft,'color',c1,'LineWidth',2);
end

for iB=1:size(HLine_add,1)
    plot(ax,HLine_add(iB,[1,3]),HLine_add(iB,[2,4])+sft,'color',c2,'LineWidth',2);
end

for iC=1:size(VLine,1)
    plot(ax,VLine(iC,[1,3]),VLine(iC,[2,4])+sft,'color',c2,'LineWidth',1);
end

%3. save figure
if ~isempty(savePath)
    savefig(f,savePath);
end

end