function [A]=prepLines(InfCell,CIDinParts,FLTS)
A=[];
for i=1:length(CIDinParts)%for each PH part:
    curCID=CIDinParts{i};
    curMsk=InfCell{curCID(end,1)}{1,1}==curCID(end,2);
    % get foci kymo R
    i90=any(curMsk,2);
    rng=[find(i90,1,"first")-3,min([find(i90,1,"last")+3,256])];
    H=diff(rng);
    cropLoc=zeros(length(curCID),2);
    cropImg=cell(1,length(curCID));
    for j=1:length(curCID)
        curMsk=InfCell{curCID(j,1)}{1,1}==curCID(j,2);
        i45=any(curMsk,2);
        ref=round((find(i45,1,"first")+find(i45,1,"last"))/2);
        cropLoc(j,:)=[ceil(ref-H/2),floor(ref+H/2)+1];
        if floor(ref+H/2)+1>256
            dix=floor(ref+H/2)-255;
            cropLoc(j,:)=[ceil(ref-H/2)-dix,256];
        end
        r=curCID(j,1);
        cropImg{j}=FLTS(cropLoc(j,1):cropLoc(j,2),((r-1)*32+1):r*32);
    end
    R=cell2mat(cellfun(@(x) sum(x,2),cropImg,'UniformOutput',false));
    R=mat2gray(R);

    % get initial guess skeleton
    tophatFiltered = imtophat(R, strel('disk', 3));
    background = imopen(tophatFiltered, strel('disk', 5));
    backgroundSubtracted = tophatFiltered - background;
    maxFiltered = imdilate(backgroundSubtracted, strel('disk', 5));
    gaussianFiltered = imgaussfilt(maxFiltered, 3);
    [~,ridges] = ridgefilt( gaussianFiltered,6,2.5,0.8 );
    B = bwskel(ridges);

    



    figure;imshow(mat2gray(B));
end




end