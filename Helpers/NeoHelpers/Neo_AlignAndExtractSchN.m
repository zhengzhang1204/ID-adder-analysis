function [uitPH4seg,uitPH,uitCIx,narrowParams,nFrm,nSchN]=Neo_AlignAndExtractSchN(cropParam,storePH,narrow,psf_c)
Cent=cellfun(@(x) (mean(x,2))', cropParam,'UniformOutput',false);
nFrm=length(cropParam);
nSchN=length(Cent{1,1});
Map=nan(length(Cent),nSchN);
uitPH4seg=uint16(zeros(256,32,1,nFrm*nSchN));%format of uitPH4seg: [4-D uint16] 256 x 32 x 1 x [nFrm x nSchN]
uitPH=uint16(zeros(256,32*nSchN,nFrm));%format of uitPH [uint16] 256 x [nSchN x 32] x nFrm
uitCIx=single(nan(2, nSchN, nFrm));%format of uitCI: [3-D uint8] 2 x nSchN x nFrm

% Core. Fill up the align map.
for iSch=1:nSchN%for each side channel
    ref=Cent{1,1}(iSch);
    Map(1,iSch)=iSch;
    for iFrm=2:length(Cent)%for each frame
        [value,idx]=min(abs(Cent{iFrm}-ref));
        if value<9
            Map(iFrm,iSch)=idx;
            ref=Cent{iFrm}(idx);
        else
            break;
        end
    end
end

% Initialization 2:
if narrow.isNarrowChn
    a=floor((narrow.width-32)/2);
    b=narrow.width-a-32;
    narrowParams=[a,b,narrow.width];
else
    narrowParams=[];
end

% Core. Fill up outputs. uitPH (tilized PH side channel images, sorted by channels), uitCI (cropping parameters for side channels)
for iSchN=1:nSchN%for each side channel
    if narrow.isNarrowChn
        for iFrm=1:nFrm%for each frame
            if ~isnan(Map(iFrm,iSchN))
                idx=Map(iFrm,iSchN);
                uitPH4seg(:,:,1,(iSchN-1)*nFrm+iFrm)=DL_stretchNarrowSchNV2(storePH{iFrm}(:,(idx-1)*32+1:idx*32),narrowParams,psf_c);
                uitPH(:,(iSchN-1)*32+1:iSchN*32,iFrm)=storePH{iFrm}(:,(idx-1)*32+1:idx*32);
                uitCIx(:,iSchN,iFrm)=(cropParam{iFrm}(idx,:))';
            end
        end
    else
        for iFrm=1:nFrm%for each frame
            if ~isnan(Map(iFrm,iSchN))
                idx=Map(iFrm,iSchN);
                uitPH4seg(:,:,1,(iSchN-1)*nFrm+iFrm)=storePH{iFrm}(:,(idx-1)*32+1:idx*32);
                uitPH(:,(iSchN-1)*32+1:iSchN*32,iFrm)=storePH{iFrm}(:,(idx-1)*32+1:idx*32);
                uitCIx(:,iSchN,iFrm)=(cropParam{iFrm}(idx,:))';
            end
        end
    end
end
end