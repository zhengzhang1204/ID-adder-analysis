function sepOnes = genSeparableOnes(sizeOfOnes)
%GENSEPARABLEONES will generate a separable version of the ones
dims = numel(sizeOfOnes);
sepOnes = cell(1,dims);
for i = 1:dims
    reshapeVec = ones(dims,1);
    reshapeVec(i) = sizeOfOnes(i);
    currGauss = reshape(ones(1,sizeOfOnes(i)),reshapeVec');
    sepOnes{i} = currGauss;
end


end

