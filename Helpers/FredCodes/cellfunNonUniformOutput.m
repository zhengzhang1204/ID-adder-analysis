function output = cellfunNonUniformOutput(varargin)
%CELLFUNNONUNIFORMOUTPUT Summary of this function goes here
%   Detailed explanation goes here

output = cellfun(varargin{:},'UniformOutput',false);
end

