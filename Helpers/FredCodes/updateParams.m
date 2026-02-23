function params = updateParams(params,userParams)
%UPDATEPARAMS takes userParams, which can be a structure or a named,value
%cell list, and updates the params structure fields.
% userParams can be 
% 1) {'name',value,...}
% 2) struct
% 3) {struct}
% 4) {'name', value, struct}

% if ~isempty(userParams)
%    if iscell(userParams)
%        % are there structs 
%    elseif isstruct(userParams)
%         params = setstructfields(params,userParams);
%    else
%       error('userParams is neither a cell or a struct');
%    end
% end

if isempty(userParams)
   return; 
end

if iscell(userParams)
    % check if there is a param struct 
    findTheStructs = find(cellfun(@isstruct,userParams)==1);
    % make sure that it is in an odd position so that it is not a
    % named-value pair
    findTheParamStructs = find(isodd(findTheStructs));
    findTheParamStructs = findTheStructs(findTheParamStructs);
    % update params
    for i = 1:numel(findTheParamStructs)
       params = setstructfields(params,userParams{findTheParamStructs(i)}); 
    end
    % delete the param structs
    userParams(findTheParamStructs) = [];
    % update the named value pairs
    if isodd(numel(userParams))
        % if it is odd
        if numel(userParams) == 1
            flattendParams = flattenCellArray(userParams);
            params = parsepropvalFC(params,flattendParams{:});
        else
           error('named value pairs is not even'); 
        end
        
    else
       params = parsepropvalFC(params,userParams{:}); 
    end
    
    
elseif isstruct(userParams)
    params = setstructfields(params,userParams);
else
   error('userParam argument neither a cell of named-value pairs or a struct'); 
end

% if ~isempty(userParams)
%    if iscell(userParams)
%        if isstruct(userParams{1})
%            params = setstructfields(params,userParams{1});
%        else
%            params = parsepropvalFC(params,userParams{:});
%        end   
%    else
%        if isstruct(userParams)
%            params = setstructfields(params,userParams);
%        else
%            error('userParams is not the type expected');
%        end
%    end
% end

end

