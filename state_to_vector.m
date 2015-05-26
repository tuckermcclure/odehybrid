function vector = state_to_vector(state)

% state_to_vector
%
% Build a state vector by recursively moving through the elements of a more
% complex state (matrix, cell array, or struct), saving the numerical 
% values along the way in a vector. This is used in odehybrid to convert a
% complex series of states into a vector to be used for continuous state
% updates. The reverse function is vector_to_state.
%
% vector = state_to_vector(state);
% 
% Inputs:
% 
% state   A numeric, cell array, or struct array type, the contents of
%         which consist of numeric, cell array, or struct array types, etc.
% 
% Outputs:
%
% vector  A vector containing the numerical values from the state (doubles)
%
% Example:
% 
% x = {[1 3; 4 2], struct('a', {5; 6}, 'bcd', {7:9; 10}), pi}
% v = state_to_vector(x)
% x2 = vector_to_state(v, x)
% isequal(x, x2)
%
% See also: vector_to_state.
%
% Online doc: http://www.anuncommonlab.com/doc/odehybrid/state_to_vector.html
%
% Copyright 2014 An Uncommon Lab

    % Ignore growth in loops. We're just going to have this problem.
    %#ok<*AGROW>
    
	% If it's numbers, convert to doubles and store.
    if isnumeric(state) || ischar(state)
        
        vector = double(state(:));
        
    % If it's a cell, convert each element.
    elseif iscell(state)
    
        vector = [];
        for k = 1:numel(state)
            vector = [vector; state_to_vector(state{k})];
        end
        
    % If it's a struct, convert each field.
    elseif isstruct(state)
        
        vector = [];
        fields = sort(fieldnames(state));
        for n = 1:length(state)
            for k = 1:length(fields)
                vector = [vector; state_to_vector(state(n).(fields{k}))];
            end
        end

    % Otherwise, we don't know what to do.
    else
        error('I don''t know to use the %s type.', class(state));
    end

end
