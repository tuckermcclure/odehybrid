function [state, count] = vector_to_state(vector, state, count)

% vector_to_state
%
% Build a state from a state vector by recursively moving through the 
% elements of a more complex example state (matrix, cell array, or struct),
% placing the relevant numerical values from the vector into the 
% appropriate places. This is used in odehybrid to convert a state vector
% into the original complex setup of states. The reverse function is
% state_to_vector.
% 
% state = vector_to_state(vector, state);
% 
% Inputs:
% 
% vector  A vector containing the numerical values from the object
% state   A numeric, cell array, or struct array type representing the
%         sizes and types for the values in vector
% count   (Internal use only)
% 
% Outputs:
%
% state   A numeric, cell array, or struct array type
% count   (Internal use only)
%
% Example:
% 
% x = {[1 3; 4 2], struct('a', {5; 6}, 'bcd', {7:9; 10}), pi}
% v = state_to_vector(x)
% x2 = vector_to_state(v, x)
% isequal(x, x2)
%
% See also: state_to_vector.
%
% Online doc: http://www.anuncommonlab.com/doc/odehybrid/vector_to_state.html
%
% Copyright 2014 An Uncommon Lab

    % If the user didn't input count, start from the beginnning.
    if nargin < 3, count = 0; end;
    
    % If it's a cell, convert each element.
    if iscell(state)
        
        for k = 1:numel(state)
            [state{k}, count] = vector_to_state(vector, state{k}, count);
        end
        
	% If it's a struct, convert and store in each field.
    elseif isstruct(state)
        
        fields = sort(fieldnames(state));
        for n = 1:length(state)
            for k = 1:length(fields)
                [state(n).(fields{k}), count] = vector_to_state( ...
                                                  vector, ...
                                                  state(n).(fields{k}), ...
                                                  count);
            end
        end
        
	% If it's a bunch of numbers or chars, convert them from doubles to the
	% appropriate type and store.
    elseif isnumeric(state) || ischar(state)

        % Since we overwrite state, it implicitly casts as necessary so
        % that the output state will be the right type, whereas vector will
        % always be a double.
        state(:) = vector(count+(1:numel(state)));
        count = count + numel(state);
        
	% Otherwise, we don't know what to do with it.
    else
        error('I don''t know to use the %s type.', class(state));
    end

end
