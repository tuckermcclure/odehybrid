function varargout = vectors_to_states(yv, varargin)
% vectors_to_states
%
% Convert a state vector history to a series of states given those states
% as inputs. If the total number of numeric states is ns and the number of
% time steps represented is nt, then yv is nt-by-ns. The remaining
% arguments are example states. The outputs will have the same form as the
% example states, but will be arrays of these states, each nt long, with
% appropriate values from yv.
% 
% When a state is a column or row vector of numeric types, the output will
% be an nt-by-x array where x is length(state)
%
% When a state is an n-dimensional matrix (neither column nor vector), the
% output will have the dimensions of the matrix concatenated on the first
% unused dimension. That is, if the state is a-by-b, the output will be
% a-by-b-by-nt.
% 
% When a state is a struct, the output will be a struct array of length nt.
% 
% Otherwise, the state will be a cell array of length nt, each containing
% the example type filled in with values from yv.
%
% [y1, y2, ... yn] = vectors_to_states(yv, y1, y2, ... yn);
% 
% Inputs:
%
% yv              State vector history (nt-by-ns)
% y1, y2, ... yn  Example states
%
% Outputs:
% 
% y1, y2, ... yn  State histories as described above
%
% Examples:
% 
% A history of 3 samples of a row vector:
% 
% yv = [state_to_vector([1 2]).'; ... % Sample 1
%       state_to_vector([7 8]).'; ... % Sample 2
%       state_to_vector([13 14]).'];  % Sample 3
% y  = zeros(1, 2);                   % Example state
% ys = vectors_to_states(yv, y)       % 3-by-2 output
%
% A history of 3 samples of 3 different states: a column vector, a cell
% array, and a struct:
% 
% % Make the example states.
% y1  = zeros(1, 3).';
% y2  = {1, 2; 3 4};
% y3  = struct('a', 0, 'b', 0);
% 
% % Make the history.
% y1v = [state_to_vector([1 2 3].').'; ...
%        state_to_vector([7 8 9].').'; ...
%        state_to_vector([13 14 15].').'];
% y2v = [state_to_vector({1 2; 7 8}).'; ...
%        state_to_vector({3 4; 9 0}).'; ...
%        state_to_vector({5 6; 1 2}).'];
% y3v = [1 2; 3 4; 5 6];
% yv  = [y1v, y2v, y3v];
% 
% % Get the state histories.
% [y1s, y2s, y3s] = vectors_to_states(yv, y1, y2, y3)
% 
% See also: vector_to_state, state_to_vector.
%
% Online doc: http://www.anuncommonlab.com/doc/odehybrid/vectors_to_states.html
%
% Copyright 2014 An Uncommon Lab
    
    % Number of samples
    ns = size(yv, 1);
    
    % Number of states
    nx = length(varargin);
    
    % For each type of state...
    count = 0;
    varargout = cell(1, nx);
    
    % If there are no states, return empties.
    if isempty(yv)
        return;
    end
    
    pre2008a = verLessThan('matlab', '8.1');
    
    for x = 1:length(varargin)
        
        % See if it's just a bunch of numbers (easy!).
        if isnumeric(varargin{x}) || ischar(varargin{x})
            
            nv  = numel(varargin{x});
            yvx = yv(:, count + (1:nv));
            szs = size(varargin{x});
            if length(szs) == 2 && any(szs(1:2) == 1)
                if pre2008a
                    varargout{x} = cast(yvx, class(varargin{x}));
                else
                    varargout{x} = cast(yvx, 'like', varargin{x});
                end
            else
                if pre2008a
                    varargout{x} = cast(reshape(yvx.', [szs ns]), ...
                                        class(varargin{x}));
                else
                    varargout{x} = cast(reshape(yvx.', [szs ns]), ...
                                        'like', varargin{x});
                end
            end
            count = count + nv;
            
        % Cell arrays are pretty easy too.
        elseif iscell(varargin{x})
            
            varargout{x} = cell(ns, 1);
            dummy = state_to_vector(varargin{x});
            yvx = yv(:, count + (1:length(dummy)));
            for k = 1:ns
                varargout{x}{k} = vector_to_state(yvx(k, :), varargin{x});
            end
            count = count + length(dummy);
            
        % Structures are wasteful but useful.
        elseif isstruct(varargin{x}) && length(varargin{x}) == 1
            
            % Make the output struct like the input struct.
            varargout{x} = varargin{x};
            dummy = state_to_vector(varargin{x});
            yvx = yv(:, count + (1:length(dummy)));
            for k = 1:ns
                varargout{x}(k,1) = vector_to_state(yvx(k,:), varargin{x});
            end
            count = count + length(dummy);
            
        else
            
            varargout{x} = cell(ns, 1);
            dummy = state_to_vector(varargin{x});
            yvx = yv(:, count + (1:length(dummy)));
            for k = 1:ns
                varargout{x}{k} = vector_to_state(yvx(k, :), varargin{x});
            end
            
        end
        
    end

end
