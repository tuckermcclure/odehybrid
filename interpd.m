function xi = interpd(t, x, ti, mode, varargin)

% interpd
%
% An interpolation function for time series which may have multiple values
% at a single timestep (e.g., discrete updates). It works like MATLAB's 
% built-in interp1.
% 
% xi = interpd(t, x, ti);
% xi = interpd(t, x, ti, mode);
% xi = interpd(t, x, ti, mode, ...);
%
% Inputs:
% 
% t     Times (n-by-1)
% x     States (n-by-m)
% ti    Output times (p-by-1)
% mode  '-' for left-most values, '+' for right-most values (default)
% ...   Any additional arguments to be based along to interp1 (e.g., 
%       'nearest').
%
% Outputs:
%
% xi  States corresponding to output times (p-by-m)
% 
% Example:
% 
% % Create some times and states. Note that there are two states at t=3.
% t = [2.76, 2.91, 3,   3,   3.12].';
% x = [0.2,  0.3,  0.4, 1.1, 1.2].';
% 
% % Create the desired output times.
% ti = (2.8:0.1:3.1).';
%
% % Interpolate (t, x) linearly at ti, keeping the right-most values.
% xi = interpd(t, x, ti).'
% 
% % Interpolate (t, x) linearly at ti, keeping the left-most values.
% xi = interpd(t, x, ti, '-').'
% 
% % Interpolate (t, x) using 'nearest' at ti, keeping the left-most values.
% xi = interpd(t, x, ti, '-', 'nearest').'
% 
% See also: examples_odehybrid.
%
% Online doc: http://www.anuncommonlab.com/doc/odehybrid/interpd.html
%
% Copyright 2014 An Uncommon Lab
    
    % Check that t and ti are columns.
    if ~any(size(t) == 1)
        error('interpd:InvalidInputs', 'Input t should be a column.');
    else
        t = t(:);
    end
    if ~any(size(ti) == 1)
        error('interpd:InvalidInputs', 'Input ti should be a column.');
    else
        ti = ti(:);
    end
    
    % Check that the rows of x match the rows of t.
    if size(x, 1) ~= length(t) && size(x, 2) == length(t)
        x = x.';
    end

    % Set the default mode.
    if nargin < 4 || isempty(mode)
        mode = '+';
    end

    % Create the output.
    xi = zeros(length(ti), size(x, 2));
    
    % Find where there are doubled steps (t(doubled(k)) and t(doubled(k)+1)
    % are the same.
    doubled = find(t(1:end-1) == t(2:end));
    
    % If nothing is doubled, just pass along to interp1.
    if isempty(doubled)
        xi = interp1(t, x, ti, varargin{:});
        return;
    end
    
    % We will break the interpolation into separate pieces -- those between
    % the doubled (or tripled or otherwise duplicated) time steps. Then
    % we'll interpolate that span. We'll keep the edges according to
    % whether we're using - or + mode.
    
    % Make sure we get the first and last spans.
    if doubled(end) ~= length(t)
        doubled(end+1, 1) = length(t);
    end
    
    % Start with the initial point.
    if strcmp(mode, '-') && doubled(1) ~= 1
        xi(1, :) = interp1(t(1:2), x(1:2, :), ti(1), varargin{:});
    else
        % We may overwrite this for +.
        xi(1, :) = x(1, :);
    end
    
    % For each doubled index...
    k_last = 1;
    for k = doubled.'
        
        % Get the span from the last index to this one.
        tk  = t(k_last:k);
        xk  = x(k_last:k, :);
        
        % Get the range for the requested outputs.
        out_indices = ti >= t(k_last) & ti <= t(k);
        
        % If we're using - mode, drop the first index; we don't want to
        % overwrite the value from the end of the last span with the value
        % from the beginning of this one.
        if strcmp(mode, '-')
            out_indices(find(out_indices, 1, 'first')) = false;
        end
        
        % If there are any outputs in this span...
        if any(out_indices)
            
            % Select the subset of output times.
            tik = ti(out_indices);

            % If there's only one point, don't interp. Otherwise, use
            % interp1.
            if k - k_last == 0
                xi(out_indices, :) = xk;
            else
                xi(out_indices, :) = interp1(tk, xk, tik, varargin{:});
            end
            
        end
        
        % Get ready to move to the next span.
        k_last = k + 1;
        
    end
    
end
