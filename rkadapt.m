function [t, x] = rkadapt(ode, ts, x0, options, a, b, c, p)

% rkadapt
% 
% Runge-Kutta adaptive-step integration using the specified weights, nodes,
% and Runge-Kutta matrix (or Dormand-Prince 5(4) by default).
% 
% Implements numerical propagation of an ordinary differential equation
% from some initial value over the desired range. This function is similar
% to MATLAB's variable-step ODE propagators (e.g., ode45), but is 
% generalized to any Butcher tableau. Unlike ode45, however, there's no
% inherent minimum number of steps that this function will take; if it can
% finish integration with a single step, it will, whereas ode45 will always
% take at least a certain number of steps, no matter how small the time
% interval and how smooth the dynamics. Therefore, this function is
% commonly more useful with odehybrid than ode45.
%
% This function is generic for all adaptive-step Runge-Kutta methods. That
% is, any adaptive-step Runge-Kutta propagator can be created by passing
% the weightes, nodes, and Runge-Kutta matrix (together, the Butcher 
% tableau) into this function. It will correctly take advantage of the 
% first-same-as-last property, eliminating one function evaluation per 
% step when the table has this property. See the examples below.
%
% [t, x] = rk4adapt(ode, ts, x0);
% [t, x] = rk4adapt(ode, ts, x0, options);
% [t, x] = rk4adapt(ode, ts, x0, options, a, b, c, p);
%
% Inputs:
% 
% ode     Ordinary differential equation function (s.t. x_dot = ode(t, x);)
% ts      Time span, [t_start, t_end]
% x0      Initial state (column vector)
% options Options structure from odeset. This function uses the
%         InitialStep, MaxStep, OutputFcn, RelTol, and AbsTol fields only.
% a       Runge-Kutta matrix (s-by-s)
% b       Weights (2-by-s), the top row containing the weights for the 
%         state update
% c       Nodes (s-by-1)
% p       Lowest order method represented by b (e.g., for
%         Runge-Kutta-Fehlberg, this is 4).
%
% Outputs:
% 
% t      Time history
% x      State history, with each row containing the state corresponding to
%        the time in the same row of t.
%
% Example:
% 
% % Simulate an undamped harmonic oscillator for 10s starting from an 
% % initial state of [1; 0] using the default a, b, c, and p
% % (Dormand-Prince 5(4)).
% [t, x] = rkadapt(@(t,x) [-x(2); x(1)], [0 10], [1; 0]);
% plot(t, x);
% 
% % Add options.
% options = odeset('MaxStep', 0.1, 'InitialStep', 0.05);
% [t, x] = rkadapt(@(t,x) [-x(2); x(1)], [0 10], [1; 0], options);
% plot(t, x);
% 
% % Now simulate with Bogacki-Shampine 3(2).
% a = [  0   0   0 0; ...
%      1/2   0   0 0; ...
%        0 3/4   0 0; ...
%      2/9 1/3 4/9 0];
% b = [ 2/9 1/3 4/9   0; ... % For the update
%      7/24 1/4 1/3 1/8];    % For the error prediction
% c = [0 1/2 3/4 1];
% p = 2;                     % The lower of the two orders
% [t, x] = rkadapt(@(t,x) [-x(2); x(1)], [0 10], [1; 0], [], a, b, c, p);
% plot(t, x);
%
% See "Runge-Kutta methods" on Wikipedia for discussion of the Butcher
% tableau (a, b, and c).
%
% Reference: 
%
% Dormand, J. R., and P. J. Prince. "A family of embedded Ruge-Kutta
% formulae." _Journal of Computational and Applied Mathematics_ 6.1 (1980):
% 19-26. Web. 22 Apr. 2014.
%
% See also: odehybrid, ode45, odeset, rkfixed.
%
% Online doc: http://www.anuncommonlab.com/doc/odehybrid/rkadapt.html
%
% Copyright 2014 An Uncommon Lab

    % Check the arguments.
    if nargin < 3
        error('rkadapt:IncorrectArguments', ...
              'Incorrect number of arguments.');
    elseif nargin >= 5 && nargin < 8
        error('rkadapt:IncorrectArguments', ...
              ['If any of a, b, c, and p are specified, all must ' ...
               'be specified.']);
    end
    
    % If there's no options structure, make one.
    if nargin == 3 || (nargin >= 4 && isempty(options))
        options = odeset();
    end
    
    % If there's no table, default to Dormand-Prince 5(4).
    if nargin < 8
        [a, b, c, p] = dptable();
    end

    % Unpack what we care about from the options.
    rel_tol = default(options.RelTol, 1e-3);
    abs_tol = default(options.AbsTol, 1e-6);
    dt_max  = default(options.MaxStep, inf);
    dt      = min(default(options.InitialStep, 1e-3), dt_max);
    
    % See if this table has the first-same-as-last property.
    fsal = all(a(end, :) == b(1, :));
    
    ex = 1/(p+1);                     % Time step correction exponent
    dt_min = 100 * eps(max(abs(ts))); % Minimum time step

    t   = ts(1);       % Time history
    x   = x0(:).';     % State history
    tk  = ts(1);       % Current time
    xk  = x0(:);       % Current state
    xkp = xk;          % Proposed state update
    
    nx = numel(x0);    % Number of states
    s  = size(b, 2);   % Length of weights
    d  = zeros(nx, s); % Matrix of derivatives
    db = b(2, :) - b(1, :);
    
    % If the user provided an OutputFcn, call it.
    if ~isempty(options.OutputFcn)
        options.OutputFcn(ts, x0, 'init');
    end
    
    % If using a table with FSAL, start the derivative here. 
    if fsal
        d(:, 1) = ode(tk, xk);
    end
    
    % Propagate, updating k to k+1 until we're done until t(k) == ts(2).
    stop = false;
    while tk < ts(2) - dt_min
        
        % Don't step past the end time.
        if tk + dt > ts(2)
            dt = ts(2) - tk;
        end
        
        % Evaluate the first bit of the derivatives.
        if ~fsal
            d(:, 1) = ode(tk, xk);
        end
        
        % Try to take a step.
        while true
            
            % Calculate derivatives.
            for z = 2:s
                dxk = sum(bsxfun(@times, dt*a(z, 1:z-1), d(:, 1:z-1)), 2);
                d(:, z) = ode(tk + c(z) * dt, xk + dxk);
            end

            % Update the state and error estimate.
            tkp = tk + dt;
            xkp = xk;
            er  = zeros(nx, 1);
            for z = 1:s
                xkp = xkp + (dt * b(1, z)) * d(:, z);
                er  = er  + (dt * db(z))   * d(:, z);
            end

            % Create the scaled error. This single number can now be
            % compared to 1 to see if the worst error satisfies both 
            % absolute and relative tolerances. The error is acceptable if
            % it is either less than the relative error or less than the
            % aboslute error.
            scaled_error = max(abs(er) ./ (abs_tol + rel_tol * abs(xkp)));

            % If the error is small, update the time step and break to the
            % next step.
            if scaled_error <= 1
                
                dt = min(0.9 * dt * scaled_error^-ex, min(3*dt, dt_max));
                break;
                
            % If the error is too large, update the time step and continue
            % the search for an appropriate time step.
            else

                % If that was already the smallest step we can take, give
                % up.
                if dt == dt_min
                    warning('Could not meet integration tolerance.');
                    break;
                end
                
                % Cut the time step down faster than we build it up.
                dt = max(0.7 * dt * scaled_error^-ex, dt_min);
                
            end
            
        end

        % Accept the update.
        tk = tkp;
        xk = xkp;
        
        % For first-same-as-last, we can copy over the derivative at k.
        if fsal
            d(:, 1) = d(:, end);
        end
        
        % Store.
        t = [t; tk];      %#ok<AGROW>
        x = [x; xk(:).']; %#ok<AGROW>
        
        % If the user provided an OutputFcn, call it.
        if ~isempty(options.OutputFcn)
            stop = options.OutputFcn(tk, xk, '');
            if stop
                break;
            end
        end
        
    end

    % Rectify the ending.
    if ~stop
        t(end) = ts(2);
    end
    
    % If the user provided an OutputFcn, call it.
    if ~isempty(options.OutputFcn)
        options.OutputFcn([], [], 'done');
    end

end

% Use the value in f or use v as the default.
function v = default(f, v)
    if ~isempty(f)
        v = f;
    end
end

% Dormand-Prince Butcher tableau
function [a, b, c, p] = dptable()
    a = [    0           0          0         0         0        0    0;...
             1/5         0          0         0         0        0    0;...
             3/40        9/40       0         0         0        0    0;...
            44/45      -56/15      32/9       0         0        0    0;...
         19372/6561 -25360/2187 64448/6561 -212/729     0        0    0;...
          9017/3168 -355/33     46732/5247   49/176 -5103/18656  0    0;...
            35/384     0          500/1113  125/192 -2187/6784  11/84 0];
    b = [35/384     0  500/1113  125/192  -2187/6784    11/84   0; ...
         5179/57600 0 7571/16695 393/640 -92097/339200 187/2100 1/40];
    c = [0 1/5 3/10 4/5 8/9 1 1];
    p = 4;
end
