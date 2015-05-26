function [t, x] = rkfixed(ode, ts, x0, dt, a, b, c)

% rkfixed
% 
% Runge-Kutta fixed-step integration using the specified weights, nodes,
% and Runge-Kutta matrix (or the Runge-Kutta 4th order "3/8" method by
% default).
% 
% Implements numerical propagation of an ordinary differential equation
% from some initial value over the desired range. This function is similar
% to MATLAB's variable-step ODE propagators (e.g., ode45), but uses a
% fixed step method. This is useful either when one knows an appropriate 
% step size or when a process is interrupted frequently (ode45 and the
% similar functions in MATLAB will always make at least a certain number of
% steps between ts(1) and ts(2), which may be very many more than are
% necessary).
%
% This function is generic for all fixed-step Runge-Kutta methods. That is,
% any fixed-step Runge-Kutta propagator can be created by passing the
% weightes, nodes, and Runge-Kutta matrix (together, the Butcher tableau) 
% into this function. See the example below.
%
% [t, x] = rk4(ode, ts, x0, dt);
% [t, x] = rk4(ode, ts, x0, dt, a, b, c);
% [t, x] = rk4(ode, ts, x0, options);
% [t, x] = rk4(ode, ts, x0, options, a, b, c);
%
% Inputs:
%
% ode     Ordinary differential equation function
% ts      Time span, [t_start, t_end]
% x0      Initial state (column vector)
% dt      Time step
% options Alternately, one can specify an options structure instead of dt
%         so that this function is compatible with ode45 and its ilk. The
%         only valid fields are MaxStep (the time step) and OutputFcn
% a       Runge-Kutta matrix
% b       Weights
% c       Nodes
%
% Outputs:
%
% t      Time history
% x      State history, with each row containing the state corresponding to
%        the time in the same row of t.
%
% Example:
% 
% % Simulate an undamped harmonic oscillator for 10s with a 0.1s time
% % step, starting from an initial state of [1; 0] using RK 4th order
% % integration (via the Butcher tableau specified by a, b, and c). This
% % is exactly the same as the rk4 function.
% a = [0   0   0 0; ...
%      0.5 0   0 0; ...
%      0   0.5 0 0; ...
%      0   0   1 0];
% b = [1 2 2 1]/6;
% c = [0 0.5 0.5 1];
% [t, x] = rkfixed(@(t,x) [-x(2); x(1)], [0 10], [1; 0], 0.1, a, b, c);
% plot(t, x);
%
% See "Runge-Kutta methods" on Wikipedia for discussion of the Butcher
% tableau (a, b, and c).
%
% See also: odehybrid, ode45, odeset, rk4.
%
% Online doc: http://www.anuncommonlab.com/doc/odehybrid/rkfixed.html
%
% Copyright 2014 An Uncommon Lab

    % Default to Runge and Kutta's 3/8 formulation (4th order).
    if nargin == 4
        a = [ 0    0 0 0; ...
              1/3  0 0 0; ...
             -1/3  1 0 0; ...
              1   -1 1 0];
        b = [1 3 3 1]/8;
        c = [0 1/3 2/3 1];
    end

    % Allow an alternate input syntax to be similar with ode45.
    if isstruct(dt)
        options = dt;
        dt = options.MaxStep;
        if isempty(dt)
            error('Specify the time step with the MaxStep option.');
        end
    else
        options.OutputFcn = [];
    end

    % Time history
    t  = (ts(1):dt:ts(end)).';
    if t(end) ~= ts(end)
        t(end+1, 1) = ts(end);
    end
    
    ns = length(t);                   % Number of samples
    nx = numel(x0);                   % Number of states
    x  = [x0(:).'; zeros(ns-1, nx)];  % State history
    xk = x0(:);                       % Current state
    
    s = length(b);                    % Length of weights
    d = zeros(nx, s);                 % Matrix of derivatives
    
    % If the user provided an OutputFcn, use it.
    if ~isempty(options.OutputFcn)
        options.OutputFcn(ts, x0, 'init');
    end
    
    % Propagate.
    for k = 1:ns-1
        
        % The last sample may be cut short.
        if k == ns-1
            dt = t(k+1) - t(k);
        end
        
        % Current time
        tk = t(k);
        
        % Calculate derivatives.
        d(:, 1) = dt * ode(tk, xk);
        for z = 2:s
            dxk = sum(bsxfun(@times, a(z, 1:z-1), d(:, 1:z-1)), 2);
            d(:, z) = dt * ode(tk + c(z) * dt, xk + dxk);
        end
        
        % Update the state.
        for z = 1:s
            xk = xk + b(z) * d(:, z);
        end
        
        % Store.
        x(k+1, :) = xk(:).';
        
        % If the user provided an OutputFcn, use it.
        if ~isempty(options.OutputFcn)
            if options.OutputFcn(t(k+1), xk, '');
                t = t(1:k+1);
                x = x(1:k+1, :);
                break;
            end
        end
        
    end

    % If the user provided an OutputFcn, use it.
    if ~isempty(options.OutputFcn)
        options.OutputFcn([], [], 'done');
    end

end
