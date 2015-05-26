function varargout = odehybrid(solver, ode, de, dt, ts, yc0, yd0, varargin)

% odehybrid
% 
% Hybrid continuous and discrete propagation.
%
% This function propagates an ordinary differential equation along with a
% discrete update equation in a manner similar to MATLAB's ode45 function
% (which only propagates an ordinary differential equation). This is useful
% for implementing discrete-time controllers or simulating processes that
% are updated at discrete intervals.
% 
% A large number of examples can be found in examples_odehybrid or by
% entering the following at the command line:
%
%   home = fileparts(which('examples_odehybrid'));
%   web(fullfile(home, 'html', 'examples_odehybrid.html'));
%
% Interfaces:
% 
% [t, yc, td, yd] = odehybrid(solver, ode, ...
%                             de, dt, ...
%                             ts, ...
%                             yc0, yd0);
% [t, yc1..m, td, yd1..n] = odehybrid(solver, ode, ...
%                                     de, dt, ...
%                                     ts, ...
%                                     {yc1, yc2..m}, {yd1, yd2..n});
% [t, ..., td, ..., te, yc1..m, yd1..n, ie] = odehybrid(solver, ode, ...
%                                                       de, dt, ...
%                                                       ts, ...
%                                                       yc0, yd0);
% [...] = odehybrid(solver, ode, ...
%                   {de1, de2, ...}, [dt1, dt2, ...], ...
%                   ts, ...
%                   yc0, yd0);
% [...] = odehybrid(..., [options], [log]);
% sol = odehybrid(...);
%
% Inputs:
%
% solver   Continuous-time propagator to use, e.g. @ode45
% ode      Ordinary differential equation to use with solver. The interface
%          should be fun(t, xc1, xc2, ..., xcm, xd1, xd2, ..., xdn) where
%          xc1, xc2, ..., xcm are the continuous states and xd1, xd2, ..., 
%          xdn are the discrete states. It should return the derivatives of
%          the continuous states (m outputs).
% de       Discrete update equation(s) (either a function handle or cell
%          array of function handles) with the same inputs as ode but
%          outputing the updated continuous and discrete states (n+m
%          outputs).
% dt       Time step(s) of discrete update equation(s). If de is a cell
%          array of multiple functions, this should be an array of the same
%          size.
% ts       Time span, [t_start, t_end]
% yc0      Cell array of initial continuous states
% yd0      Cell array of initial discrete states
% options  (Optional) options structure from odeset
% log      (Optional) TimeSeriesLogger for logging in the ode and de. If a
%          log is passed in, both ode and de *must* be able to accomodate
%          having a log input or not as the final argument. E.g., |ode|
%          will be called as: ode(t, xc1, ..., xd1, ..., xdn) and
%          ode(t, xc1, ..., xd1, ..., xdn, log).
%
% Outputs:
% 
% t        Times corresponding to continuous states (nc-by-1)
% yc1..m   Continuous state outputs (m outputs, each nc-by-1)
% td       Times corresponding to discrete state updates (nd-by-1)
% yd1..n   Discrete state outputs (n outputs, each nd-by-1)
% te       Times corresponding to events (ne-by-1)
% yce1..m  Continuous states at events (m outputs, each ne-by-1)
% yde1..n  Discrete states at events (n outputs, each ne-by-1)
% ie       Index of triggering event. See documentation in odeset for more 
%          on events.
% sol      If a single output is requested, it will be a structure with
%          fields for each of the individual outputs. The various
%          continuous states will be grouped in 'yc', the discrete into
%          'yd', the continuous states at events into 'yce', and the
%          discrete states at events into 'yde'.
%
% Example:
% 
% This is a quick example of simulating an unstable continuous system with
% a stabilizing discrete-time controller.
% 
% ode = @(t, x, u) [0 1; 2 0] * x + [0; 1] * u;  % Differential equation
% de  = @(t, x, u) deal(x, -[8 4] * x);          % Discrete update equation
% dt  = 0.1;                                     % Discrete eq. time step
% ts  = [0 5];                                   % From 0 to 5s
% x0  = [1; 0];                                  % Initial continuous state
% u0  = 0;                                       % Initial discrete state
% [t, x, tu, u] = odehybrid(@rkadapt, ode, de, dt, ts, x0, u0); % Simulate!
% plot(t, x, tu, u, '.'); xlabel('Time');                       % Plot 'em.
% legend('x_1', 'x_2', 'u');                                    % Label 'em.
%
% See also: examples_odehybrid, ode45, rkadapt, rkfixed.
%
% Online doc: http://www.anuncommonlab.com/doc/odehybrid/odehybrid.html
%
% Copyright 2014 An Uncommon Lab

    % Check inputs
    if nargin < 7
        error('odehybrid:TooFewArgs', 'Too few input arguments.');
    elseif length(ts) ~= 2
        error('odehybrid:InvalidArguments', ...
              ['Propagation window should specify both start and stop '...
               'times and only start and stop times.']);
    elseif nargin >= 8 && ~isempty(varargin{1}) && ~isstruct(varargin{1})
        error('odehybrid:InvalidArguments', ...
              'Expected odeset for argument 8.');
    elseif nargin >= 9 && ~isa(varargin{2}, 'TimeSeriesLogger')
        error('odehybrid:InvalidArguments', ...
              'Expected TimeSeriesLogger for argument 9.');
    end

    % We have to intercept the output_fcn or it will get called with 'init'
    % and 'done' for every short step between the discrete steps. That's
    % probably undesirable.
    if nargin >= 8 && ~isempty(varargin{1})
        if ~isempty(varargin{1}.OutputFcn)
            f = varargin{1}.OutputFcn;
            varargin{1}.OutputFcn = @(varargin) ...
                       output_fcn(f, varargin{:}, yc0);
        end
        if ~isempty(varargin{1}.Events)
            f = varargin{1}.Events;
            varargin{1}.Events = @(varargin) ...
                       event_fcn(f, varargin{:}, yc0);
        end
    end
    
	% See if we're passing in separated states (the "full" version).
    if ~isnumeric(yc0) || ~isnumeric(yd0)

        % Get the full inputs.
        [varargout{1:nargout}] = odehybridfull(solver, ode, de, dt, ...
                                               ts, yc0, yd0, varargin{:});

	% Otherwise, just pass everything on directly to the odehybridcore.
    else
        [varargout{1:nargout}] = odehybridcore(solver, ode, de, dt, ...
                                               ts, yc0, yd0, ...
                                               varargin{:});
    end

end

% odehybridcore(solver, ode, de, dt, ts, yc0, yd0, [options], [log])
function [t, yc, td, yd, te, yce, yde, ie] = odehybridcore(solver, ...
                                                           ode, ...
                                                           de, dt, ...
                                                           ts, yc0, yd0,...
                                                           varargin)
    
    % Ignore growth in loops. We're just going to have this problem.
    %#ok<*AGROW>
    
    % For multiple sample rates, we'll expect the discrete updates to be
    % stored in a cell array. If the user is just passing in a single
    % discrete update, it might not be in a cell array. Put it into one.
    if ~iscell(de)
        de = {de};
    end
    
    % Make sure dt is a row vector.
    if size(dt, 1) > 1
        dt = dt.';
    end
    
    % Set the solver's options if the user provided them.
    if nargin >= 8 && ~isempty(varargin{1})
        options = varargin{1};
    else
        options = odeset();
    end

    % Check for the things we don't support.
    if ~isempty(options.Vectorized)
        error('odehybrid:options', ...
              'odehybrid doesn''t support vectorized functions.');
    elseif ~isempty(options.OutputSel)
        error('odehybrid:options', ...
              ['odehybrid doesn''t support the OutputSel option ' ...
               'because states need not be vectors.']);
    end
    
    % Figure out the time steps. First, calculate the resolution with which 
    % we will be able to differentiate steps. Then, make arrays for each 
    % separate time step.
    epsilon = 2 * max((diff(ts) ./ dt) .* eps(dt) + eps(ts(2)));
    tds = cell(1, length(dt));
    for k = 1:length(dt)
        tds{k} = ts(1):dt(k):ts(end);
    end
    
    % Sort the necessary discrete times into a single list and remove
    % doubled steps. We now have a list of all times at which to break for
    % a discrete step of one or more discrete update functions.
    td = sort([tds{:}]);
    td([false diff(td) < epsilon]) = [];
    td = td(:);
    
    % We're going to overwrite the output_fcn, so store the original.
    orig_outputfcn = options.OutputFcn;
    
    % We're going to overwrite the output_fcn, so store the original.
    orig_eventfcn = options.Events;
    
    % If the user passed in a log, use it.
    add_logging = nargin >= 9 && ~isempty(varargin{2});
    if add_logging
        
        % Add logging to the output (preserving the user's outputfcn if
        % it's provided).
        log = varargin{2};
        for k = 1:length(dt)
            de{k} = @(t, yc, yd) de{k}(t, yc, yd, log);
        end
        
    else
        log = [];
    end
    
    % Set the maximum time step to be the smaller of what's already
    % specified or dt.
    if isempty(options.MaxStep)
        options.MaxStep = min(dt);
    else
        options.MaxStep = min(options.MaxStep, min(dt));
    end
    
    % Set the initial step.
    if isempty(options.InitialStep)
        options.InitialStep = min(0.1 * min(dt), 0.01 * diff(ts));
    end
    
    % Initialize outputs.
    t  = ts(1); yc = yc0.';
    te = []; yce = []; ie = []; % Events
    if iscell(yd0)
        yd = cell(length(td), numel(yd0));
        [yd{1, :}] = yd0{:};
        yde = {};
    else
        yd = [yd0.'; zeros(length(td)-1, length(yd0))];
        yde = [];
    end
    
    % Call the output function if there is one.
    if ~isempty(orig_outputfcn)
        orig_outputfcn(ts, yc0, yd0, 'init');
    end
    
    % Loop through all k, updating from k to k+1.
    counts = zeros(size(dt));
    for k = 1:length(td)
        
        %%%%%%%%%%%%%%%%%%%%
        % Discrete Updates %
        %%%%%%%%%%%%%%%%%%%%
        
        if k <= length(td)

            % Find out what should discrete functions should be called at 
            % this time.
            to_tick = find(abs(counts .* dt + ts(1) - td(k)) <= epsilon);

            % Call the relevant discrete updates.
            for z = to_tick
                [yc0, yd0] = de{z}(td(k), yc0(:), yd0);
            end

            % We can now store the updated state for k. 
            if iscell(yd0)
                [yd{k, :}] = yd0{:};
            else
                yd(k, :) = yd0.';
            end

            % Tick them.
            counts(to_tick) = counts(to_tick) + 1;
    
            % Call the output function with the udpated discrete states.
            if    ~isempty(orig_outputfcn) ...
               && orig_outputfcn(td(k), yc0, yd0, '');
                td = td(1:k);
                if iscell(yd)
                    yd = yd(1:k);
                else
                    yd = yd(1:k, :);
                end
                break;
            end
            
        end
        
        %%%%%%%%%%%%%%%%%%%%%%
        % Continuous Updates %
        %%%%%%%%%%%%%%%%%%%%%%
        
        % If we're at the end of the discrete steps, see if there's still a
        % little more left to simulation in continuous time. If so, create
        % the correct time span. If there's no continuous time left, just
        % break (we're done). And if we're not at the end of the discrete
        % steps, just create the right time window.
        if k == length(td)
            if ts(2) - td(k) > epsilon
                tsk = [td(k), ts(2)];
            else
                break;
            end
        else
            tsk = [td(k), td(k+1)];
        end
        
        % Make the output function with the new discrete states.
        if add_logging || ~isempty(options.OutputFcn)
            options.OutputFcn = @(t, y, flag) log_frame(t, y, flag, ...
                                                        ode, yd0, log, ...
                                                        orig_outputfcn);
        end
        if ~isempty(options.Events)
            options.Events = @(t, y) orig_eventfcn(t, y, yd0);
        end
        
        % If there are events, output one way. If not, output the basic
        % way.
        if isempty(options.Events)
            [tk, yck] = solver(@(t, y) ode(t, y, yd0), ...
                               tsk, ...
                               yc0, ...
                               options);
        else
            [tk, yck, tek, yek, iek] = solver(@(t, y) ode(t, y, yd0), ...
                                              tsk, ...
                                              yc0, ...
                                              options);
            if nargout == 1 || nargout >= 3, te  = [te; tek];  end;
            if nargout == 1 || nargout >= 4, yce = [yce; yek]; end;
            if nargout == 1 || nargout >= 6, ie  = [ie; iek];  end;
            if nargout == 1 || nargout >= 5
                for z = 1:length(tek)
                    if iscell(yd0)
                        [yde{end+1, 1:numel(yd0)}] = yd0{:};
                    else
                        yde = [yde; yd0(:).'];
                    end
                end
            end;
        end
        
        % Add the outputs to the list.
        t  = [t;  tk];
        yc = [yc; yck];
        
        % See if the function terminated early.
        if abs(tk(end) - tsk(end)) > 2*eps(tsk(end))
            td = td(1:k);
            if iscell(yd)
                yd = yd(1:k);
            else
                yd = yd(1:k, :);
            end
            break;
        end
        
        % Set the initial step for next time to be the largest step we took
        % during this iteration.
        options.InitialStep = max(diff(tk));
        
        % Get ready for the next step.
        yc0 = yc(end, :)';
        
    end
    
    % Call the output function if there is one.
    if ~isempty(orig_outputfcn)
        orig_outputfcn(ts, yc0, yd0, 'done');
    end

    % If there's only one output, use the structure format.
    if nargout == 1
        sol = struct('t',   t, ...
                     'yc',  yc, ...
                     'td',  td, ...
                     'yd',  yd, ...
                     'te',  te, ...
                     'yce', yce, ...
                     'yde', yde, ...
                     'ie',  ie);
        t = sol;
    end
    
end

% We customize the output function to ignore 'init' and 'done'.
function status = output_fcn(f, t, ycv, yd, flag, yc0)

    % We're going to rely on cell arrays to pass these to the original
    % output function, so make sure they're cell arrays.
    if ~iscell(yc0), yc0 = {yc0}; end;
    if ~iscell(yd),  yd  = {yd};  end;
    
    switch flag
        
        % For 'done', all of the states should be empty.
        case 'done'

            % Pull out the continuous state first (this is how we store it).
            yc = vector_to_state(ycv, yc0);

            for k = 1:length(yc)
                yc{k} = [];
            end
            for k = 1:length(yd)
                yd{k} = [];
            end
            status = f([], yc{:}, yd{:}, flag);

        % For init, we pass the time span along with the initial state.
        case 'init'

            % Pull out the continuous state first.
            yc = vector_to_state(ycv, yc0);
            status = f(t, yc{:}, yd{:}, flag);
    
        % Otherwise, pass along the states and the flag.
        otherwise

            for k = 1:length(t)

                % Pull out the continuous state first.
                yc = vector_to_state(ycv(:, k), yc0);
                status = f(t(k), yc{:}, yd{:}, flag);

            end
            
    end
        
end

% We customize the output function to ignore 'init' and 'done'.
function varargout = event_fcn(f, t, ycv, yd, yc0)

    % We're going to rely on cell arrays to pass these to the original
    % output function, so make sure they're cell arrays.
    if ~iscell(yc0), yc0 = {yc0}; end;
    if ~iscell(yd),  yd  = {yd};  end;
    
    % Pull out the continuous state first.
    yc = vector_to_state(ycv, yc0);
    [varargout{1:nargout}] = f(t, yc{:}, yd{:});

end

% Run the ODE again, this time using the log.
function status = log_frame(t, y, flag, ode, yd, log, f)

    % We never tell integration to stop (but the user's function can
    % override this).
    status = 0;

    % Only call the ODE with the logger on "normal" steps (steps with 
    % neither 'init' nor 'done' flags).
    if isempty(flag)

        % If there's a TimeSeriesLogger, use it.
        if ~isempty(log)
            for k = 1:length(t)
                ode(t(k), y(:, k), yd, log);
            end
        end

        % Call the user's OutputFcn if provided.
        if ~isempty(f)
            status = f(t, y, yd, flag);
        end
        
    end

end

% Run the continuous-discrete-input version of odehybrid.
function varargout = odehybridfull(solver, ode, de, dt, ...
                                   ts, yc0, yd0, varargin)

    % Let the user pass in cells or anything else, but always make sure we
    % work with cells. This simplifies life when we have to dynamically
    % pass these into functions.
    if ~iscell(yc0), yc0 = {yc0}; end;
    if ~iscell(yd0), yd0 = {yd0}; end;
    if ~iscell(de),  de  = {de};  end;
    
    % Create the initial continuous state vector.
    yc0v = state_to_vector(yc0);

    % Create a function to expand the vector into the state, pass the state
    % to the ODE, and turn the result back into a vector.
    ode_v = @(t,yc,yd,varargin) run_ode(ode, t, yc, yd, yc0, varargin{:});
    for k = 1:length(dt)
        de_v{k} = @(t, yc, yd, varargin) run_de(de{k}, t, yc, yd, yc0, ...
                                                varargin{:});
    end

    % Determine what outputs we need.
    if nargout == 1                                     % Structure
        n_outputs = 8;
    elseif nargout <= 1 + numel(yc0)                    % Continuous
        n_outputs = 2;
    elseif nargout <=   1 + numel(yc0) ...              % Continuous
                      + 1 + numel(yd0)                  % Discrete
        n_outputs = 4;
    elseif nargout <=   1 + numel(yc0) ...              % Continuous
                      + 1 + numel(yd0) ...              % Discrete
                      + 1 + numel(yc0) + numel(yd0) + 1 % Events
        n_outputs = 8;
    else
        error('Too many outputs are requested from odoehybrid.');
    end
    
    % Propagate the vector versions of the ODE and DE.
    outputs = cell(1, n_outputs);
    [outputs{:}] = odehybridcore(solver, ...
                                 ode_v, de_v, dt, ts, yc0v, yd0, ...
                                 varargin{:});

	% For states at events, convert to the appropriate type from the big
	% cell array of stored values.
    if n_outputs >= 8
        
        % For the continuous part...
        continuous_states = cell(1, numel(yc0));
        [continuous_states{:}] = vectors_to_states(outputs{6}, yc0{:});

        % Make the output.
        outputs = [outputs(1:5), ...
                   continuous_states, ...
                   output_discrete_states(outputs{7}, yd0), ...
                   outputs(8:end)];
        
    end

	% For discrete states, convert to the appropriate type from the big
	% cell array of stored values.
    if n_outputs >= 4
        outputs = [outputs(1:3), ...
                   output_discrete_states(outputs{4}, yd0), ...
                   outputs(5:end)];
    end
    
    % For continuous state, convert back to original types from the state
    % vectors.
    if n_outputs >= 2
        continuous_states = cell(1, numel(yc0));
        [continuous_states{:}] = vectors_to_states(outputs{2}, yc0{:});
        outputs = [outputs(1), ...
                   continuous_states, ...
                   outputs(3:end)];
    end
    
	% Return separately the states as cell arrays.
    varargout = outputs;

    % If there's only one output, use the structure format.
    if nargout == 1
        
        sol.t   = outputs{1};
        sol.yc  = outputs(1+(1:numel(yc0)));
        sol.td  = outputs{1+numel(yc0)+1};
        sol.yd  = outputs(1+numel(yc0)+1+(1:numel(yd0)));
        sol.te  = outputs{1+numel(yc0)+1+numel(yd0)+1};
        sol.yce = outputs(1+numel(yc0)+1+numel(yd0)+1+(1:numel(yc0)));
        sol.yde = outputs(  1+numel(yc0)+1+numel(yd0)+1+numel(yc0) ...
                          + (1:numel(yd0)));
        sol.ie  = outputs{end};
        varargout{1} = sol;
    end

end

% Convert a state vector to the appropriate state, run the ODE, and
% turn the result back into a vector.
function dyvdt = run_ode(ode, t, ycv, yd, yc0, varargin)

    % Pull out the continuous state first (this is how we store it).
    yc = vector_to_state(ycv, yc0);
    
    % Get the state differences.
    [state_difference{1:numel(yc0)}] = ode(t, yc{:}, yd{:}, varargin{:});
    
    % Convert the derivative to a vector.
    dyvdt = state_to_vector(state_difference);
    
end

% Convert a state vector to the appropriate state, run the ODE, and
% turn the result back into a vector.
function [ycv, yd] = run_de(de, t, ycv, yd, yc0, varargin)

    % Pull out the continuous state first (this is how we store it).
    yc = vector_to_state(ycv, yc0);

    % Run the update.
    [yc{:}, yd{:}] = de(t, yc{:}, yd{:}, varargin{:});
    
    % Convert to vectors.
    ycv = state_to_vector(yc);
    
end

% Turn an n-by-m cell array of discrete states into a 1-by-m cell array of
% discrete state lists each with n entries.
function discrete_states = output_discrete_states(states, yd0)

    % Make space for the individual outputs.
    discrete_states = cell(1, numel(yd0));

    % For each output...
    for k = 1:numel(yd0)

        % Convert to a matrix with rows containing vectors.
        if isnumeric(yd0{k}) || ischar(yd0{k})

            dim = find([size(yd0{k}), 1] == 1, 1, 'first');
            if dim == 2
                discrete_states{k} = cat(dim, states{:, k}).';
            else
                discrete_states{k} = cat(dim, states{:, k});
            end

        % Convert to an array of structs.
        elseif isstruct(yd0{k}) && length(yd0{k}) == 1

            discrete_states{k} = cellfun(@(v) v, states(:, k));

        % Give up and just use the cell array of states.
        else

            discrete_states{k} = states(:, k);

        end

    end
        
end
