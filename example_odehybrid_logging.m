% A simple example of using logging with odehybrid.
function example_odehybrid_logging()

    dt  = 0.1;                    % Discrete eq. time step
    ts  = [0 5];                  % Simulation time
    x0  = struct('p', 1, 'v', 0); % Initial continuous states
    d0  = struct('u', 0, 'i', 0); % Initial discrete states
    log = TimeSeriesLogger();     % Create the logger.
    
    % Simulate.
    [t, sig, td, ctrl] = odehybrid(@ode45, @ode, @de, dt, ts, x0, d0, ...
                                   [], log);
    
    % Plot.
    plot(t, [sig.p], t, [sig.v], td, [ctrl.u], '.', td, [ctrl.i], '.');
    xlabel('Time');
    legend('p', 'v', 'u', '\int p(t)/2 dt');
    
    % Add log output.
    log.plot();
    xlabel('Time');
    
end

% Continuous differential equation
function dsdt = ode(t, signal, controller, log)

    % Calculate the derivatives.
    dsdt.p = signal.v;
    dsdt.v = 2 * signal.p + controller.u + 1;
    
    % Log the acceleration. We *must* check to see if the log is passed in;
    % it won't always be passed in.
    if nargin >= 4
        log.add('acceleration', t, dsdt.v);
    end
    
end

% Discrete update equation
function [signal, controller] = de(t, signal, controller, log)

    % Update the discrete state.
    controller.u = -8 * signal.p - 4*signal.v - controller.i;
    controller.i = controller.i + 0.5 * signal.p;
    
    % Log the velocity as it was sampled. Logs are always passed to the
    % discrete update functions, so we don't explicitly need to check.
    log.add('sampled v', t, signal.v);
    
end
