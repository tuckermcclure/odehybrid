% A simple example of using structs with odehybrid.
function example_odehybrid_structs()

    dt = 0.1;                    % Discrete eq. time step
    ts = [0 5];                  % Simulation time
    x0 = struct('p', 1, 'v', 0); % Initial continuous states
    d0 = struct('u', 0, 'i', 0); % Initial discrete states
    
    % Simulate.
    [t, sig, td, ctrl] = odehybrid(@ode45, @ode, @de, dt, ts, x0, d0);
    
    % Plot.
    plot(t, [sig.p], t, [sig.v], td, [ctrl.u], '.', td, [ctrl.i], '.');
    xlabel('Time');
    legend('p', 'v', 'u', '\int p(t)/2 dt');
    
end

% Continuous differential equation
function dsdt = ode(t, signal, controller) %#ok<INUSL>
    dsdt.p = signal.v;
    dsdt.v = 2 * signal.p + controller.u + 1;
end

% Discrete update equation
function [signal, controller] = de(t, signal, controller) %#ok<INUSL>
    controller.u = -8 * signal.p - 4*signal.v - controller.i;
    controller.i = controller.i + 0.5 * signal.p;
end
