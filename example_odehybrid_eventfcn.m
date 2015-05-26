% A simple example Event function used with examples_odehybrid.
function [h, t, d] = example_odehybrid_eventfcn(t, x, p, v)
    h = x(2); % An "event" occurs when the velocity is 0.
    t = true; % Terminate on event; stop the simulation when h=0.
    d = 1;    % Trigger when going positive from negative
end
