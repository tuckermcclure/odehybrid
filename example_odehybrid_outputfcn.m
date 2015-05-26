% A simple example OutputFcn used with examples_odehybrid.
function status = example_odehybrid_outputfcn(t, x, p, v, flag)

    % Return true to terminate the propagation.
    status = 0;
    
    switch flag
        
        % On init, we receive the time span and init. states.
        case 'init'
            fprintf('Simulating from %fs to %fs.\n', t);
            
        % When done, times and states are all empty.
        case 'done'
            fprintf('Finished the simulation.\n');
            
        % Otherwise, we receive 1 or more samples.
        otherwise
            status = x(1) < 0; % End the simulation with x(1) < 0.
            
    end
    
end
