classdef TimeSeriesLogger < handle

% TimeSeriesLogger class
% 
% A class for keeping track of numerous pieces of data over time throughout
% a process. It's designed to be easy to use and relatively quickly. It
% logs numeric types of any size (scalars, vectors, or matrices), as long
% as the size is consistent from sample to sample.
% 
% Example:
% 
%   log = TimeSeriesLogger(); % Make a new logger.
%   x   = [0; 0];
%   for t = 1:100             % Make a random walk.
%     x = x + randn(2, 1);    % Take a single step.
%     log.add('walk', t, x);  % Log a single sample, x, at time t.
%   end
%   log.plot();               % Plot everything.
% 
% We can also access specific logs by their names. In the above, we only
% have one log ('walk'). Let's get that log from the logger.
% 
%   [t_log, x_log] = log.get_log('walk'); % Get a specific log.
%   plot(x_log(:, 1), x_log(:, 2));       % Do something with it.
%
% We can make sure a log exists before trying to do anything with it:
% 
%   if log.contains('walk')
%     x_log = log.get_log('walk');
%     plot(x_log(:, 1), x_log(:, 2));
%   end
%
% If we want to log something but don't want it plotted when we call
% |plot|, then we can pass in false to the logger when we add the signal.
% 
%   log.add('var1', t, x, false);
% 
% To group items into different figures when they're plotted, give them
% "groups":
% 
%   log.add('var1', t, x, true, 'group1');
%   log.add('var2', t, y, true, 'group2');
%   log.add('foo',  t, bar, true, 'group1');
% 
% All of the items with common group names will be added to the same
% figures as subplots.
%
% Finally, we can clear out the logger with |initialize|. This deletes all
% data and returns the logger to its initial state.
%
%   log.initialize();
% 
% For more information on the various methods, see the individual methods
% with |doc TimeSeriesLogger|:
% 
% <a href="matlab:doc TimeSeriesLogger;">doc TimeSeriesLogger</a>
% 
% Methods
% -------
% 
% initialize  Clear out all data.
% add         Add a single data point to a time series.
% contains    See if a time series exists (by name).
% get_log     Return a log by name.
% plot        Plot the logged series.
%
% Notes
% -----
%
% TimeSeriesLogger was created as the logging mechanism for odehybrid.
% However, it's not dependent on that project and so can be used anywhere.
%  
% TimeSeriesLogger is not related to the 'timeseries' class in MATLAB.
%
% Since this is such a simple class, all properties are public, preventing
% the need for getters and setters. However, under normal use cases, there
% would be no reason to access any properties of the class.
% 
% Since this class never knows how much more data is going to be logged, it
% can't preallocate the appropriate amount of space. However, increasing
% the size of its data stores by one row every time a new sample is added
% is very slow. To combat this, the data store starts off small and doubles
% whenever the store is at capacity. Say we're logging a 4x1 signal. When
% the first sample is added (say it's log k), data{k}{2} will be 1x4. When
% the second signal is added, it becomes 2x4. For the third, it's 4x4, then
% 8x4, 16x4, etc. A separate counter stores how much of this allocated 
% space is currently used. This reduces the number of allocations from n to
% log2(n). Practically, it saves a little time during logging without too
% much complexity.
%
% Online doc: http://www.anuncommonlab.com/doc/odehybrid/TimeSeriesLogger.html
%
% Copyright 2014 An Uncommon Lab

    properties
        names = {}; % Unique identifier (within the group)
        data  = {}; % Data, stored as {t, x}
        show  = {}; % True iff the data should be plotted
        sizes = {}; % Amount of allocated space
        count = {}; % Amount of allocated space currently used
        group = {}; % Plot group identifier
    end
    
    methods
        
        function initialize(this)
        % Clear everything out and start fresh.
            this.names = {};
            this.data  = {};
            this.show  = {};
            this.sizes = {};
            this.count = {};
            this.group = {};
        end
        
        function added_to_list = add(this, name, t, x, show_it, group)
        % Add a new log or append to an existing log by name.
        % 
        % log.add('var1', t, x);
        % log.add('var1', t, x, true);  % Show log in plots (default)
        % log.add('var1', t, x, false); % Don't show log in plots
        % 
        % % Signals can be grouped together into figures by given them
        % a common group argument. Here, both var1 and var2 logs will
        % be plotted together.
        % log.add('var1', t,  x, true, 'group1');
        % log.add('var2', t2, y, true, 'group1');
        % 
        % Returns true iff a new log was created.
        %
        % The signals are stored as time-data pairs, {t, x}, where each
        % row of x corresponds to each row of t. So data{k}{1} contains
        % the logged times for the kth signal, and data{k}{2} contains
        % the data for the logged signal.
            
            % Defaults
            if nargin < 5, show_it = true; end; % Plot by default.
            if nargin < 6, group = '';     end; % Default group is ''.
            
            % See if we've already added this variable.
            index = find(strcmp(this.names, name));
            
            % We'll add it to the list if it's not there (and we return
            % this status).
            added_to_list = isempty(index);
            
            % If not, add it to the list.
            if added_to_list
                
                this.data{end+1}  = {t, x(:).'};
                this.names{end+1} = name;
                this.show{end+1}  = show_it;
                this.sizes{end+1} = 1;
                this.count{end+1} = 1;
                this.group{end+1} = group;
                
            % Otherwise, append.
            else
                
                % If we've used up all reserved space, double it.
                c = this.count{index} + 1;
                if c > this.sizes{index}
                    this.data{index}{1} = [this.data{index}{1}; ...
                                           t; ...
                                           zeros(this.sizes{index}-1, 1)];
                    this.data{index}{2} = [this.data{index}{2}; ...
                                           x(:).'; ...
                                           zeros(this.sizes{index}-1, ...
                                                 length(x(:).'))];
                    this.sizes{index} = 2 * this.sizes{index};
                else
                    this.data{index}{1}(c)    = t;
                    this.data{index}{2}(c, :) = x(:).';
                end
                this.count{index} = c;
                
                % Without doubling:
                % this.data{index}{1} = [this.data{index}{1}; t];
                % this.data{index}{2} = [this.data{index}{2}; x(:).'];
                
                % We should show this data if the user has ever asked to
                % plot it.
                this.show{index} = show_it || any(this.show{index});
                this.group{index} = group;
                
            end
            
        end
        
        function result = contains(this, name)
        % Return true iff the log contains this name.
            result = any(strcmp(this.names, name));
        end
        
        function varargout = get_log(this, name)
        % Return a specific log by name.
        % 
        % If one output is request, it returns the data. If two are 
        % requested, it returns the time and the data. Returns empty if
        % the logs don't exist (use the 'contains' function to test for
        % this).
        % 
        % x = log.get_log('signal1');
        % [t, x] = log.get_log('signal1');
        %
        % Note that |t| will be ns-by-1 and x will be ns-by-nx, where
        % nx is the number of elements in a single sample.
            
            index = find(strcmp(this.names, name));
            if isempty(index)
                varargout = cell(1, nargout);
                return;
            else
                t = this.data{index}{1}(1:this.count{index});
                x = this.data{index}{2}(1:this.count{index}, :);
                % Without doubling:
                % t = this.data{index}{1};
                % x = this.data{index}{2};
            end
            if nargout == 1
                varargout = {x};
            elseif nargout == 2
                varargout = {t, x};
            end
        end
        
        function plot(this, x_label)
        % Plot all of the signals, grouped appropriately into figures.
        %
        % A custom x label can be added to figures as well (second input).
            
            % Get the unique groups.
            groups = unique(this.group);
            
            % For each group...
            for g = 1:length(groups)

                % See what is to be shown and count them.
                to_show = find(  strcmp(this.group, groups{g}) ...
                               & [this.show{:}]);
                n = length(to_show);

                % Bail if there's no need to show anything.
                if n == 0, continue; end;

                % Set up the figure.
                unique_figure(sprintf('LoggerDefaultPlot%d', g), ...
                              'Name', [groups{g} ' Log']);
                clf();
                for k = 1:n
                    subplot(n, 1, k);
                    plot(this.data{to_show(k)}{1}(1:this.count{to_show(k)}), ...
                         this.data{to_show(k)}{2}(1:this.count{to_show(k)}, :));
                    % plot(this.data{to_show(k)}{1}, ...
                    %      this.data{to_show(k)}{2});
                    ylabel(this.names{to_show(k)});
                end
                if nargin >= 2, xlabel(x_label); end;
            
            end
                
        end
        
    end
    
end

function h = unique_figure(id, varargin)
% h = unique_figure(id, varargin)
% 
% Creates a figure with the given ID (tag) or selects it if it already
% exists, allowing one to easily reuse the same figure window identified
% with text instead of handles. This is useful when, e.g., running a script
% many times after clearing between runs without having to hard-code figure
% numbers (which can become hard to keep track of).
% 
% Any additional arguments are passed along to the figure's |set| method.
%
% Example:
% 
%   h = unique_figure('trajectory');
%
% Example with extra arguments:
% 
%   h = unique_figure('trajectory', 'Name', 'Trajectory');

    % See if there's a figure with this ID.
    h = findobj('Type', 'figure', 'Tag', id);
    
    % If not, make it, passing along any options for the figure.
    if isempty(h)
        
        h = figure('Tag', id, varargin{:});
        
    else
        
        % If we found it, select it.
        figure(h);
        
        % If there were any additional arguments, pass them along. Note
        % that if there weren't and we called set(h), then it would print h
        % to the command window -- probably not what we want, hence the
        % condition.
        if nargin > 1
            set(h, varargin{:});
        end
        
    end

end
