% =========================================================================
% -- Optimization Hyperparameters
% -------------------------------------------------------------------------
% Description:
%   Set the iteration count, step-size configuration, random
%   initializations, and regularization weight for the selected scenario
%   and objective.
%
% Input:
%   var: Simulation state structure containing the channel, noise, and
%   scheduler variables.
%   par: Simulation parameter structure.
%
% Output:
%   var: Updated simulation state structure.
%
% -- (c) 2026 Victoria Palhares, Christoph Studer
% =========================================================================

function var = param_optimization(var,par)

var.grad_iterations = 100; % Number of gradient-descent iterations
var.adaptive_stepsize = 'on';
var.tau = 1e-2;
switch(par.system_type)
    case 'mmwave'
        switch(par.sim_scenario)
            case 'a' % B = 16, U = 16, T=2
                % Number of random initializations
                var.num_initializations = 80;
                switch(var.scheduling_method)
                    case 'MMSE'
                        var.lambda = 1e-1;
                    case 'Capacity'
                        var.lambda = 10;
                    case 'Rate'
                        var.lambda = 1;
                    otherwise
                        error('Unsupported scheduling method.');
                end
            case 'b' % B = 32, U = 32, T=2
                % Number of random initializations
                var.num_initializations = 10;
                switch(var.scheduling_method)
                    case 'MMSE'
                        var.lambda = 1e-1;
                    case 'Capacity'
                        var.lambda = 5;
                    case 'Rate'
                        var.lambda = 1;
                    otherwise
                        error('Unsupported scheduling method.');
                end
            case 'c' % B = 32, U = 64, T=4
                % Number of random initializations
                var.num_initializations = 3;
                switch(var.scheduling_method)
                    case 'MMSE'
                        var.lambda = 1;
                    case 'Capacity'
                        var.lambda = 1;
                    case 'Rate'
                        var.lambda = 1e-1;
                    otherwise
                        error('Unsupported scheduling method.');
                end
            case 'd' % B = 32, U = 64, T=8, T_max=4
                % Number of random initializations
                var.num_initializations = 3;
                switch(var.scheduling_method)
                    case 'MMSE'
                        var.lambda = 1;
                    case 'Capacity'
                        var.lambda = 1;
                    case 'Rate'
                        var.lambda = 1e-1;
                    otherwise
                        error('Unsupported scheduling method.');
                end
            otherwise
                error('Unsupported simulation scenario.');
        end
    case 'cell_free'
        switch(par.sim_scenario)
            case 'a' % B = 16, U = 16, T=2
                % Number of random initializations
                var.num_initializations = 80;
                switch(var.scheduling_method)
                    case 'MMSE'
                        var.lambda = 10;
                    case 'Capacity'
                        var.lambda = 5;
                    case 'Rate'
                        var.lambda = 5;
                    otherwise
                        error('Unsupported scheduling method.');
                end
            case 'b' % B = 32, U = 32, T=2
                % Number of random initializations
                var.num_initializations = 10;
                switch(var.scheduling_method)
                    case 'MMSE'
                        var.lambda = 1e-5;
                    case 'Capacity'
                        var.lambda = 20;
                    case 'Rate'
                        var.lambda = 1;
                    otherwise
                        error('Unsupported scheduling method.');
                end
            case 'c' % B = 32, U = 64, T=4
                % Number of random initializations
                var.num_initializations = 3;
                switch(var.scheduling_method)
                    case 'MMSE'
                        var.lambda = 1e-3;
                    case 'Capacity'
                        var.lambda = 1;
                    case 'Rate'
                        var.lambda = 1e-2;
                    otherwise
                        error('Unsupported scheduling method.');
                end
            otherwise
                error('Unsupported simulation scenario.');
        end
    otherwise
        error('Unsupported system type.');

end
end
