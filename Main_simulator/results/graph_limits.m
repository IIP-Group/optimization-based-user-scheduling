% =========================================================================
% -- Result Plot Limits
% -------------------------------------------------------------------------
% Description:
%   Return the plotting grids and axis limits used for each system type and
%   simulation scenario.
%
% Input:
%   par: Simulation parameter structure.
%
% Output:
%   limits: Structure containing plotting grids and axis limits.
%
% -- (c) 2026 Victoria Palhares, Christoph Studer
% =========================================================================

function [limits] = graph_limits(par)

switch(par.system_type)
    case 'mmwave'
        switch(par.sim_scenario)
            case 'a'
                limits.ber.xgrid = logspace(-3,0,100);
                limits.ber.axis = [1e-3 1 0 1];
                limits.mi.xgrid = linspace(0,2.5,100);
                limits.mi.axis = [0 2.5 0 1];
                limits.mse.xgrid = linspace(-10,5,100);
                limits.mse.axis = [-10 5 0 1];
                limits.rmsse.xgrid = linspace(0,100,100);
                limits.rmsse.axis = [0 100 0 1];
                limits.sinr.xgrid = linspace(0,15,100);
                limits.sinr.axis = [0 15 0 1];
                limits.rate.xgrid = linspace(0,3,100);
                limits.rate.axis = [0 3 0 1];
            case 'b'
                limits.ber.xgrid = logspace(-3,0,100);
                limits.ber.axis = [1e-3 1 0 1];
                limits.mi.xgrid = linspace(0,2.5,100);
                limits.mi.axis = [0 2.5 0 1];
                limits.mse.xgrid = linspace(-10,5,100);
                limits.mse.axis = [-10 5 0 1];
                limits.rmsse.xgrid = linspace(0,100,100);
                limits.rmsse.axis = [0 100 0 1];
                limits.sinr.xgrid = linspace(0,15,100);
                limits.sinr.axis = [0 15 0 1];
                limits.rate.xgrid = linspace(0,3,100);
                limits.rate.axis = [0 3 0 1];
            case 'c'
                limits.ber.xgrid = logspace(-2,0,100);
                limits.ber.axis = [1e-2 1 0 1];
                limits.mi.xgrid = linspace(0,2,100);
                limits.mi.axis = [0 2 0 1];
                limits.mse.xgrid = linspace(-10,5,100);
                limits.mse.axis = [-10 5 0 1];
                limits.rmsse.xgrid = linspace(0,100,100);
                limits.rmsse.axis = [0 100 0 1];
                limits.sinr.xgrid = linspace(0,15,100);
                limits.sinr.axis = [0 15 0 1];
                limits.rate.xgrid = linspace(0,1,100);
                limits.rate.axis = [0 1 0 1];
            case 'd'
                limits.ber.xgrid = logspace(-2,0,100);
                limits.ber.axis = [1e-2 1 0 1];
                limits.mi.xgrid = linspace(0,2,100);
                limits.mi.axis = [0 2 0 1];
                limits.mse.xgrid = linspace(-10,5,100);
                limits.mse.axis = [-10 5 0 1];
                limits.rmsse.xgrid = linspace(0,100,100);
                limits.rmsse.axis = [0 100 0 1];
                limits.sinr.xgrid = linspace(0,15,100);
                limits.sinr.axis = [0 15 0 1];
                limits.rate.xgrid = linspace(0,1,100);
                limits.rate.axis = [0 1 0 1];
            otherwise
                error('Unsupported simulation scenario.');
        end
    case 'cell_free'
        switch(par.sim_scenario)
            case 'a'
                limits.ber.xgrid = logspace(-6,0,100);
                limits.ber.axis = [1e-6 1 0 1];
                limits.mi.xgrid = linspace(0,7,100);
                limits.mi.axis = [0 7 0 1];
                limits.mse.xgrid = linspace(-40,5,100);
                limits.mse.axis = [-40 5 0 1];
                limits.rmsse.xgrid = linspace(0,100,100);
                limits.rmsse.axis = [0 100 0 1];
                limits.sinr.xgrid = linspace(0,15,100);
                limits.sinr.axis = [0 15 0 1];
                limits.rate.xgrid = linspace(0,10,100);
                limits.rate.axis = [0 10 0 1];
            case 'b'
                limits.ber.xgrid = logspace(-6,0,100);
                limits.ber.axis = [1e-6 1 0 1];
                limits.mi.xgrid = linspace(0,7,100);
                limits.mi.axis = [0 7 0 1];
                limits.mse.xgrid = linspace(-40,0,100);
                limits.mse.axis = [-40 0 0 1];
                limits.rmsse.xgrid = linspace(0,100,100);
                limits.rmsse.axis = [0 100 0 1];
                limits.sinr.xgrid = linspace(0,15,100);
                limits.sinr.axis = [0 15 0 1];
                limits.rate.xgrid = linspace(0,10,100);
                limits.rate.axis = [0 10 0 1];
            case 'c'
                limits.ber.xgrid = logspace(-6,0,100);
                limits.ber.axis = [1e-6 1 0 1];
                limits.mi.xgrid = linspace(0,7,100);
                limits.mi.axis = [0 7 0 1];
                limits.mse.xgrid = linspace(-40,0,100);
                limits.mse.axis = [-40 0 0 1];
                limits.rmsse.xgrid = linspace(0,100,100);
                limits.rmsse.axis = [0 100 0 1];
                limits.sinr.xgrid = linspace(0,15,100);
                limits.sinr.axis = [0 15 0 1];
                limits.rate.xgrid = linspace(0,6,100);
                limits.rate.axis = [0 6 0 1];
            otherwise
                error('Unsupported simulation scenario.');
        end
end
