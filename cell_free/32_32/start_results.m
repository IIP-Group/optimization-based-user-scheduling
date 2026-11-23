% =========================================================================
% -- Scenario S5 - Cell-Free Result Post-Processing
% -------------------------------------------------------------------------
% Description:
%   Post-process saved results for the cell-free scenario with B=32 receive
%   antennas, U=32 UEs, and T=2 time slots.
%
% Input:
%   Configured Main_simulator path and directory containing saved
%   simulation MAT files.
%
% Output:
%   CDF/CCDF figures and 90th-percentile summary MAT files.
%
% -- (c) 2026 Victoria Palhares, Christoph Studer
% =========================================================================

clear all;clc;close all;

par.plot_type = 'CDF'; % CDF/CCDF post-processing
par.system_type = 'cell_free'; % 'mmwave' or 'cell_free'
par.sim_scenario = 'b';
par.main_path = '/path/to/ue-scheduling/Main_simulator/';
addpath(par.main_path);
par = param_config(par);
addpath([par.main_path,'results/']);

plot_figures_cdf(par,'/path/to/results/cell_free_32_32/');
