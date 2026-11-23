% =========================================================================
% -- T=8 mmWave Result Post-Processing
% -------------------------------------------------------------------------
% Description:
%   Post-process saved results for the mmWave experiment with B=32 receive
%   antennas, U=64 UEs, T=8 time slots, and Tmax=4.
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
par.system_type = 'mmwave'; % 'mmwave' or 'cell_free'
par.sim_scenario = 'd';
par.main_path = '/path/to/ue-scheduling/Main_simulator/';
addpath(par.main_path);
par = param_config(par);
addpath([par.main_path,'results/']);

plot_figures_cdf(par,'/path/to/results/mmwave_32_64_T_8_Tmax_4/');
