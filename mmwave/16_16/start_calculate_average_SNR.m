% =========================================================================
% -- Scenario S1 - mmWave Average SNR Post-Processing
% -------------------------------------------------------------------------
% Description:
%   Aggregate saved results for the mmWave scenario with B=16 receive
%   antennas, U=16 UEs, and T=2 time slots.
%
% Input:
%   Configured Main_simulator path and directory containing saved
%   simulation MAT files.
%
% Output:
%   MAT file containing average SNR statistics.
%
% -- (c) 2026 Victoria Palhares, Christoph Studer
% =========================================================================

clear all;clc;close all;

par.plot_type = 'CDF'; % CDF/CCDF simulation configuration
par.system_type = 'mmwave'; % 'mmwave' or 'cell_free'
par.sim_scenario = 'a';
par.main_path = '/path/to/ue-scheduling/Main_simulator/';
addpath(par.main_path);
par = param_config(par);
addpath([par.main_path,'results/']);

calculate_average_SNR(par,'/path/to/results/mmwave_16_16/');

close all;

