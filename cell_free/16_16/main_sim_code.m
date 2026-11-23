% =========================================================================
% -- Scenario S4 - Cell-Free Simulation Runner
% -------------------------------------------------------------------------
% Description:
%   Run the cell-free scenario with B=16 receive antennas, U=16 UEs, and
%   T=2 time slots.
%
% Input:
%   proc_ID: Zero-based process identifier, supplied as a numeric value
%   or string.
%
% Output:
%   None. Simulation jobs save their results to MAT files.
%
% -- (c) 2026 Victoria Palhares, Christoph Studer
% =========================================================================

function main_sim_code(proc_ID)
if ischar(proc_ID) || isstring(proc_ID)
    proc_ID = str2double(proc_ID);
end
if ~isscalar(proc_ID) || isnan(proc_ID) || proc_ID < 0 || proc_ID ~= floor(proc_ID)
    error('proc_ID must be a nonnegative integer.');
end
tic;

par.plot_type = 'CDF'; % CDF/CCDF simulation configuration
par.system_type = 'cell_free'; % 'mmwave' or 'cell_free'
par.sim_scenario = 'a';
par = param_config(par);

%% Compute the total number of jobs
par.parameters = [];
% Enumerate simulation jobs
for type=1:length(par.scheduling_chosen)
    for c = 1:par.channels
        for int_snr=1:length(par.SNRdB_list)
            par.parameters = [par.parameters;c,par.SNRdB_list(int_snr),string(par.scheduling_chosen{type})];
        end
    end
end

%% Run the simulation job
main_function(par,proc_ID);
fprintf('Process %d finished at %s.\n', ...
    proc_ID,datestr(now,'HH:MM:SS.FFF'));

end
