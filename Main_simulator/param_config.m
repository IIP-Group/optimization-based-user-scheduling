% =========================================================================
% -- Simulation Parameter Configuration
% -------------------------------------------------------------------------
% Description:
%   Configure system dimensions, scheduling methods, channel estimation,
%   modulation, noise, power control, and job indexing for the selected
%   scenario.
%
% Input:
%   par: Simulation parameter structure.
%
% Output:
%   par: Updated simulation parameter structure.
%
% -- (c) 2026 Victoria Palhares, Christoph Studer
% =========================================================================

function par = param_config(par)


par.trials = 1e5; % Number of Monte Carlo trials
switch (par.plot_type)
    case 'CDF'
        par.SNRdB_list = 0; % SNR values in dB
        par.channels = 1e3; % Number of channel realizations
    case 'SNR'
        par.SNRdB_list = 5:5:25; % SNR values in dB
        par.channels = 1e2; % Number of channel realizations
end
par.scheduling_options = {'Random','ES-MSE','ES-Capacity','ES-Rate','MMSE','Capacity','Rate',...
    'SUS','CSS','Greedy','No Scheduling','LoFi (K=1)','LoFi (K=4)','LoFi++ (K=1)','LoFi++ (K=4)'};
switch(par.sim_scenario)
    case 'a'
        par.scheduling_chosen = {'Random','ES-MSE','ES-Capacity','ES-Rate','MMSE','Capacity','Rate',...
            'SUS','CSS','Greedy','No Scheduling','LoFi (K=1)','LoFi (K=4)','LoFi++ (K=1)','LoFi++ (K=4)'};
        % System parameters
        par.B = 16; % Number of receive antennas
        par.U = 16; % Number of UEs
        par.timeslots = 2; % Number of time slots
        par.min_sum_columns = 8; % U_min
        par.max_sum_columns = 10; % U_max
        par.min_sum_rows = 1; % T_min
        par.max_sum_rows = 2; % T_max
    case 'b'
        par.scheduling_chosen = {'Random','MMSE','Capacity','Rate',...
    'SUS','CSS','Greedy','No Scheduling','LoFi (K=1)','LoFi (K=4)','LoFi++ (K=1)','LoFi++ (K=4)'};
        % System parameters
        par.B = 32; % Number of receive antennas
        par.U = 32; % Number of UEs
        par.timeslots = 2; % Number of time slots
        par.min_sum_columns = 16; % U_min
        par.max_sum_columns = 20; % U_max
        par.min_sum_rows = 1; % T_min
        par.max_sum_rows = 2; % T_max
    case 'c'
        par.scheduling_chosen = {'Random','MMSE','Capacity','Rate',...
    'SUS','CSS','Greedy','No Scheduling','LoFi (K=1)','LoFi (K=4)'};
        % System parameters
        par.B = 32; % Number of receive antennas
        par.U = 64; % Number of UEs
        par.timeslots = 4; % Number of time slots
        par.min_sum_columns = 16; % U_min
        par.max_sum_columns = 20; % U_max
        par.min_sum_rows = 1; % T_min
        par.max_sum_rows = 2; % T_max
    case 'd'
        par.scheduling_chosen = {'Random','MMSE','Capacity','Rate',...
            'SUS','CSS','Greedy','No Scheduling','LoFi (K=1)','LoFi (K=4)'};
        % System parameters
        par.B = 32; % Number of receive antennas
        par.U = 64; % Number of UEs
        par.timeslots = 8; % Number of time slots
        par.min_sum_columns = 8; % U_min
        par.max_sum_columns = 10; % U_max
        par.min_sum_rows = 1; % T_min
        par.max_sum_rows = 4; % T_max
    otherwise
        error('Unsupported simulation scenario.');
end
par.total_jobs = par.channels*length(par.SNRdB_list)*length(par.scheduling_chosen);

% Enable exhaustive search only when selected
if(par.U <= 24)
    par.status_es = 'on';
else
    par.status_es = 'off';
end

par.Us = par.min_sum_columns;

% Power control
par.pwr_cntrl = 'on';
par.P_db = 6; % Maximum receive-power dynamic range [dB]
% Configure the post-equalization evaluation metric
par.evaluation = 'mmse';
% Configure channel estimation
switch(par.system_type)
    case 'mmwave'
        par.channel_estimation = 'BEACHES'; % Perfect, BEACHES or LS
    case 'cell_free'
        par.channel_estimation = 'LS'; % Perfect, BEACHES or LS
end

% Configure modulation and noise power
switch (par.plot_type)
    case 'CDF'
        switch(par.system_type)
            case 'mmwave'
                par.bandwidth = 100e6; % Bandwidth [Hz]
                % Supported modulation: BPSK, QPSK, 16QAM, or 64QAM
                par.modulation = 'QPSK';
                % Configure transmit power
                par.transmit_power = 100e-3; % Transmit power: 100 mW
            case 'cell_free'
                par.bandwidth = 100e6; % Bandwidth [Hz]
                % Supported modulation: BPSK, QPSK, 16QAM, or 64QAM
                par.modulation = '64QAM';
                par.transmit_power = 10e-3; % Transmit power: 10 mW
        end
        par.T0 = 290; % Noise temperature [K]
        par.noise_figure_db = 9; % Noise figure [dB]
        par.noise_figure = 10^(par.noise_figure_db/10); % Noise figure
        % Thermal-noise power [W]
        par.Np = par.bandwidth*1.381e-23*par.T0*par.noise_figure;
        par.rho_u = par.transmit_power/par.Np;
        par.N0 = 1; % Unit noise variance after channel/SNR normalization
end

% Save the parameter configuration
FileName=['parameters.mat'];
save (FileName);
