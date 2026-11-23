% =========================================================================
% -- Average SNR Aggregation
% -------------------------------------------------------------------------
% Description:
%   Load saved simulation results and compute the average SNR over the
%   selected scheduling methods and time slots.
%
% Input:
%   par: Simulation parameter structure.
%   path: Directory containing saved simulation result files.
%
% Output:
%   None. Aggregated SNR statistics are saved to a MAT file.
%
% -- (c) 2026 Victoria Palhares, Christoph Studer
% =========================================================================

function calculate_average_SNR(par,path)

my_path = pwd;

index_process = [];
for type=1:length(par.scheduling_chosen)
    scheduler_idx = find(strcmp(par.scheduling_options,par.scheduling_chosen{type}));
    for c = 1:par.channels
        for int_snr=1:length(par.SNRdB_list)
            index = (c-1)*length(par.SNRdB_list)+int_snr;
            index_process = [index_process;((scheduler_idx-1)*(length(par.SNRdB_list)*par.channels)+index)-1];
        end
    end
end

% Locate result files for all simulation jobs
for i=1:length(index_process)
    file_name = [path,['USER_SCHEDULING_',par.system_type,'_',num2str(par.B),'_',num2str(par.U),'_T_',num2str(par.timeslots),'_',num2str(index_process(i)),'.mat']];
    file{i} = file_name;
end

res.SNR = cell(length(par.scheduling_options),1);
res.average_SNR_per_timeslot = cell(length(par.scheduling_options),1);
res.average_SNR = cell(length(par.scheduling_options),1);

% Aggregate SNR data across methods and jobs
for n_proc=1:length(index_process)
    load(file{n_proc},'results');

    ind = find(strcmp(par.scheduling_options,results.scheduling_method));

    res.SNR{ind} = [res.SNR{ind};results.SNR_dB.'];
end

for i=1:length(par.scheduling_chosen)

    ind = find(strcmp(par.scheduling_options,par.scheduling_chosen{i}));

    res.average_SNR_per_timeslot{ind} = mean(res.SNR{ind},1);
    res.average_SNR{ind} = mean(res.average_SNR_per_timeslot{ind});
end

FileName=[datestr(now, 'yyyy_mmm_dd'),'_USER_SCHEDULING_',num2str(par.B),'_',num2str(par.U),'_',num2str(par.timeslots),'_SNR.mat'];
save (FileName);
