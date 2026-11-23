% =========================================================================
% -- CDF and CCDF Result Generation
% -------------------------------------------------------------------------
% Description:
%   Aggregate simulation-job result files, generate the paper CDF/CCDF
%   figures, and compute percentile summaries.
%
% Input:
%   par: Simulation parameter structure.
%   path: Directory containing saved simulation result files.
%
% Output:
%   None. Figures and summary MAT files are written to disk.
%
% -- (c) 2026 Victoria Palhares, Christoph Studer
% =========================================================================

function plot_figures_cdf(par,path)

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
    file_name = [path,['USER_SCHEDULING_',par.system_type,'_',num2str(par.B),'_',num2str(par.U),'_T_',num2str(par.frequency_slots),'_',num2str(index_process(i)),'.mat']];
    file{i} = file_name;
end

res.H = zeros(par.B,par.U,par.total_jobs);
res.C = zeros(par.U,par.frequency_slots,par.total_jobs);

res.BER_per_user = cell(length(par.scheduling_options),1);
res.MI_per_user = cell(length(par.scheduling_options),1);
res.MSE_per_user = cell(length(par.scheduling_options),1);
res.RMSSE_per_user = cell(length(par.scheduling_options),1);
res.rate_per_user = cell(length(par.scheduling_options),1);
res.SINR_per_user = cell(length(par.scheduling_options),1);

res.real_H_norm = cell(length(par.scheduling_options),1);
res.status = cell(length(par.scheduling_options),1);
res.final_lambda = cell(length(par.scheduling_options),1);

mmse_idx = find(strcmp(par.scheduling_options,'MMSE'));
rate_idx = find(strcmp(par.scheduling_options,'Rate'));

% Aggregate per-UE metrics across methods and jobs
for n_proc=1:length(index_process)
    load(file{n_proc},'results');

    ind = find(strcmp(par.scheduling_options,results.scheduling_method));

    % BER samples
    res.BER_per_user{ind} = [res.BER_per_user{ind};results.ber_per_user];
    % MI samples
    res.MI_per_user{ind} = [res.MI_per_user{ind};results.mi_per_user];
    % MSE samples
    res.MSE_per_user{ind} = [res.MSE_per_user{ind};results.mse_per_user];
    % RMSSE samples
    res.RMSSE_per_user{ind} = [res.RMSSE_per_user{ind};results.rmsse_per_user];
    % Rate samples
    res.rate_per_user{ind} = [res.rate_per_user{ind};results.rate_per_user];
    % SINR samples
    res.SINR_per_user{ind} = [res.SINR_per_user{ind};results.sinr_per_user];

    scheduler_idx = find(strcmp(par.scheduling_options,results.scheduling_method));

    if(scheduler_idx >= mmse_idx && scheduler_idx <= rate_idx)
        load(file{n_proc},'status');
        res.status{ind} = [res.status{ind};status];
        res.final_lambda{ind} = [res.final_lambda{ind};results.lambda];
    end
end

for i=1:length(res.rate_per_user)
    res.rate_per_user{i} = res.rate_per_user{i}./par.frequency_slots;
    res.SINR_per_user{i} = res.SINR_per_user{i}./par.frequency_slots;
end

FileName=['Intermediate_file.mat'];
save (FileName);

FileName=['Metrics.mat'];
save (FileName,"res");

%% Generate figures
% Define plotting styles
colors = define_colors();
for_publication = true;
[new_order,legends,style_color,style_line,style_marker] = plot_specifications(par,for_publication,colors);

switch(par.system_type)
    case 'mmwave'
        system_type = 'mmWave';
    case 'cell_free'
        system_type = 'cell_free';
end
[limits] = graph_limits(par);

p=1;
% BER CDF
h = figure(p);
for i=1:length(new_order)
    ind = find(strcmp(par.scheduling_options,new_order{i}));
    temp = res.BER_per_user{ind};
    semilogx(limits.ber.xgrid,cdfgen(limits.ber.xgrid,temp),'Color',style_color{i},'LineStyle',style_line{i},'Marker',style_marker{i},'MarkerFaceColor',style_color{i},'LineWidth',2,'MarkerSize',5,'DisplayName',legends{i});
    hold on;
end
legend('Location','northeast')
legend('show', 'FontSize',12);
xlabel('x');
ylabel('Pr[BER <= x]');
grid;
axis(limits.ber.axis);
plot_title = ([system_type,', B=',num2str(par.B),', U=',num2str(par.U),', Us=',num2str(par.Us),', Ch=',num2str(par.channels)]);
title(plot_title);
file_title = ([system_type,'_B=',num2str(par.B),'_U=',num2str(par.U),'_Us=',num2str(par.Us),'_Ch=',num2str(par.channels)]);
savefig(h,[my_path '/figures/' 'BER_CDF_' file_title])
saveas(h,[my_path '/figures/' 'BER_CDF_' file_title,'.pdf'])
clear('h');

% MI CDF
p=p+1;
h = figure(p);
for i=1:length(new_order)
    ind = find(strcmp(par.scheduling_options,new_order{i}));
    temp = res.MI_per_user{ind};
    plot(limits.mi.xgrid ,cdfgen(limits.mi.xgrid ,temp),'Color',style_color{i},'LineStyle',style_line{i},'Marker',style_marker{i},'MarkerFaceColor',style_color{i},'LineWidth',2,'MarkerSize',5,'DisplayName',legends{i});
    hold on;
end
legend('Location','northeast')
legend('show', 'FontSize',12);
xlabel('x [bits/symbol]');
ylabel('Pr[MI <= x]');
grid;
axis(limits.mi.axis);
plot_title = ([system_type,', B=',num2str(par.B),', U=',num2str(par.U),', Us=',num2str(par.Us),', Ch=',num2str(par.channels)]);
title(plot_title);
file_title = ([system_type,'_B=',num2str(par.B),'_U=',num2str(par.U),'_Us=',num2str(par.Us),'_Ch=',num2str(par.channels)]);
savefig(h,[my_path '/figures/' 'MI_CDF_' file_title])
saveas(h,[my_path '/figures/' 'MI_CDF_' file_title,'.pdf'])
clear('h');

% MI CCDF
p=p+1;
h = figure(p);
for i=1:length(new_order)
    ind = find(strcmp(par.scheduling_options,new_order{i}));
    temp = res.MI_per_user{ind};
    plot(limits.mi.xgrid ,1-cdfgen(limits.mi.xgrid ,temp),'Color',style_color{i},'LineStyle',style_line{i},'Marker',style_marker{i},'MarkerFaceColor',style_color{i},'LineWidth',2,'MarkerSize',5,'DisplayName',legends{i});
    hold on;
end
legend('Location','northeast')
legend('show', 'FontSize',12);
xlabel('x [bits/symbol]');
ylabel('Pr[MI > x]');
grid;
axis(limits.mi.axis);
plot_title = ([system_type,', B=',num2str(par.B),', U=',num2str(par.U),', Us=',num2str(par.Us),', Ch=',num2str(par.channels)]);
title(plot_title);
file_title = ([system_type,'_B=',num2str(par.B),'_U=',num2str(par.U),'_Us=',num2str(par.Us),'_Ch=',num2str(par.channels)]);
savefig(h,[my_path '/figures/' 'MI_CCDF_' file_title])
saveas(h,[my_path '/figures/' 'MI_CCDF_' file_title,'.pdf'])
clear('h');

% MSE CDF
p=p+1;
h = figure(p);
for i=1:length(new_order)
    ind = find(strcmp(par.scheduling_options,new_order{i}));
    temp = res.MSE_per_user{ind};
    plot(limits.mse.xgrid,cdfgen(limits.mse.xgrid,10*log10(temp)),'Color',style_color{i},'LineStyle',style_line{i},'Marker',style_marker{i},'MarkerFaceColor',style_color{i},'LineWidth',2,'MarkerSize',5,'DisplayName',legends{i});
    hold on;
end
legend('Location','northeast')
legend('show', 'FontSize',12);
xlabel('x [dB]');
ylabel('Pr[MSE <= x]');
grid;
axis(limits.mse.axis);
plot_title = ([system_type,', B=',num2str(par.B),', U=',num2str(par.U),', Us=',num2str(par.Us),', Ch=',num2str(par.channels)]);
title(plot_title);
file_title = ([system_type,'_B=',num2str(par.B),'_U=',num2str(par.U),'_Us=',num2str(par.Us),'_Ch=',num2str(par.channels)]);
savefig(h,[my_path '/figures/' 'MSE_CDF_' file_title])
saveas(h,[my_path '/figures/' 'MSE_CDF_' file_title,'.pdf'])
clear('h');

% Rate CDF
p=p+1;
h = figure(p);
for i=1:length(new_order)
    ind = find(strcmp(par.scheduling_options,new_order{i}));
    temp = res.rate_per_user{ind};
    plot(limits.rate.xgrid,cdfgen(limits.rate.xgrid,temp),'Color',style_color{i},'LineStyle',style_line{i},'Marker',style_marker{i},'MarkerFaceColor',style_color{i},'LineWidth',2,'MarkerSize',5,'DisplayName',legends{i});
    hold on;
end
legend('Location','northeast')
legend('show', 'FontSize',12);
xlabel('x [bits/channel use]');
ylabel('Pr[Rate_u <= x]');
grid;
axis(limits.rate.axis);
plot_title = ([system_type,', B=',num2str(par.B),', U=',num2str(par.U),', Us=',num2str(par.Us),', Ch=',num2str(par.channels)]);
title(plot_title);
file_title = ([system_type,'_B=',num2str(par.B),'_U=',num2str(par.U),'_Us=',num2str(par.Us),'_Ch=',num2str(par.channels)]);
savefig(h,[my_path '/figures/' 'Rate_CDF_' file_title])
saveas(h,[my_path '/figures/' 'Rate_CDF_' file_title,'.pdf'])
clear('h');

% Rate CCDF
p=p+1;
h = figure(p);
for i=1:length(new_order)
    ind = find(strcmp(par.scheduling_options,new_order{i}));
    temp = res.rate_per_user{ind};
    plot(limits.rate.xgrid,1-cdfgen(limits.rate.xgrid,temp),'Color',style_color{i},'LineStyle',style_line{i},'Marker',style_marker{i},'MarkerFaceColor',style_color{i},'LineWidth',2,'MarkerSize',5,'DisplayName',legends{i});
    hold on;
end
legend('Location','northeast')
legend('show', 'FontSize',12);
xlabel('x [bits/channel use]');
ylabel('Pr[Rate_u > x]');
grid;
axis(limits.rate.axis);
plot_title = ([system_type,', B=',num2str(par.B),', U=',num2str(par.U),', Us=',num2str(par.Us),', Ch=',num2str(par.channels)]);
title(plot_title);
file_title = ([system_type,'_B=',num2str(par.B),'_U=',num2str(par.U),'_Us=',num2str(par.Us),'_Ch=',num2str(par.channels)]);
savefig(h,[my_path '/figures/' 'Rate_CCDF_' file_title])
saveas(h,[my_path '/figures/' 'Rate_CCDF_' file_title,'.pdf'])
clear('h');

FileName=[datestr(now, 'yyyy_mmm_dd'),'_USER_SCHEDULING_',num2str(par.B),'_',num2str(par.U),'_',num2str(par.frequency_slots),'.mat'];
save (FileName);

get_90_percent_point(my_path,par,file_title);

