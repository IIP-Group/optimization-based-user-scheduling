% =========================================================================
% -- Single Simulation Job
% -------------------------------------------------------------------------
% Description:
%   Run one simulation job for a selected channel realization, SNR point,
%   and scheduling method and save the resulting performance metrics.
%
% Input:
%   par: Simulation parameter structure.
%   proc_ID: Zero-based simulation process identifier.
%
% Output:
%   None. Simulation metrics are saved to a MAT file.
%
% -- (c) 2026 Victoria Palhares, Christoph Studer
% =========================================================================

function main_function(par,proc_ID)

% Decode the channel ID, SNR, and scheduling method for this job
var.id = proc_ID+1;
parameters = par.parameters(var.id,:);
var.ch = str2num(parameters(1)); % Channel ID
var.SNR_dB = str2num(parameters(2)); % SNR in dB
var.scheduling_method = parameters(3); % Selected scheduling method

snr_idx = find(par.SNRdB_list == var.SNR_dB);
scheduler_idx = find(par.scheduling_options == var.scheduling_method);

index = (var.ch-1)*length(par.SNRdB_list)+snr_idx;
rng(1000*index + 0, 'twister');

% Configure the modulation scheme
% Define the Gray-mapped constellation alphabet according to IEEE 802.11
switch (par.modulation)
    case 'BPSK'
        par.symbols = [ -1 1 ];
    case 'QPSK'
        par.symbols = [ -1-1i,-1+1i, ...
            +1-1i,+1+1i ];
    case '16QAM'
        par.symbols = [ -3-3i,-3-1i,-3+3i,-3+1i, ...
            -1-3i,-1-1i,-1+3i,-1+1i, ...
            +3-3i,+3-1i,+3+3i,+3+1i, ...
            +1-3i,+1-1i,+1+3i,+1+1i ];
    case '64QAM'
        par.symbols = [ -7-7i,-7-5i,-7-1i,-7-3i,-7+7i,-7+5i,-7+1i,-7+3i, ...
            -5-7i,-5-5i,-5-1i,-5-3i,-5+7i,-5+5i,-5+1i,-5+3i, ...
            -1-7i,-1-5i,-1-1i,-1-3i,-1+7i,-1+5i,-1+1i,-1+3i, ...
            -3-7i,-3-5i,-3-1i,-3-3i,-3+7i,-3+5i,-3+1i,-3+3i, ...
            +7-7i,+7-5i,+7-1i,+7-3i,+7+7i,+7+5i,+7+1i,+7+3i, ...
            +5-7i,+5-5i,+5-1i,+5-3i,+5+7i,+5+5i,+5+1i,+5+3i, ...
            +1-7i,+1-5i,+1-1i,+1-3i,+1+7i,+1+5i,+1+1i,+1+3i, ...
            +3-7i,+3-5i,+3-1i,+3-3i,+3+7i,+3+5i,+3+1i,+3+3i ];
    otherwise
        error('Unsupported modulation scheme.')
end
% Normalize the constellation and compute the average symbol energy
% Normalize to unit average energy
par.symbols = par.symbols/sqrt(mean(abs(par.symbols).^2));
var.Es = mean(abs(par.symbols).^2);

% Precompute bit labels
var.Q = log2(length(par.symbols)); % Number of bits per symbol
var.bits = de2bi(0:length(par.symbols)-1,var.Q,'left-msb');

% Load the channel realization
[var] = gen_channel(par,var);
var.H_full_CSI = var.H;

% Configure the training-noise variance
switch (par.plot_type)
    case 'SNR'
        var.H_norm = norm(var.H,'fro')^2;
        var.SNR = 10^(var.SNR_dB/10);
        var.N0_training = ((var.H_norm*var.Es)/(var.SNR*par.B));
        % Compute the resulting SNR
        SNR = (var.H_norm*var.Es)/(par.B*var.N0_training);
        SNR_dB = 10*log10(SNR);
    case 'CDF'
        var.N0_training = par.N0;
end
% Estimate the channel
switch(par.channel_estimation)
    case 'BEACHES'
        N0_est = var.N0_training/(par.U*var.Es);
        channel_noise = (randn(par.B,par.U)+1i*randn(par.B,par.U));
        Hnoisy = var.H + channel_noise*sqrt(N0_est/2);
        % Estimate the channel using BEACHES
        [Hest,~,var.MSE] = BEACHES(par,Hnoisy,N0_est,'antenna_domain');
        % Use the estimated channel for scheduling and detection
        var.H = Hest;
        var.H_error = var.H-var.H_full_CSI; % Channel-estimation error
        var.covariance_channel_error = var.MSE;
    case 'LS' % Least-squares channel estimation
        par.pilot_length = par.U;
        S = hadamard(par.pilot_length);
        S = S/sqrt(mean(abs(S(:)).^2)); % Normalize to unit average energy
        pilots_variance = mean(abs(S(:)).^2);
        N = sqrt(0.5*var.N0_training)*(randn(par.B,par.pilot_length)+1i*randn(par.B,par.pilot_length));
        Y = var.H_full_CSI*S  + N;
        Hest = Y*S'/(S*S');
        var.H = Hest;
        var.H_error = var.H-var.H_full_CSI; % Channel-estimation error
        % Channel-estimation error covariance
        var.covariance_channel_error = (var.N0_training/par.U)*ones(par.U,1);
end

var.N0_scheduling = var.N0_training;

% Reinitialize the random-number generator
rng(1000*index + 0, 'twister');
% Select the scheduling algorithm
switch (var.scheduling_method)
    case 'Random'
        % Random scheduling
        [C_rand,best_result] = random_assignment(par,var);
        C = C_rand;
    case {'ES-MSE','ES-Capacity','ES-Rate'}
        % Exhaustive-search scheduling
        switch (par.status_es)
            case 'on'
            [C_es,best_result] = exhaustive_search(par,var);
            C = C_es(:,:,1);
            otherwise
            C = zeros(par.U,par.timeslots);
            best_result = 0;
        end
    case {'MMSE','Capacity','Rate'}
        % Optimization-based scheduling
        var.regularization = 'L2 Norm';
        [par,var,C_not_quant,C_quant,best_result,best_initialization,status,f_vec,results] = ...
            GD_scheduling(par,var);
        C = C_quant;
    case 'SUS'
        % SUS scheduling
        [C_sus] = SUS(par,var);
        C = C_sus;
    case 'CSS'
        % CSS scheduling
        [C_css] = CSS(par,var);
        C = C_css;
    case 'Greedy'
        % Greedy scheduling
        [C_greedy] = greedy(par,var);
        C = C_greedy;
    case 'No Scheduling'
        % AUA: all-user assignment without scheduling
        C = ones(par.U,par.timeslots);
    case "LoFi (K=1)"
        % LoFi (K=1)
        [C_lofi_K_1] = lofi_K_1(par,var);
        C = C_lofi_K_1;
    case "LoFi (K=4)"
        % LoFi (K=4)
        num_restarts=8;
        [C_lofi_K_N] = lofi_K_N(par,var,num_restarts);
        C = C_lofi_K_N;
    case "LoFi++ (K=1)"
        % LoFi++ (K=1)
        [C_lofi_2_K_1] = lofi_2_K_1(par,var);
        C = C_lofi_2_K_1;
    case "LoFi++ (K=4)"
        % LoFi++ (K=4)
        num_restarts = 4;
        [C_lofi_2_K_N] = lofi_2_K_N(par,var,num_restarts);
        C = C_lofi_2_K_N;
    otherwise
        error('Unsupported scheduling method.');
end

% Reinitialize the random-number generator so all scheduling methods
% use identical transmitted bits and noise realizations.
rng(1000*index + 0, 'twister');

% Generate random transmitted bits
% One symbol vector per trial and time slot
var.fixed_bits = randi([0 1],par.U,var.Q,par.timeslots,par.trials);

% Initialize result arrays
SNR_dB = zeros(par.timeslots,1);
ber = zeros(par.U,par.timeslots,par.trials);
mse = zeros(par.U,par.timeslots,par.trials);
SINR = zeros(par.U,par.timeslots);
achievable_rate = zeros(par.U,par.timeslots);
capacity = zeros(par.timeslots,1);

s_schedule = zeros(par.U,par.timeslots,par.trials);
possible_values = [0,1]; % Binary alphabet
count_b_b_hat = zeros(length(possible_values),length(possible_values),var.Q,par.U);

for t=1:par.timeslots
    Dc = diag(C(:,t));
    Hc_est = var.H*Dc; % Scheduled estimated channel
    Hc_real = var.H_full_CSI*Dc; % Scheduled true channel
    Hc_error = var.H_error*Dc;

    % Configure the detection-noise variance
    switch (par.plot_type)
        case 'SNR'
            var.real_H_norm(t) = norm(Hc_real,'fro')^2;
            var.SNR = 10^(var.SNR_dB/10);
            var.real_N0(t) = ((var.real_H_norm(t)*var.Es)/(var.SNR*par.B));
            % Compute the resulting SNR
            SNR = (var.real_H_norm(t)*var.Es)/(par.B*var.real_N0(t));
            SNR_dB(t) = 10*log10(SNR);
        case 'CDF'
            var.real_H_norm(t) = norm(Hc_real,'fro')^2;
            var.real_N0(t) = par.N0;
            % Compute the resulting SNR
            SNR = (var.real_H_norm(t)*var.Es)/(par.B*var.real_N0(t));
            SNR_dB(t) = 10*log10(SNR);
    end

    % LMMSE equalizer for the scheduled channel
    Wc = Hc_est/(Hc_est'*Hc_est+(var.real_N0(t)/var.Es)*eye(par.U));

    for tt=1:par.trials
        % Generate the transmit symbol vector
        idx = bi2de(var.fixed_bits(:,:,t,tt),'left-msb')+1;
        s = par.symbols(idx).';

        % Generate additive noise
        n = sqrt(0.5*var.real_N0(t))*(randn(par.B,1)+1i*randn(par.B,1));

        yc = Hc_real*s+n;
        s_schedule(:,t,tt) = Dc*s;

        switch(par.evaluation)
            case('mmse')
                s_hat_schedule = Dc*Wc'*yc;
                var.vector = real(diag(Dc*Wc'*Hc_real));
            otherwise
                error('Unsupported evaluation metric.')
        end

        % Accumulate bit errors
        [~,idxhat] = min(abs(s_hat_schedule*ones(1,length(par.symbols))-var.vector*par.symbols).^2,[],2);
        bithat = var.bits(idxhat,:);
        bits = C(:,t).*var.fixed_bits(:,:,t,tt);
        ber(:,t,tt) = sum(bits ~= bithat,2);

        % Accumulate symbol MSE
        mse(:,t,tt) = abs(s_schedule(:,t,tt)-s_hat_schedule).^2;

        % Accumulate the joint PMF used for conditional entropy
        for u=1:par.U
            if(C(u,t) == 1) % Include active UEs in the statistics
                for uu=1:var.Q
                    idx_row = find(possible_values == bits(u,uu));
                    idx_column = find(possible_values == bithat(u,uu));
                    % Count each transmitted/detected bit pair
                    count_b_b_hat(idx_row,idx_column,uu,u) = count_b_b_hat(idx_row,idx_column,uu,u) + 1;
                end
            end
        end

    end

    % Compute achievable rates
    for u=1:par.U
        % Desired-signal power
        des_sig = (var.Es*abs(Wc(:,u)'*Hc_est(:,u))^2);
        % Interference power
        interf = (var.Es*((abs(Wc(:,u)'*Hc_est).^2*ones(par.U,1))- (abs(Wc(:,u)'*Hc_est(:,u))^2)));
        % Channel-estimation error power
        ch_error = var.Es*(abs(Wc(:,u)'*Hc_error).^2*ones(par.U,1));
        % Noise power
        noise = (var.real_N0(t)*norm(Wc(:,u))^2);

        SINR(u,t) = des_sig/(interf+ch_error+noise);

        if(des_sig == 0 && interf == 0 && ch_error==0 && noise == 0)
            SINR(u,t) = 0;
        end
        achievable_rate(u,t) = log2(1+SINR(u,t));
    end

    % Compute channel capacity for the scheduled channel
    % Capacity for the scheduled channel
    capacity(t) = log2(det(eye(par.B)+(var.Es/var.real_N0(t))*(Hc_real*Hc_real')));

end

% Store the average noise variance
switch (par.plot_type)
    case 'SNR'
        results.avg_real_N0 = mean(var.real_N0);
    case 'CDF'
        results.avg_real_N0 = par.N0;
end

results.SNR_dB = SNR_dB;

results.C = C;
% BER
results.timeslots_per_user = sum(C,2);
sum_all_C = sum(sum(C));
% Average over scheduled time slots and trials
results.ber_per_user = sum(sum(ber,2),3)./(var.Q*par.trials.*results.timeslots_per_user);
% Average BER over scheduled UE-slots
results.ber_all_users = sum(ber(:))/(var.Q*par.trials*sum_all_C);

% MSE
% Average over scheduled time slots and trials
results.mse_per_user = sum(sum(mse,2),3)./(par.trials.*results.timeslots_per_user);
results.mse_all_users = sum(mse(:))/(par.trials*sum_all_C);

results.rmsse_per_user = zeros(par.U,1);
for u=1:par.U
    results.rmsse_per_user(u) = sqrt(sum(mse(u,:,:),'all')/sum(abs(s_schedule(u,:,:)).^2,'all'));
end

% Achievable rate
% Sum over scheduled time slots
results.rate_per_user = sum(achievable_rate,2);
% Average rate per UE
results.rate_all_users = sum(achievable_rate(:))/par.U;

% SINR
results.sinr_per_user = sum(SINR,2); % Sum over scheduled time slots
results.sinr_all_users = sum(SINR(:))/par.U; % Average SINR per UE

% Capacity
results.capacity = sum(capacity)/par.U; % Average capacity per UE

% Mutual information
% Compute the source entropy
p_b = 1/2; % Equiprobable binary source
entropy = (-2*p_b*log2(p_b));

n_counts = par.trials.*results.timeslots_per_user;
% Compute the conditional entropy
conditional_entropy = zeros(par.U,var.Q);
probability_b_b_hat = zeros(length(possible_values),length(possible_values),var.Q,par.U);
for u=1:par.U
    for uu=1:var.Q
        % Joint PMF for each bit and UE
        probability_b_b_hat(:,:,uu,u) = count_b_b_hat(:,:,uu,u)./n_counts(u);
        % Sum over transmitted bits to compute p(b_hat)
        p_b_hat = sum(probability_b_b_hat(:,:,uu,u),1);
        for i=1:length(possible_values)
            for j=1:length(possible_values)
                if(probability_b_b_hat(i,j,uu,u)~= 0)
                    conditional_entropy(u,uu) = conditional_entropy(u,uu)+probability_b_b_hat(i,j,uu,u)*(log2(probability_b_b_hat(i,j,uu,u)/p_b_hat(j)));
                end
            end
        end
    end
end

% Conditional entropy per UE and bit position
conditional_entropy = -conditional_entropy;

% Mutual information per UE and bit position
mutual_information = entropy - conditional_entropy;

% Sum mutual information across bit positions
results.mi_per_user = sum(mutual_information,2);
results.mi_all_users = sum(results.mi_per_user);

% Store scheduling-specific output variables
results.scheduling_method = var.scheduling_method;
switch (var.scheduling_method)
    case {'MMSE','Capacity','Rate'}
        results.lambda = var.lambda;
end

index_process = ((scheduler_idx-1)*(length(par.SNRdB_list)*par.channels)+index)-1;

path_to_save = (['/path/to/results/']);

FileName=[path_to_save,'USER_SCHEDULING_',par.system_type,'_',num2str(par.B),'_',num2str(par.U),'_T_',num2str(par.timeslots),'_',num2str(index_process),'.mat'];
save (FileName,'results');

switch (var.scheduling_method)
    case {'MMSE','Capacity','Rate'}
        save(FileName,'status','-append')
end

toc;

end

