% =========================================================================
% -- LoFi++ Scheduling (K>1)
% -------------------------------------------------------------------------
% Description:
%   Generate a LoFi++ schedule over multiple random restarts and worst-user
%   swap candidates.
%
% Input:
%   par: Simulation parameter structure.
%   var: Simulation state structure containing the channel, noise, and
%   scheduler variables.
%   num_restarts: Number of random scheduling restarts.
%
% Output:
%   C_chosen: Selected scheduling matrix.
%
% -- (c) 2023 Reinhard Wiesmayr, Christoph Studer
% =========================================================================

function [C_chosen] = lofi_2_K_N(par,var,num_restarts)

snr_min_ar = zeros(num_restarts*2,1);

C_matrix = zeros(par.U,par.frequency_slots,2*num_restarts);

for i=1:num_restarts
    % Generate one candidate per random restart
    chosen_ues_1 = randsample(par.U, par.Us).';
    complement_1 = 1:par.U;
    complement_1(chosen_ues_1) = [];

    C_matrix(chosen_ues_1,1,2*i-1) = 1;
    C_matrix(complement_1,2,2*i-1) = 1;

    t=1;
    Hc = var.H*diag(C_matrix(:,t,2*i-1));
    L = inv(Hc*Hc'+(var.N0_scheduling/var.Es)*eye(par.B));

    L_H = zeros(par.B,par.U);
    H_L = zeros(par.U,par.B);
    for n = 1:par.U
        L_H(:,n) = L*var.H(:,n);
        H_L(n,:) = var.H(:,n)'*L;
    end

    SINR_u = zeros(par.U,2);
    for u=chosen_ues_1
        % Equalizer vector for UE u
        w_c_u = C_matrix(u,t,2*i-1)*(H_L(u,:));
        % Channel vector for UE u
        h_c_u = C_matrix(u,t,2*i-1)*(var.H(:,u));

        d_c = (var.Es*abs(w_c_u*h_c_u)^2); % Desired-signal power
        % Interference power
        i_c = (var.Es*((abs(w_c_u*Hc).^2*ones(par.U,1))- ...
            (abs(w_c_u*h_c_u)^2)));
        n_c = (var.N0_scheduling*norm(w_c_u)^2); % Noise power

        SINR_u(u,1) = d_c/(i_c+ n_c);
    end

    t=2;
    Hc = var.H*diag(C_matrix(:,t,2*i-1));
    L = inv(Hc*Hc'+(var.N0_scheduling/var.Es)*eye(par.B));

    L_H = zeros(par.B,par.U);
    H_L = zeros(par.U,par.B);
    for n = 1:par.U
        L_H(:,n) = L*var.H(:,n);
        H_L(n,:) = var.H(:,n)'*L;
    end

    for u=complement_1
        % Equalizer vector for UE u
        w_c_u = C_matrix(u,t,2*i-1)*(H_L(u,:));
        % Channel vector for UE u
        h_c_u = C_matrix(u,t,2*i-1)*(var.H(:,u));

        d_c = (var.Es*abs(w_c_u*h_c_u)^2); % Desired-signal power
        i_c = (var.Es*((abs(w_c_u*Hc).^2*ones(par.U,1))- ...
            (abs(w_c_u*h_c_u)^2))); % Interference power
        n_c = (var.N0_scheduling*norm(w_c_u)^2); % Noise power

        SINR_u(u,1) = d_c/(i_c+ n_c);
    end

    % Identify the weakest UE in each time slot
    [min_candidate_1_1, idx_1] = min(SINR_u(chosen_ues_1,1)); % Time slot 1
    [min_candidate_1_2, idx_2] = min(SINR_u(complement_1,1)); % Time slot 2
    idx_1 = chosen_ues_1(idx_1); % UE index
    idx_2 = complement_1(idx_2);

    % Swap the weakest UEs between the two time slots
    chosen_ues_2 = chosen_ues_1;
    chosen_ues_2(chosen_ues_1==idx_1) = idx_2;
    complement_2 = complement_1;
    complement_2(complement_1==idx_2) = idx_1;

    % Evaluate the swapped schedule
    C_matrix(chosen_ues_2,1,2*i) = 1;
    C_matrix(complement_2,2,2*i) = 1;

    t=1;
    Hc = var.H*diag(C_matrix(:,t,2*i));
    L = inv(Hc*Hc'+(var.N0_scheduling/var.Es)*eye(par.B));

    L_H = zeros(par.B,par.U);
    H_L = zeros(par.U,par.B);
    for n = 1:par.U
        L_H(:,n) = L*var.H(:,n);
        H_L(n,:) = var.H(:,n)'*L;
    end

    for u=chosen_ues_2
        % Equalizer vector for UE u
        w_c_u = C_matrix(u,t,2*i)*(H_L(u,:));
        % Channel vector for UE u
        h_c_u = C_matrix(u,t,2*i)*(var.H(:,u));

        d_c = (var.Es*abs(w_c_u*h_c_u)^2); % Desired-signal power
        i_c = (var.Es*((abs(w_c_u*Hc).^2*ones(par.U,1))- ...
            (abs(w_c_u*h_c_u)^2))); % Interference power
        n_c = (var.N0_scheduling*norm(w_c_u)^2); % Noise power

        SINR_u(u,2) = d_c/(i_c+ n_c);
    end

    t=2;
    Hc = var.H*diag(C_matrix(:,t,2*i));
    L = inv(Hc*Hc'+(var.N0_scheduling/var.Es)*eye(par.B));

    L_H = zeros(par.B,par.U);
    H_L = zeros(par.U,par.B);
    for n = 1:par.U
        L_H(:,n) = L*var.H(:,n);
        H_L(n,:) = var.H(:,n)'*L;
    end

    for u=complement_2
        % Equalizer vector for UE u
        w_c_u = C_matrix(u,t,2*i)*(H_L(u,:));
        % Channel vector for UE u
        h_c_u = C_matrix(u,t,2*i)*(var.H(:,u));

        d_c = (var.Es*abs(w_c_u*h_c_u)^2); % Desired-signal power
        i_c = (var.Es*((abs(w_c_u*Hc).^2*ones(par.U,1))- ...
            (abs(w_c_u*h_c_u)^2))); % Interference power
        n_c = (var.N0_scheduling*norm(w_c_u)^2); % Noise power

        SINR_u(u,2) = d_c/(i_c+ n_c);
    end

    % Minimum SINR for each candidate schedule
    SINR_mins = min(SINR_u, [], 1);
    snr_min_ar(2*i-1:2*i) = SINR_mins(:);
end

[~, idx_max] = max(snr_min_ar);

C_chosen = C_matrix(:,:,idx_max); % Selected candidate schedule

end

