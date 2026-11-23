% =========================================================================
% -- LoFi Scheduling (K>1)
% -------------------------------------------------------------------------
% Description:
%   Generate a LoFi schedule over multiple random restarts and select the
%   candidate with the best minimum SINR.
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

function [C_chosen] = lofi_K_N(par,var,num_restarts)

% Generate candidate schedules over all time slots
C_matrix = zeros(par.U,par.frequency_slots,num_restarts);
SINR_u = zeros(par.U,num_restarts);
min_candidate = zeros(num_restarts,1);

for s=1:num_restarts

    all_users = 1:par.U;
    remaining_users_ = all_users;

    for t = 1:par.frequency_slots
        chosen_ues = randsample(remaining_users_, par.Us);
        C_matrix(chosen_ues, t, s) = 1;
        remaining_users_(ismember(remaining_users_, chosen_ues)) = [];
    end

    for t=1:par.frequency_slots
        Hc = var.H(:,:,t)*diag(C_matrix(:,t,s));
        L = inv(Hc*Hc'+(var.N0_scheduling/var.Es)*eye(par.B));

        L_H = zeros(par.B,par.U);
        H_L = zeros(par.U,par.B);
        for n = 1:par.U
            L_H(:,n) = L*var.H(:,n,t);
            H_L(n,:) = var.H(:,n,t)'*L;
        end

        idx = find(C_matrix(:,t,s) == 1);
        for u=1:length(idx)
            % Equalizer vector for UE u
            w_c_u = C_matrix(idx(u),t,s)*(H_L(idx(u),:));
            % Channel vector for UE u
            h_c_u = C_matrix(idx(u),t,s)*(var.H(:,idx(u),t));

            d_c = (var.Es*abs(w_c_u*h_c_u)^2); % Desired-signal power
            % Interference power
            i_c = (var.Es*((abs(w_c_u*Hc).^2*ones(par.U,1))- (abs(w_c_u*h_c_u)^2)));
            n_c = (var.N0_scheduling*norm(w_c_u)^2); % Noise power

            SINR_u(idx(u),s) = d_c/(i_c+ n_c);
        end
    end

    min_candidate(s) = min(SINR_u(:,s));
end

[~,idx_chosen] = max(min_candidate);

C_chosen = C_matrix(:,:,idx_chosen);

end

