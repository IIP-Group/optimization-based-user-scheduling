% =========================================================================
% -- Exhaustive-Search Scheduling
% -------------------------------------------------------------------------
% Description:
%   Enumerate feasible scheduling matrices and select the schedule that
%   optimizes the configured exhaustive-search metric.
%
% Input:
%   par: Simulation parameter structure.
%   var: Simulation state structure containing the channel, noise, and
%   scheduler variables.
%
% Output:
%   best_C: Best scheduling matrix found by exhaustive search.
%   best_measure: Objective value of the best exhaustive-search schedule.
%
% -- (c) 2026 Victoria Palhares, Christoph Studer
% =========================================================================

function [best_C,best_measure] = exhaustive_search(par,var)

combos = nchoosek(1:par.U,par.Us);
vec_k = 1:par.U;
[var.num_poss,~] = size(combos);
C_matrix = zeros(par.U,par.timeslots,var.num_poss);

for i=1:par.U
    for j=1:var.num_poss
        C_matrix(i,1,j) = ismember(vec_k(i),combos(j,:));
        C_matrix(i,2,j) = (C_matrix(i,1,j) == 0);
    end
end

measure = zeros(par.timeslots,var.num_poss);
sum_measure = zeros(var.num_poss,1);

for n=1:var.num_poss
    for t=1:par.timeslots
        C = C_matrix(:,t,n);
        Hc = var.H*diag(C);
        switch(var.scheduling_method)
            case 'ES-MSE'
                % Evaluate the MSE objective
                L = eye(par.B)/(Hc*Hc'+(var.N0_scheduling/var.Es)*eye(par.B));
                M = Hc*diag(C.^2)*Hc';
                measure(t,n) = var.Es*trace(diag(C.^2))-var.Es*trace(L*M);
            case 'ES-Capacity'
                % Evaluate the capacity objective
                measure(t,n) = log2(det(eye(par.B)+(var.Es/var.N0_scheduling)*(Hc*Hc')));
            case 'ES-Rate'
                L = eye(par.B)/(Hc*Hc'+(var.N0_scheduling/var.Es)*eye(par.B));
                rate = zeros(par.U,1);
                SINR = zeros(par.U,1);
                for u=1:par.U
                    % Equalizer vector for UE u
                    w_c_u = C(u)*(var.H(:,u)'*L);
                    h_c_u = C(u)*(var.H(:,u)); % Channel vector for UE u

                    % Desired-signal power
                    d_c = (var.Es*abs(w_c_u*h_c_u)^2);
                    % Interference power
                    i_c = (var.Es*((abs(w_c_u*Hc).^2*ones(par.U,1))- (abs(w_c_u*h_c_u)^2)));
                    n_c = (var.N0_scheduling*norm(w_c_u)^2); % Noise power
                    % Channel-estimation error power
                    e_c_1 = var.Es*(norm(w_c_u)^2);
                    e_c_2 = ones(1,par.U)*(diag(C.^2)*diag(var.covariance_channel_error))*ones(par.U,1);
                    e_c = e_c_1*e_c_2;

                    SINR(u) = d_c/(i_c+n_c+e_c);

                    if(d_c == 0 && i_c == 0 && n_c == 0 && e_c == 0)
                        SINR(u) = 0;
                    end

                    rate(u) = log2(1+SINR(u));
                end
                measure(t,n) = sum(rate);
            otherwise
                error('Unsupported exhaustive-search objective. Choose ES-MSE, ES-Capacity, or ES-Rate.');
        end
    end
    sum_measure(n) = sum(measure(:,n));
end

switch(var.scheduling_method)
    case 'ES-MSE'
        best_measure = min(sum_measure);
    case {'ES-Capacity','ES-Rate'}
        best_measure = max(sum_measure);
    otherwise
        error('Unsupported exhaustive-search objective. Choose ES-MSE, ES-Capacity, or ES-Rate.');
end
best_poss = find(sum_measure == best_measure);
best_C = C_matrix(:,:,best_poss);

end

