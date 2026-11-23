% =========================================================================
% -- Scheduling Matrix Quantization
% -------------------------------------------------------------------------
% Description:
%   Quantize a continuous scheduling matrix while enforcing the configured
%   lower and upper resource constraints.
%
% Input:
%   par: Simulation parameter structure.
%   var: Simulation state structure containing the channel, noise, and
%   scheduler variables.
%   C: Scheduling matrix.
%
% Output:
%   status: Constraint-satisfaction status after quantization.
%   C_final_quant: Quantized scheduling matrix.
%
% -- (c) 2026 Victoria Palhares, Christoph Studer
% =========================================================================

function [status,C_final_quant] = calculate_quantization(par,var,C)

% Elementwise threshold quantization
C_quant = C;

[status,~] = satisfy_constraints(par,C_quant);

% Greedily enforce the scheduling constraints
if(status ~= 1)
    f_vec = [];

    %% Enforce lower-bound constraints
    C_quant = zeros(par.U,par.frequency_slots);
    % Greedily satisfy the lower-bound constraints
    U_set = 1:par.U;
    for t = 1:par.frequency_slots
        % Select the largest U_min entries among the available UEs.
        [~,idx_max] = maxk(C(U_set,t),par.min_sum_columns);
        C_quant(U_set(idx_max),t) = 1; % Activate the selected UEs
        % Remove UEs that already satisfy the minimum slot count
        % from the candidate set.
        idx_sum = intersect(U_set.',find(sum(C_quant,2) >= par.min_sum_rows));
        for ii=1:length(idx_sum)
            U_set(find(U_set == idx_sum(ii))) = [];
        end
    end

    % Evaluate the current objective value
    var.calc_grad = 'no';
    [~,f_no_reg,~,~,~,~] = calculate_objective_and_gradient(par,var,C_quant);
    f = real(f_no_reg);
    f_vec = [f_vec;f];

    %% Evaluate additional time-slot assignments
    extra_users_per_timeslot = par.max_sum_columns-par.min_sum_columns;

    U_set = 1:par.U;
    for t = 1:par.frequency_slots
        % UEs not yet assigned to time slot t
        users_not_chosen = find(C_quant(:,t) == 0);
        idx_intersection = intersect(U_set.',users_not_chosen);
        [~,idx_max] = maxk(C(idx_intersection,t),extra_users_per_timeslot);
        for ii=1:length(idx_max) % Evaluate candidate UEs for time slot t
            C_temp = C_quant; % Candidate scheduling matrix
            % Add the candidate UE to time slot t
            C_temp(idx_intersection(idx_max(ii)),t) = 1;
            var.calc_grad = 'no'; % Evaluate the objective value only
            [~,f_no_reg,~,~,~,~] = calculate_objective_and_gradient(par,var,C_temp);
            if(real(f_no_reg)<f) % Accept an improving objective value
                f = real(f_no_reg);
                f_vec = [f_vec;f]; % Store the accepted objective value
                C_quant = C_temp; % Accept the UE assignment
                vec_sum_rows = sum(C_quant,2);
                % UE has reached the maximum number of assigned time slots
                if(vec_sum_rows(idx_intersection(idx_max(ii))) >= par.max_sum_rows)
                    % Remove the UE from the candidate set.
                    U_set(find(U_set == idx_intersection(idx_max(ii)))) = [];
                end
            end
        end
    end

end

C_final_quant = C_quant;

end
