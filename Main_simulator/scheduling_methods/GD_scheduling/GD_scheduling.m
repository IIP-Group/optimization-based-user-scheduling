% =========================================================================
% -- Optimization-Based Scheduling
% -------------------------------------------------------------------------
% Description:
%   Run the optimization-based scheduling procedure over multiple random
%   initializations and select the best quantized schedule.
%
% Input:
%   par: Simulation parameter structure.
%   var: Simulation state structure containing the channel, noise, and
%   scheduler variables.
%
% Output:
%   par: Updated simulation parameter structure.
%   var: Updated simulation state structure.
%   chosen_C_final: Selected continuous scheduling matrix.
%   chosen_C_final_quantized: Selected quantized scheduling matrix.
%   best_result: Objective value of the selected schedule.
%   best_initialization: Index of the selected random initialization.
%   status_C: Constraint-satisfaction status of the selected schedule.
%   f_vec: Objective-function history for the selected schedule.
%   results: Optimization result structure.
%
% -- (c) 2026 Victoria Palhares, Christoph Studer
% =========================================================================

function [par,var,chosen_C_final,chosen_C_final_quantized,best_result,best_initialization,status_C,f_vec,results] = GD_scheduling(par,var)

var = param_optimization(var,par);


% Initialize arrays for all random initializations
% Continuous scheduling matrices
C = zeros(par.U,par.timeslots,var.grad_iterations,var.num_initializations);
% Final continuous schedules
C_final = zeros(par.U,par.timeslots,var.num_initializations);
% Final quantized schedules
C_final_quantized = zeros(par.U,par.timeslots,var.num_initializations);
% Regularized objective values
f_with_regularizer = zeros(var.grad_iterations,var.num_initializations);
% Unregularized objective values
f_without_regularizer = zeros(var.grad_iterations,var.num_initializations);
% Final regularized objective values
f_with_regularizer_quant_final = zeros(var.num_initializations,1);
% Final unregularized objective values
f_without_regularizer_quant_final = zeros(var.num_initializations,1);
gradient_f_with_regularizer = zeros(par.U,par.timeslots,var.grad_iterations,var.num_initializations);


status_C = zeros(var.num_initializations,1);
tau = zeros(var.grad_iterations,var.num_initializations);

while(sum(status_C) == 0)
    for l=1:var.num_initializations % Iterate over random initializations
        if(strcmp(var.adaptive_stepsize,'on'))
            tau(1,l) = var.tau;
        else
            tau(:,l) = var.tau*ones(var.grad_iterations,1);
        end
        % Random schedule initialization
        C(:,:,1,l) = rand(par.U,par.timeslots);
        % Project onto the scheduling constraints
        [C(:,:,1,l)] = ProjMatrix(par,C(:,:,1,l));

        % Evaluate the objective at the current schedule
        var.calc_grad = 'yes';
        [f_with_regularizer(1,l),f_without_regularizer(1,l),~,...
            gradient_f_with_regularizer(:,:,1,l),~,~] = calculate_objective_and_gradient(par,var,C(:,:,1,l));

        for i=1:var.grad_iterations
            %% Compute the objective gradient

            % Gradient-descent step
            Z = C(:,:,i,l) - tau(i,l).*gradient_f_with_regularizer(:,:,i,l);

            % Project onto the scheduling constraints
            [C(:,:,i+1,l)] = ProjMatrix(par,Z);

            % Evaluate the objective at the current schedule
            var.calc_grad = 'yes';
            [f_with_regularizer(i+1,l),f_without_regularizer(i+1,l),~,...
                gradient_f_with_regularizer(:,:,i+1,l),~,~] = calculate_objective_and_gradient(par,var,C(:,:,i+1,l));

            % Adaptive step size
            if(strcmp(var.adaptive_stepsize,'on'))
                Dx = C(:,:,i+1,l)-C(:,:,i,l);
                %% Compute the BB/spectral step size for the next iteration
                Dg = gradient_f_with_regularizer(:,:,i+1,l)-gradient_f_with_regularizer(:,:,i,l);
                dotprod = real(dot(Dx(:),Dg(:)));
                tau_s = norm(Dx(:))^2/ dotprod;  % First BB step-size rule
                % Alternate BB step-size rule
                tau_m = dotprod / norm(Dg(:))^2;
                tau_m = max(tau_m,0);
                if 2*tau_m > tau_s   % Select the adaptive BB step size
                    tau(i+1,l) = tau_m;
                else
                    tau(i+1,l) = tau_s - .5*tau_m;  % Combined BB step size
                end
                % Enforce a finite positive step size
                if tau(i+1,l) <=0 || isinf(tau(i+1,l)) || isnan(tau(i+1,l))
                    % Increase tau; backtracking limits overly large steps
                    tau(i+1,l) = tau(i,l)*1.5;
                end
            end
        end
        % Final schedule for this initialization
        C_final(:,:,l) = C(:,:,i+1,l);
        % Quantize the final scheduling matrix.
        [status_C,C_final_quantized(:,:,l)] = calculate_quantization(par,var,C_final(:,:,l));
        % Evaluate the objective for the final quantized schedules
        var.calc_grad = 'no';
        [f_with_regularizer_quant_final(l),f_without_regularizer_quant_final(l),~,~,~,~] = calculate_objective_and_gradient(par,var,C_final_quantized(:,:,l));
    end

    switch(var.scheduling_method)
        case{'MMSE'}
            best_result = inf;
        case{'Capacity','Rate'}
            best_result = -inf;
        otherwise
            error('Unsupported optimization objective. Choose MMSE, Capacity, or Rate.');
    end
    for l=1:var.num_initializations
        % Check scheduling constraints
        [status_C(l),~] = satisfy_constraints(par,C_final_quantized(:,:,l));
        if(status_C(l) == 1)
            switch(var.scheduling_method)
                case {'MMSE'}
                    if(f_without_regularizer_quant_final(l) < best_result)
                        best_result = f_without_regularizer_quant_final(l);
                        best_initialization = l;
                    end
                case {'Capacity','Rate'}
                    if(f_without_regularizer_quant_final(l) > best_result)
                        best_result = f_without_regularizer_quant_final(l);
                        best_initialization = l;
                    end
                otherwise
                    error('Unsupported optimization objective. Choose MMSE, Capacity, or Rate.');
            end
        end
    end
    if(sum(status_C) == 0)
        var.lambda = 2*var.lambda;
        warning('Scheduling constraints were not satisfied.');
    else
        break;
    end
end

% Select the best feasible random initialization
chosen_C_final = C_final(:,:,best_initialization);
chosen_C_final_quantized = C_final_quantized(:,:,best_initialization);

results.f_with_regularizer_quant_final = f_with_regularizer_quant_final;
results.f_without_regularizer_quant_final = f_without_regularizer_quant_final;

% Objective history for the selected initialization
f_vec = real(f_without_regularizer(:,best_initialization));
end
