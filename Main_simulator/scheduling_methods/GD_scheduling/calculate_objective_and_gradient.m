% =========================================================================
% -- Objective and Gradient Evaluation
% -------------------------------------------------------------------------
% Description:
%   Evaluate the selected scheduling objective, regularizer, and their
%   gradients for the current scheduling matrix.
%
% Input:
%   par: Simulation parameter structure.
%   var: Simulation state structure containing the channel, noise, and
%   scheduler variables.
%   C: Scheduling matrix.
%
% Output:
%   f_with_regularizer: Objective value including regularization.
%   f_without_regularizer: Objective value before regularization.
%   regularization: Regularization value.
%   gradient_f_with_regularizer: Gradient of the regularized objective.
%   grad_f_without_regularizer: Gradient of the unregularized objective.
%   grad_regularization: Gradient of the regularizer.
%
% -- (c) 2026 Victoria Palhares, Christoph Studer
% =========================================================================

function [f_with_regularizer,f_without_regularizer,regularization,gradient_f_with_regularizer,...
    grad_f_without_regularizer,grad_regularization] = calculate_objective_and_gradient(par,var,C)

factor = 1;

switch(var.regularization)
    case('L2 Norm')
        [regularization,grad_regularization] = L2_norm(par,var,C);
    otherwise
        error('Unsupported regularization setting.');
end
switch(var.scheduling_method)
    case 'MMSE'
        [f_without_regularizer,grad_f_without_regularizer] = MMSE(par,var,C);
    case 'Capacity'
        [f_without_regularizer,grad_f_without_regularizer] = Capacity(par,var,C);
    case 'Rate'
        [f_without_regularizer,grad_f_without_regularizer] = Rate(par,var,C);
    otherwise
        error('Unsupported optimization objective. Choose MMSE, Capacity, or Rate.');
end
switch(var.scheduling_method)
    case {'Capacity','Rate'}
        factor = -1;
end

f_with_regularizer = factor*f_without_regularizer - regularization;
gradient_f_with_regularizer = factor*grad_f_without_regularizer - grad_regularization;

end
