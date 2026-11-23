% =========================================================================
% -- Binarization Regularizer
% -------------------------------------------------------------------------
% Description:
%   Evaluate the regularizer that promotes binary scheduling values and its
%   gradient.
%
% Input:
%   par: Simulation parameter structure.
%   var: Simulation state structure containing the channel, noise, and
%   scheduler variables.
%   C: Scheduling matrix.
%
% Output:
%   regularization: Regularization value.
%   grad_regularization: Gradient of the regularizer.
%
% -- (c) 2026 Victoria Palhares, Christoph Studer
% =========================================================================

function [regularization,grad_regularization] = L2_norm(par,var,C)

regularization = 0;
grad_regularization = zeros(par.U,par.frequency_slots);
for t = 1:par.frequency_slots
    regularization = regularization+(var.lambda*norm(C(:,t)-(0.5.*ones(par.U,1)),'fro')^2);

    switch(var.calc_grad)
        case 'yes'
            grad_regularization(:,t) = (var.lambda*(2*C(:,t)-ones(par.U,1)));
    end
end

end
