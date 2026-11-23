% =========================================================================
% -- Post-LMMSE MSE Objective
% -------------------------------------------------------------------------
% Description:
%   Evaluate the post-LMMSE mean-squared-error scheduling objective and its
%   gradient.
%
% Input:
%   par: Simulation parameter structure.
%   var: Simulation state structure containing the channel, noise, and
%   scheduler variables.
%   C: Scheduling matrix.
%
% Output:
%   f_without_regularizer: Objective value before regularization.
%   grad_f_without_regularizer: Gradient of the unregularized objective.
%
% -- (c) 2026 Victoria Palhares, Christoph Studer
% =========================================================================

function [f_without_regularizer,grad_f_without_regularizer] = MMSE(par,var,C)

f_without_regularizer = 0;
grad_f_without_regularizer = zeros(par.U,par.timeslots);
for t = 1:par.timeslots
    Hc = var.H*diag(C(:,t));
    L = eye(par.B)/(Hc*Hc'+(var.N0_scheduling/var.Es)*eye(par.B));
    M = Hc*diag(C(:,t).^2)*Hc';

    f_without_regularizer = f_without_regularizer+var.Es*trace(diag(C(:,t).^2))-var.Es*trace(L*M);

    switch(var.calc_grad)
        case 'yes'
            for u=1:par.U
                L_H = L*var.H(:,u);
                H_L = var.H(:,u)'*L;
                H_L_M = H_L*M;
                % Gradient of the unregularized objective
                grad_f_without_regularizer(u,t) = var.Es*(2*C(u,t)+2*C(u,t)*H_L_M*L_H-4*C(u,t)^3*var.H(:,u)'*L_H);
            end
    end

end

end
