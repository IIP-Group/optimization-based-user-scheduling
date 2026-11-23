% =========================================================================
% -- Channel Capacity Objective
% -------------------------------------------------------------------------
% Description:
%   Evaluate the channel-capacity scheduling objective and its gradient.
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

function [f_without_regularizer,grad_f_without_regularizer] = Capacity(par,var,C)

f_without_regularizer = 0;
grad_f_without_regularizer = zeros(par.U,par.frequency_slots);
rho = (var.Es/var.N0_scheduling);

for t = 1:par.frequency_slots
    Hc = var.H(:,:,t)*diag(C(:,t));
    L = eye(par.B)+rho*(Hc*Hc'); % L
    inv_L = inv(L);

    f_without_regularizer = f_without_regularizer+log2(det(L));

    switch(var.calc_grad)
        case 'yes'
            for u=1:par.U
                L_H = inv_L*var.H(:,u,t);
                grad_f_without_regularizer(u,t) = 2*rho*C(u)*log2(exp(1))*var.H(:,u,t)'*L_H;
            end
    end
end

end

