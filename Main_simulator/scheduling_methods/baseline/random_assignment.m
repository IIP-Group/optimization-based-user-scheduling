% =========================================================================
% -- Random Scheduling Baseline
% -------------------------------------------------------------------------
% Description:
%   Generate a random scheduling assignment and evaluate its
%   post-equalization MSE.
%
% Input:
%   par: Simulation parameter structure.
%   var: Simulation state structure containing the channel, noise, and
%   scheduler variables.
%
% Output:
%   C_matrix: Generated scheduling matrix.
%   sum_MSE: Sum-MSE value associated with the random assignment.
%
% -- (c) 2026 Victoria Palhares, Christoph Studer
% =========================================================================

function [C_matrix,sum_MSE] = random_assignment(par,var)

C_matrix = zeros(par.U,par.frequency_slots);

for i=1:par.frequency_slots
    for j=1:par.Us
        C_matrix((i-1)*par.Us+j,i) = 1;
    end
end

MSE = zeros(par.frequency_slots,1);
for t=1:par.frequency_slots
    C = C_matrix(:,t);
    Hc = var.H(:,:,t)*diag(C);
    L = eye(par.B)/(Hc*Hc'+(var.N0_scheduling/var.Es)*eye(par.B));
    M = Hc*diag(C.^2)*Hc';
    MSE(t) = var.Es*trace(diag(C.^2))-var.Es*trace(L*M);

end
sum_MSE = real(sum(MSE));

end

