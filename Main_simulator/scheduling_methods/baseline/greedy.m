% =========================================================================
% -- Greedy Scheduling
% -------------------------------------------------------------------------
% Description:
%   Generate a scheduling matrix using the greedy sum-rate scheduling
%   baseline.
%
% Input:
%   par: Simulation parameter structure.
%   var: Simulation state structure containing the channel, noise, and
%   scheduler variables.
%
% Output:
%   C_matrix: Generated scheduling matrix.
%
% -- (c) 2026 Victoria Palhares, Christoph Studer
% =========================================================================

function [C_matrix] = greedy(par,var)

% Transform the channel to the DFT beam domain
A_0 = dftmtx(par.B); % DFT matrix
A = zeros(par.B); % Normalize each DFT column
for i=1:par.B
    A(:,i) = A_0(:,i)/norm(A_0(:,i));
end
% Initialize the candidate UE set
U_group = 1:par.U;
order = zeros(par.Us,par.frequency_slots);
order_asc = zeros(par.Us,par.frequency_slots);
C_matrix = zeros(par.U,par.frequency_slots);

for t=1:par.frequency_slots
    H_b = A*var.H(:,:,t); % Beam-domain channel matrix
    i=1;
    H_gre = [];
    while (i<=par.Us)
        % Evaluate the sum rate for each candidate UE
        sum_rate = zeros(1,length(U_group));
        for j=1:length(U_group)
            H_test = [H_gre H_b(:,U_group(j))];
            [sum_rate(j)] = calculate_sum_rate(par,var,H_test);
        end
        [~,sel] = max(sum_rate);
        order(i,t) = U_group(sel);

        % Add the best candidate UE to the scheduled set
        H_gre = [H_gre H_b(:,U_group(sel))];

        % Update the candidate UE set
        U_group(sel) = [];

        i=i+1;
    end
    order_asc(:,t) = sort(order(:,t),'ascend');
    C_matrix(order_asc(:,t),t) = 1;
end

end
