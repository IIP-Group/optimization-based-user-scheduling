% =========================================================================
% -- Baseline Sum-Rate Evaluation
% -------------------------------------------------------------------------
% Description:
%   Compute the post-equalization sum rate used by baseline scheduling
%   methods for a candidate user set.
%
% Input:
%   par: Simulation parameter structure.
%   var: Simulation state structure containing the channel, noise, and
%   scheduler variables.
%   H_part: Channel matrix for the candidate scheduled user set.
%
% Output:
%   sum_rate: Post-equalization sum rate for the candidate user set.
%
% -- (c) 2026 Victoria Palhares, Christoph Studer
% =========================================================================

function [sum_rate] = calculate_sum_rate(par,var,H_part)

[~,Us] = size(H_part);

W = H_part/(H_part'*H_part+(var.N0_scheduling/var.Es)); % LMMSE equalizer

SINR_u = zeros(Us,1);
R_u = zeros(Us,1);

for u=1:Us
    des_sig = (var.Es*abs(W(:,u)'*H_part(:,u))^2); % Desired-signal power
    % Interference power
    interf = (var.Es*((abs(W(:,u)'*H_part).^2*ones(Us,1))- (abs(W(:,u)'*H_part(:,u))^2)));
    noise = (var.N0_scheduling*norm(W(:,u))^2); % Noise power

    SINR_u(u) = des_sig/(interf+ noise);
    R_u(u) = log2(1+real(SINR_u(u))); % Achievable rate for the current UE
end

sum_rate = sum(R_u);

end
