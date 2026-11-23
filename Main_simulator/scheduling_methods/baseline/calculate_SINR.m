% =========================================================================
% -- Baseline SINR Evaluation
% -------------------------------------------------------------------------
% Description:
%   Compute the minimum post-equalization SINR used by baseline scheduling
%   methods for a candidate user set.
%
% Input:
%   var: Simulation state structure containing the channel, noise, and
%   scheduler variables.
%   H_part: Channel matrix for the candidate scheduled user set.
%
% Output:
%   SINR_u_final: Minimum post-equalization SINR for the candidate user
%   set.
%
% -- (c) 2026 Victoria Palhares, Christoph Studer
% =========================================================================

function [SINR_u_final] = calculate_SINR(var,H_part)

[~,Us] = size(H_part);

W = H_part/(H_part'*H_part+(var.N0_scheduling/var.Es)); % MMSE Equalizer

SINR_u = zeros(Us,1);

for u=1:Us
    des_sig = (var.Es*abs(W(:,u)'*H_part(:,u))^2); % Desired-signal power
    % Interference power
    interf = (var.Es*((abs(W(:,u)'*H_part).^2*ones(Us,1))- (abs(W(:,u)'*H_part(:,u))^2)));
    noise = (var.N0_scheduling*norm(W(:,u))^2); % Noise power

    SINR_u(u) = des_sig/(interf+ noise);
end

SINR_u_final = SINR_u(Us);

end
