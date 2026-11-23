% =========================================================================
% -- Channel Realization Loader
% -------------------------------------------------------------------------
% Description:
%   Load a channel realization, apply power control, and form the effective
%   channel used by the simulator.
%
% Input:
%   par: Simulation parameter structure.
%   var: Simulation state structure containing the channel, noise, and
%   scheduler variables.
%
% Output:
%   var: Updated simulation state structure.
%
% -- (c) 2026 Victoria Palhares, Christoph Studer
% =========================================================================

function [var] = gen_channel(par,var)

%% Load channel realization
FileLoad=['/path/to/channel/',par.system_type,'_',num2str(par.B),'_',num2str(par.U),'/channels_',...
    par.system_type,'_',num2str(par.B),'_',num2str(par.U),'_',num2str(var.ch),'.mat'];
load (FileLoad,'H');

var.G = H(:,1:par.U);

% Apply power control
switch(par.pwr_cntrl)
    case 'on' % Power control enabled
        % Squared L2 norm of each channel vector
        gpwr= sum(abs(var.G).^2,1);
        gmin = min(gpwr); % Minimum channel-vector power
        p = 10^(par.P_db/10);
        lambda = min([gpwr;gmin*p*ones(size(gpwr))])./gpwr;
        Lambda_vector = sqrt(lambda);
    otherwise
        Lambda_vector = ones(par.U,1);
end

switch (par.plot_type)
    case 'CDF'
        % Form the effective channel after power control
        var.H = sqrt(par.rho_u)*var.G*diag(Lambda_vector);
    case 'SNR'
        var.H = var.G*diag(Lambda_vector);
end

% Compute the receive-power dynamic range
gpwr= sum(abs(var.H).^2,1);
gpwr_db = 10.*log10(gpwr); % [dB]
dynamic_range = max(gpwr_db)-min(gpwr_db); % [dB]

end
