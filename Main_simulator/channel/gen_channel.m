% =========================================================================
% -- Channel Realization Loader
% -------------------------------------------------------------------------
% Description:
%   Load a frequency-selective channel realization, apply power control
%   independently to each frequency slot, and form the effective channel.
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

%% Load frequency-selective channel realization
FileLoad=['/path/to/channel/',par.system_type,'_',num2str(par.B),'_',...
    num2str(par.U),'_QuadMMLoS_100MHz/channels_',par.system_type,'_',...
    num2str(par.B),'_',num2str(par.U),'_freq_sch_',num2str(var.ch),'.mat'];
load (FileLoad,'H');

for s=1:par.frequency_slots
    var.G(:,:,s) = H(:,1:par.U,s);

    % Apply power control independently to each frequency slot
    switch(par.pwr_cntrl)
        case 'on'
            gpwr = sum(abs(var.G(:,:,s)).^2,1);
            gmin = min(gpwr);
            pwr_ratio = 10^(par.P_db/10);
            lambda = min([gpwr;gmin*pwr_ratio*ones(size(gpwr))])./gpwr;
            Lambda_vector = sqrt(lambda);
        otherwise
            Lambda_vector = ones(par.U,1);
    end

    switch (par.plot_type)
        case 'CDF'
            var.H(:,:,s) = sqrt(par.rho_u)*var.G(:,:,s)*diag(Lambda_vector);
        case 'SNR'
            var.H(:,:,s) = var.G(:,:,s)*diag(Lambda_vector);
    end

    % Compute the receive-power dynamic range
    gpwr = sum(abs(var.H(:,:,s)).^2,1);
    gpwr_db = 10.*log10(gpwr);
    dynamic_range = max(gpwr_db)-min(gpwr_db);
end

end
