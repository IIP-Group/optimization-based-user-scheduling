% =========================================================================
% -- BEACHES Channel Denoiser
% -------------------------------------------------------------------------
% Description:
%   Denoise noisy channel vectors using SURE-based thresholding in the
%   antenna or beamspace domain.
%
% Input:
%   par: Simulation parameter structure.
%   Hn: Noisy channel matrix.
%   N0: Noise variance.
%   sim_mode: Domain selector for channel denoising.
%
% Output:
%   Hdenoised: Denoised channel matrix.
%   mse_avg: Average estimated denoising MSE.
%   mse: Per-user estimated denoising MSE.
%
% -- 2018 (c) rg548@cornell.edu, studer@cornell.edu, sm2675@cornell.edu
% =========================================================================

function [Hdenoised, mse_avg,mse] = BEACHES(par,Hn, N0, sim_mode)
scalemax = 0.8;
hdenoised = zeros(size(Hn));
SURE = zeros(par.B,1);
mse = zeros(par.U,1);
for uu=1:par.U

    if strcmp(sim_mode, 'antenna_domain')
        hnoisy = fft(Hn(:,uu))/sqrt(par.B);
    else
        hnoisy = Hn(:,uu);
    end

    tau_max = scalemax*max(abs(hnoisy));
    N = length(hnoisy);
    sorth = sort(abs(hnoisy),'ascend');
    recip = 1./sorth;
    recip(sorth < 2^(-10)) = 0;
    cumsum2 = 0;
    cumsuminv = sum(recip);
    tau_opt = inf;
    suremin = inf;
    for bb = 1:N
        tau = sorth(bb);
        SURE(bb) = cumsum2 + (N-bb+1)*tau^2 + N*N0 - 2*N0*(bb-1)-tau*N0*cumsuminv;
        cumsum2 = cumsum2 + sorth(bb).^2;
        cumsuminv = cumsuminv - recip(bb);
        if SURE(bb)<suremin
            suremin = SURE(bb);
            tau_opt = tau;
        end
    end
    tau_opt = min(tau_opt, tau_max);
    hdenoised(:,uu) = exp(1i*angle(hnoisy)).*max(abs(hnoisy)-tau_opt,0);
    mse(uu) = suremin/N;
end
if strcmp(sim_mode, 'antenna_domain')
        Hdenoised = ifft(hdenoised)*sqrt(par.B);
    else
        Hdenoised = hdenoised;
end
mse_avg = mean(mse);
end
