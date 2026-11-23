% =========================================================================
% -- Achievable-Rate Objective
% -------------------------------------------------------------------------
% Description:
%   Evaluate the post-equalization achievable-rate scheduling objective and
%   its gradient.
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

function [f_without_regularizer,grad_f_without_regularizer] = Rate(par,var,C)
f_without_regularizer = 0;
grad_f_without_regularizer = zeros(par.U,par.frequency_slots);
for t = 1:par.frequency_slots
    Hc = var.H(:,:,t)*diag(C(:,t));
    L = eye(par.B)/(Hc*Hc'+(var.N0_scheduling/var.Es)*eye(par.B));

    L_H = zeros(par.B,par.U);
    H_L = zeros(par.U,par.B);
    for n = 1:par.U
        L_H(:,n) = L*var.H(:,n,t);
        H_L(n,:) = var.H(:,n,t)'*L;
    end

    % Evaluate the achievable-rate objective
    dF_rate = zeros(par.U,par.U);
    for u=1:par.U
        w_c_u = C(u,t)*(H_L(u,:)); % Equalizer vector for UE u
        h_c_u = C(u,t)*(var.H(:,u,t)); % Channel vector for UE u

        d_c = (var.Es*abs(w_c_u*h_c_u)^2); % Desired-signal power
        % Interference power
        i_c = (var.Es*((abs(w_c_u*Hc).^2*ones(par.U,1))- (abs(w_c_u*h_c_u)^2)));
        n_c = (var.N0_scheduling*norm(w_c_u)^2); % Noise power
        % Channel-estimation error power
        e_c_1 = var.Es*(norm(w_c_u)^2);
        e_c_2 = ones(1,par.U)*(diag(C(:,t).^2)*diag(var.covariance_channel_error(:,t)))*ones(par.U,1);
        e_c = e_c_1*e_c_2;

        SINR_u = d_c/(i_c+n_c+e_c);

        if(d_c == 0 && i_c == 0 && n_c == 0 && e_c == 0)
            SINR_u = 0;
            dF_SINR = zeros(par.U,1);
        else
            % Compute the objective gradient
            switch(var.calc_grad)
                case 'yes'
                    h_L_h = (var.H(:,u,t)'*L_H(:,u));

                    grad_interference = zeros(par.U,1);
                    grad_desired_signal = zeros(par.U,1);
                    grad_noise = zeros(par.U,1);
                    % First channel-estimation-error gradient term
                    grad_est_1 = zeros(par.U,1);

                    d_1_u = 4*var.Es*C(u,t)^3*(h_L_h)^2*(1-C(u,t)^2*(h_L_h));
                    factor_1 = 2*C(u,t)*(H_L(u,:)*L_H(:,u)-2*C(u,t)^2*H_L(u,:)*L_H(:,u)*h_L_h);
                    n_1_u = var.N0_scheduling*factor_1;
                    e_1_u = var.Es*factor_1;
                    for j=1:par.U
                        if j == u
                            grad_desired_signal(j) = d_1_u;
                            grad_noise(j) = n_1_u;
                            % First channel-estimation-error gradient term
                            grad_est_1(j) = e_1_u;
                        else
                            grad_desired_signal(j) = -4*var.Es*C(u,t)^4*C(j,t)*h_L_h*(var.H(:,j,t)'*L_H(:,u)*(var.H(:,u,t)'*L_H(:,j)));
                            factor_2 = -2*C(u,t)^2*C(j,t)*(H_L(j,:)*L_H(:,u)*(var.H(:,u,t)'*L_H(:,j))...
                                +H_L(u,:)*L_H(:,j)*(var.H(:,j,t)'*L_H(:,u)));
                            grad_noise(j) = var.N0_scheduling*factor_2;
                            % First channel-estimation-error gradient term
                            grad_est_1(j) = var.Es*factor_2;
                        end
                        interf = 0;
                        for i=1:par.U
                            if(i ~= u)
                                if(j == u)
                                    interf = interf+2*C(u,t)*C(i,t)^2*(var.H(:,i,t)'*L_H(:,u)*(var.H(:,u,t)'*L_H(:,i))...
                                        -2*C(u,t)^2*var.H(:,i,t)'*L_H(:,u)*h_L_h*(var.H(:,u,t)'*L_H(:,i)));
                                elseif(j == i)
                                    interf = interf+2*C(i,t)*C(u,t)^2*(var.H(:,i,t)'*L_H(:,u)*(var.H(:,u,t)'*L_H(:,i))...
                                        -2*C(i,t)^2*(var.H(:,u,t)'*L_H(:,i))*(var.H(:,i,t)'*L_H(:,i))*var.H(:,i,t)'*L_H(:,u));
                                else
                                    interf = interf-2*C(u,t)^2*C(i,t)^2*C(j,t)*(var.H(:,i,t)'*L_H(:,j)*(var.H(:,j,t)'*L_H(:,u))*(var.H(:,u,t)'*L_H(:,i))...
                                        +var.H(:,i,t)'*L_H(:,u)*(var.H(:,u,t)'*L_H(:,j))*(var.H(:,j,t)'*L_H(:,i)));
                                end
                            end
                        end
                        grad_interference(j) = var.Es*interf;
                    end
                    % Second channel-estimation-error gradient term
                    grad_est_2 = 2.*C(:,t).*var.covariance_channel_error(:,t);
                    % Apply the product rule
                    % Channel-estimation-error gradient
                    grad_ch_error= grad_est_1.*e_c_2+e_c_1.*grad_est_2;
                    grad_ch_error = real(grad_ch_error);

                    % Compute the SINR derivative
                    dF_SINR = ((grad_desired_signal.*(i_c+n_c+e_c)) - (d_c.*(grad_interference+grad_ch_error+grad_noise)))./((i_c+n_c+e_c)^2);
            end
        end

        switch(var.calc_grad)
            case 'yes'
                dF_rate(:,u) = log2(exp(1))*(1/(1+SINR_u))*dF_SINR;
        end

        f_without_regularizer = f_without_regularizer + log2(1+SINR_u);

    end
    switch(var.calc_grad)
        case 'yes'
            % Gradient of the sum rate over all UEs with respect to
            % each scheduling coefficient c_i.
            grad_f_without_regularizer(:,t) = sum(dF_rate,2);
    end
end
end
