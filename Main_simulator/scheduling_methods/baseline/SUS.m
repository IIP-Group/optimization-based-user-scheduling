% =========================================================================
% -- Semi-Orthogonal User Selection
% -------------------------------------------------------------------------
% Description:
%   Generate a scheduling matrix using the semi-orthogonal user selection
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

function  [C_matrix] = SUS(par,var)

Ns = par.Us;

% Initialize the candidate UE set
U_group = 1:par.U;

order = zeros(par.Us,par.timeslots);
order_asc = zeros(par.Us,par.timeslots);
C_matrix = zeros(par.U,par.timeslots);

for t=1:par.timeslots
    epsilon = 0.5;

    while (nnz(order(:,t)) < Ns)
        epsilon = epsilon+0.1;
        U_group_temp = U_group;
        H_sus = [];
        g_0 = [];
        i = 1;
        Tlen = [];
        while(i<=Ns && isempty(U_group_temp)~=1)

            Tlen(i) = length(U_group_temp);
            g = zeros(par.B,length(U_group_temp));
            g_norm = zeros(length(U_group_temp),1);
            for k = 1:length(U_group_temp)

                % Compute components orthogonal to previously selected UEs
                temp = zeros(par.B);
                for j=1:i-1
                    temp = temp + (g_0(:,j)*g_0(:,j)')/norm(g_0(:,j))^2;
                end

                g(:,k) = (eye(par.B)-temp) * var.H(:,U_group_temp(k));
                g_norm(k) = norm(g(:,k));

            end
            % Select the UE with the largest orthogonal component
            [~,sel] = max(g_norm);
            order(i,t) = U_group_temp(sel);

            H_sus = [H_sus var.H(:,U_group_temp(sel))];
            g_0 = [g_0 g(:,sel)];

            % Prune UEs that are insufficiently semi-orthogonal
            U_group_temp(sel) = [];
            temp = [];
            for j = 1:length(U_group_temp)
                if abs(var.H(:,U_group_temp(j))'*g_0(:,i))/(norm(var.H(:,U_group_temp(j)))*norm(g_0(:,i))) < epsilon
                    temp = [temp U_group_temp(j)];
                end

            end
            U_group_temp = temp;

            i=i+1;
        end
    end
    order_asc(:,t) = sort(order(:,t),'ascend');
    C_matrix(order_asc(:,t),t) = 1;
    idx = zeros(length(order),1);
    for i=1:length(order)
        [~,idx(i)] = find(U_group == order_asc(i,t));
    end
    U_group(idx) = [];

end

end
