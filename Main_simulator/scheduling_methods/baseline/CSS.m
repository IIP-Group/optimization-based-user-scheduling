% =========================================================================
% -- Channel Structure-Based Scheduling
% -------------------------------------------------------------------------
% Description:
%   Generate a scheduling matrix using the channel structure-based
%   scheduling baseline.
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

function [C_matrix] = CSS(par,var)

% Transform the channel to the DFT beam domain
A_0 = dftmtx(par.B); % DFT matrix
A = zeros(par.B); % Normalize each DFT column
for i=1:par.B
    A(:,i) = A_0(:,i)/norm(A_0(:,i));
end
L = 4;
Nb = 2*L;
Ns = par.Us;

Hbeam = [];

% Initialize the candidate UE set
U_group = 1:par.U;

order = zeros(par.Us,par.frequency_slots);
order_asc = zeros(par.Us,par.frequency_slots);
C_matrix = zeros(par.U,par.frequency_slots);

for t=1:par.frequency_slots
    H_b = A*var.H(:,:,t); % Beam-domain channel matrix
    [~, temp] = sort(H_b,'descend');
    Bk = Nb+1 > temp; % Select the Nb strongest beam indices
    idx = Bk == 1; % Record the selected beam support
    Bk = rand(size(H_b));
    Bk(idx) = 0;

    epsilon = 0.3;
    Nol = 2;

    while (nnz(order(:,t)) < Ns)
        Nol = Nol+1;
        epsilon = epsilon+0.1;

        U_group_temp = U_group;
        H_css = [];
        f_0 = [];

        i = 1;

        while (i<=Ns && (sum(U_group_temp) ~= 0))
            % Select the candidate with the highest SINR
            R_u = zeros(1,length(U_group_temp));
            SINR_u = zeros(1,length(U_group_temp));
            f = zeros(par.B,length(U_group_temp));
            for j=1:length(U_group_temp)
                H_test = [H_css H_b(:,U_group_temp(j))];
                [SINR_u(j)] = calculate_SINR(var,H_test);
            end
            [~,sel] = max(SINR_u);
            order(i,t) = U_group_temp(sel);
            H_css = [H_css H_b(:,U_group_temp(sel))];

            % Compute the orthogonal channel component
            % Compute the component orthogonal to the selected subspace
            temp = zeros(par.B);
            for j=1:i-1
                temp = temp + (f_0(:,j)*f_0(:,j)')/norm(f_0(:,j))^2;
            end
            f(:,sel) = (eye(par.B)-temp) * H_b(:,U_group_temp(sel));
            f_0 = [f_0 f(:,sel)];

            % Apply semi-orthogonality and beam-overlap criteria
            % Store the beam support of the selected UE
            Hbeam = [Hbeam, Bk(:,U_group_temp(sel))];

            % Update the candidate UE set
            U_group_temp(sel) = [];
            temp = [];
            for j = 1:length(U_group_temp)
                if abs(f_0(:,i)'*H_b(:,U_group_temp(j)))/(norm(f_0(:,i))*norm(H_b(:,U_group_temp(j)))) < epsilon
                    % Limit beam-support overlap
                    if sum(Hbeam(:,i) == Bk(:,U_group_temp(j)) ) <= Nol
                        temp = [temp U_group_temp(j)];
                    end
                end
            end
            U_group_temp = temp;

            % Continue selecting UEs for the current time slot
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

