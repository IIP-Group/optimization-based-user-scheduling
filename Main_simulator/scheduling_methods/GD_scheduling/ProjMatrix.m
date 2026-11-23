% =========================================================================
% -- Projection onto Scheduling Constraints
% -------------------------------------------------------------------------
% Description:
%   Project a scheduling matrix onto the intersection of row and column
%   resource constraints using Douglas-Rachford splitting.
%
% Input:
%   par: Simulation parameter structure.
%   Z: Matrix to be projected.
%
% Output:
%   final_X: Projected scheduling matrix.
%
% -- (c) 2026 Victoria Palhares, Christoph Studer
% =========================================================================

function final_X = ProjMatrix(par,Z)

beta = 1; % Douglas-Rachford tuning parameter
tol = 1e-10;

% First indicator function: column-sum constraints
% Second indicator function: row-sum constraints

err_2 = inf;

maxiter = 2000;

size_Z = size(Z);

Yk = zeros(size_Z);

for ii=1:maxiter
    % Apply the column-sum proximal operator
    % Proximal update for the column constraints
    Xk = Proximal_Algorithm(((beta*Z+Yk)/(beta+1)),par.min_sum_columns,par.max_sum_columns);
    sum_rows = sum(Xk,2);
    % Check the row-sum bounds
    temp_1 = find(sum_rows < par.min_sum_rows);
    temp_2 = find(sum_rows > par.max_sum_rows);
    if(~isempty(temp_1)) % At least one sum is below the lower bound
        % Maximum absolute bound violation
        err_1 = norm(sum_rows(temp_1)-(ones(size(temp_1,1),1)*par.min_sum_rows),'inf');
    elseif(~isempty(temp_2)) % At least one sum is above the upper bound
        % Maximum absolute bound violation
        err_1 = norm(sum_rows(temp_2)-(ones(size(temp_2,1),1)*par.max_sum_rows),'inf');
    else % Row-sum bounds are satisfied
        err_1 = 0;
    end
    if (err_1<tol && err_2<tol)
        break;
    end

    % Apply the row-sum proximal operator
    Yk = (Yk-Xk+Proximal_Algorithm(((beta*Z+(2*Xk-Yk))/(beta+1)).',par.min_sum_rows,par.max_sum_rows).');
    sum_columns = sum(Xk,1);
    % Check the column-sum bounds
    temp_1 = find(sum_columns < par.min_sum_columns);
    temp_2 = find(sum_columns > par.max_sum_columns);
    if(~isempty(temp_1)) % At least one sum is below the lower bound
        % Maximum absolute bound violation
        err_2 = norm(sum_columns(temp_1)-(ones(1,size(temp_1,2))*par.min_sum_columns),'inf');
    elseif(~isempty(temp_2)) % At least one sum is above the upper bound
        % Maximum absolute bound violation
        err_2 = norm(sum_columns(temp_2)-(ones(1,size(temp_2,2))*par.max_sum_columns),'inf');
    else % Column-sum bounds are satisfied
        err_2 = 0;
    end
    if (err_1<tol && err_2<tol)
        break;
    end

end
if ii==maxiter
    warning('Maximum projection iterations (%i) reached.',maxiter);
end

final_X = Xk;

end

