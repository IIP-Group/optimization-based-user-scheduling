% =========================================================================
% -- Scheduling Constraint Check
% -------------------------------------------------------------------------
% Description:
%   Check whether a scheduling matrix satisfies the configured per-slot and
%   per-user resource constraints.
%
% Input:
%   par: Simulation parameter structure.
%   C: Scheduling matrix.
%
% Output:
%   status_C: Constraint-satisfaction status of the selected schedule.
%   constraints: Structure containing evaluated scheduling constraints.
%
% -- (c) 2026 Victoria Palhares, Christoph Studer
% =========================================================================

function [status_C,constraints] = satisfy_constraints(par,C)

[U] = size(C);

tol = 1e-10;

sum_columns = sum(C,1);
sum_rows = sum(C,2);

for i=1:length(sum_columns)
    if(sum_columns(i) < par.min_sum_columns)
        if(abs(sum_columns(i) - par.min_sum_columns) < tol)
            sum_columns(i) = par.min_sum_columns;
        end
    elseif(sum_columns(i) > par.max_sum_columns)
        if(abs(sum_columns(i) - par.max_sum_columns) < tol)
            sum_columns(i) = par.max_sum_columns;
        end
    end
end

for i=1:length(sum_rows)
    if(sum_rows(i) < par.min_sum_rows)
        if(abs(sum_rows(i) - par.min_sum_rows) < tol)
            sum_rows(i) = par.min_sum_rows;
        end
    elseif(sum_rows(i) > par.max_sum_rows)
        if(abs(sum_rows(i) - par.max_sum_rows) < tol)
            sum_rows(i) = par.max_sum_rows;
        end
    end
end

constraints = zeros(3,1);

% Check binary-valued entries
constraints(1) = (sum(sum((abs(C-0)<tol) | (abs(C-1)<tol))) == (U(1)*U(2)));
% Check row-sum constraints
constraints(2) = (sum(sum_rows >= par.min_sum_rows & sum_rows <= par.max_sum_rows)== U(1));
% Check column-sum constraints
constraints(3) = (sum(sum_columns >= par.min_sum_columns & sum_columns <= par.max_sum_columns)== U(2));

if (sum(constraints) == length(constraints))
    status_C = 1;
else
    status_C = 0;
end

end
