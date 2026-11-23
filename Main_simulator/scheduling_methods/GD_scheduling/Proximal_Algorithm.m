% =========================================================================
% -- Simplex Projection with Inequality Bounds
% -------------------------------------------------------------------------
% Description:
%   Project each matrix column onto the box-constrained simplex defined by
%   lower and upper sum bounds.
%
% Input:
%   Z: Matrix to be projected.
%   a: Lower bound on the projected column sum.
%   b: Upper bound on the projected column sum.
%
% Output:
%   X: Projected matrix.
%
% -- (c) 2026 Victoria Palhares, Christoph Studer
% =========================================================================

function [X] = Proximal_Algorithm(Z,a,b)

[size_Z] = size(Z);

X = zeros(size_Z);

for t=1:size_Z(2)

%% Simplex projection
z = real(Z(:,t)); % Current column
n = length(z);
z_sorted = sort(z,'descend'); % Sorted column
z_shift_lambda = z_sorted-z_sorted(1); % Shifted sorted column

%% Linear search over breakpoints
% Enumerate clipping breakpoints
% Sorted clipping breakpoints
[cumulative_shifts,idx] = sort([-z_shift_lambda;(1-z_shift_lambda)]);
see_shifts = zeros(n,length(cumulative_shifts));
for i=1:length(cumulative_shifts)
    see_shifts(:,i) = z_shift_lambda+cumulative_shifts(i);
end
shifts = (idx > n); % Identify lower- and upper-clipping breakpoints
% Candidate lambda offsets
difference = zeros(length(cumulative_shifts),1);
difference(1) = -z_sorted(1)+cumulative_shifts(1); % First lambda offset
intervals = zeros(length(cumulative_shifts),2);
intervals(1,:) = [-inf,-z_sorted(1)];

x = z_shift_lambda; % Shifted column
store_x = zeros(n,length(cumulative_shifts));
store_x(:,1) = x;
sum_x = zeros(length(cumulative_shifts),1);
sum_x(1) = sum(x);
final_MSE = inf;

valid_lambda = [];

m = sum(max(min(x,1),0)==1); % Number of upper-clipped entries
p = sum(max(min(x,1),0)==0);

% Construct intervals between consecutive breakpoints
for i=2:length(cumulative_shifts)
    current_valid_lambda = [];

    if(sum_x(i-1,1) > b)
        break;
    end
    % Within each interval, the clipping pattern remains fixed
    % Absolute lambda value
    difference(i) = -z_sorted(1)+cumulative_shifts(i);
    % Interval with a fixed clipping pattern
    intervals(i,:) = [difference(i-1),difference(i)];

    x = z+difference(i); % Evaluate x at the interval upper endpoint
    store_x(:,i) = x;
    sum_x(i,1) = sum(max(min(x,1),0)); % Compute the projected-column sum

    if(shifts(i-1) == 1)
        m = m+1; % Update the number of upper-clipped entries
    end
    if(shifts(i-1) == 0)
        p = p-1; % Update the number of lower-clipped entries
    end

    if(sum_x(i) >= a)
        % Current interval endpoints
        endpoints = [intervals(i,1),intervals(i,2)];

        % Sum entries that remain inside the box constraints
        z_j = sum(z_sorted(m+1:n-p));
        % Candidate lambda values
        current_lambdas = [(1/(n-p-m))*(a-m-z_j),(1/(n-p-m))*(b-m-z_j),0];
        % Tighten the lower endpoint when the previous sum is below a
        if(sum_x(i-1) < a)
            % Lambda satisfying sum(x) = a
            endpoints(1) = current_lambdas(1);
        end
        % Tighten the upper endpoint when the current sum exceeds b
        if(sum_x(i) > b)
            % Lambda satisfying sum(x) = b
            endpoints(2) = current_lambdas(2);
        end

        % Check the first candidate lambda
        check_bound_1 = current_lambdas(1) >= endpoints(1) && current_lambdas(1) <= endpoints(2);
        % Check the second candidate lambda
        check_bound_2 = current_lambdas(2) >= endpoints(1) && current_lambdas(2) <= endpoints(2);
        % Check the third candidate lambda
        check_bound_3 = current_lambdas(3) >= endpoints(1) && current_lambdas(3) <= endpoints(2);

        % At least one candidate lies in the current interval
        if(check_bound_1 == 1 || check_bound_2 == 1 || check_bound_3 == 1)
            if(check_bound_1 == 1) % First candidate is feasible
                % Add the feasible lambda candidate
                current_valid_lambda = [current_valid_lambda;current_lambdas(1)];
            end
            if(check_bound_2 == 1) % Second candidate is feasible
                % Add the feasible lambda candidate
                current_valid_lambda = [current_valid_lambda;current_lambdas(2)];
            end
            if(check_bound_3 == 1) % Third candidate is feasible
                % Add the feasible lambda candidate
                current_valid_lambda = [current_valid_lambda;current_lambdas(3)];
            end
        % If no candidate is feasible, use the lower interval endpoint
        else
            current_valid_lambda = [current_valid_lambda;intervals(i,1)];
        end
        % Store all feasible lambda candidates
        valid_lambda = [valid_lambda;current_valid_lambda];
        % Evaluate the feasible lambda candidates
        for j=1:length(current_valid_lambda)
            % Form the projected candidate
            possible_x = max(min(z+current_valid_lambda(j),1),0);
            % Evaluate the projection objective
            MSE = 1/2*(norm(possible_x-z))^2;
            % Retain the candidate with the lowest projection error
            if(MSE < final_MSE)
                final_MSE = MSE; % Update the best projection error
                % Update the best lambda
                final_lambda = current_valid_lambda(j);
                X(:,t) = possible_x; % Store the current best projection
                final_m = m;
                final_p = p;
            end
        end
    end
end
end

end

