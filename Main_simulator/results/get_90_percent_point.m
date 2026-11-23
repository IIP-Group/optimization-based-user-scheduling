% =========================================================================
% -- 90th-Percentile Extraction
% -------------------------------------------------------------------------
% Description:
%   Extract 90th-percentile operating points from the generated CDF and
%   CCDF figure files and save the values.
%
% Input:
%   path: Directory containing saved simulation result files.
%   par: Simulation parameter structure.
%   file_title: Figure filename suffix identifying the simulation setup.
%
% Output:
%   None. Percentile values and a formatted table are saved to MAT files.
%
% -- (c) 2026 Victoria Palhares, Christoph Studer
% =========================================================================

function get_90_percent_point(path,par,file_title)
close all;

types = {'BER_CDF_','MI_CCDF_','MSE_CDF_','Rate_CCDF_'};
x_90_percent = cell(length(types),1);
for i=1:length(types)
    close;
    openfig([path,'/figures/',types{i},file_title]);
    h = findobj(gca,'Type','line');
    x=get(h,'Xdata');
    y=get(h,'Ydata');
    x_interest = zeros(length(x),1);
    x_interest_transform = zeros(length(x),1);
    a = zeros(length(x),1);
    b = zeros(length(x),1);
    x_transform = cell(length(x),1);
    for j=1:length(x)
        y_interest = [0.90];
        switch(types{i})
            case {'BER_CDF_','MSE_CDF_'}
                x_transform{j} = log10(x{j});
            case {'MI_CCDF_','Rate_CCDF_'}
                x_transform{j} = x{j};
        end
        % Find the samples nearest to the target probability
        temp = abs(y{j}-y_interest);
        idx_1 = find(temp == min(temp));
        switch(types{i})
            case {'BER_CDF_','MSE_CDF_'}
                idx_1 = idx_1(1);
                closest_y_1 = y{j}(idx_1);
                if(y{j}(idx_1) > y_interest)
                    idx_2 = idx_1-1;
                    closest_y_2 = y{j}(idx_2);
                elseif(y{j}(idx_1) < y_interest)
                    idx_2 = idx_1+1;
                    closest_y_2 = y{j}(idx_2);
                end
            case {'MI_CCDF_','Rate_CCDF_'}
                idx_1 = idx_1(length(idx_1));
                closest_y_1 = y{j}(idx_1);
                if(y{j}(idx_1) > y_interest)
                    idx_2 = idx_1+1;
                    closest_y_2 = y{j}(idx_2);
                elseif(y{j}(idx_1) < y_interest)
                    idx_2 = idx_1-1;
                    closest_y_2 = y{j}(idx_2);
                end
        end
        % Linearly interpolate between the two neighboring samples
            a(j) = (closest_y_1-closest_y_2)/(x_transform{j}(idx_1)-x_transform{j}(idx_2));
            b(j) = closest_y_1-a(j)*x_transform{j}(idx_1);
        % Evaluate the interpolant at the target probability
        x_interest(j) = (y_interest-b(j))/a(j);
        switch(types{i})
            case {'BER_CDF_','MSE_CDF_'}
                x_interest_transform(j) = 10^(x_interest(j));
            case {'MI_CCDF_','Rate_CCDF_'}
                x_interest_transform(j) = x_interest(j);
        end

    end
    clear h;
    x_90_percent{i} = flip(real(x_interest_transform));
end

FileName=[datestr(now, 'yyyy_mmm_dd'),'_USER_SCHEDULING_',num2str(par.B),'_',num2str(par.U),'_',num2str(par.timeslots),'_90_percentile.mat'];
save (FileName,'x_90_percent');
close all;

formatting_table(FileName,types,par);

end
