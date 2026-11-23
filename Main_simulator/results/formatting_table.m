% =========================================================================
% -- 90th-Percentile Table Formatter
% -------------------------------------------------------------------------
% Description:
%   Format saved 90th-percentile metric values into the table string used
%   for reporting simulation results.
%
% Input:
%   FileName: MAT-file name containing saved percentile data.
%   types: Metric labels used for table formatting.
%   par: Simulation parameter structure.
%
% Output:
%   None. The formatted table string is saved to a MAT file.
%
% -- (c) 2026 Victoria Palhares, Christoph Studer
% =========================================================================

function formatting_table(FileName,types,par)

load(FileName,'x_90_percent');

switch(par.sim_scenario)
    case 'a'
        new_order = {'AUA','random','SUS [16]','CSS [5]','greedy [5]','LoFi (K=1) [34]',...
            'LoFi (K=4) [34]','LoFi++ (K=1) [34]','LoFi++ (K=4) [34]','MSE-LMMSE','channel capacity',...
            'achievable rate','ES-MSE','ES-capacity','ES-rate'};
    case {'b'}
        new_order = {'AUA','random','SUS [16]','CSS [5]','greedy [5]','LoFi (K=1) [34]',...
            'LoFi (K=4) [34]','LoFi++ (K=1) [34]','LoFi++ (K=4) [34]','MSE-LMMSE','channel capacity',...
            'achievable rate'};
    case {'c'}
        new_order = {'AUA','random','SUS [16]','CSS [5]','greedy [5]','LoFi (K=1) [34]',...
            'LoFi (K=4) [34]','MSE-LMMSE','channel capacity',...
            'achievable rate'};
    case {'d'}
        new_order = {'AUA','random','SUS [16]','CSS [5]','greedy [5]','LoFi (K=1) [34]',...
            'LoFi (K=4) [34]','MSE-LMMSE','channel capacity',...
            'achievable rate'};
end

new_string = '';
for i=1:length(new_order)
    new_string = [new_string,new_order{i}];
    for j=1:length(types)
        new_string = ([new_string,' & $' num2str(round(x_90_percent{j}(i),3)),'$']);
    end
    new_string = [new_string,'\\'];
end

FileName=[datestr(now, 'yyyy_mmm_dd'),'_table_90_percentile.mat'];
save (FileName,'new_string');

end

