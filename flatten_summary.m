function S_flat = flatten_summary(S_metrics)
%%fatten_summary reorganizes data in structrual array for stat analysis
%
%   input :
%       S_metrics : A structural array containing data from one or more
%       tracking experimental groups.
%
%   output :
%       S_flat : A structural array containing the same data but in array
%       format for faster statistical comparisions

fnames = fieldnames(S_metrics);
for f = 1:numel(fnames)
    try
        S_flat.(fnames{f}).par_mat = flatten_diag_cell(S_metrics.(fnames{f}).par_cell);
        S_flat.(fnames{f}).perp_mat = flatten_diag_cell(S_metrics.(fnames{f}).perp_cell);
        S_flat.(fnames{f}).overlap_mat = flatten_diag_cell(S_metrics.(fnames{f}).overlap_cell);
        S_flat.(fnames{f}).dist_array = flatten_diag_cell(S_metrics.(fnames{f}).dist_mat');
        S_flat.(fnames{f}).msd_mat = vertcat(S_metrics.(fnames{f}).msd_cell{:});
        S_metrics = rmfield(S_metrics, (fnames{f}));
    catch
        continue;
    end
end
