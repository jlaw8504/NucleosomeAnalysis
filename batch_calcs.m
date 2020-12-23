function S_metrics = batch_calcs(S, tpc_taus)
%%batch_calcs Appends metrics to input structural array
%
%   input :
%       S :
%
%   output :
%       S_metrics :
%
%% Iterate over fieldnames
fnames = fieldnames(S);
fwait = waitbar(0, 'Iterating over fields...');
for f = 1:numel(fnames)
    %% Try to calculate MSD matrix and ensemble MSD per field/condition
    try
        x_cell = S.(fnames{f}).x_cell;
        y_cell = S.(fnames{f}).y_cell;
    catch
        warning('Unable to parse x_cell or y_cell from %s', fnames{f});
        waitbar(f/numel(fnames), fwait);
        continue;
    end
    %pre-allocate cell arrays for each "cell" in set
    msd_cell = cell([numel(S.(fnames{f}).x_cell, 1), 1]);
    dist_cell = cell([numel(S.(fnames{f}).x_cell, 1), 1]);
    %add second dimension to correspond to tau
    par_cell = cell([numel(S.(fnames{f}).x_cell, 1), ...
        numel(tpc_taus)]);
    perp_cell = cell([numel(S.(fnames{f}).x_cell, 1), ...
        numel(tpc_taus)]);
    overlap_cell = cell([numel(S.(fnames{f}).x_cell, 1), ...
        numel(tpc_taus)]);
    % iterate over "cells"
    iwait = waitbar(0, 'Iterating over cells...');
    for i = 1:numel(x_cell)
        msd_cell{i} = calc_msd( ...
            x_cell{i}, ...
            y_cell{i} ...
            );
        dist_cell{i} = calc_dist( ...
            x_cell{i}, ...
            y_cell{i} ...
            );
        tau_idx = 1;
        for tau = tpc_taus
            [par_cell{i,tau_idx}, perp_cell{i,tau_idx}, overlap_cell{i,tau_idx}] = ...
                calc_tpc( ...
                x_cell{i}, ...
                y_cell{i}, ...
                tau, ...
                0 ...
                );
            tau_idx = tau_idx + 1;
        end
        waitbar(i/numel(x_cell), iwait);
    end
    %% Append data to S
    S.tpc_taus = tpc_taus;
    S.(fnames{f}).msd_cell = msd_cell;
    S.(fnames{f}).dist_mat = dist_cell;
    S.(fnames{f}).par_cell = par_cell;
    S.(fnames{f}).perp_cell = perp_cell;
    S.(fnames{f}).overlap_cell = overlap_cell;
    waitbar(f/numel(fnames), fwait);
    close(iwait);
end
close(fwait);
S_metrics = S;