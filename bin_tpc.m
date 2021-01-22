function S_bin = bin_tpc(S, pixel_size_um, min_overlap)
%%binned_plot plot a metric versus binned distances of 1 micron
%
%   inputs :
%       metric_array : A column vector of metrics that has the same order
%       as the dist_array.
%
%       dist_array : A column vector of foci to foci distances in pixels
%
%       pixel_size_um : A scalar variable specifying the pixel size in
%       microns
%
%       min_overlap : An integer variable specifying the minimum number of
%       steps two foci must overlap
%% Filter par, perp and dist metrics by min_overlap
S.par_mat(S.overlap_mat < min_overlap) = nan;
S.perp_mat(S.overlap_mat < min_overlap) = nan;
%% Create 1 um bins from 0 to max(dist)
S.dist_array = S.dist_array * pixel_size_um;
for n = 0:ceil(max(S.dist_array))
    idx_mat = repmat(S.dist_array > n & S.dist_array <= (n+1), [1, size(S.par_mat,2)]);
    par_mat = S.par_mat;
    perp_mat = S.perp_mat;
    par_mat(~idx_mat) = nan;
    perp_mat(~idx_mat) = nan;
    S_bin.par_mean_mat(n+1, :) = nanmean(par_mat);
    S_bin.perp_mean_mat(n+1, :) = nanmean(perp_mat);
    S_bin.par_std_mat(n+1, :) = nanstd(par_mat);
    S_bin.perp_std_mat(n+1, :) = nanstd(perp_mat);
    S_bin.par_num_mat(n+1, :) = sum(~isnan(par_mat));
    S_bin.perp_num_mat(n+1, :) = sum(~isnan(perp_mat));
end
S_bin.par_sem_mat = S_bin.par_std_mat./sqrt(S_bin.par_num_mat);
S_bin.perp_sem_mat = S_bin.perp_std_mat./sqrt(S_bin.perp_num_mat);
    