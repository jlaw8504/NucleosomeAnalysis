function [par_mat, perp_mat, dist_mat] = calc_tpc(x_mat, y_mat, tau, overlap, keep_duplicates)
%%calc_tpc calculated the parallel and perpendicular correlation values and
%%foci-foci distances
%
%   input :
%       x_mat : A 2D matrix where each row represents a foci,
%       every column represents a frame, and each value is an x position in
%       pixels.
%
%       y_mat : A 2D matrix where each row represents a foci,
%       every column represents a frame, and each value is an y position in
%       pixels.
%
%   outputs :
%       par_mat : A symmetrical matrix  containing the parallel two-point
%       correlation values for each foci pair.
%
%       perp_mat : A symmetrical matrix  containing the perpendicular
%       two-point correlation values for each foci pair.
%
%       dist_mat : A symmetrical matrix  containing the foci-foci distance
%% Pre-allocate parallel and perpendicular rho matrices
par_mat = nan(size(x_mat,1));
perp_mat = par_mat;
%% Iterate over sorted data to calculate LOC point correlation
for i = 1:size(x_mat, 1)
    for j = 1:size(x_mat, 1)
        % evaluate if indices have overlapping timepoints
        time_idx = ~isnan(x_mat(i,:)) & ~isnan(x_mat(j,:));
        % if there are no overlapping timepoints assign rho of nan
        % This is a very STRINGENT filter, needs to have all timepoints
        % present!
        if sum(time_idx) < overlap
            par_mat(i, j) = nan;
            perp_mat(i, j) = nan;
            continue;
        end
        xD = x_mat(i, 1:end-tau) - x_mat(j, 1:end-tau);
        yD = y_mat(i, 1:end-tau) - y_mat(j, 1:end-tau);
        r = sqrt((xD.^2)+(yD.^2));
        LOC = [xD./r; yD./r];
        iMoves = [x_mat(i,tau+1:end) - x_mat(i,1:end-tau);...
            y_mat(i,tau+1:end) - y_mat(i,1:end-tau)];
        jMoves = [x_mat(j,tau+1:end) - x_mat(j,1:end-tau);...
            y_mat(j,tau+1:end) - y_mat(j,1:end-tau)];
        rhos = (dot(iMoves,LOC)).*(dot(jMoves,LOC));
        delts = [dot(iMoves,LOC).^2; dot(jMoves,LOC).^2];
        normz = (nanmean(delts(1,:)))*(nanmean(delts(2,:)));
        par_mat(i,j) = nanmean(rhos)/sqrt(normz);
        % Find correlation perpendicular to LOC
        LOCperp = [-yD./r; xD./r];
        rhos = (dot(iMoves,LOCperp)).*(dot(jMoves,LOCperp));
        delts = [dot(iMoves,LOCperp).^2; dot(jMoves,LOCperp).^2];
        normz = (nanmean(delts(1,:)))*(nanmean(delts(2,:)));
        perp_mat(i,j) = nanmean(rhos)/sqrt(normz);
    end
end
%% Calcuate the pair-wise Euclidean distance between all points
 dist_mat = squareform(pdist([nanmean(x_mat,2), nanmean(y_mat,2)]));
 dist_mat(dist_mat == 0) = nan;
 %remove distances that do not have associated TPC values
 dist_mat(isnan(par_mat)) = nan;
%% Remove duplicates
if strcmp(keep_duplicates, 'off')
    idx = tril(ones(size(par_mat)));
    par_mat = par_mat .* idx;
    par_mat(par_mat == 0) = nan;
    perp_mat = perp_mat .* idx;
    perp_mat(perp_mat == 0) = nan;
    % only keep distances with par and perp correlations
    dist_mat(isnan(par_mat)) = nan;
end