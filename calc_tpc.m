function [par_mat, perp_mat, overlap_mat] = calc_tpc(x_mat, y_mat, tau, keep_duplicates)
%%calc_tpc calculated the parallel and perpendicular correlation values and
%%foci-foci distances
%
%   inputs :
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
%       overlap_mat : A symmetrical matrix containing the number of
%       timepoints two foci overlap at a given tau.
%% Pre-allocate parallel and perpendicular rho matrices
par_mat = nan(size(x_mat,1));
perp_mat = nan(size(x_mat,1));
overlap_mat = nan(size(x_mat,1));
%% Iterate over sorted data to calculate LOC point correlation
for i = 1:size(x_mat, 1)
    for j = 1:size(x_mat, 1)
        if i == j
            continue;
        end
        %% TPC calculations
        xD = x_mat(i, 1:end-tau) - x_mat(j, 1:end-tau); %x-distance
        yD = y_mat(i, 1:end-tau) - y_mat(j, 1:end-tau); %y-distance
        r = sqrt((xD.^2)+(yD.^2)); %p-distance
        LOC = [xD./r; yD./r]; %normalize x and y distances by p-distance, creating a unit vector
        iMoves = [x_mat(i,tau+1:end) - x_mat(i,1:end-tau);...
            y_mat(i,tau+1:end) - y_mat(i,1:end-tau)];
        jMoves = [x_mat(j,tau+1:end) - x_mat(j,1:end-tau);...
            y_mat(j,tau+1:end) - y_mat(j,1:end-tau)];
        % old code that doesn't make sense to me
%         rhos = (dot(iMoves,LOC)).*(dot(jMoves,LOC));
%         delts = [dot(iMoves,LOC).^2; dot(jMoves,LOC).^2];
%         normz = (nanmean(delts(1,:)))*(nanmean(delts(2,:)));
%         par_mat(i,j) = nanmean(rhos)/sqrt(normz);
        %Dot product of Moves and LOC, it same as ||a||*cos(theta) or
        %magnitude of motion adjacent to
        iMovesLOC = dot(iMoves,LOC)';
        jMovesLOC = dot(jMoves,LOC)';
        overlap_mat(i,j) = sum(~isnan(iMovesLOC) & ~isnan(jMovesLOC));
        par_mat(i,j) = corr(iMovesLOC, jMovesLOC, 'rows', 'complete');
        % Find correlation perpendicular to LOC
        LOCperp = [-yD./r; xD./r];
        % old code that doesn't make sense to me
%         rhos = (dot(iMoves,LOCperp)).*(dot(jMoves,LOCperp));
%         delts = [dot(iMoves,LOCperp).^2; dot(jMoves,LOCperp).^2];
%         normz = (nanmean(delts(1,:)))*(nanmean(delts(2,:)));
%         perp_mat(i,j) = nanmean(rhos)/sqrt(normz);
        iMovesLOCperp = dot(iMoves,LOCperp)';
        jMovesLOCperp = dot(jMoves,LOCperp)';
        perp_mat(i,j) = corr(iMovesLOCperp, jMovesLOCperp, 'rows', 'complete');
    end
end
%% Remove duplicates
if strcmp(keep_duplicates, 'off')
    idx = tril(ones(size(par_mat)));
    par_mat = par_mat .* idx;
    par_mat(par_mat == 0) = nan;
    perp_mat = perp_mat .* idx;
    perp_mat(perp_mat == 0) = nan;
    % only keep distances with par and perp correlations
    overlap_mat(isnan(par_mat)) = nan;
end