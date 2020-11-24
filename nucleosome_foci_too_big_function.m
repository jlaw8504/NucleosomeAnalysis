function tau_S = nucleosome_foci_too_big_function(S, tau, overlap, keep_duplicates)
%% Iterate over all S fields/Experimental Conditions
fnames = fieldnames(S);
for f = 1:numel(fnames)
    x_mat = S.(fnames{f}).x_matrix;
    y_mat = S.(fnames{f}).y_matrix;
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
    %% Aggregate and filter matrices
    if strcmp(keep_duplicates, 'off')
        idx_mat = find(tril(ones(size(par_mat))));
        par_mat(idx_mat) = nan;
        perp_mat(idx_mat) = nan;
        dist_mat(idx_mat) = nan;
    end
    tau_S.(fnames{f}).par_mat = par_mat;
    tau_S.(fnames{f}).perp_mat = perp_mat;
    tau_S.(fnames{f}).dist_mat = dist_mat;
    tau_S.(fnames{f}).par_array = par_mat(~isnan(par_mat));
    tau_S.(fnames{f}).perp_array = perp_mat(~isnan(perp_mat));
    tau_S.(fnames{f}).dist_array = dist_mat(~isnan(par_mat));
end