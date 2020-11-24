function S_msd = batch_msd(S)

fnames = fieldnames(S);
for f = 1:numel(fnames)
    try
        num_cells = numel(S.(fnames{f}).x_cell);
        msd_mat = zeros( ...
            [ ...
            size(S.(fnames{f}).x_cell{1}, 2)-1, ...
            num_cells ...
            ] ...
            );
        for c = 1:num_cells
            msd_mat(:,c) = calc_msd( ...
                S.(fnames{f}).x_cell{c}, ...
                S.(fnames{f}).y_cell{c} ...
                );
        end
        S.(fnames{f}).msd_mat = msd_mat;
        S.(fnames{f}).ens_msd = nanmean(msd_mat,2);
    catch
        continue;
    end
end
S_msd = S;