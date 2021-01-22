load('flat_metrics.mat');

fnames = fieldnames(S_flat);

for f = 1:numel(fnames)
    figure;
    scatter(S_flat.(fnames{f}).overlap_mat(:,1), S_flat.(fnames{f}).par_mat(:,1));
    title(fnames{f});
    xlabel('Number of Overlapping Moves');
    ylabel('Parallel TPC');
end