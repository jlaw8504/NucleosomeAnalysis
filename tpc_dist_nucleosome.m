function [mean_array, dist_cell, sem_array, bins ] = tpc_dist_nucleosome(S_cell,label_cell, metric, bin_size, pixel_size, max_dist)

%% extract distance and metric
data_cell = cellfun(@(x) x.(metric), S_cell, 'UniformOutput', false);
dist_cell = cellfun(@(x) x.dist_array * pixel_size, S_cell, 'UniformOutput', false);
%% Create distance bins
bins = 1:bin_size:max_dist;
%% Iterate over bins to create series of idx_cells
% preallocate data
mean_array = zeros([numel(bins)-1, numel(S_cell)]);
sem_array = mean_array;
vals = cell([numel(bins)-1, numel(S_cell)]);
for n = 1:numel(bins) - 1
    min_dist = bins(n);
    max_dist = bins(n+1);
    %% filter data by distance criteria
    min_idx_cell = cellfun(@(x) x>min_dist, dist_cell, 'UniformOutput', false);
    max_idx_cell = cellfun(@(x) x<=max_dist, dist_cell, 'UniformOutput', false);
    idx_cell = cellfun(@(x,y) x & y, min_idx_cell, max_idx_cell, 'UniformOutput', false);
    vals(n,:) = cellfun(@(x,y) x(y), data_cell, idx_cell, 'UniformOutput', false);
    mean_array(n,:) = cellfun(@(x,y) nanmean(x(y)), data_cell, idx_cell);
    sem_array(n,:) = cellfun(@(x,y) nanstd(x(y))/sqrt(sum(y)), data_cell, idx_cell);
end
for i = 1:numel(label_cell)
    if i == 1
        figure;
        errorbar(mean_array(:,i), sem_array(:,i));
        hold on;
    else
        errorbar(mean_array(:,i), sem_array(:,i));
    end
end
hold off;
legend(label_cell);
xlabel('Binned Distance (µm)');
ylabel(metric);
xticklabel_cell = arrayfun(@(x,y) sprintf('%d-%d', x, y), bins(1:end-1), bins(2:end), 'UniformOutput', false);
xticks(bins);
xticklabels(xticklabel_cell);