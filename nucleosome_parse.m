function S = nucleosome_parse(csv_filename)
%%nucleosome_parse parse the nucleosome foci tracks CSV-file
%
%   input :
%       csv_filename : A string variable pointing to a csv file of foci
%       trajectories to parse. First column is x position, second is y
%       position (in pixels), third column is frame/timepoint, and fourth
%       column is foci index.
%   output :
%       S : A structural array where each field is an experimental
%       condition. Each experimental condition is a structural array
%       contianing the following fields :
%
%           x_cell : A cell array containing 2D matrices where each row represents a foci,
%           every column represents a frame, and each value is an x position in
%           pixels. Each cell represents a single timelapse
%
%           y_cell : A cell array containing 2D matrices where each row represents a foci,
%           every column represents a frame, and each value is an y position in
%           pixels. Each cell represents a single timelpase.
%
%           border_cell : A cell array containing 2D matrices where each
%           row represents a foci, each column represents a frame, and each
%           value is the distance of the foci from the nearest nuclear
%           border. Each cell represents a single timelapse

%% Store metadata in the structural array
S.in_file = csv_filename;
%% read csv file
T = import_nucleosome_all_data(csv_filename);
sub_data_all = table2array(T(:, {'x', 'y', 'frame', 'particle', 'dist_to_boundary'}));
exp_labels = unique(T.exp_label);
for i = 1:numel(exp_labels)
    %% Create fieldname from exp_label
    label = char((exp_labels(i)));
    label = label(5:end);
    %% Determine number of cells per exp_label and pre-allocate
    cell_labels = table2cell(unique(...
        T(T.exp_label == exp_labels(i), 'raw_data')));
    S.(label).x_cell = cell([numel(cell_labels), 1]);
    S.(label).y_cell = cell([numel(cell_labels), 1]);
    S.(label).border_cell = cell([numel(cell_labels), 1]);
    for j = 1:numel(cell_labels)
        cell_idx = T.raw_data == cell_labels{j};
        subset_exp = T(cell_idx, 'exp_label');
        if numel(unique(subset_exp)) ~= 1
            error('Experiment %s has error in raw_data labeling!', cell_labels{j});
        end
        sub_data = sub_data_all(cell_idx, :);
        %% pre-allocate cell arrays of x,y, and dist matrices
        % since foci indices are missing values, we need to re-index them
        foci_idxs = unique(sub_data(:,4));
        S.(label).x_cell{j} = zeros([numel(foci_idxs), max(sub_data(:,3)) + 1]);
        S.(label).y_cell{j} = S.(label).x_cell{j};
        S.(label).border_cell{j} = S.(label).x_cell{j};
        %% iterate over sub_data matrix
        for n = 1:size(sub_data,1)
            foci_idx = find(foci_idxs == sub_data(n,4));
            S.(label).x_cell{j}(foci_idx, sub_data(n,3) + 1) = sub_data(n,1);
            S.(label).y_cell{j}(foci_idx, sub_data(n,3) + 1) = sub_data(n,2);
            S.(label).border_cell{j}(foci_idx, sub_data(n,3) + 1) = sub_data(n,5);
        end
        %% replace missing entries with NaNs
        S.(label).x_cell{j}(S.(label).x_cell{j} == 0) = nan;
        S.(label).y_cell{j}(S.(label).y_cell{j} == 0) = nan;
        S.(label).border_cell{j}(S.(label).border_cell{j} == 0) = nan;       
    end
end