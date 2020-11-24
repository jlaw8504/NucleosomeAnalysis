function S = nucleosome_foci_sort(csv_filename, tau, overlap, run_length)
%%foci_sort parses a csv file of foci positions and places them in matrices
%%where every row is a timepoint and every column is a foci
%
%   inputs :
%       csv_filename : A string variable pointing to a csv file of foci
%       trajectories to parse. First column is x position, second is y
%       position (in pixels), third column is frame/timepoint, and fourth
%       column is foci index.
%
%       tau : A scalar variable that specifies the multiple of the base
%       imaging rate.
%
%       overlap : A scalar variable that specifies the amount of overlap
%       two traces need to have to be measured.
%
%
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
%
%           par_cell : A cell array containing matrices
%           of the parallel two-point correlation values for 
%           each foci pair. Each matrix had upper-right values set to NAN
%           to remove duplicate values. Each cell represents a single
%           timelpase.
%
%           perp_cell : A cell array containing matrices
%           of the perpendicular two-point correlation values for 
%           each foci pair. Each matrix had upper-right values set to NAN
%           to remove duplicate values. Each cell represents a single
%           timelpase.
%
%           dist_cell : A cell array containing matrices
%           of the foci-foci distances for 
%           each foci pair. Each matrix had upper-right values set to NAN
%           to remove duplicate values. Each cell represents a single
%           timelpase.
%
%           par_array, perp_array, and dist_array are column vectors of the
%           values in par_cell, perp_cell, and dist_cell. These vectors
%           combine all the non-NAN values from all timelapses.
%% Store metadata in the structural array
S.in_file = csv_filename;
S.tau = tau;
S.overlap = overlap;
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
    S.(label).par_cell = cell([numel(cell_labels), 1]);
    S.(label).perp_cell = cell([numel(cell_labels), 1]);
    S.(label).dist_cell = cell([numel(cell_labels), 1]);
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
        %% filter out all but longest run per foci trace
        for foci = 1:size(S.(label).x_cell{j}, 1)
            row = S.(label).x_cell{j}(foci,:);
            row_cc = bwconncomp(~(row == 0));
            run_sizes = cellfun(@length, row_cc.PixelIdxList);
            [max_run_size, max_run_idx] = max(run_sizes);
            filt_row = zeros(size(row));
            if max_run_size >= run_length
                filt_row(row_cc.PixelIdxList{max_run_idx}) = 1;
            end
            S.(label).x_cell{j}(foci,:) = S.(label).x_cell{j}(foci,:) .* filt_row;
            S.(label).y_cell{j}(foci,:) = S.(label).y_cell{j}(foci,:) .* filt_row;
            S.(label).border_cell{j}(foci,:) = S.(label).border_cell{j}(foci,:) .* filt_row;
        end
        %% replace missing entries with NaNs
        S.(label).x_cell{j}(S.(label).x_cell{j} == 0) = nan;
        S.(label).y_cell{j}(S.(label).y_cell{j} == 0) = nan;
        S.(label).border_cell{j}(S.(label).border_cell{j} == 0) = nan;
        %% Calculate TPC and distance matrices
        [ ...
            S.(label).par_cell{j}, ...
            S.(label).perp_cell{j}, ...
            S.(label).dist_cell{j} ...
            ] = calc_tpc( S.(label).x_cell{j},  S.(label).y_cell{j}, tau, overlap, 'off');
        
    end
    %% Gather up all TPC values and distances into arrays
    par_filt_cell = ...
        cellfun(@(x) x(~isnan(x)), ...
        S.(label).par_cell, 'UniformOutput', false);
    S.(label).par_array = vertcat(par_filt_cell{:});
    
    perp_filt_cell = ...
        cellfun(@(x) x(~isnan(x)), ...
        S.(label).perp_cell, 'UniformOutput', false);
    S.(label).perp_array = vertcat(perp_filt_cell{:});
    
    dist_filt_cell = ...
        cellfun(@(x) x(~isnan(x)), ...
        S.(label).dist_cell, 'UniformOutput', false);
    S.(label).dist_array = vertcat(dist_filt_cell{:});
end