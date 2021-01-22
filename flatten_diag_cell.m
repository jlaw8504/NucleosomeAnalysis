function mat = flatten_diag_cell(cell_array)
%%flatten_cell parses cell array of diagonal matrices and returns matrix
%
%   input :
%       cell_array : A cell array containing diagonal matrices where each
%       row represents a biological cell and each column represents a time
%       increment (tau).
%
%   output :
%       mat : A 2D matrix where each row represents a biologcial cell and
%       each column represents a time increment (tau).

flat_cell = cellfun(@(x) x(:), cell_array, 'UniformOutput', false);
mat = zeros([numel(vertcat(flat_cell{:,1})), size(cell_array, 2)]);
for n = 1:size(cell_array,2)
    mat(:,n) = vertcat(flat_cell{:,n});
end