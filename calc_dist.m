function dist_mat = calc_dist(x_mat, y_mat)
%%calc_dist calculate the pair-wise distance of each foci to all other foci
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
%   output :
%       dist_mat : A symmetrical matrix  containing the foci-foci distance

%% Calcuate the pair-wise Euclidean distance between all points
 dist_mat = squareform(pdist([nanmean(x_mat,2), nanmean(y_mat,2)]));
 dist_mat(dist_mat == 0) = nan;