function msd_mat = calc_msd(x_mat, y_mat)
%%calc_msd calculates mean squared displacement for every possible tau
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
%       ens_msd_mat : A 2D matrix where each row represents a foci trace
%       and each column represents a tau
%% Pre-allocate
msd_mat = zeros([size(x_mat, 1), size(x_mat, 2)-1]);
for tau = 1:size(x_mat, 2)-1
    % X-coords
    X_t1 = x_mat(:,1:end-tau);
    X_t2 = x_mat(:,tau+1:end);
    X_sub = X_t2 - X_t1;
    % Y-coords
    Y_t1 = y_mat(:,1:end-tau);
    Y_t2 = y_mat(:,tau+1:end);
    Y_sub = Y_t2 - Y_t1;
    % Create squared displacement
    Disp_sq = X_sub.^2 + Y_sub.^2;
    msd_mat(:,tau) = nanmean(Disp_sq,2);
end