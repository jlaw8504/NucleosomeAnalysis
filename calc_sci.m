function sci = calc_sci(coords, L)
%%calc_sci calculates the speed correlation index for a given time window
%%of size L
%
%   inputs :
%       coords : A 2D matrix variable where each row represents a time
%       point and the columns are X and Y coordinates respectively.
%
%       L : A scalar varible specifying the averaging time window in number
%       of time points, i.e. 2.
%
%   output :
%       sci : The mean linear correlation of different time windows over
%       time.

%% Define terms
N = size(coords,1); %number time points
i = 1:N; %set of time points
filter_i = (i >= L) & (i <= (N -L)); %index of time points to include in sci calc
C_mat = zeros([N,L]);
C_mat(C_mat == 0) = nan;
for k = 1:L %number of sets of time window segments;
    shift = k - 1;
    start = coords(1+shift:L:end,:);
    start_i = i(1+shift:L:end)';
    finish = coords(L+shift:L:end,:);
    finish_i = i(L+shift:L:end)';
    try
        disp = start - finish;
    catch
        start = start(1:end-1, :);
        start_i = start_i(1:end-1);
        disp = start - finish;
    end
    for j = 1:size(disp, 1)-1
        C{j,1} = dot(disp(j,:)', disp(j+1,:)')/...
            (norm(disp(j,:)')*norm(disp(j+1,:)'));
        C{j,2} = start_i(j):finish_i(j); 
    end
    for i2 = 1:N
        for j2 = 1:size(C,1)
            if ismember(i2, C{j2,2})
                C_mat(i2,k) = C{j2,1};
            end
        end
    end
    clear C;
end
sci = mean(C_mat, 2);
