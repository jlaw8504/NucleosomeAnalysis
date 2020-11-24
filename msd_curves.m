fnames = fieldnames(S_msd);
idx = 1;
for f = 1:numel(fnames)
    try
        ens_msd = S_msd.(fnames{f}).ens_msd;
        msd_mat = S_msd.(fnames{f}).msd_mat;
        if strcmp(fnames{f}, 'LivingB') || ...
                strcmp(fnames{f}, 'LivingM')
            errorbar( ...
                (1/20)*(1:20), ...
                ens_msd(1:20) * 163^2, ...
                nanstd(msd_mat(1:20, :), 0, 2)/sqrt(size(msd_mat,2)) * 163^2 ...
                );
        else
            errorbar( ...
                (1/20)*(1:40), ...
                ens_msd(1:40) * 163^2, ...
                nanstd(msd_mat(1:40, :), 0, 2)/sqrt(size(msd_mat,2)) * 163^2 ...
                );
        end
        if idx == 1
            hold on;
        end
        legend_cell{idx} = fnames{f};
        idx = idx + 1;
    catch
        continue;
    end
end
hold off;
xlabel('Tau');
ylabel('Ensemble MSD');
legend(legend_cell);