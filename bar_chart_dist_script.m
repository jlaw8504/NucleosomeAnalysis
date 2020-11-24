load 0_minTrace_80_overlap_tau1.mat;
mba = {'MOCK', 'BLM', 'Fixed'};
for n = 1:numel(mba)
    par = S.(mba{n}).par_array;
    perp = S.(mba{n}).perp_array;
    idx = S.(mba{n}).dist_array > 0/0.163 & S.(mba{n}).dist_array < 1.7/0.163;
    Sdist.(mba{n}).par_mean = mean(par(idx));
    Sdist.(mba{n}).par_sem = std(par(idx))/sqrt(sum(idx));
    Sdist.(mba{n}).perp_mean = mean(perp(idx));
    Sdist.(mba{n}).perp_sem = std(perp(idx))/sqrt(sum(idx));
    Sdist.(mba{n}).par = par(idx);
    Sdist.(mba{n}).perp = perp(idx);
end

bar([Sdist.MOCK.par_mean, Sdist.BLM.par_mean, Sdist.Fixed.par_mean]);
set(gca, 'XTickLabel', {'Mock', 'Bleomycin', 'Fixed'});
hold on;
errorbar([Sdist.MOCK.par_mean, Sdist.BLM.par_mean, Sdist.Fixed.par_mean], ...
    [Sdist.MOCK.par_sem, Sdist.BLM.par_sem, Sdist.Fixed.par_sem],...
    'LineStyle', 'none', 'Marker', 'none', 'Color', 'black');
ylabel('Parallel Correlation');
ylim([0, 0.15]);
xlim([0.5, 3.5]);
figure;
bar([Sdist.MOCK.perp_mean, Sdist.BLM.perp_mean, Sdist.Fixed.perp_mean]);
set(gca, 'XTickLabel', {'Mock', 'Bleomycin', 'Fixed'});
hold on;
errorbar([Sdist.MOCK.perp_mean, Sdist.BLM.perp_mean, Sdist.Fixed.perp_mean], ...
    [Sdist.MOCK.perp_sem, Sdist.BLM.perp_sem, Sdist.Fixed.perp_sem],...
    'LineStyle', 'none', 'Marker', 'none', 'Color', 'black');
ylabel('Perpendicular Correlation');
ylim([0, 0.15]);
xlim([0.5, 3.5]);
figure;
bar([Sdist.MOCK.par_mean, Sdist.BLM.par_mean]);
set(gca, 'XTickLabel', {'Mock', 'Bleomycin'});
hold on;
errorbar([Sdist.MOCK.par_mean, Sdist.BLM.par_mean], ...
    [Sdist.MOCK.par_sem, Sdist.BLM.par_sem],...
    'LineStyle', 'none', 'Marker', 'none', 'Color', 'black');
ylabel('Parallel Correlation');
ylim([0, 0.1]);
xlim([0.5, 2.5]);
figure;
bar([Sdist.MOCK.perp_mean, Sdist.BLM.perp_mean]);
set(gca, 'XTickLabel', {'Mock', 'Bleomycin'});
hold on;
errorbar([Sdist.MOCK.perp_mean, Sdist.BLM.perp_mean], ...
    [Sdist.MOCK.perp_sem, Sdist.BLM.perp_sem],...
    'LineStyle', 'none', 'Marker', 'none', 'Color', 'black');
ylabel('Perpendicular Correlation');
ylim([0, 0.1]);
xlim([0.5, 2.5]);