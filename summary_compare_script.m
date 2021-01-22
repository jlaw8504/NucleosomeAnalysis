%% Mock vs BLM vs ATP vs Fixed
tau = 50;
[c_drugs_par, m_drugs_par, g_drugs_par] = jgl_anova1( ...
    { ...
    S_flat.MOCK.par_mat(:,tau), ...
    S_flat.BLM.par_mat(:,tau), ...
    S_flat.ATP.par_mat(:,tau), ...
    S_flat.Fixed.par_mat(:,tau) ...
    }, ...
    { ...
    'Mock', 'Bleomycin', 'ATP-depeleted', 'Fixed'});
[c_drugs_perp, m_drugs_perp, g_drugs_perp] = jgl_anova1( ...
    { ...
    S_flat.MOCK.perp_mat(:,tau), ...
    S_flat.BLM.perp_mat(:,tau), ...
    S_flat.ATP.perp_mat(:,tau), ...
    S_flat.Fixed.perp_mat(:,tau) ...
    }, ...
    { ...
    'Mock', 'Bleomycin', 'ATP-depeleted', 'Fixed'});
%% Plot Parallel
figure;
bar(m_drugs_par(:,1));
hold on;
errorbar(m_drugs_par(:,1), m_drugs_par(:,2), 'Marker', 'none', 'LineStyle', 'none', 'Color', 'k');
hold off;
xticklabels({'Mock', 'Bleomycin', 'ATP-depeleted', 'Fixed'});
ylabel('TPC Parallel');
xlim([0.5, 4.5]);
ylim([0, 0.1]);
if tau > 1
    txt = ' timesteps';
else
    txt = ' timestep';
end
title(['Tau is ', num2str(tau), txt]);

%% Plot Perpindicular
figure;
bar(m_drugs_perp(:,1));
hold on;
errorbar(m_drugs_perp(:,1), m_drugs_perp(:,2), 'Marker', 'none', 'LineStyle', 'none', 'Color', 'k');
hold off;
xticklabels({'Mock', 'Bleomycin', 'ATP-depeleted', 'Fixed'});
ylabel('TPC Perpindicular');
xlim([0.5, 4.5]);
ylim([0, 0.1]);
if tau > 1
    txt = ' timesteps';
else
    txt = ' timestep';
end
title(['Tau is ', num2str(tau), txt]);