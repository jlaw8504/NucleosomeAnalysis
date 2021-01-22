load S_tpc_summary.mat;
tau = 10;

figure;
errorbar(S_tpc_summary.MOCK.par_mean_mat(:,tau), S_tpc_summary.MOCK.par_sem_mat(:,tau)); hold on;
errorbar(S_tpc_summary.Fixed.par_mean_mat(:,tau), S_tpc_summary.Fixed.par_sem_mat(:,tau));
errorbar(S_tpc_summary.BLM.par_mean_mat(:,tau), S_tpc_summary.BLM.par_sem_mat(:,tau));
errorbar(S_tpc_summary.ATP.par_mean_mat(:,tau), S_tpc_summary.ATP.par_sem_mat(:,tau));
hold off;
xlim([0,10]);
legend({'Mock', 'Fixed', 'BLM', 'ATP'});

figure;
errorbar(S_tpc_summary.MOCK.perp_mean_mat(:,tau), S_tpc_summary.MOCK.perp_sem_mat(:,tau)); hold on;
errorbar(S_tpc_summary.Fixed.perp_mean_mat(:,tau), S_tpc_summary.Fixed.perp_sem_mat(:,tau));
errorbar(S_tpc_summary.BLM.perp_mean_mat(:,tau), S_tpc_summary.BLM.perp_sem_mat(:,tau));
errorbar(S_tpc_summary.ATP.perp_mean_mat(:,tau), S_tpc_summary.ATP.perp_sem_mat(:,tau));
hold off;
xlim([0,10]);
legend({'Mock', 'Fixed', 'BLM', 'ATP'});