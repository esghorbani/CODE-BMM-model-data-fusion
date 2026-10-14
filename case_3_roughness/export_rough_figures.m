% EXPORT_ROUGH_FIGURES  CSV data for the pgfplots versions of the roughness figures of the letter.
load('roughness_results.mat');
F = '../../PAPER-2026-BMM-FUSION/JSV_Submission/rebuttal/2026-09-11-round-1/figures/';
xg = linspace(0, par.L, 1201)';
T = table(xg, 1e3*eval_roughness(results(1).prof, xg), 1e3*eval_roughness(results(2).prof, xg), 1e3*eval_roughness(results(3).prof, xg), ...
    'VariableNames', {'x','rA','rB','rC'}); writetable(T, [F 'rough_profiles.csv']);
t = results(1).t; k = find(t <= 15); k = k(1:5:end);
T = table(t(k)', results(1).M_est(k)', results(2).M_est(k)', results(3).M_est(k)', results(1).v_est(k)', results(2).v_est(k)', results(3).v_est(k)', ...
    'VariableNames', {'t','MA','MB','MC','vA','vB','vC'}); writetable(T, [F 'rough_params.csv']);
k = find(t >= 12 & t <= 16); k = k(1:2:end);
fC_lp = lowpass_bw(results(3).fc_true, F_LP, dt); fCu_lp = lowpass_bw(results(3).fc_est, F_LP, dt); fCk_lp = lowpass_bw(known_data(3).fc_est, F_LP, dt);
T = table(t(k)', 1e-3*results(1).fc_true(k)', 1e-3*results(3).fc_true(k)', 1e-3*fC_lp(k)', 1e-3*fCu_lp(k)', 1e-3*fCk_lp(k)', ...
    'VariableNames', {'t','fA','fC','fC_lp','fC_unk_lp','fC_known_lp'}); writetable(T, [F 'rough_force.csv']);
k = find(t >= 13 & t <= 13.4);
T = table(t(k)', 1e-3*results(1).fc_true(k)', 1e-3*results(1).fc_est(k)', 1e-3*known_data(1).fc_est(k)', 1e-3*results(3).fc_true(k)', 1e-3*results(3).fc_est(k)', 1e-3*known_data(3).fc_est(k)', ...
    'VariableNames', {'t','fA_true','fA_unk','fA_known','fC_true','fC_unk','fC_known'}); writetable(T, [F 'rough_zoom.csv']);
disp('csv ok')
