% EXPORT_SENSITIVITY_FIGURE  CSV data for the pgfplots version of the identifiability figure (letter, R2-C5).
load('reviewer2_results.mat');
F = '../../PAPER-2026-BMM-FUSION/JSV_Submission/rebuttal/2026-09-11-round-1/figures/';
t = D.t; k = find(t <= 8); k = k(1:2:end);
writetable(table(t(k)', 1e3*D.S_M(k)', 1e3*D.S_v(k)', 'VariableNames', {'t','SM','Sv'}), [F 'sens_sensitivity.csv']);
k = find(t >= 1 & t <= 15); k = k(1:5:end);
writetable(table(t(k)', D.crlb_M(k)', D.crlb_v(k)', 'VariableNames', {'t','cM','cv'}), [F 'sens_crlb.csv']);
disp('csv ok')
