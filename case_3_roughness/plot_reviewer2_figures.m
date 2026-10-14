% PLOT_REVIEWER2_FIGURES  Figures of BMM_UKF_Reviewer2_Studies (standalone).
% Loads reviewer2_results.mat when the study variables are not in the workspace.
if ~exist('A', 'var') || ~exist('C', 'var') || ~exist('D', 'var')
    load('reviewer2_results.mat');
end
fs = 10; fn = 'Times New Roman';
n = par.n_modes; t = D.t; S_M = D.S_M; S_v = D.S_v; crlb_M = D.crlb_M; crlb_v = D.crlb_v;
Qs = unique([A.Q]); als = unique([B.alfa]); vs = [C.v];
d_static = M1*par.g*par.L^3/(48*par.EI);
cols = lines(max(5, numel(vs)));
figure('Position', [100 100 500 230])
subplot(1,2,1); hold on
Ad = A([A.init] == 2); for k = 1:numel(Ad), plot(t, Ad(k).M, 'Color', cols(k,:), 'LineWidth', 1.1); end
plot(t, M1*ones(size(t)), 'k--', 'LineWidth', 1.2)
legend([arrayfun(@(q) sprintf('Q = 10^{%.1f}', log10(q)), Qs, 'UniformOutput', false), {'Ground Truth'}], 'Location', 'northoutside', 'NumColumns', 2, 'FontSize', 7)
ylabel('Moving mass (kg)', 'FontName', fn, 'FontSize', fs); xlabel('Time (sec)', 'FontName', fn, 'FontSize', fs)
ylim([0.5*M1 1.2*M1]); xlim([0 15]); set(gca, 'FontName', fn, 'FontSize', fs)
subplot(1,2,2); hold on
Bd = B([B.init] == 2); for k = 1:numel(Bd), plot(t, Bd(k).M, 'Color', cols(k,:), 'LineWidth', 1.1); end
plot(t, M1*ones(size(t)), 'k--', 'LineWidth', 1.2)
legend([arrayfun(@(a) sprintf('\\alpha = %.2f', a), als, 'UniformOutput', false), {'Ground Truth'}], 'Location', 'northoutside', 'NumColumns', 2, 'FontSize', 7)
ylabel('Moving mass (kg)', 'FontName', fn, 'FontSize', fs); xlabel('Time (sec)', 'FontName', fn, 'FontSize', fs)
ylim([0.5*M1 1.2*M1]); xlim([0 15]); set(gca, 'FontName', fn, 'FontSize', fs)
exportgraphics(gcf, 'Q_alpha_sensitivity.pdf', 'ContentType', 'vector', 'BackgroundColor', 'none');

figure('Position', [100 100 750 170])
subplot(1,3,1); hold on
for k = 1:numel(C), plot(linspace(0,1,round(par.L/vs(k)/dt)+1), 1e3*C(k).disp(1:round(par.L/vs(k)/dt)+1), 'Color', cols(k,:), 'LineWidth', 1.1); end
plot([0 1], 1e3*d_static*[1 1], 'k--', 'LineWidth', 1)
legend([arrayfun(@(v) sprintf('%d m/s', v), vs, 'UniformOutput', false), {'static'}], 'Location', 'best', 'FontSize', 7)
ylabel('Mid-span displacement (mm)', 'FontName', fn, 'FontSize', fs); xlabel('vt / L (first crossing)', 'FontName', fn, 'FontSize', fs)
set(gca, 'FontName', fn, 'FontSize', fs)
subplot(1,3,2); hold on
plot(vs, [C.dM2], 'k-o', 'LineWidth', 1, 'MarkerFaceColor', 'k'); plot(vs, [C.dv2], 'k:s', 'LineWidth', 1)
legend({'\Delta M', '\Delta v'}, 'Location', 'best', 'FontSize', 7)
ylabel('Converged error (%)', 'FontName', fn, 'FontSize', fs); xlabel('Speed (m/s)', 'FontName', fn, 'FontSize', fs)
set(gca, 'FontName', fn, 'FontSize', fs); xlim([8 35])
subplot(1,3,3); hold on
plot(vs, [C.DAF], 'k-o', 'LineWidth', 1, 'MarkerFaceColor', 'k')
ylabel('DAF of mid-span displacement', 'FontName', fn, 'FontSize', fs); xlabel('Speed (m/s)', 'FontName', fn, 'FontSize', fs)
yyaxis right; plot(vs, [C.fdyn_max], 'r-s', 'LineWidth', 1, 'MarkerFaceColor', 'r'); ylabel('max |f_c - Mg| / Mg (%)', 'FontSize', fs)
set(gca, 'FontName', fn, 'FontSize', fs); xlim([8 35])
exportgraphics(gcf, 'Speed_study.pdf', 'ContentType', 'vector', 'BackgroundColor', 'none');

figure('Position', [100 100 500 170])
subplot(1,2,1); hold on
plot(t, 1e3*S_M, 'LineWidth', 1.1); plot(t, 1e3*S_v, 'LineWidth', 1.1)
legend({'M \partial y/\partial M', 'v \partial y/\partial v'}, 'Location', 'best', 'FontSize', 7)
ylabel('Sensitivity (mm per 100 %)', 'FontName', fn, 'FontSize', fs); xlabel('Time (sec)', 'FontName', fn, 'FontSize', fs)
xlim([0 8]); set(gca, 'FontName', fn, 'FontSize', fs)
subplot(1,2,2)
semilogy(t, crlb_M, 'LineWidth', 1.1); hold on; semilogy(t, crlb_v, 'LineWidth', 1.1)
legend({'M', 'v'}, 'Location', 'best', 'FontSize', 7)
ylabel('Cramer-Rao bound, rel. std (%)', 'FontName', fn, 'FontSize', fs); xlabel('Time (sec)', 'FontName', fn, 'FontSize', fs)
xlim([0 15]); ylim([1e-2 1e2]); set(gca, 'FontName', fn, 'FontSize', fs)
exportgraphics(gcf, 'Velocity_sensitivity.pdf', 'ContentType', 'vector', 'BackgroundColor', 'none');
