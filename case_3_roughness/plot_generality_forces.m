% PLOT_GENERALITY_FORCES  Figures 9 (dyn_force_vs_statics.pdf) and 10 (two_cases.pdf) of the
% manuscript from the corrected-filter generality run (BMM_UKF_Generality_Forces, 4 parts).
function plot_generality_forces(mode)
if nargin < 1, mode = 'clean'; end
sfx = ''; if strcmp(mode,'error'), sfx = '_error'; end
close all
% ---- load the 4 parts of the generality run ----
R = [];
for p = 1:4, S = load(sprintf('generality_forces%s_part%d_results.mat', sfx, p)); R = [R, S.results]; end
% order as in the manuscript dataset: mass outer (4000:100:6000), speed inner (13:3:31)
[~, order] = sortrows([[R.M]', [R.v]']); R = R(order); nC = numel(R);
fprintf('%d cases; converged |dM|<=5%%: %d; mean dM %.2f%%, mean |dv| %.3f%%\n', nC, sum(abs([R.dM])<=5), mean([R.dM]), mean(abs([R.dv])));
Fdyn = [R.Fdyn_mean3]; Fstat = [R.Fstat_mean3]; DAF = Fdyn./Fstat;   % mean dynamic / static force over the last 3 s, all cases
fprintf('Fdyn/Fstat (last 3 s): mean %.4f  std %.4f  min %.4f  max %.4f\n', mean(DAF), std(DAF), min(DAF), max(DAF));
fn = 'Times New Roman'; fs = 13.4;   % figure font
% ---------------- Figure 9: representative case (5000 kg, 28 m/s) + all-case scatter
k = find([R.M]==5000 & [R.v]==28, 1); r = R(k);   % representative case for panel (a)
figure('Position',[100 100 850 300])
subplot(1,2,1)
i = r.t_full >= 2;   % skip the initial transient
fc_s = movmean(r.fc_est, 31);   % same 31-sample (0.12 s) moving average as the submitted figure
fprintf('smoothed force 2-15 s: min %.3e max %.3e N; static last 3 s %.3e N; max/static %.2f min/static %.2f\n', min(fc_s(i)), max(fc_s(i)), r.Fstat_mean3, max(fc_s(i))/r.Fstat_mean3, min(fc_s(i))/r.Fstat_mean3);
plot(r.t_full(i), fc_s(i), 'b', 'LineWidth', 1.0); hold on
plot(r.t_full(i), r.F_static(i), '--r', 'LineWidth', 1.6)
legend('dynamic force', 'static force', 'Location', 'northwest'); grid on; xlim([0 15]); ylim([0 10e4])
xlabel('Time (s)'); ylabel('Force (N)'); set(gca, 'FontName', fn, 'FontSize', fs)
% panel (b): mean dynamic force vs static weight, one dot per case, identity line
subplot(1,2,2)
scatter(Fstat, Fdyn, 40, 'filled'); hold on
plot([min(Fstat) max(Fstat)], [min(Fstat) max(Fstat)], 'k--', 'LineWidth', 1.5); grid on; axis equal
xlabel('Average static force (N)'); ylabel('Average dynamic force (N)'); set(gca, 'FontName', fn, 'FontSize', fs)
axs_ = findobj(gcf,'Type','axes'); [~,ord_] = sort(arrayfun(@(a) a.Position(1), axs_)); labs_ = {'(a)','(b)','(c)','(d)'};
for ii_ = 1:numel(axs_), text(axs_(ord_(ii_)), -0.12, -0.24, labs_{ii_}, 'Units','normalized','FontName','Times New Roman','FontSize',fs); end   % panel labels
exportgraphics(gcf, ['dyn_force_vs_statics' sfx '.pdf'], 'ContentType', 'vector', 'BackgroundColor', 'none');
% ---------------- Figure 10: two representative cases (41: 4500 kg / 28 m/s, 128: 5800 kg / 16 m/s)
ks = [find([R.M]==4500 & [R.v]==28, 1), find([R.M]==5800 & [R.v]==16, 1)];   % cases A and B of the paper
col = [0 0.447 0.741; 0.85 0.325 0.098]; lab = {'A','B'};   % case 41 -> A, case 128 -> B
figure('Position',[100 100 800 320])
% mass panel (a) and speed panel (b): 2-sigma band, estimate, truth (dashed)
for j = 1:2
    r = R(ks(j)); t = r.t;
    subplot(1,2,1); hold on
    fill([t fliplr(t)], [r.M_est+2*r.sig_M, fliplr(r.M_est-2*r.sig_M)], col(j,:), 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility', 'off')
    plot(t, r.M_est, 'Color', col(j,:), 'LineWidth', 1.6, 'DisplayName', sprintf('Case %s (M=%d kg)', lab{j}, r.M))
    plot(t, r.M*ones(size(t)), '--', 'Color', col(j,:), 'LineWidth', 1.4, 'HandleVisibility', 'off')
    subplot(1,2,2); hold on
    fill([t fliplr(t)], [r.v_est+2*r.sig_v, fliplr(r.v_est-2*r.sig_v)], col(j,:), 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility', 'off')
    plot(t, r.v_est, 'Color', col(j,:), 'LineWidth', 1.6, 'DisplayName', sprintf('Case %s (v=%.1f m/s)', lab{j}, r.v))
    plot(t, r.v*ones(size(t)), '--', 'Color', col(j,:), 'LineWidth', 1.4, 'HandleVisibility', 'off')
end
subplot(1,2,1); xlim([0 15]); ylim([2500 6000]); grid on; xlabel('Time (s)'); ylabel('Mass (kg)'); legend('Location', 'southeast'); set(gca, 'FontName', fn, 'FontSize', fs)
subplot(1,2,2); xlim([0 15]); ylim([0 40]); grid on; xlabel('Time (s)'); ylabel('Velocity (m/s)'); legend('Location', 'southeast'); set(gca, 'FontName', fn, 'FontSize', fs)
axs_ = findobj(gcf,'Type','axes'); [~,ord_] = sort(arrayfun(@(a) a.Position(1), axs_)); labs_ = {'(a)','(b)','(c)','(d)'};
for ii_ = 1:numel(axs_), text(axs_(ord_(ii_)), -0.12, -0.24, labs_{ii_}, 'Units','normalized','FontName','Times New Roman','FontSize',fs); end   % panel labels
exportgraphics(gcf, ['two_cases' sfx '.pdf'], 'ContentType', 'vector', 'BackgroundColor', 'none');
for j = 1:2, r = R(ks(j)); fprintf('Case %d: M=%d v=%g start %g dM %.2f%% dv %.2f%%\n', ks(j), r.M, r.v, r.start, r.dM, r.dv); end
disp('figures written')
% Section 6.1 style data (case number, estimated and true M and v) for a generality_cases counterpart
T = table((1:nC)', [R.M]', [R.v]', arrayfun(@(r) mean(r.M_est(r.t>=10)), R)', arrayfun(@(r) mean(r.v_est(r.t>=10)), R)', [R.dM]', [R.dv]', [R.start]', 'VariableNames', {'case','M_true','v_true','M_est','v_est','dM_pct','dv_pct','start'});
writetable(T, ['generality_cases' sfx '_from_forces.csv']);
end
