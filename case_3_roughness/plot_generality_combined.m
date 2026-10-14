% PLOT_GENERALITY_COMBINED  Clean (Section 5 model) and model-error (quarter-car + class C + 5 % damping,
% filter uninformed) results of the 147-case study in one figure each:
%   dyn_force_vs_statics_combined.pdf, two_cases_combined.pdf
close all
load_all = @(sfx) cell2mat(arrayfun(@(p) getfield(load(sprintf('generality_forces%s_part%d_results.mat', sfx, p)), 'results'), 1:4, 'UniformOutput', false));
Rc = load_all(''); Re = load_all('_error');
[~, o] = sortrows([[Rc.M]', [Rc.v]']); Rc = Rc(o); [~, o] = sortrows([[Re.M]', [Re.v]']); Re = Re(o);
fn = 'Times New Roman'; fs = 10; cb = [0 0.447 0.741]; co = [0.85 0.325 0.098];
% ---------------- Figure 9 combined
k = find([Rc.M]==5000 & [Rc.v]==28, 1);
figure('Position',[100 100 900 300])
subplot(1,2,1); hold on
i = Rc(k).t_full >= 2;
plot(Rc(k).t_full(i), subsref(movmean(Rc(k).fc_est,31), struct('type','()','subs',{{i}})), 'Color', cb, 'LineWidth', 0.9)
plot(Re(k).t_full(i), subsref(movmean(Re(k).fc_est,31), struct('type','()','subs',{{i}})), 'Color', co, 'LineWidth', 0.9)
plot(Rc(k).t_full(i), Rc(k).F_static(i), '--', 'Color', cb, 'LineWidth', 1.6)
plot(Re(k).t_full(i), Re(k).F_static(i), '--', 'Color', co, 'LineWidth', 1.6)
legend({'dynamic force, no model error','dynamic force, model error','static force, no model error','static force, model error'}, 'Location', 'northwest', 'FontSize', 8)
grid on; xlim([0 15]); ylim([0 10e4]); xlabel('Time (s)'); ylabel('Force (N)'); set(gca, 'FontName', fn, 'FontSize', fs)
subplot(1,2,2); hold on
scatter([Rc.Fstat_mean3], [Rc.Fdyn_mean3], 36, cb, 'filled')
scatter([Re.Fstat_mean3], [Re.Fdyn_mean3], 36, co, '^', 'filled')
lim = [min([Rc.Fstat_mean3, Re.Fstat_mean3]) max([Rc.Fstat_mean3, Re.Fstat_mean3])];
plot(lim, lim, 'k--', 'LineWidth', 1.5)
legend({'no model error','model error','1:1'}, 'Location', 'northwest', 'FontSize', 8)
grid on; axis equal; xlabel('Average static force (N)'); ylabel('Average dynamic force (N)'); set(gca, 'FontName', fn, 'FontSize', fs)
exportgraphics(gcf, 'dyn_force_vs_statics_combined.pdf', 'ContentType', 'vector', 'BackgroundColor', 'none');
% ---------------- Figure 10 combined: two cases, clean solid / model error dash-dot, bands for both
ks = [find([Rc.M]==4500 & [Rc.v]==28, 1), find([Rc.M]==5800 & [Rc.v]==16, 1)]; col = [cb; co];
% Re-run the two clean cases with the speed start selected by the bank in the model-error run,
% so that both curves of a case begin from the same initial condition.
par.n_modes = 4; par.m = 200; par.EI = 5.3e9; par.L = 60; par.g = 9.81;
par.mu_w = 1; par.k_v = 0; par.c_v = 0; par.prof = []; par.zeta = 0;
dt = 0.004; tspan = [0 30]; n = par.n_modes; phi_center = mode_shapesPP(par.L/2, n, par.L)';
amplitude_factor = 0.01; offset = 0.01; frequency = 5; dec = 5;
for j = 1:2
    M1 = Rc(ks(j)).M; v1 = Rc(ks(j)).v; v0 = Re(ks(j)).start;
    [t, X] = ode4(@(tt,X) BMM_SS_rough(tt, X, M1, v1, par), tspan, zeros(2*n+2,1), dt);
    t = t(:)'; d = phi_center*X(:,1:n)'; rng(1234); Nt = numel(t); tau = linspace(0,2*pi,Nt);
    y = d + (offset + amplitude_factor*abs(d) + amplitude_factor*sin(frequency*tau)).*randn(1,Nt);
    [X_est, P_saved] = run_UKF_rough(t, y, M1, v1, par, dt, struct('v_init', v0, 'Q', 1e-6));
    r = Rc(ks(j)); r.start = v0; r.M_est = X_est(2*n+3,1:dec:end); r.v_est = X_est(2*n+4,1:dec:end);
    r.sig_M = sqrt(squeeze(P_saved(2*n+3,2*n+3,1:dec:end)))'; r.sig_v = sqrt(squeeze(P_saved(2*n+4,2*n+4,1:dec:end)))';
    idx = t >= 10; r.dM = 100*(mean(X_est(2*n+3,idx))-M1)/M1; r.dv = 100*(mean(X_est(2*n+4,idx))-v1)/v1;
    fprintf('clean case %d rerun with start %g m/s: dM %.2f%%  dv %.2f%%\n', ks(j), v0, r.dM, r.dv);
    Rc(ks(j)) = r;
end
figure('Position',[100 100 850 330])
for j = 1:2
    for s = 1:2
        if s == 1, r = Rc(ks(j)); ls = '-'; tag = 'no model error'; else, r = Re(ks(j)); ls = '-.'; tag = 'model error'; end
        t = r.t;
        subplot(1,2,1); hold on
        fill([t fliplr(t)], [r.M_est+2*r.sig_M, fliplr(r.M_est-2*r.sig_M)], col(j,:), 'FaceAlpha', 0.12, 'EdgeColor', 'none', 'HandleVisibility', 'off')
        plot(t, r.M_est, ls, 'Color', col(j,:), 'LineWidth', 1.5, 'DisplayName', sprintf('Case %d (M=%d kg), %s', ks(j), r.M, tag))
        subplot(1,2,2); hold on
        fill([t fliplr(t)], [r.v_est+2*r.sig_v, fliplr(r.v_est-2*r.sig_v)], col(j,:), 'FaceAlpha', 0.12, 'EdgeColor', 'none', 'HandleVisibility', 'off')
        plot(t, r.v_est, ls, 'Color', col(j,:), 'LineWidth', 1.5, 'DisplayName', sprintf('Case %d (v=%.0f m/s), %s', ks(j), r.v, tag))
    end
    subplot(1,2,1); plot([0 15], [1 1]*Rc(ks(j)).M, ':', 'Color', col(j,:), 'LineWidth', 1.3, 'HandleVisibility', 'off')
    subplot(1,2,2); plot([0 15], [1 1]*Rc(ks(j)).v, ':', 'Color', col(j,:), 'LineWidth', 1.3, 'HandleVisibility', 'off')
end
subplot(1,2,1); xlim([0 15]); ylim([2500 6000]); grid on; xlabel('Time (s)'); ylabel('Mass (kg)'); legend('Location', 'southeast', 'FontSize', 7); set(gca, 'FontName', fn, 'FontSize', fs)
subplot(1,2,2); xlim([0 15]); ylim([0 40]); grid on; xlabel('Time (s)'); ylabel('Velocity (m/s)'); legend('Location', 'southeast', 'FontSize', 7); set(gca, 'FontName', fn, 'FontSize', fs)
exportgraphics(gcf, 'two_cases_combined.pdf', 'ContentType', 'vector', 'BackgroundColor', 'none');
disp('combined figures written')
