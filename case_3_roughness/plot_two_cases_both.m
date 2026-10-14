% PLOT_TWO_CASES_BOTH  Convergence of cases A (4500 kg, 28 m/s) and B (5800 kg, 16 m/s) of the
% 147-case study, without and with model error, in one mass panel and one speed panel
% (two_cases_both.pdf). For each case both scenarios start from the same speed, the value selected
% by the filter bank in the model-error run (26 m/s for A, 18 m/s for B); the clean case is re-run
% from that start with the Section 5 model and the time-varying noise of the generality study.
close all
% ---- load the clean and model-error runs (4 parts each) ----
load_all = @(sfx) cell2mat(arrayfun(@(p) getfield(load(sprintf('generality_forces%s_part%d_results.mat', sfx, p)), 'results'), 1:4, 'UniformOutput', false));
Rc = load_all(''); Re = load_all('_error');
[~, o] = sortrows([[Rc.M]', [Rc.v]']); Rc = Rc(o); [~, o] = sortrows([[Re.M]', [Re.v]']); Re = Re(o);   % mass outer, speed inner
ks = [find([Rc.M]==4500 & [Rc.v]==28, 1), find([Rc.M]==5800 & [Rc.v]==16, 1)]; lab = {'A','B'};   % cases A and B
col = [0 0.447 0.741; 0.85 0.325 0.098]; fn = 'Times New Roman'; fs = 13.4;
% re-run the clean cases from the model-error bank start
par.n_modes = 4; par.m = 200; par.EI = 5.3e9; par.L = 60; par.g = 9.81;   % Section 2 model: truth and filter model of the clean case
par.mu_w = 1; par.k_v = 0; par.c_v = 0; par.prof = []; par.zeta = 0;
dt = 0.004; tspan = [0 30]; n = par.n_modes; phi_center = mode_shapesPP(par.L/2, n, par.L)';
amplitude_factor = 0.01; offset = 0.01; frequency = 5; dec = 5;
for j = 1:2
    M1 = Rc(ks(j)).M; v1 = Rc(ks(j)).v; v0 = Re(ks(j)).start;   % v0: start speed chosen by the bank in the model-error run
    % truth, noisy measurement, filter started at v0
    [t, X] = ode4(@(tt,X) BMM_SS_rough(tt, X, M1, v1, par), tspan, zeros(2*n+2,1), dt);
    t = t(:)'; d = phi_center*X(:,1:n)'; rng(1234); Nt = numel(t); tau = linspace(0,2*pi,Nt);
    y = d + (offset + amplitude_factor*abs(d) + amplitude_factor*sin(frequency*tau)).*randn(1,Nt);
    [X_est, P_saved] = run_UKF_rough(t, y, M1, v1, par, dt, struct('v_init', v0, 'Q', 1e-6));
    % replace the stored clean entry by the re-run (decimated like the stored histories)
    r = Rc(ks(j)); r.start = v0; r.M_est = X_est(2*n+3,1:dec:end); r.v_est = X_est(2*n+4,1:dec:end);
    r.sig_M = sqrt(squeeze(P_saved(2*n+3,2*n+3,1:dec:end)))'; r.sig_v = sqrt(squeeze(P_saved(2*n+4,2*n+4,1:dec:end)))';
    idx = t >= 10; r.dM = 100*(mean(X_est(2*n+3,idx))-M1)/M1; r.dv = 100*(mean(X_est(2*n+4,idx))-v1)/v1;
    Rc(ks(j)) = r;
end
% ---- figure: solid = no model error, dashed = model error, dotted = truth ----
figure('Position',[100 100 850 330])
for j = 1:2
  for s = 1:2
    if s == 1, r = Rc(ks(j)); ls = '-'; tag = 'no model error'; else, r = Re(ks(j)); ls = '--'; tag = 'model error'; end
    t = r.t;
    fprintf('case %s, %s: M=%d v=%g start %g dM %.2f%% dv %.2f%%\n', lab{j}, tag, r.M, r.v, r.start, r.dM, r.dv);
    subplot(1,2,1); hold on
    fill([t fliplr(t)], [r.M_est+2*r.sig_M, fliplr(r.M_est-2*r.sig_M)], col(j,:), 'FaceAlpha', 0.12, 'EdgeColor', 'none', 'HandleVisibility', 'off')
    plot(t, r.M_est, ls, 'Color', col(j,:), 'LineWidth', 1.5, 'DisplayName', sprintf('Case %s, %s', lab{j}, tag))
    subplot(1,2,2); hold on
    fill([t fliplr(t)], [r.v_est+2*r.sig_v, fliplr(r.v_est-2*r.sig_v)], col(j,:), 'FaceAlpha', 0.12, 'EdgeColor', 'none', 'HandleVisibility', 'off')
    plot(t, r.v_est, ls, 'Color', col(j,:), 'LineWidth', 1.5, 'DisplayName', sprintf('Case %s, %s', lab{j}, tag))
  end
  subplot(1,2,1); plot([0 15], [1 1]*Rc(ks(j)).M, ':', 'Color', 'k', 'LineWidth', 2.0, 'HandleVisibility', 'off')
  subplot(1,2,2); plot([0 15], [1 1]*Rc(ks(j)).v, ':', 'Color', 'k', 'LineWidth', 2.0, 'HandleVisibility', 'off')
end
subplot(1,2,1); xlim([0 15]); ylim([2500 6000]); grid on; xlabel('Time (s)'); ylabel('Mass (kg)'); set(gca, 'FontName', fn, 'FontSize', fs); legend('Location', 'southeast', 'FontSize', 9.5)
subplot(1,2,2); xlim([0 15]); ylim([0 40]); grid on; xlabel('Time (s)'); ylabel('Velocity (m/s)'); set(gca, 'FontName', fn, 'FontSize', fs); legend('Location', 'southeast', 'FontSize', 9.5)
axs_ = findobj(gcf,'Type','axes'); [~,ord_] = sort(arrayfun(@(a) a.Position(1), axs_)); labs_ = {'(a)','(b)'};
for ii_ = 1:numel(axs_), text(axs_(ord_(ii_)), -0.12, -0.24, labs_{ii_}, 'Units','normalized','FontName',fn,'FontSize',fs); end
exportgraphics(gcf, 'two_cases_both.pdf', 'ContentType', 'vector', 'BackgroundColor', 'none');
disp('two_cases_both written')
