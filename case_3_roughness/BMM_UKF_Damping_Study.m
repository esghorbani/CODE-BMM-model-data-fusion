% BMM_UKF_DAMPING_STUDY
% Reviewer 1, comments C3 and C5: sensitivity to structural damping.
%
% Truth: quarter-car (mu_w = 0.1) on a smooth deck, beam with modal damping
%        ratio zeta_true in {0, 1, 2, 5} %  (same ratio for all modes).
% Filter: (a) undamped model, as in the paper  -> neglected damping
%         (b) damping ratio known              -> reference
% Reports mass / velocity errors, mean reconstructed force error and the
% converged adaptive measurement-noise variance R, to see whether the
% neglected damping is absorbed by M, v or by R.
% Outputs: Damping_params.pdf, Damping_R.pdf, damping_results.mat, damping_log.txt

clc; close all; clear;
rng('default'); rng(1234);

par.n_modes = 4; par.m = 200; par.EI = 5.3e9; par.L = 60; par.g = 9.81;
par.prof = []; par.mu_w = 0.1;
M1 = 5000; v1 = 15; dt = 0.004; tspan = [0 30];
f_v = 1.5; zeta_v = 0.20; m_s = (1-par.mu_w)*M1;
par.k_v = m_s*(2*pi*f_v)^2; par.c_v = 2*zeta_v*m_s*(2*pi*f_v);
n = par.n_modes; phi_center = mode_shapesPP(par.L/2, n, par.L)';
amplitude_factor = 0.01; offset = 0.01; frequency = 5;

zetas = [0 0.01 0.02 0.05];
omega1 = sqrt(par.EI/par.m)*(pi/par.L)^2;
fprintf('first beam frequency %.2f Hz\n', omega1/2/pi);

results = struct([]);
for c = 1:numel(zetas)
    z = zetas(c);
    par_true = par; par_true.zeta = z;
    [t, X] = ode4(@(tt,X) BMM_SS_rough(tt, X, M1, v1, par_true), tspan, zeros(2*n+2,1), dt);
    t = t(:)'; disp_true = phi_center*X(:,1:n)';
    fc_true = contact_force_rough(t, X', M1, v1, par_true);
    rng(1234); Nt = numel(t); tau = linspace(0,2*pi,Nt);
    nl = offset + amplitude_factor*abs(disp_true) + amplitude_factor*sin(frequency*tau);
    y = disp_true + nl.*randn(1,Nt);
    idx = t >= 10;

    for known = [false true]
        par_est = par; par_est.zeta = z*known;
        [X_est, P_saved, R0] = run_UKF_rough(t, y, M1, v1, par_est, dt);
        M_est = X_est(2*n+3,:); v_est = X_est(2*n+4,:);
        fc_est = contact_force_rough(t, X_est(1:2*n+2,:), M_est, v_est, par_est);
        tag = 'neglected'; if known, tag = 'known    '; end
        fprintf('zeta_true = %4.1f %%  filter damping %s : dM %+.2f %%  dv %+.3f %%  mean f_c err %+.2f %%  R(t>10) %.2e  max disp %.1f mm\n', ...
            100*z, tag, 100*(mean(M_est(idx))-M1)/M1, 100*(mean(v_est(idx))-v1)/v1, ...
            100*(mean(fc_est(idx))-mean(fc_true(idx)))/(M1*par.g), mean(R0(idx)), 1e3*max(disp_true(idx)));
        r.zeta = z; r.known = known; r.t = t; r.M_est = M_est; r.v_est = v_est; r.R0 = R0;
        r.disp_true = disp_true; r.fc_true = fc_true; r.fc_est = fc_est;
        if isempty(results), results = r; else, results(end+1) = r; end
    end
end
save('damping_results.mat', 'results', 'par', 'M1', 'v1', 'dt', '-v7.3');

%% figures
fs = 10; fn = 'Times New Roman'; cols = lines(numel(zetas));
neg = results(~[results.known]);
figure('Position', [100 100 500 170])
subplot(1,2,1); hold on
for c = 1:numel(neg), plot(t, neg(c).M_est, 'Color', cols(c,:), 'LineWidth', 1.2); end
plot(t, M1*ones(size(t)), 'k--', 'LineWidth', 1.2)
legend([arrayfun(@(z) sprintf('\\zeta = %g%%', 100*z), zetas, 'UniformOutput', false), {'Ground Truth'}], 'Location', 'southeast', 'FontSize', 7)
ylabel('Moving mass (kg)', 'FontName', fn, 'FontSize', fs); xlabel('Time (sec)', 'FontName', fn, 'FontSize', fs)
ylim([0.5*M1 1.2*M1]); xlim([0 15]); set(gca, 'FontName', fn, 'FontSize', fs)
subplot(1,2,2); hold on
for c = 1:numel(neg), plot(t, neg(c).v_est, 'Color', cols(c,:), 'LineWidth', 1.2); end
plot(t, v1*ones(size(t)), 'k--', 'LineWidth', 1.2)
ylabel('Velocity (m/s)', 'FontName', fn, 'FontSize', fs); xlabel('Time (sec)', 'FontName', fn, 'FontSize', fs)
ylim([0 1.4*v1]); xlim([0 15]); set(gca, 'FontName', fn, 'FontSize', fs)
exportgraphics(gcf, 'Damping_params.pdf', 'ContentType', 'vector', 'BackgroundColor', 'none');

figure('Position', [100 100 500 170])
subplot(1,2,1); hold on
for c = 1:numel(neg), plot(t, neg(c).disp_true, 'Color', cols(c,:), 'LineWidth', 1); end
legend(arrayfun(@(z) sprintf('\\zeta = %g%%', 100*z), zetas, 'UniformOutput', false), 'Location', 'best', 'FontSize', 7)
ylabel('Mid-span displacement (m)', 'FontName', fn, 'FontSize', fs); xlabel('Time (sec)', 'FontName', fn, 'FontSize', fs)
xlim([12 16]); set(gca, 'FontName', fn, 'FontSize', fs)
subplot(1,2,2); hold on
for c = 1:numel(neg), plot(t, sqrt(neg(c).R0), 'Color', cols(c,:), 'LineWidth', 1); end
ylabel('$\sqrt{R}$ (m)', 'Interpreter', 'latex', 'FontSize', fs); xlabel('Time (sec)', 'FontName', fn, 'FontSize', fs)
xlim([5 20]); set(gca, 'FontName', fn, 'FontSize', fs)
exportgraphics(gcf, 'Damping_R.pdf', 'ContentType', 'vector', 'BackgroundColor', 'none');
