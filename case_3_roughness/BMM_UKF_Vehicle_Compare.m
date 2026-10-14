% BMM_UKF_VEHICLE_COMPARE
% Reviewer 1, comment C2: single moving mass vs. quarter-car (SDOF) vehicle.
% Smooth deck (no roughness) so that only the vehicle model differs.
%
% Cases (truth model -> filter model):
%   1  rigid mass      -> rigid mass         (original paper, mu_w = 1)
%   2  quarter-car     -> quarter-car        (mu_w = 0.1, suspension known)
%   3  pure SDOF       -> pure SDOF          (mu_w = 0, suspension known)
%   4  quarter-car     -> rigid mass         (paper's estimator applied to
%                                             data with suspension dynamics:
%                                             model-form error)
% Outputs: Vehicle_compare_params.pdf, Vehicle_compare_force.pdf,
%          vehicle_compare_results.mat, compare_log.txt

clc; close all; clear;
rng('default'); rng(1234);

par.n_modes = 4; par.m = 200; par.EI = 5.3e9; par.L = 60; par.g = 9.81;
par.prof = [];                           % smooth deck
M1 = 5000; v1 = 15; dt = 0.004; tspan = [0 30];
f_v = 1.5; zeta_v = 0.20;
n = par.n_modes;
phi_center = mode_shapesPP(par.L/2, n, par.L)';
amplitude_factor = 0.01; offset = 0.01; frequency = 5;

cases = struct( ...
    'name',   {'Rigid mass', 'Quarter-car', 'Pure SDOF', 'Quarter-car data, rigid-mass filter'}, ...
    'mu_true',{1, 0.1, 0, 0.1}, ...
    'mu_filt',{1, 0.1, 0, 1});

results = struct([]);
for c = 1:numel(cases)
    fprintf('\n=== %s ===\n', cases(c).name);

    % suspension properties defined from the TRUE sprung mass
    m_s = (1 - cases(c).mu_true) * M1;
    par_true = par;  par_true.mu_w = cases(c).mu_true;
    par_true.k_v = m_s*(2*pi*f_v)^2;  par_true.c_v = 2*zeta_v*m_s*(2*pi*f_v);

    [t, X] = ode4(@(tt,X) BMM_SS_rough(tt, X, M1, v1, par_true), tspan, zeros(2*n+2,1), dt);
    t = t(:)';
    disp_true = phi_center * X(:,1:n)';
    fc_true   = contact_force_rough(t, X', M1, v1, par_true);

    rng(1234);
    Nt = numel(t);  tau = linspace(0, 2*pi, Nt);
    noise_level = offset + amplitude_factor*abs(disp_true) + amplitude_factor*sin(frequency*tau);
    y = disp_true + noise_level .* randn(1, Nt);

    par_est = par_true;  par_est.mu_w = cases(c).mu_filt;
    tic
    [X_est, P_saved] = run_UKF_rough(t, y, M1, v1, par_est, dt);
    fprintf('AUKF time %.1f s\n', toc);

    M_est = X_est(2*n+3,:);  v_est = X_est(2*n+4,:);
    sig_M = sqrt(squeeze(P_saved(2*n+3,2*n+3,:)))';
    sig_v = sqrt(squeeze(P_saved(2*n+4,2*n+4,:)))';
    fc_est = contact_force_rough(t, X_est(1:2*n+2,:), M_est, v_est, par_est);

    idx = t >= 10;
    fprintf('M  : mean %.0f kg, rel. error %.2f %%\n', mean(M_est(idx)), 100*(mean(M_est(idx))-M1)/M1);
    fprintf('v  : mean %.3f m/s, rel. error %.2f %%\n', mean(v_est(idx)), 100*(mean(v_est(idx))-v1)/v1);
    fprintf('true f_c: dynamic amplification (max/Mg) %.3f, std %.1f %% of Mg\n', ...
        max(fc_true(idx))/(M1*par.g), 100*std(fc_true(idx))/(M1*par.g));
    fprintf('f_c: RMS error %.1f %% of Mg; mean-force error %.2f %%\n', ...
        100*rms(fc_est(idx)-fc_true(idx))/(M1*par.g), 100*(mean(fc_est(idx))-mean(fc_true(idx)))/(M1*par.g));
    fprintf('midspan disp: max %.4f m\n', max(disp_true(idx)));

    results(c).name = cases(c).name; results(c).t = t; results(c).disp_true = disp_true;
    results(c).M_est = M_est; results(c).v_est = v_est; results(c).sig_M = sig_M; results(c).sig_v = sig_v;
    results(c).fc_true = fc_true; results(c).fc_est = fc_est; results(c).X_est = X_est;
end
save('vehicle_compare_results.mat', 'results', 'par', 'M1', 'v1', 'dt', '-v7.3');

%% figures
fs = 10; fn = 'Times New Roman'; cols = lines(4);
figure('Position', [100 100 500 170])
subplot(1,2,1); hold on
for c = 1:4, plot(t, results(c).M_est, 'Color', cols(c,:), 'LineWidth', 1.2); end
plot(t, M1*ones(size(t)), 'k--', 'LineWidth', 1.2)
legend({'Rigid mass', 'Quarter-car', 'Pure SDOF', 'QC data / RM filter', 'Ground Truth'}, 'Location', 'southeast', 'FontSize', 7)
ylabel('Moving mass (kg)', 'FontName', fn, 'FontSize', fs); xlabel('Time (sec)', 'FontName', fn, 'FontSize', fs)
ylim([0.5*M1 1.2*M1]); xlim([0 15]); set(gca, 'FontName', fn, 'FontSize', fs)
subplot(1,2,2); hold on
for c = 1:4, plot(t, results(c).v_est, 'Color', cols(c,:), 'LineWidth', 1.2); end
plot(t, v1*ones(size(t)), 'k--', 'LineWidth', 1.2)
ylabel('Velocity (m/s)', 'FontName', fn, 'FontSize', fs); xlabel('Time (sec)', 'FontName', fn, 'FontSize', fs)
ylim([0 1.4*v1]); xlim([0 15]); set(gca, 'FontName', fn, 'FontSize', fs)
exportgraphics(gcf, 'Vehicle_compare_params.pdf', 'ContentType', 'vector', 'BackgroundColor', 'none');

figure('Position', [100 100 500 170])
subplot(1,2,1); hold on   % true contact force, one passage
for c = 1:3, plot(t, 1e-3*results(c).fc_true, 'Color', cols(c,:), 'LineWidth', 1.1); end
plot(t, 1e-3*M1*par.g*ones(size(t)), 'k--', 'LineWidth', 1)
legend({'Rigid mass', 'Quarter-car', 'Pure SDOF', 'M g'}, 'Location', 'best', 'FontSize', 7)
ylabel('Contact force (kN)', 'FontName', fn, 'FontSize', fs); xlabel('Time (sec)', 'FontName', fn, 'FontSize', fs)
xlim([12 16]); set(gca, 'FontName', fn, 'FontSize', fs)
subplot(1,2,2); hold on   % reconstructed vs true, quarter-car
plot(t, 1e-3*results(2).fc_est, 'Color', cols(2,:), 'LineWidth', 1.2)
plot(t, 1e-3*results(2).fc_true, 'k--', 'LineWidth', 1)
legend({'UKF reconstructed', 'Ground Truth'}, 'Location', 'best', 'FontSize', 7)
ylabel('Contact force (kN)', 'FontName', fn, 'FontSize', fs); xlabel('Time (sec)', 'FontName', fn, 'FontSize', fs)
xlim([12 16]); set(gca, 'FontName', fn, 'FontSize', fs)
exportgraphics(gcf, 'Vehicle_compare_force.pdf', 'ContentType', 'vector', 'BackgroundColor', 'none');
