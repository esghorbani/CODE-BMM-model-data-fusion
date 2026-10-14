% BMM_UKF_ROBUSTNESS_COMBINED
% Robustness analysis for the new manuscript section.
% TRUTH : quarter-car (sprung mass + suspension, mu_w = 0.1) crossing a beam
%         with modal damping zeta_b on an ISO 8608 deck profile, Case-1 noise.
% FILTER: the ORIGINAL model of the paper -> rigid moving mass (mu_w = 1),
%         undamped beam, smooth deck, default tuning and 0.7 initialisation.
% Grid  : roughness {C, B, A, none} x beam damping {5, 2, 0} %.
% Also quantifies the model-form error: mid-span displacement of the filter
% model driven with the TRUE M, v minus the true displacement.
% Outputs: robustness_combined_results.mat, robustness_combined_log.txt
clc; close all; clear;
rng('default'); rng(1234);
diary('robustness_combined_log.txt'); diary on
par.n_modes = 4; par.m = 200; par.EI = 5.3e9; par.L = 60; par.g = 9.81;
M1 = 5000; v1 = 15; dt = 0.004; tspan = [0 30];
n = par.n_modes; phi_center = mode_shapesPP(par.L/2, n, par.L)';
amplitude_factor = 0.01; offset = 0.01; frequency = 5;
n_min = 0.01; n_max = 3.0; N_harm = 300;
f_v = 1.5; zeta_v = 0.20; mu_true = 0.10;

rough_list = {'C','B','A','none'};
zeta_list  = [0.05 0.02 0];
Q_default  = 10^-6.2;

% ---- filter model: original paper model (rigid mass, undamped, smooth) ----
par_est = par; par_est.mu_w = 1; par_est.k_v = 0; par_est.c_v = 0; par_est.prof = []; par_est.zeta = 0;

% ---- original-model displacement with TRUE parameters (for model error) ----
[t0, X0] = ode4(@(tt,X) BMM_SS_rough(tt, X, M1, v1, par_est), tspan, zeros(2*n+2,1), dt);
disp_model = phi_center * X0(:,1:n)';

results = struct([]); k = 0;
for ir = 1:numel(rough_list)
  for iz = 1:numel(zeta_list)
    k = k + 1; rc = rough_list{ir}; zb = zeta_list(iz);
    fprintf('\n=== roughness %s, beam damping %.0f %% ===\n', rc, 100*zb);
    par_true = par; par_true.mu_w = mu_true; par_true.zeta = zb;
    m_s = (1-mu_true)*M1; par_true.k_v = m_s*(2*pi*f_v)^2; par_true.c_v = 2*zeta_v*m_s*(2*pi*f_v);
    if strcmp(rc,'none'), par_true.prof = []; else, par_true.prof = iso8608_profile(rc, n_min, n_max, N_harm, 42); end

    [t, X] = ode4(@(tt,X) BMM_SS_rough(tt, X, M1, v1, par_true), tspan, zeros(2*n+2,1), dt);
    t = t(:)'; disp_true = phi_center * X(:,1:n)';
    fc_true = contact_force_rough(t, X', M1, v1, par_true);
    rng(1234); Nt = numel(t); tau = linspace(0,2*pi,Nt);
    nl = offset + amplitude_factor*abs(disp_true) + amplitude_factor*sin(frequency*tau);
    y = disp_true + nl.*randn(1,Nt);
    idx = t >= 10;

    % model-form error in the measurement
    e_model = disp_model - disp_true;
    rms_model = rms(e_model(idx)); rms_noise = rms(nl(idx).*1);  % noise std ~ nl
    fprintf('model error: RMS %.3f mm (var %.2e m^2); noise RMS %.3f mm; ratio %.2f; var/Q %.1f\n', ...
        1e3*rms_model, rms_model^2, 1e3*rms(y(idx)-disp_true(idx)), rms_model/rms(y(idx)-disp_true(idx)), rms_model^2/Q_default);

    tic; [X_est, P_saved, R0] = run_UKF_rough(t, y, M1, v1, par_est, dt); fprintf('AUKF %.0f s\n', toc);
    M_est = X_est(2*n+3,:); v_est = X_est(2*n+4,:);
    dM = 100*(mean(M_est(idx))-M1)/M1; dv = 100*(mean(v_est(idx))-v1)/v1;
    sM = 100*std(M_est(t>=20))/M1;
    fc_est = contact_force_rough(t, X_est(1:2*n+2,:), M_est, v_est, par_est);
    dfbar = 100*(mean(fc_est(idx))-mean(fc_true(idx)))/(M1*par.g);
    fprintf('M  : %.0f kg, error %.2f %% (std over t>=20 s: %.2f %%)\n', mean(M_est(idx)), dM, sM);
    fprintf('v  : %.3f m/s, error %.2f %%\n', mean(v_est(idx)), dv);
    fprintf('mean force error %.2f %% of Mg; converged R %.2e\n', dfbar, R0(end));

    results(k).rough = rc; results(k).zeta_b = zb; results(k).t = t;
    results(k).M_est = M_est; results(k).v_est = v_est; results(k).R0 = R0;
    results(k).dM = dM; results(k).dv = dv; results(k).sM = sM; results(k).dfbar = dfbar;
    results(k).rms_model = rms_model; results(k).rms_noise = rms(y(idx)-disp_true(idx));
    results(k).disp_true = disp_true; results(k).y = y; results(k).e_model = e_model;
  end
end
save('robustness_combined_results.mat', 'results', 'par', 'M1', 'v1', 'dt', '-v7.3');
fprintf('\nSUMMARY rough/zeta_b : dM %% | dv %% | model-error RMS mm | var/Q\n');
for k = 1:numel(results), r = results(k);
  fprintf('%-4s %4.0f%% : %6.2f | %6.2f | %6.3f | %6.1f\n', r.rough, 100*r.zeta_b, r.dM, r.dv, 1e3*r.rms_model, r.rms_model^2/Q_default); end
diary off
