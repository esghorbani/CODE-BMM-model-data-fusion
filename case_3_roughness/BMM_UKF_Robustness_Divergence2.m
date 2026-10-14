% BMM_UKF_ROBUSTNESS_DIVERGENCE
% Continuation of BMM_UKF_Robustness_Combined: increase the roughness beyond
% class C (G_d(n0) = 256e-6 .. 8192e-6 m^3, i.e. classes C..F) with the
% quarter-car truth and 5 % beam damping, filter = original model, to find
% where the model-form error reaches the Q = 1e-4 level and where the
% estimates diverge.  Outputs: robustness_divergence2_results.mat / _log.txt
clc; close all; clear; rng('default'); rng(1234);
diary('robustness_divergence2_log.txt'); diary on
par.n_modes = 4; par.m = 200; par.EI = 5.3e9; par.L = 60; par.g = 9.81;   % beam: 4 modes, m [kg/m], EI [N m^2], L [m]
M1 = 5000; v1 = 15; dt = 0.004; tspan = [0 30];   % true mass [kg], speed [m/s], RK4 step [s], simulated time [s]
n = par.n_modes; phi_center = mode_shapesPP(par.L/2, n, par.L)';
amplitude_factor = 0.01; offset = 0.01; frequency = 5;   % Case-1 noise model: std = offset + 1 % of |disp| + sinusoidal term
n_min = 0.01; n_max = 3.0; N_harm = 300;   % ISO 8608 band [cycles/m] and number of harmonics
f_v = 1.5; zeta_v = 0.20; mu_true = 0.10; zb = 0.05;   % quarter-car: f_v [Hz], zeta_v, m_w/M; beam damping 5 %
Gd_list = [16384 32768 65536]*1e-6;  names = {'F','F/G','G'};   % roughness coefficient G_d(n0) [m^3], beyond the ISO classes
Q_default = 10^-6.2;
par_est = par; par_est.mu_w = 1; par_est.k_v = 0; par_est.c_v = 0; par_est.prof = []; par_est.zeta = 0;   % filter model: rigid mass, undamped, smooth
% filter model driven with the true M, v: reference for the model-form error
[t0, X0] = ode4(@(tt,X) BMM_SS_rough(tt, X, M1, v1, par_est), tspan, zeros(2*n+2,1), dt);
disp_model = phi_center * X0(:,1:n)';
results = struct([]);
for k = 1:numel(Gd_list)
    Gd = Gd_list(k); fprintf('\n=== G_d(n0) = %.0fe-6 m^3 (%s), beam damping %.0f %% ===\n', 1e6*Gd, names{k}, 100*zb);
    par_true = par; par_true.mu_w = mu_true; par_true.zeta = zb;   % data-generating model
    m_s = (1-mu_true)*M1; par_true.k_v = m_s*(2*pi*f_v)^2; par_true.c_v = 2*zeta_v*m_s*(2*pi*f_v);
    par_true.prof = iso8608_profile(Gd, n_min, n_max, N_harm, 42);   % numeric G_d(n0) instead of a class letter
    fprintf('profile RMS %.1f mm\n', 1e3*par_true.prof.sigma_r);
    % synthesise the truth and the noisy mid-span measurement
    [t, X] = ode4(@(tt,X) BMM_SS_rough(tt, X, M1, v1, par_true), tspan, zeros(2*n+2,1), dt);
    t = t(:)'; disp_true = phi_center * X(:,1:n)';
    rng(1234); Nt = numel(t); tau = linspace(0,2*pi,Nt);
    nl = offset + amplitude_factor*abs(disp_true) + amplitude_factor*sin(frequency*tau);
    noise = nl.*randn(1,Nt); y = disp_true + noise; idx = t >= 10;   % idx: converged window
    e_model = disp_model - disp_true; var_model = var(e_model(idx)); var_noise = var(noise(idx));   % model-form error vs measurement noise
    fprintf('model error var %.2e m^2 (RMS %.1f mm), /Q %.0f ; noise var %.2e m^2\n', var_model, 1e3*rms(e_model(idx)), var_model/Q_default, var_noise);
    % AUKF with the original filter model
    tic; [X_est, P_saved, R0] = run_UKF_rough(t, y, M1, v1, par_est, dt); fprintf('AUKF %.0f s\n', toc);
    M_est = X_est(2*n+3,:); v_est = X_est(2*n+4,:);
    dM = 100*(mean(M_est(idx))-M1)/M1; dv = 100*(mean(v_est(idx))-v1)/v1; sM = 100*std(M_est(t>=20))/M1;   % converged errors [%], mass scatter after 20 s
    fprintf('M error %.2f %% (std %.2f %%), v error %.2f %%, converged R %.2e\n', dM, sM, dv, R0(end));
    % store for the robustness_errors figure and the summary
    results(k).Gd = Gd; results(k).name = names{k}; results(k).sigma_r = par_true.prof.sigma_r;
    results(k).t = t; results(k).M_est = M_est; results(k).v_est = v_est; results(k).R0 = R0;
    results(k).dM = dM; results(k).dv = dv; results(k).sM = sM;
    results(k).var_model = var_model; results(k).var_noise = var_noise; results(k).e_model = e_model; results(k).noise = noise;
end
save('robustness_divergence2_results.mat', 'results', 'M1', 'v1', 'dt', '-v7.3');
fprintf('\nSUMMARY Gd(e-6) | RMS r mm | model var | /Q | dM %% | dv %% | R_end\n');
for k = 1:numel(results), r = results(k); fprintf('%5.0f | %5.1f | %.2e | %5.0f | %7.2f | %6.2f | %.2e\n', 1e6*r.Gd, 1e3*r.sigma_r, r.var_model, r.var_model/Q_default, r.dM, r.dv, r.R0(end)); end
diary off
