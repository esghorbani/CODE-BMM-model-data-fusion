% BMM_UKF_ROBUSTNESS_DYNAMICS
% Model-form errors that change the DYNAMICS (not just the measurement):
% heavy beam damping, beam stiffness error, beam mass error. Truth = quarter-car
% on class C deck (as in the combined study) plus the listed error; filter =
% original model (rigid, undamped, smooth, nominal EI and m).
% Question: when does R plateau while the estimates degrade (i.e. when would Q
% need to adapt)?  Outputs: robustness_dynamics_results.mat / _log.txt
clc; close all; clear; rng('default'); rng(1234);
diary('robustness_dynamics_log.txt'); diary on
par.n_modes = 4; par.m = 200; par.EI = 5.3e9; par.L = 60; par.g = 9.81;
M1 = 5000; v1 = 15; dt = 0.004; tspan = [0 30];
n = par.n_modes; phi_center = mode_shapesPP(par.L/2, n, par.L)';
amplitude_factor = 0.01; offset = 0.01; frequency = 5;
n_min = 0.01; n_max = 3.0; N_harm = 300; f_v = 1.5; zeta_v = 0.20; mu_true = 0.10;
Q_default = 10^-6.2;
cases = struct('name',{},'zeta',{},'EIf',{},'mf',{});
for z = [0.05 0.10 0.20 0.40], cases(end+1) = struct('name',sprintf('zeta_b %.0f%%',100*z),'zeta',z,'EIf',1,'mf',1); end
for f = [0.8 0.9 1.1 1.2],     cases(end+1) = struct('name',sprintf('EI x%.1f',f),'zeta',0.05,'EIf',f,'mf',1); end
for f = [0.9 1.1],             cases(end+1) = struct('name',sprintf('m x%.1f',f),'zeta',0.05,'EIf',1,'mf',f); end
par_est = par; par_est.mu_w = 1; par_est.k_v = 0; par_est.c_v = 0; par_est.prof = []; par_est.zeta = 0;
[t0, X0] = ode4(@(tt,X) BMM_SS_rough(tt, X, M1, v1, par_est), tspan, zeros(2*n+2,1), dt);
disp_model = phi_center * X0(:,1:n)';
prof = iso8608_profile('C', n_min, n_max, N_harm, 42);
results = struct([]);
for k = 1:numel(cases)
    c = cases(k); fprintf('\n=== %s (class C deck, quarter-car) ===\n', c.name);
    par_true = par; par_true.mu_w = mu_true; par_true.zeta = c.zeta; par_true.EI = par.EI*c.EIf; par_true.m = par.m*c.mf;
    m_s = (1-mu_true)*M1; par_true.k_v = m_s*(2*pi*f_v)^2; par_true.c_v = 2*zeta_v*m_s*(2*pi*f_v); par_true.prof = prof;
    [t, X] = ode4(@(tt,X) BMM_SS_rough(tt, X, M1, v1, par_true), tspan, zeros(2*n+2,1), dt);
    t = t(:)'; disp_true = phi_center * X(:,1:n)';
    rng(1234); Nt = numel(t); tau = linspace(0,2*pi,Nt);
    nl = offset + amplitude_factor*abs(disp_true) + amplitude_factor*sin(frequency*tau);
    noise = nl.*randn(1,Nt); y = disp_true + noise; idx = t >= 10;
    e_model = disp_model - disp_true; var_model = var(e_model(idx));
    fprintf('model error var %.2e m^2 (RMS %.1f mm), /Q %.0f; max |disp| true %.1f mm\n', var_model, 1e3*rms(e_model(idx)), var_model/Q_default, 1e3*max(abs(disp_true)));
    tic; [X_est, P_saved, R0] = run_UKF_rough(t, y, M1, v1, par_est, dt); fprintf('AUKF %.0f s\n', toc);
    M_est = X_est(2*n+3,:); v_est = X_est(2*n+4,:);
    dM = 100*(mean(M_est(idx))-M1)/M1; dv = 100*(mean(v_est(idx))-v1)/v1; sM = 100*std(M_est(t>=20))/M1;
    fprintf('M error %.2f %% (std %.2f %%), v error %.2f %%, converged R %.2e\n', dM, sM, dv, R0(end));
    results(k).name = c.name; results(k).zeta = c.zeta; results(k).EIf = c.EIf; results(k).mf = c.mf;
    results(k).t = t; results(k).M_est = M_est; results(k).v_est = v_est; results(k).R0 = R0;
    results(k).dM = dM; results(k).dv = dv; results(k).sM = sM; results(k).var_model = var_model; results(k).var_noise = var(noise(idx));
end
save('robustness_dynamics_results.mat', 'results', 'M1', 'v1', 'dt', '-v7.3');
fprintf('\nSUMMARY case | model var | /Q | dM %% | dv %% | R_end\n');
for k = 1:numel(results), r = results(k); fprintf('%-12s | %.2e | %5.0f | %7.2f | %6.2f | %.2e\n', r.name, r.var_model, r.var_model/Q_default, r.dM, r.dv, r.R0(end)); end
diary off
