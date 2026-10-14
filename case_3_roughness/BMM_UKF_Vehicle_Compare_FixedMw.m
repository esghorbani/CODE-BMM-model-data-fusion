% BMM_UKF_VEHICLE_COMPARE_FIXEDMW
% Same four cases as BMM_UKF_Vehicle_Compare.m, but the filter treats the
% wheel (unsprung) mass m_w, the spring k_v and the dashpot c_v as KNOWN
% vehicle properties (fixed numbers), and estimates the sprung mass m_s
% through the total mass state M = m_w + m_s.  Truth data identical to the
% ratio-based run.  Outputs: vehicle_compare_fixedmw_results.mat,
% compare_fixedmw_log.txt
clc; close all; clear;
rng('default'); rng(1234);
diary('compare_fixedmw_log.txt'); diary on
par.n_modes = 4; par.m = 200; par.EI = 5.3e9; par.L = 60; par.g = 9.81;
par.prof = [];
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
    fprintf('\n=== %s (fixed m_w in filter) ===\n', cases(c).name);
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
    if cases(c).mu_filt < 1          % oscillator in the filter: wheel mass is a known number
        par_est.m_w_fixed = cases(c).mu_filt * M1;   % 500 kg (quarter-car) or 0 (SDOF)
    end                              % rigid-mass filter: m_w = M, no sprung mass
    fprintf('filter: m_w fixed = %s, k_v = %.3e, c_v = %.3e\n', ...
        mat2str(cases(c).mu_filt*M1*(cases(c).mu_filt<1)), par_est.k_v, par_est.c_v);
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
    fprintf('f_c: RMS error %.1f %% of Mg; mean-force error %.2f %%\n', ...
        100*rms(fc_est(idx)-fc_true(idx))/(M1*par.g), 100*(mean(fc_est(idx))-mean(fc_true(idx)))/(M1*par.g));

    results(c).name = cases(c).name; results(c).t = t;
    results(c).M_est = M_est; results(c).v_est = v_est; results(c).sig_M = sig_M; results(c).sig_v = sig_v;
    results(c).fc_true = fc_true; results(c).fc_est = fc_est;
    results(c).dM = 100*(mean(M_est(idx))-M1)/M1; results(c).dv = 100*(mean(v_est(idx))-v1)/v1;
end
save('vehicle_compare_fixedmw_results.mat', 'results', 'par', 'M1', 'v1', 'dt', '-v7.3');
fprintf('\nSUMMARY (fixed m_w): '); for c=1:4, fprintf('%s dM=%.2f%% dv=%.2f%% | ', results(c).name, results(c).dM, results(c).dv); end; fprintf('\n');
diary off
