% BMM_UKF_INNOVATION_TEST  Combined 12-case robustness study (quarter-car truth,
% roughness none/A/B/C x damping 0/2/5 %, original filter model) with the
% INNOVATION form of the adaptive R instead of the residual (posterior) form.
clc; clear; rng('default'); rng(1234); diary('innovation_test_log.txt'); diary on
par.n_modes = 4; par.m = 200; par.EI = 5.3e9; par.L = 60; par.g = 9.81;
M1 = 5000; v1 = 15; dt = 0.004; tspan = [0 30]; n = par.n_modes; phi_center = mode_shapesPP(par.L/2, n, par.L)';
amplitude_factor = 0.01; offset = 0.01; frequency = 5; f_v = 1.5; zeta_v = 0.20; mu_true = 0.10;
par_est = par; par_est.mu_w = 1; par_est.k_v = 0; par_est.c_v = 0; par_est.prof = []; par_est.zeta = 0;
fprintf('%-6s %4s | %8s %8s %10s | %8s %8s %10s\n','rough','zeta','res dM','res dv','res R','inn dM','inn dv','inn R');
for rc = {'none','A','B','C'}
  for zb = [0 0.02 0.05]
    par_true = par; par_true.mu_w = mu_true; par_true.zeta = zb;
    m_s = (1-mu_true)*M1; par_true.k_v = m_s*(2*pi*f_v)^2; par_true.c_v = 2*zeta_v*m_s*(2*pi*f_v);
    if strcmp(rc{1},'none'), par_true.prof = []; else, par_true.prof = iso8608_profile(rc{1}, 0.01, 3, 300, 42); end
    [t, X] = ode4(@(tt,X) BMM_SS_rough(tt, X, M1, v1, par_true), tspan, zeros(2*n+2,1), dt); t = t(:)';
    d = phi_center*X(:,1:n)'; rng(1234); Nt = numel(t); tau = linspace(0,2*pi,Nt);
    y = d + (offset + amplitude_factor*abs(d) + amplitude_factor*sin(frequency*tau)).*randn(1,Nt); idx = t >= 10;
    out = zeros(2,3); forms = {'posterior','innovation'};
    for f = 1:2
        [Xe,~,R0] = run_UKF_rough(t, y, M1, v1, par_est, dt, struct('R_covariance', forms{f}));
        out(f,:) = [100*(mean(Xe(2*n+3,idx))-M1)/M1, 100*(mean(Xe(2*n+4,idx))-v1)/v1, mean(R0(idx))];
    end
    fprintf('%-6s %4.0f | %8.2f %8.2f %10.2e | %8.2f %8.2f %10.2e\n', rc{1}, 100*zb, out(1,:), out(2,:));
  end
end
fprintf('DONE\n'); diary off
