% BMM_UKF_VELOCITYSTATE  Reviewer 2, C5: quality of the mid-span VELOCITY STATE
% estimate when only the displacement is measured.  Rigid mass, v = 15 m/s,
% Case 1 noise model scaled by a factor (1 = manuscript level, ~10-20 mm).
% Reports displacement/velocity state errors, bias, and the error after a
% 5-Hz low-pass; exports CSV for the letter figure.
clc; clear;
par.n_modes = 4; par.m = 200; par.EI = 5.3e9; par.L = 60; par.g = 9.81;
par.prof = []; par.mu_w = 1; par.k_v = 0; par.c_v = 0;
M1 = 5000; v1 = 15; dt = 0.004; tspan = [0 30]; n = 4;
phi_center = mode_shapesPP(par.L/2, n, par.L)';
F = '../../PAPER-2026-BMM-FUSION/JSV_Submission/rebuttal/2026-09-11-round-1/figures/';
[t, X] = ode4(@(tt,X) BMM_SS_rough(tt, X, M1, v1, par), tspan, zeros(2*n+2,1), dt); t = t(:)';
d_true = phi_center*X(:,1:n)'; vel_true = phi_center*X(:,n+1:2*n)';
idx = t >= 10; Nt = numel(t); tau = linspace(0,2*pi,Nt);
fprintf('true mid-span velocity: RMS %.4f m/s, peak %.3f m/s\n', rms(vel_true(idx)), max(abs(vel_true(idx))));
scales = [1 0.3 0.1]; out = [];
for s = scales
    rng(1234); nl = s*(0.01 + 0.01*abs(d_true) + 0.01*sin(5*tau)); y = d_true + nl.*randn(1,Nt);
    v0 = initial_speed_estimate(t, y, par.L);
    Xe = run_UKF_rough(t, y, M1, v1, par, dt, struct('init',[0.7 v0/v1],'P0_v',(0.05*v0)^2));
    d_est = phi_center*Xe(1:n,:); vel_est = phi_center*Xe(n+1:2*n,:);
    e = vel_est - vel_true; vl = lowpass_bw(vel_est, 5, dt); vt = lowpass_bw(vel_true, 5, dt); e_lp = vl - vt; corr_lp = corr(vl(idx)', vt(idx)');
    fprintf('noise std %5.1f mm : disp RMS err %5.2f mm | vel RMS err %.4f m/s (%.0f %% of signal RMS), bias %+.5f m/s, corr %.2f | after 5-Hz low-pass: RMS err %.4f m/s, corr %.2f | dM %+.2f %%, dv %+.3f %%\n', ...
        1e3*mean(nl), 1e3*rms(d_est(idx)-d_true(idx)), rms(e(idx)), 100*rms(e(idx))/rms(vel_true(idx)), mean(e(idx)), corr(vel_est(idx)',vel_true(idx)'), ...
        rms(e_lp(idx)), corr_lp, 100*(mean(Xe(2*n+3,idx))-M1)/M1, 100*(mean(Xe(2*n+4,idx))-v1)/v1);
    if s == 1, vel_est_paper = vel_est; vel_lp_paper = lowpass_bw(vel_est, 5, dt); end
end
k = find(t >= 10 & t <= 15); k = k(1:2:end);
writetable(table(t(k)', vel_true(k)', vel_est_paper(k)', vel_lp_paper(k)', 'VariableNames', {'t','v_true','v_est','v_est_lp'}), [F 'velstate.csv']);
disp('csv ok')
