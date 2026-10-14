% Velocity-state error of the ORIGINAL Case 1 filter (10 states, init 0.7, run_UKF_displacement)
% on the same data as BMM_UKF_VelocityState, plus export of the 1-mm-noise estimate for the letter figure.
clc; clear; addpath('../case_1');
global n_modes m EI L; n_modes = 4; m = 200; EI = 5.3e9; L = 60; n = 4;
par.n_modes = 4; par.m = 200; par.EI = 5.3e9; par.L = 60; par.g = 9.81; par.prof = []; par.mu_w = 1; par.k_v = 0; par.c_v = 0;
M1 = 5000; v1 = 15; dt = 0.004; tspan = [0 30]; phi_center = mode_shapesPP(L/2, n, L)';
[t, X] = ode4(@(tt,X) BMM_SS(tt, X, n, M1, m, EI, L, v1), tspan, zeros(2*n,1), dt); t = t(:)';
d_true = phi_center*X(:,1:n)'; vel_true = phi_center*X(:,n+1:2*n)'; Nt = numel(t); tau = linspace(0,2*pi,Nt); idx = t >= 10;
rng('default'); rng(1234); nl = 0.01 + 0.01*abs(d_true) + 0.01*sin(5*tau); y = d_true + nl.*randn(1,Nt);
param.dt = dt; Xa = run_UKF_displacement(t, y, M1, v1, param);
va = phi_center*Xa(n+1:2*n,:);
fprintf('ORIGINAL filter: vel state RMS err %.3f m/s, peak |est| %.2f m/s, peak |true| %.3f m/s, corr %.2f\n', rms(va(idx)-vel_true(idx)), max(abs(va(idx))), max(abs(vel_true(idx))), corr(va(idx)',vel_true(idx)'));
% revised filter at 1 mm noise -> export for figure
rng(1234); nl1 = 0.1*nl; y1 = d_true + nl1.*randn(1,Nt); v0 = initial_speed_estimate(t, y1, L);
Xe = run_UKF_rough(t, y1, M1, v1, par, dt, struct('init',[0.7 v0/v1],'P0_v',(0.05*v0)^2)); ve1 = phi_center*Xe(n+1:2*n,:);
fprintf('REVISED filter, 1 mm noise: vel state RMS err %.4f m/s, peak |est| %.3f m/s, corr %.2f\n', rms(ve1(idx)-vel_true(idx)), max(abs(ve1(idx))), corr(ve1(idx)',vel_true(idx)'));
F = '../../PAPER-2026-BMM-FUSION/JSV_Submission/rebuttal/2026-09-11-round-1/figures/';
T = readtable([F 'velstate.csv']); k = find(t >= 10 & t <= 15); k = k(1:2:end);
assert(numel(k) == height(T)); T.v_est_lownoise = ve1(k)'; writetable(T, [F 'velstate.csv']); disp('csv ok')
