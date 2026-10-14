% Compare filter tunings on rigid-mass data (Case 1 setting, smooth deck).
clc; clear; addpath('../case_1');
global n_modes m EI L
n_modes = 4; m = 200; EI = 5.3e9; L = 60;
par.n_modes = 4; par.m = 200; par.EI = 5.3e9; par.L = 60; par.g = 9.81; par.prof = [];
par.mu_w = 1; par.k_v = 0; par.c_v = 0;
M1 = 5000; v1 = 15; dt = 0.004; tspan = [0 30]; n = 4;
phi_center = mode_shapesPP(L/2, n, L)';
[t, X] = ode4(@(tt,X) BMM_SS(tt, X, n, M1, m, EI, L, v1), tspan, zeros(2*n,1), dt); t = t(:)';
disp_true = phi_center*X(:,1:n)';
rng('default'); rng(1234);
Nt = numel(t); tau = linspace(0,2*pi,Nt);
nl = 0.01 + 0.01*abs(disp_true) + 0.01*sin(5*tau);
y = disp_true + nl.*randn(1,Nt);
idx = t >= 10;
param.dt = dt;
Xa = run_UKF_displacement(t, y, M1, v1, param);
fprintf('Case1 run_UKF_displacement (alfa .7, Q 1e-6, 10 states): M err %.2f %%\n', 100*(mean(Xa(2*n+1,idx))-M1)/M1);
Xb = run_UKF_rough(t, y, M1, v1, par, dt, struct('alfa',0.7,'Q',1e-6));
fprintf('rough, mu_w=1, alfa .7, Q 1e-6, 12 states           : M err %.2f %%\n', 100*(mean(Xb(2*n+3,idx))-M1)/M1);
Xc = run_UKF_rough(t, y, M1, v1, par, dt, struct('alfa',0.9,'Q',10^-6.2));
fprintf('rough, mu_w=1, alfa .9, Q 10^-6.2, 12 states        : M err %.2f %%\n', 100*(mean(Xc(2*n+3,idx))-M1)/M1);
Xd = run_UKF_rough(t, y, M1, v1, par, dt, struct('alfa',0.7,'Q',1e-6,'P0_veh',1e-12));
fprintf('rough, mu_w=1, alfa .7, Q 1e-6, veh P0 ~0           : M err %.2f %%\n', 100*(mean(Xd(2*n+3,idx))-M1)/M1);
par2 = par; par2.mu_w = 0.1; m_s = 0.9*M1; par2.k_v = m_s*(2*pi*1.5)^2; par2.c_v = 2*0.2*m_s*(2*pi*1.5);
[t2, X2] = ode4(@(tt,X) BMM_SS_rough(tt, X, M1, v1, par2), tspan, zeros(2*n+2,1), dt);
d2 = phi_center*X2(:,1:n)'; rng('default'); rng(1234); y2 = d2 + (0.01 + 0.01*abs(d2) + 0.01*sin(5*tau)).*randn(1,Nt);
Xe = run_UKF_rough(t, y2, M1, v1, par2, dt, struct('alfa',0.9,'Q',10^-6.2));
fprintf('quarter-car, alfa .9, Q 10^-6.2                     : M err %.2f %%\n', 100*(mean(Xe(2*n+3,idx))-M1)/M1);
