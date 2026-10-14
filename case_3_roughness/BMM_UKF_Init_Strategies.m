% BMM_UKF_INIT_STRATEGIES  Alternatives to the data-driven initial speed, keeping
% the paper's 0.7 initialisation.  Truth = paper model (rigid, smooth, undamped).
%  A: 0.7 M, 0.7 v, initial R0 = 1e-4 instead of 1e-1
%  B: filter bank: 0.7 M and v_init in {6,10,...,34} m/s (never equal to the
%     true speed), pick the filter with the largest innovation log-likelihood
%     over the first 10 s
clc; close all; clear; rng('default'); rng(1234);
diary('init_strategies_log.txt'); diary on
par.n_modes = 4; par.m = 200; par.EI = 5.3e9; par.L = 60; par.g = 9.81;
par.mu_w = 1; par.k_v = 0; par.c_v = 0; par.prof = []; par.zeta = 0;
dt = 0.004; tspan = [0 30]; n = par.n_modes; phi_center = mode_shapesPP(par.L/2, n, par.L)';
amplitude_factor = 0.01; offset = 0.01; frequency = 5;
pts = [3000 8; 5000 8; 3000 12; 5000 16; 6000 20; 5000 28; 5000 32; 6000 32];
bank = 6:4:34;
fprintf('%6s %4s | %8s %8s | %8s %8s | %8s %8s %6s\n','M','v','A dM','A dv','B dM','B dv','0.7 dM','0.7 dv','vsel');
out = struct([]);
for i = 1:size(pts,1)
    M1 = pts(i,1); v1 = pts(i,2);
    [t, X] = ode4(@(tt,X) BMM_SS_rough(tt, X, M1, v1, par), tspan, zeros(2*n+2,1), dt);
    t = t(:)'; d = phi_center*X(:,1:n)';
    rng(1234); Nt = numel(t); tau = linspace(0,2*pi,Nt);
    y = d + (offset + amplitude_factor*abs(d) + amplitude_factor*sin(frequency*tau)).*randn(1,Nt); idx = t >= 10;
    err = @(Xe) deal(100*(mean(Xe(2*n+3,idx))-M1)/M1, 100*(mean(Xe(2*n+4,idx))-v1)/v1);
    Xe = run_UKF_rough(t, y, M1, v1, par, dt); [dM0, dv0] = err(Xe);
    Xe = run_UKF_rough(t, y, M1, v1, par, dt, struct('R0_init',1e-4)); [dMA, dvA] = err(Xe);
    LLs = zeros(size(bank)); dMs = LLs; dvs = LLs;
    for b = 1:numel(bank)
        [Xe,~,~,LL] = run_UKF_rough(t, y, M1, v1, par, dt, struct('v_init',bank(b)));
        LLs(b) = sum(LL(t<=10 & t>0)); [dMs(b), dvs(b)] = err(Xe);
    end
    [~,bs] = max(LLs);
    fprintf('%6d %4d | %8.2f %8.2f | %8.2f %8.2f | %8.2f %8.2f %6d\n', M1, v1, dMA, dvA, dMs(bs), dvs(bs), dM0, dv0, bank(bs));
    fprintf('        bank dM: %s\n        bank LL: %s\n', mat2str(round(dMs,1)), mat2str(round(LLs)));
    out(i).M = M1; out(i).v = v1; out(i).dMA = dMA; out(i).dvA = dvA; out(i).dMs = dMs; out(i).dvs = dvs; out(i).LLs = LLs; out(i).sel = bank(bs);
end
save('init_strategies_results.mat','out','bank');
fprintf('DONE\n'); diary off
