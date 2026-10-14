% Data-driven v0 with a prior variance consistent with its accuracy, all speeds.
clc; clear;
par.n_modes = 4; par.m = 200; par.EI = 5.3e9; par.L = 60; par.g = 9.81;
par.prof = []; par.mu_w = 1; par.k_v = 0; par.c_v = 0;
M1 = 5000; dt = 0.004; tspan = [0 30]; n = 4; phi_center = mode_shapesPP(par.L/2, n, par.L)';
for v = [8 10 15 20 25 30 33]
    [t, X] = ode4(@(tt,X) BMM_SS_rough(tt, X, M1, v, par), tspan, zeros(2*n+2,1), dt); t = t(:)';
    d = phi_center*X(:,1:n)'; rng(1234); Nt = numel(t); tau = linspace(0,2*pi,Nt);
    y = d + (0.01 + 0.01*abs(d) + 0.01*sin(5*tau)).*randn(1,Nt); idx = t >= 10;
    v0 = initial_speed_estimate(t, y, par.L);
    out = sprintf('v = %2d (v0 %+.1f %%):', v, 100*(v0-v)/v);
    for P0v = [10, (0.05*v0)^2, (0.02*v0)^2, 0.1]
        Xe = run_UKF_rough(t, y, M1, v, par, dt, struct('init', [0.7 v0/v], 'P0_v', P0v));
        out = [out sprintf('  P0v=%6.3f: dM %+6.2f dv %+6.3f |', P0v, 100*(mean(Xe(2*n+3,idx))-M1)/M1, 100*(mean(Xe(2*n+4,idx))-v)/v)];
    end
    disp(out)
end
