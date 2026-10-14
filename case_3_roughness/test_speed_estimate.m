clc; clear;
par.n_modes = 4; par.m = 200; par.EI = 5.3e9; par.L = 60; par.g = 9.81;
par.prof = []; par.mu_w = 1; par.k_v = 0; par.c_v = 0;
M1 = 5000; dt = 0.004; tspan = [0 30]; n = 4; phi_center = mode_shapesPP(par.L/2, n, par.L)';
for v = [8 15 20 30 45 60 80]
    [t, X] = ode4(@(tt,X) BMM_SS_rough(tt, X, M1, v, par), tspan, zeros(2*n+2,1), dt); t = t(:)';
    d = phi_center*X(:,1:n)'; rng(1234); Nt = numel(t); tau = linspace(0,2*pi,Nt);
    y = d + (0.01 + 0.01*abs(d) + 0.01*sin(5*tau)).*randn(1,Nt);
    v0 = initial_speed_estimate(t, y, par.L);
    v0s = initial_speed_estimate(t(t<=12), y(t<=12), par.L);      % shorter record
    fprintf('v = %2d : v0 = %6.2f (%+5.1f %%)   12-s record: %6.2f (%+5.1f %%)\n', v, v0, 100*(v0-v)/v, v0s, 100*(v0s-v)/v);
end
