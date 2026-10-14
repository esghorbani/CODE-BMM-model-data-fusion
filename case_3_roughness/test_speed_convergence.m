% Basin of convergence of the AUKF at higher speeds (rigid mass, smooth deck).
clc; clear;
par.n_modes = 4; par.m = 200; par.EI = 5.3e9; par.L = 60; par.g = 9.81;
par.prof = []; par.mu_w = 1; par.k_v = 0; par.c_v = 0;
M1 = 5000; dt = 0.004; tspan = [0 30]; n = 4;
phi_center = mode_shapesPP(par.L/2, n, par.L)';
variants = {
  'paper: init 0.7, P0v 10, Q 10^-6.2', struct();
  'init 0.9',                            struct('init',[0.7 0.9]);
  'init 1.2',                            struct('init',[0.7 1.2]);
  'P0v 100',                             struct('P0_v',100);
  'Q 1e-8',                              struct('Q',1e-8);
  'Q 1e-8, P0v 100',                     struct('Q',1e-8,'P0_v',100);
  'init 0.9, Q 1e-8',                    struct('init',[0.7 0.9],'Q',1e-8);
};
for v = [30 45 60]
    [t, X] = ode4(@(tt,X) BMM_SS_rough(tt, X, M1, v, par), tspan, zeros(2*n+2,1), dt); t = t(:)';
    d = phi_center*X(:,1:n)'; rng(1234); Nt = numel(t); tau = linspace(0,2*pi,Nt);
    y = d + (0.01 + 0.01*abs(d) + 0.01*sin(5*tau)).*randn(1,Nt); idx = t >= 10;
    fprintf('\n--- v = %d m/s ---\n', v);
    for k = 1:size(variants,1)
        Xe = run_UKF_rough(t, y, M1, v, par, dt, variants{k,2});
        fprintf('%-36s : dM %+7.2f %%   dv %+7.2f %%\n', variants{k,1}, 100*(mean(Xe(2*n+3,idx))-M1)/M1, 100*(mean(Xe(2*n+4,idx))-v)/v);
    end
end
