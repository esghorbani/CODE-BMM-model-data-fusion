% Does the ORIGINAL 10-state Case 1 filter (run_UKF_displacement, init 0.7) converge at 30-32 m/s?
clc; clear; addpath('../case_1');
global n_modes m EI L; n_modes = 4; m = 200; EI = 5.3e9; L = 60; n = 4;
M1 = 5000; dt = 0.004; tspan = [0 30]; phi_center = mode_shapesPP(L/2, n, L)'; param.dt = dt;
for v = [15 25 30 32]
    [t, X] = ode4(@(tt,X) BMM_SS(tt, X, n, M1, m, EI, L, v), tspan, zeros(2*n,1), dt); t = t(:)';
    d = phi_center*X(:,1:n)'; rng('default'); rng(1234); Nt = numel(t); tau = linspace(0,2*pi,Nt);
    y = d + (0.01 + 0.01*abs(d) + 0.01*sin(5*tau)).*randn(1,Nt); idx = t >= 10;
    Xa = run_UKF_displacement(t, y, M1, v, param);
    fprintf('Case1 filter, v = %2d m/s : dM %+7.2f %%   dv %+7.2f %%\n', v, 100*(mean(Xa(2*n+1,idx))-M1)/M1, 100*(mean(Xa(2*n+2,idx))-v)/v);
end
