% BMM_UKF_ALPHA_SWEEP  Forgetting-factor sensitivity (R2-C3), rigid mass, smooth
% deck, manuscript initialisation (0.7 M, 0.7 v), Q = 10^-6.2, corrected filter
% (innovation-form adaptive R).  Exports alpha_sensitivity_init07.csv for the letter.
clc; clear; rng('default'); rng(1234); diary('alpha_sweep_log.txt'); diary on
par.n_modes = 4; par.m = 200; par.EI = 5.3e9; par.L = 60; par.g = 9.81;
par.mu_w = 1; par.k_v = 0; par.c_v = 0; par.prof = []; par.zeta = 0;
M1 = 5000; v1 = 15; dt = 0.004; tspan = [0 30]; n = par.n_modes; phi_center = mode_shapesPP(par.L/2, n, par.L)';
amplitude_factor = 0.01; offset = 0.01; frequency = 5;
[t, X] = ode4(@(tt,X) BMM_SS_rough(tt, X, M1, v1, par), tspan, zeros(2*n+2,1), dt); t = t(:)';
d = phi_center*X(:,1:n)'; rng(1234); Nt = numel(t); tau = linspace(0,2*pi,Nt);
y = d + (offset + amplitude_factor*abs(d) + amplitude_factor*sin(frequency*tau)).*randn(1,Nt); idx = t >= 10;
als = [0.5 0.7 0.75 0.8 0.85 0.9 0.95]; k = find(t <= 15); k = k(1:5:end); Mx = t(k)'; Vx = t(k)'; names = {'t'};
for a = als
    [Xe,~,R0] = run_UKF_rough(t, y, M1, v1, par, dt, struct('alfa', a));
    Me = Xe(2*n+3,:); ve = Xe(2*n+4,:);
    fprintf('alpha = %.2f : dM %7.2f %%  dv %7.2f %%  R(t>=10) %.2e\n', a, 100*(mean(Me(idx))-M1)/M1, 100*(mean(ve(idx))-v1)/v1, mean(R0(idx)));
    Mx = [Mx, Me(k)']; Vx = [Vx, ve(k)']; names{end+1} = sprintf('a%03d', round(100*a));
end
F = '../../PAPER-2026-BMM-FUSION/JSV_Submission/rebuttal/2026-09-11-round-1/figures/';
T = array2table([Mx, Vx(:,2:end)], 'VariableNames', [names, strrep(names(2:end), 'a', 'v')]);
writetable(T, [F 'alpha_sensitivity_init07.csv']); fprintf('DONE\n'); diary off
