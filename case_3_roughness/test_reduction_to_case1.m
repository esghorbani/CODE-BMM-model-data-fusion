% TEST_REDUCTION_TO_CASE1
% Sanity check: with mu_w = 1 (all mass unsprung) and a smooth deck, the
% extended model VBI_ROUGH_RHS must reproduce case_1/BMM_SS.m exactly.
% Also checks the analytic roughness derivatives against finite differences.
clc; clear;
addpath('../case_1');

par.n_modes = 4; par.m = 200; par.EI = 5.3e9; par.L = 60; par.g = 9.81;
par.mu_w = 1; par.k_v = 0; par.c_v = 0; par.prof = [];
M1 = 5000; v1 = 15; dt = 0.004; tspan = [0 8];
n = par.n_modes;

[t, Xa] = ode4(@(tt,X) BMM_SS(tt, X, n, M1, par.m, par.EI, par.L, v1), tspan, zeros(2*n,1), dt);
[~, Xb] = ode4(@(tt,X) BMM_SS_rough(tt, X, M1, v1, par), tspan, zeros(2*n+2,1), dt);

err = max(abs(Xa(:) - reshape(Xb(:,1:2*n),[],1)));
fprintf('max |BMM_SS - BMM_SS_rough(mu_w=1, r=0)| = %.3e\n', err);
assert(err < 1e-12, 'Reduction to Case 1 failed');

% roughness derivative check
prof = iso8608_profile('C', 0.01, 3, 300, 7);
x = linspace(0, 60, 2001); hx = 1e-4;
[r, dr, ddr] = eval_roughness(prof, x);
[rp] = eval_roughness(prof, x + hx); [rm] = eval_roughness(prof, x - hx);
fprintf('max |dr  - FD| = %.3e\n', max(abs(dr  - (rp - rm)/(2*hx))));
fprintf('max |ddr - FD| = %.3e\n', max(abs(ddr - (rp - 2*r + rm)/hx^2)));
fprintf('class C: theoretical RMS = %.2f mm, sample RMS = %.2f mm\n', 1e3*prof.sigma_r, 1e3*std(r));
disp('OK');
