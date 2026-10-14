% BMM_UKF_REVIEWER2_STUDIES
% Numerical studies for Reviewer 2:
%   (A) R2-C3(1): sensitivity to the (constant) process-noise level Q
%   (B) R2-C3(2): sensitivity to the forgetting factor alpha of the adaptive R
%   (C) R2-C4   : speeds 10-33 m/s (36-120 km/h), dynamic amplification and interaction force
%   (D) R2-C5   : sensitivity of the mid-span displacement to M and v and the
%                 resulting Cramer-Rao bounds (why v converges more slowly)
% All runs: rigid moving mass of the paper (mu_w = 1), smooth deck, Case 1
% noise model, 12-state filter of case_3 (run_UKF_rough).
% Outputs: reviewer2_results.mat, reviewer2_log.txt and, via PLOT_REVIEWER2_FIGURES,
%          Q_alpha_sensitivity.pdf, Speed_study.pdf, Velocity_sensitivity.pdf

clc; close all; clear;
rng('default'); rng(1234);
par.n_modes = 4; par.m = 200; par.EI = 5.3e9; par.L = 60; par.g = 9.81;
par.prof = []; par.mu_w = 1; par.k_v = 0; par.c_v = 0;
M1 = 5000; v1 = 15; dt = 0.004; tspan = [0 30]; n = par.n_modes;
phi_center = mode_shapesPP(par.L/2, n, par.L)';
amplitude_factor = 0.01; offset = 0.01; frequency = 5;
fs = 10; fn = 'Times New Roman';

simulate = @(M, v, p) ode4(@(tt,X) BMM_SS_rough(tt, X, M, v, p), tspan, zeros(2*n+2,1), dt);
function [y, nl] = add_noise(d, offset, af, freq)
    rng(1234); Nt = numel(d); tau = linspace(0, 2*pi, Nt);
    nl = offset + af*abs(d) + af*sin(freq*tau);
    y = d + nl.*randn(1, Nt);
end
conv_time = @(x, xt, tol, t) t(min(find(abs(x - xt)/xt > tol, 1, 'last') + 1, numel(t)));  % first time after which |err| < tol (= t_end if never)

[t, X] = simulate(M1, v1, par); t = t(:)';
d0 = phi_center*X(:,1:n)'; [y0, nl0] = add_noise(d0, offset, amplitude_factor, frequency);
idx = t >= 10;
v0 = initial_speed_estimate(t, y0, par.L);
fprintf('data-driven initial speed: %.2f m/s (true %d, error %+.1f %%)\n', v0, v1, 100*(v0-v1)/v1);
% data-driven initialisation: v0 from the signal, prior std = SIG_V0 * v0
SIG_V0 = 0.05;
inits = {'paper init (0.7 v)', struct(); 'data-driven v0', struct('init', [0.7 v0/v1], 'P0_v', (SIG_V0*v0)^2)};

% ---------------------------------------------------------------- (A) Q ---
fprintf('\n(A) process-noise sensitivity (alpha = 0.9)\n');
Qs = [1e-8 1e-7 10^-6.2 1e-5 1e-4];
A = struct([]);
for ii = 1:2, fprintf('  -- %s --\n', inits{ii,1});
for k = 1:numel(Qs)
    o = inits{ii,2}; o.Q = Qs(k); [Xe, Pe, R0] = run_UKF_rough(t, y0, M1, v1, par, dt, o);
    Me = Xe(2*n+3,:); ve = Xe(2*n+4,:);
    r = struct('init', ii, 'Q', Qs(k), 'M', Me, 'v', ve, 'R0', R0, ...
        'dM', 100*(mean(Me(idx))-M1)/M1, 'dv', 100*(mean(ve(idx))-v1)/v1, ...
        'sM', 100*std(Me(idx))/M1, 'sv', 100*std(ve(idx))/v1, ...
        'tM', conv_time(Me, M1, 0.05, t), 'tv', conv_time(ve, v1, 0.02, t));
    if isempty(A), A = r; else, A(end+1) = r; end
    fprintf('Q = %8.1e : dM %+6.2f %%  std %5.2f %%   dv %+6.3f %%  std %5.3f %%   t(|dM|<5%%) %5.1f s  t(|dv|<2%%) %5.1f s\n', ...
        Qs(k), r.dM, r.sM, r.dv, r.sv, r.tM, r.tv);
end
end

% ------------------------------------------------------------ (B) alpha ---
fprintf('\n(B) forgetting-factor sensitivity (Q = 10^-6.2)\n');
als = [0.5 0.7 0.9 0.95 0.99];
B = struct([]);
for ii = 1:2, fprintf('  -- %s --\n', inits{ii,1});
for k = 1:numel(als)
    o = inits{ii,2}; o.alfa = als(k); [Xe, Pe, R0] = run_UKF_rough(t, y0, M1, v1, par, dt, o);
    Me = Xe(2*n+3,:); ve = Xe(2*n+4,:);
    trackR = 100*rms(sqrt(R0(idx)) - nl0(idx))/mean(nl0(idx));   % how well sqrt(R) follows the true noise std
    r = struct('init', ii, 'alfa', als(k), 'M', Me, 'v', ve, 'R0', R0, ...
        'dM', 100*(mean(Me(idx))-M1)/M1, 'dv', 100*(mean(ve(idx))-v1)/v1, ...
        'sM', 100*std(Me(idx))/M1, 'trackR', trackR, 'tM', conv_time(Me, M1, 0.05, t));
    if isempty(B), B = r; else, B(end+1) = r; end
    fprintf('alpha = %4.2f : dM %+6.2f %%  std %5.2f %%   dv %+6.3f %%   t(|dM|<5%%) %5.1f s   RMS(sqrt(R)-sigma_n)/mean(sigma_n) %5.1f %%\n', ...
        als(k), r.dM, r.sM, r.dv, r.tM, trackR);
end
end

% ------------------------------------------------------------ (C) speed ---
fprintf('\n(C) speed study (rigid mass)\n');
vs = [10 15 20 25 30 33];     % 36-120 km/h: posted road speed limits (R2-C4)
w1 = sqrt(par.EI/par.m)*(pi/par.L)^2; vcr = w1*par.L/pi;      % critical speed 2 f1 L
d_static = M1*par.g*par.L^3/(48*par.EI);
fprintf('f1 = %.2f Hz, critical speed %.0f m/s, static mid-span deflection %.1f mm\n', w1/2/pi, vcr, 1e3*d_static);
C = struct([]);
for k = 1:numel(vs)
    v = vs(k);
    [tv, Xv] = simulate(M1, v, par); tv = tv(:)';
    dv_ = phi_center*Xv(:,1:n)'; yv = add_noise(dv_, offset, amplitude_factor, frequency);
    fcv = contact_force_rough(tv, Xv', M1, v, par);
    v0v = initial_speed_estimate(tv, yv, par.L);
    [Xe, Pe]   = run_UKF_rough(tv, yv, M1, v, par, dt);                                 % paper init 0.7 v
    [Xe2, Pe2] = run_UKF_rough(tv, yv, M1, v, par, dt, struct('init', [0.7 v0v/v], 'P0_v', (SIG_V0*v0v)^2));   % data-driven v0, matched prior
    Me = Xe(2*n+3,:); ve = Xe(2*n+4,:); Me2 = Xe2(2*n+3,:); ve2 = Xe2(2*n+4,:); idv = tv >= 10;
    r = struct('v', v, 'v0', v0v, 'M', Me, 'vest', ve, 'M2', Me2, 'vest2', ve2, 'disp', dv_, 'fc', fcv, 't', tv, ...
        'DAF', max(dv_(idv))/d_static, 'fdyn_max', 100*max(abs(fcv(idv)-M1*par.g))/(M1*par.g), ...
        'fdyn_std', 100*std(fcv(idv))/(M1*par.g), ...
        'dM', 100*(mean(Me(idv))-M1)/M1, 'dv', 100*(mean(ve(idv))-v)/v, ...
        'dM2', 100*(mean(Me2(idv))-M1)/M1, 'dv2', 100*(mean(ve2(idv))-v)/v, 'sv2', 100*std(ve2(idv))/v);
    if isempty(C), C = r; else, C(end+1) = r; end
    fprintf('v = %2d m/s (%3.0f km/h, v/vcr %.2f): DAF %.3f  dyn. force max %5.1f %% std %4.1f %% of Mg | v0 %.1f (%+.1f %%) | paper init: dM %+6.2f %% dv %+6.2f %% | data init: dM %+6.2f %% dv %+6.3f %% (std %.3f %%)\n', ...
        v, 3.6*v, v/vcr, r.DAF, r.fdyn_max, r.fdyn_std, v0v, 100*(v0v-v)/v, r.dM, r.dv, r.dM2, r.dv2, r.sv2);
end

% ------------------------------------------------- (D) sensitivity to M, v ---
fprintf('\n(D) sensitivity of the mid-span displacement to M and v\n');
eps_rel = 0.002;
[~, Xp] = simulate(M1*(1+eps_rel), v1, par); [~, Xm] = simulate(M1*(1-eps_rel), v1, par);
S_M = (phi_center*Xp(:,1:n)' - phi_center*Xm(:,1:n)')/(2*eps_rel);        % dy/dM * M   [m per 100 %]
[~, Xp] = simulate(M1, v1*(1+eps_rel), par); [~, Xm] = simulate(M1, v1*(1-eps_rel), par);
S_v = (phi_center*Xp(:,1:n)' - phi_center*Xm(:,1:n)')/(2*eps_rel);        % dy/dv * v
sig = 0.01;                                                                 % nominal noise std [m]
I_M = cumsum(S_M.^2)/sig^2; I_v = cumsum(S_v.^2)/sig^2; I_Mv = cumsum(S_M.*S_v)/sig^2;
crlb_M = zeros(size(t)); crlb_v = zeros(size(t)); rho = zeros(size(t));
for k = 2:numel(t)
    F = [I_M(k) I_Mv(k); I_Mv(k) I_v(k)]; Ci = inv(F);
    crlb_M(k) = 100*sqrt(Ci(1,1)); crlb_v(k) = 100*sqrt(Ci(2,2)); rho(k) = I_Mv(k)/sqrt(I_M(k)*I_v(k));
end
for T = [1 2 4 8 15]
    k = find(t >= T, 1);
    fprintf('t = %4.1f s : CRLB rel. std  M %6.2f %%   v %6.2f %%   corr(S_M,S_v) %+.2f\n', T, crlb_M(k), crlb_v(k), rho(k));
end
fprintf('RMS normalised sensitivity over first crossing: |S_M| %.4f m, |S_v| %.4f m\n', rms(S_M(t<4)), rms(S_v(t<4)));
D = struct('t', t, 'S_M', S_M, 'S_v', S_v, 'crlb_M', crlb_M, 'crlb_v', crlb_v, 'rho', rho);
save('reviewer2_results.mat', 'A', 'B', 'C', 'D', 'par', 'M1', 'v1', 'dt', 'SIG_V0', '-v7.3');

%% ---------------------------------------------------------------- figures
plot_reviewer2_figures
