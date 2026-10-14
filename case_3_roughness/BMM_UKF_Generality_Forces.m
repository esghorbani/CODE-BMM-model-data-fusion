function BMM_UKF_Generality_Forces(part, nparts, mode)
% BMM_UKF_GENERALITY_FORCES  Section 6 generality study (147 cases: M = 4000:100:6000 kg,
% v = 13:3:31 m/s, rigid mass, smooth deck, undamped, Case-1 noise) with the CORRECTED
% filter and the four-filter speed bank, storing what Figures 9-10 of the manuscript need:
% decimated mass/speed histories with their standard deviations, and the reconstructed
% contact force (point estimate from the estimated states, mass and speed) at full rate.
% part/nparts split the mass list over parallel MATLAB processes.
if nargin < 2, nparts = 1; end
if nargin < 3, mode = 'clean'; end
% mode = 'clean' : truth = rigid mass, smooth deck, undamped (as in Section 5)
% mode = 'error' : truth = quarter-car (mu_w = 0.1, f_v = 1.5 Hz, zeta_v = 0.2) + ISO 8608 class C
%                  roughness + 5 %% beam damping; the filter keeps the rigid/smooth/undamped model
sfx = ''; if strcmp(mode,'error'), sfx = '_error'; end
rng('default'); rng(1234);
diary(sprintf('generality_forces%s_part%d_log.txt', sfx, part)); diary on
par.n_modes = 4; par.m = 200; par.EI = 5.3e9; par.L = 60; par.g = 9.81;
par.mu_w = 1; par.k_v = 0; par.c_v = 0; par.prof = []; par.zeta = 0;   % filter model (rigid, smooth, undamped)
n_min = 0.01; n_max = 3.0; N_harm = 300; f_v = 1.5; zeta_v = 0.20; mu_true = 0.10;
if strcmp(mode,'error'), prof = iso8608_profile('C', n_min, n_max, N_harm, 42); end
dt = 0.004; tspan = [0 30]; n = par.n_modes; phi_center = mode_shapesPP(par.L/2, n, par.L)';
amplitude_factor = 0.01; offset = 0.01; frequency = 5;
M_list = 4000:100:6000; v_list = 13:3:31; M_list = M_list(part:nparts:end);
results = struct([]); k = 0; dec = 5;
for M1 = M_list
  for v1 = v_list
    k = k + 1; t0 = tic;
    par_true = par;
    if strcmp(mode,'error')
        par_true.mu_w = mu_true; par_true.zeta = 0.05; par_true.prof = prof;
        m_s = (1-mu_true)*M1; par_true.k_v = m_s*(2*pi*f_v)^2; par_true.c_v = 2*zeta_v*m_s*(2*pi*f_v);
    end
    [t, X] = ode4(@(tt,X) BMM_SS_rough(tt, X, M1, v1, par_true), tspan, zeros(2*n+2,1), dt);
    t = t(:)'; d = phi_center*X(:,1:n)';
    rng(1234); Nt = numel(t); tau = linspace(0,2*pi,Nt);
    y = d + (offset + amplitude_factor*abs(d) + amplitude_factor*sin(frequency*tau)).*randn(1,Nt); idx = t >= 10;
    bank = 10:8:34; LLs = zeros(size(bank)); Xb = cell(size(bank)); Pb = cell(size(bank));
    for b = 1:numel(bank)
        [Xb{b},Pb{b},~,LL] = run_UKF_rough(t, y, M1, v1, par, dt, struct('v_init', bank(b), 'Q', 1e-6));
        LLs(b) = sum(LL(t > 0 & t <= 10));
    end
    [~,bs] = max(LLs); X_est = Xb{bs}; P_saved = Pb{bs}; sel = bank(bs);
    M_est = X_est(2*n+3,:); v_est = X_est(2*n+4,:);
    sig_M = sqrt(squeeze(P_saved(2*n+3,2*n+3,:)))'; sig_v = sqrt(squeeze(P_saved(2*n+4,2*n+4,:)))';
    fc_est = contact_force_rough(t, X_est(1:2*n+2,:), M_est, v_est, par);
    fc_true = contact_force_rough(t, X', M1, v1, par_true);
    F_static = M_est*par.g;
    i3 = t >= t(end) - 3;
    dM = 100*(mean(M_est(idx))-M1)/M1; dv = 100*(mean(v_est(idx))-v1)/v1;
    fprintf('M=%5d v=%4.1f start %5.2f : dM %7.2f %%  dv %7.2f %%  Fdyn/Fstat(last 3 s) %.4f  (%.0f s)\n', ...
        M1, v1, sel, dM, dv, mean(fc_est(i3))/mean(F_static(i3)), toc(t0));
    r.M = M1; r.v = v1; r.dM = dM; r.dv = dv; r.start = sel;
    r.t = t(1:dec:end); r.M_est = M_est(1:dec:end); r.v_est = v_est(1:dec:end);
    r.sig_M = sig_M(1:dec:end); r.sig_v = sig_v(1:dec:end);
    r.t_full = t; r.fc_est = fc_est; r.fc_true = fc_true; r.F_static = F_static;
    r.Fdyn_mean3 = mean(fc_est(i3)); r.Fstat_mean3 = mean(F_static(i3)); r.Ftrue_mean3 = mean(fc_true(i3));
    if isempty(results), results = r; else, results(k) = r; end
  end
end
save(sprintf('generality_forces%s_part%d_results.mat', sfx, part), 'results', '-v7.3');
fprintf('DONE part %d\n', part); diary off
end
