function BMM_UKF_Robustness_Grid(mode, init)
% init = 'paper' (0.7 M, 0.7 v) or 'data' (0.7 M, data-driven v0 with 5 % prior)
if nargin < 2, init = 'paper'; end
opts = struct();
% BMM_UKF_ROBUSTNESS_GRID(mode)  Mass/speed grid of the paper (M 3000-6000 kg,
% v 8-32 m/s) with the ORIGINAL filter model and the paper's 0.7 initialisation.
%   mode = 'baseline' : truth = rigid mass, undamped, smooth (as in the paper)
%   mode = 'error'    : truth = quarter-car + class C roughness + 5 % damping
% Output: robustness_grid_<mode>_results.mat / _log.txt
rng('default'); rng(1234);
diary(['robustness_grid_' mode '_' init '_log.txt']); diary on
par.n_modes = 4; par.m = 200; par.EI = 5.3e9; par.L = 60; par.g = 9.81;
dt = 0.004; tspan = [0 30]; n = par.n_modes; phi_center = mode_shapesPP(par.L/2, n, par.L)';
amplitude_factor = 0.01; offset = 0.01; frequency = 5;
n_min = 0.01; n_max = 3.0; N_harm = 300; f_v = 1.5; zeta_v = 0.20; mu_true = 0.10;
M_list = [3000 4000 5000 6000]; v_list = [8 12 16 20 24 28 32];
par_est = par; par_est.mu_w = 1; par_est.k_v = 0; par_est.c_v = 0; par_est.prof = []; par_est.zeta = 0;
prof = iso8608_profile('C', n_min, n_max, N_harm, 42);
results = struct([]); k = 0;
for M1 = M_list
  for v1 = v_list
    k = k + 1; fprintf('\n=== %s: M = %d kg, v = %d m/s ===\n', mode, M1, v1);
    par_true = par;
    if strcmp(mode,'error')
        par_true.mu_w = mu_true; par_true.zeta = 0.05; par_true.prof = prof;
        m_s = (1-mu_true)*M1; par_true.k_v = m_s*(2*pi*f_v)^2; par_true.c_v = 2*zeta_v*m_s*(2*pi*f_v);
    else
        par_true.mu_w = 1; par_true.zeta = 0; par_true.prof = []; par_true.k_v = 0; par_true.c_v = 0;
    end
    [t, X] = ode4(@(tt,X) BMM_SS_rough(tt, X, M1, v1, par_true), tspan, zeros(2*n+2,1), dt);
    t = t(:)'; disp_true = phi_center * X(:,1:n)';
    rng(1234); Nt = numel(t); tau = linspace(0,2*pi,Nt);
    nl = offset + amplitude_factor*abs(disp_true) + amplitude_factor*sin(frequency*tau);
    y = disp_true + nl.*randn(1,Nt); idx = t >= 10;
    if strcmp(init,'data')
        v0 = initial_speed_estimate(t, y, par.L);
        opts = struct('init', [0.7 v0/v1], 'P0_v', (0.05*v0)^2);
        fprintf('data-driven v0 = %.2f m/s (%.1f %%)\n', v0, 100*(v0-v1)/v1);
    elseif strcmp(init,'bank')
        % filter bank over the admissible speed range, 0.7 M start, selection
        % by innovation log-likelihood over the first 10 s (no data pre-processing)
        bank = 6:4:34; LLs = zeros(size(bank)); Xb = cell(size(bank)); Pb = Xb; Rb = Xb;
        for b = 1:numel(bank)
            [Xb{b}, Pb{b}, Rb{b}, LL] = run_UKF_rough(t, y, M1, v1, par_est, dt, struct('v_init', bank(b)));
            LLs(b) = sum(LL(t > 0 & t <= 10));
        end
        [~, bs] = max(LLs); X_est = Xb{bs}; P_saved = Pb{bs}; R0 = Rb{bs};
        fprintf('bank: selected start %d m/s\n', bank(bs)); opts = [];
    elseif strncmp(init,'mult',4)
        kmul = str2double(init(5:end));          % e.g. 'mult1.5' -> v0 = 1.5 v, default prior
        opts = struct('init', [0.7 kmul]);
    elseif strcmp(init,'v085')
        opts = struct('init', [0.7 0.85]);   % 0.85 v, default prior P0_v = 10
    else
        opts = struct();
    end
    if ~isempty(opts) || ~strcmp(init,'bank'), [X_est, P_saved, R0] = run_UKF_rough(t, y, M1, v1, par_est, dt, opts); end
    M_est = X_est(2*n+3,:); v_est = X_est(2*n+4,:);
    dM = 100*(mean(M_est(idx))-M1)/M1; dv = 100*(mean(v_est(idx))-v1)/v1;
    fprintf('dM %.2f %%  dv %.2f %%  R_end %.2e\n', dM, dv, R0(end));
    results(k).M = M1; results(k).v = v1; results(k).dM = dM; results(k).dv = dv; results(k).R_end = R0(end);
    results(k).M_est = M_est(1:5:end); results(k).v_est = v_est(1:5:end); results(k).t = t(1:5:end);
  end
end
save(['robustness_grid_' mode '_' init '_results.mat'], 'results', 'M_list', 'v_list', '-v7.3');
fprintf('\nSUMMARY %s  (rows M, cols v = %s)\n dM %%:\n', mode, mat2str(v_list));
for i = 1:numel(M_list), fprintf('%5d :', M_list(i)); for j = 1:numel(v_list), fprintf(' %7.2f', results((i-1)*numel(v_list)+j).dM); end; fprintf('\n'); end
fprintf(' dv %%:\n');
for i = 1:numel(M_list), fprintf('%5d :', M_list(i)); for j = 1:numel(v_list), fprintf(' %7.2f', results((i-1)*numel(v_list)+j).dv); end; fprintf('\n'); end
diary off
end
