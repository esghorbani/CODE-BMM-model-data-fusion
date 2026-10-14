function BMM_UKF_Generality_Grid(part, nparts, init)
% BMM_UKF_GENERALITY_GRID  Re-run of the manuscript's Section 6.1 generality
% study (147 cases: M = 4000:100:6000 kg, v = 13:3:33 m/s, rigid mass, smooth,
% undamped, Case-1 noise) with the CORRECTED filter.
%   init = 'paper' : 0.7 M, 0.7 v, Q = 1e-6, alpha = 0.9 (as in the dataset script)
%   init = 'bank'  : 0.7 M, starts 10:8:34 m/s, likelihood selection over 10 s
% part/nparts split the mass list over parallel MATLAB processes.
if nargin < 3, init = 'paper'; end
rng('default'); rng(1234);
diary(sprintf('generality_%s_part%d_log.txt', init, part)); diary on
par.n_modes = 4; par.m = 200; par.EI = 5.3e9; par.L = 60; par.g = 9.81;
par.mu_w = 1; par.k_v = 0; par.c_v = 0; par.prof = []; par.zeta = 0;
dt = 0.004; tspan = [0 30]; n = par.n_modes; phi_center = mode_shapesPP(par.L/2, n, par.L)';
amplitude_factor = 0.01; offset = 0.01; frequency = 5;
M_list = 4000:100:6000; v_list = 13:3:33; M_list = M_list(part:nparts:end);
results = struct([]); k = 0;
for M1 = M_list
  for v1 = v_list
    k = k + 1;
    [t, X] = ode4(@(tt,X) BMM_SS_rough(tt, X, M1, v1, par), tspan, zeros(2*n+2,1), dt);
    t = t(:)'; d = phi_center*X(:,1:n)';
    rng(1234); Nt = numel(t); tau = linspace(0,2*pi,Nt);
    y = d + (offset + amplitude_factor*abs(d) + amplitude_factor*sin(frequency*tau)).*randn(1,Nt); idx = t >= 10;
    if strcmp(init,'bank')
        bank = 10:8:34; LLs = zeros(size(bank)); Xb = cell(size(bank));
        for b = 1:numel(bank)
            [Xb{b},~,~,LL] = run_UKF_rough(t, y, M1, v1, par, dt, struct('v_init', bank(b), 'Q', 1e-6));
            LLs(b) = sum(LL(t > 0 & t <= 10));
        end
        [~,bs] = max(LLs); X_est = Xb{bs}; sel = bank(bs);
    else
        X_est = run_UKF_rough(t, y, M1, v1, par, dt, struct('Q', 1e-6)); sel = 0.7*v1;
    end
    M_est = X_est(2*n+3,:); v_est = X_est(2*n+4,:);
    dM = 100*(mean(M_est(idx))-M1)/M1; dv = 100*(mean(v_est(idx))-v1)/v1;
    fprintf('M=%5d v=%4.1f start %5.2f : dM %7.2f %%  dv %7.2f %%\n', M1, v1, sel, dM, dv);
    results(k).M = M1; results(k).v = v1; results(k).dM = dM; results(k).dv = dv; results(k).start = sel;
    results(k).M_est = M_est(1:5:end); results(k).v_est = v_est(1:5:end); results(k).t = t(1:5:end);
  end
end
save(sprintf('generality_%s_part%d_results.mat', init, part), 'results', '-v7.3');
fprintf('DONE part %d\n', part); diary off
end
