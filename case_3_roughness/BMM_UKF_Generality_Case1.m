function BMM_UKF_Generality_Case1(part, nparts, alfa)
% BMM_UKF_GENERALITY_CASE1  The manuscript's 147-case study run with the
% manuscript's own 10-state filter (case_1/run_UKF_displacement.m, corrected:
% time alignment + exact residual-form adaptive R), 0.7 M, 0.7 v, Q = 1e-6.
% alfa: forgetting factor (dataset used 0.9; helper default 0.7).
global n_modes m EI L
if nargin < 3, alfa = 0.9; end
addpath('../case_1');
rng('default'); rng(1234);
diary(sprintf('generality_case1_a%.1f_part%d_log.txt', alfa, part)); diary on
n_modes = 4; m = 200; EI = 5.3e9; L = 60; dt = 0.004; tspan = [0 30];
phi_center = mode_shapesPP(L/2, n_modes, L)';
amplitude_factor = 0.01; offset = 0.01; frequency = 5;
param.dt = dt; param.Q = 1e-6; param.alfa = alfa;
M_list = 4000:100:6000; v_list = 13:3:33; M_list = M_list(part:nparts:end);
results = struct([]); k = 0;
for M1 = M_list
  for v1 = v_list
    k = k + 1;
    [t, X] = ode4(@(tt,X) BMM_SS(tt, X, n_modes, M1, m, EI, L, v1), tspan, zeros(2*n_modes,1), dt);
    t = t(:)'; d = phi_center*X(:,1:n_modes)';
    rng(1234); Nt = numel(t); tau = linspace(0,2*pi,Nt);
    y = d + (offset + amplitude_factor*abs(d) + amplitude_factor*sin(frequency*tau)).*randn(1,Nt); idx = t >= 10;
    X_est = run_UKF_displacement(t, y, M1, v1, param);
    M_est = X_est(2*n_modes+1,:); v_est = X_est(2*n_modes+2,:);
    dM = 100*(mean(M_est(idx))-M1)/M1; dv = 100*(mean(v_est(idx))-v1)/v1;
    fprintf('M=%5d v=%4.1f : dM %7.2f %%  dv %7.2f %%\n', M1, v1, dM, dv);
    results(k).M = M1; results(k).v = v1; results(k).dM = dM; results(k).dv = dv;
    results(k).M_est = M_est(1:5:end); results(k).v_est = v_est(1:5:end); results(k).t = t(1:5:end);
  end
end
save(sprintf('generality_case1_a%.1f_part%d_results.mat', alfa, part), 'results', '-v7.3');
fprintf('DONE part %d\n', part); diary off
end
