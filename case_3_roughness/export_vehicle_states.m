% EXPORT_VEHICLE_STATES  Sprung-mass displacement u_v and velocity u_v' (states) for the
% letter, class C, roughness not included / included in the filter model.
load('roughness_results.mat');
F = '../../PAPER-2026-BMM-FUSION/JSV_Submission/rebuttal/2026-09-11-round-1/figures/';
c = numel(results); r = results(c); n = par.n_modes; t = r.t;
uv_true = r.X_true(:, 2*n+1)'; uvd_true = r.X_true(:, 2*n+2)';
uv_unk  = r.X_est(2*n+1, :);   uvd_unk  = r.X_est(2*n+2, :);
% profile known to the filter: same setting as known_data (data-driven v0, 5 % prior)
par_k = par; par_k.prof = r.prof; v0 = initial_speed_estimate(t, r.y, par.L);
X_k = run_UKF_rough(t, r.y, M1, v1, par_k, dt, struct('init', [0.7 v0/v1], 'P0_v', (0.05*v0)^2));
uv_kn = X_k(2*n+1, :); uvd_kn = X_k(2*n+2, :);
idx = t >= 10;
fprintf('class %s: u_v true RMS %.2f mm | RMS err unknown %.2f mm (corr %.2f), known %.2f mm (corr %.2f)\n', r.class, ...
    1e3*rms(uv_true(idx)), 1e3*rms(uv_unk(idx)-uv_true(idx)), corr(uv_unk(idx)',uv_true(idx)'), 1e3*rms(uv_kn(idx)-uv_true(idx)), corr(uv_kn(idx)',uv_true(idx)'));
fprintf('         u_v'' true RMS %.3f m/s | RMS err unknown %.3f (corr %.2f), known %.3f (corr %.2f)\n', ...
    rms(uvd_true(idx)), rms(uvd_unk(idx)-uvd_true(idx)), corr(uvd_unk(idx)',uvd_true(idx)'), rms(uvd_kn(idx)-uvd_true(idx)), corr(uvd_kn(idx)',uvd_true(idx)'));
k = find(t >= 10 & t <= 15); k = k(1:2:end);
writetable(table(t(k)', 1e3*uv_true(k)', 1e3*uv_unk(k)', 1e3*uv_kn(k)', uvd_true(k)', uvd_unk(k)', uvd_kn(k)', ...
    'VariableNames', {'t','uv_true','uv_unk','uv_known','uvd_true','uvd_unk','uvd_known'}), [F 'vehicle_states.csv']);
disp('csv ok')
