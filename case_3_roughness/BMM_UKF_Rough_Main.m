% BMM_UKF_ROUGH_MAIN
% Case 3: robustness of the AUKF mass/velocity estimator to deck-surface
% roughness (Reviewer 1, comment C4).
%
% Ground truth : Euler-Bernoulli beam + moving oscillator (sprung mass on a
%                spring/dashpot over an unsprung wheel mass) crossing a
%                stochastic ISO 8608 deck profile.     -> BMM_SS_rough.m
% Filter model : same oscillator, but the deck roughness is UNKNOWN to the
%                filter (smooth deck assumed) unless ROUGH_KNOWN = true.
%                                                      -> BMM_UKF_rough.m
% Measurement  : noisy mid-span displacement, same noise model as Case 1.
%
% Outputs (PDF + roughness_results.mat):
%   Roughness_profiles.pdf   deck profiles for the ISO classes studied
%   Parameters_rough.pdf     estimated M and v vs. ground truth, all classes
%   Force_rough.pdf          true force (classes A/C); reconstructed force with
%                            roughness unknown and known to the filter (class C)
%   Displacement_rough.pdf   truth vs. noisy measurement (roughest class)

clc; close all; clear;
rng('default'); rng(1234);

% ---------------------------------------------------------------- beam ---
par.n_modes = 4;
par.m  = 200;              % beam mass per unit length [kg/m]
par.EI = 5.3e9;            % flexural rigidity [N m^2]
par.L  = 60;               % span [m]
par.g  = 9.81;

% ------------------------------------------------------------- vehicle ---
M1 = 5000;                 % total vehicle mass [kg]  (estimated)
v1 = 15;                   % speed [m/s]              (estimated)
par.mu_w = 0.10;           % unsprung mass fraction m_w/M (0 -> pure SDOF)
f_v    = 1.5;              % body bounce frequency [Hz]   (typical 1.5-3 Hz)
zeta_v = 0.20;             % suspension damping ratio
m_s      = (1 - par.mu_w) * M1;
par.k_v  = m_s * (2*pi*f_v)^2;            % [N/m]   assumed known to filter
par.c_v  = 2 * zeta_v * m_s * (2*pi*f_v); % [N s/m] assumed known to filter

% ------------------------------------------------------------ roughness ---
classes     = {'A', 'B', 'C'};   % ISO 8608 classes (very good / good / average)
n_min = 0.01;  n_max = 3.0;  N_harm = 300;   % band [cycles/m], harmonics
ROUGH_KNOWN = false;             % true -> filter is given the true profile

% ----------------------------------------------------------- simulation ---
dt    = 0.004;
tspan = [0 30];
n     = par.n_modes;
phi_center = mode_shapesPP(par.L/2, n, par.L)';

% noise model identical to Case 1 (BMM_UKF_Single_Mass_PP.m)
amplitude_factor = 0.01;  offset = 0.01;  frequency = 5;

results = struct([]);
for c = 1:numel(classes)
    fprintf('\n=== ISO 8608 class %s ===\n', classes{c});

    % same random phases for every class -> profiles differ only in amplitude
    prof = iso8608_profile(classes{c}, n_min, n_max, N_harm, 42);
    fprintf('profile RMS = %.2f mm\n', 1e3*prof.sigma_r);

    % --- ground truth ---------------------------------------------------
    par_true = par;  par_true.prof = prof;
    X0 = zeros(2*n+2, 1);
    [t, X] = ode4(@(tt,X) BMM_SS_rough(tt, X, M1, v1, par_true), tspan, X0, dt);
    t = t(:)';
    disp_true = phi_center * X(:,1:n)';
    fc_true   = contact_force_rough(t, X', M1, v1, par_true);

    % --- noisy measurement (Case 1 noise model) --------------------------
    rng(1234);
    Nt = numel(t);
    tau = linspace(0, 2*pi, Nt);
    noise_level = offset + amplitude_factor*abs(disp_true) + amplitude_factor*sin(frequency*tau);
    y = disp_true + noise_level .* randn(1, Nt);

    % --- AUKF -----------------------------------------------------------
    par_est = par;
    if ROUGH_KNOWN, par_est.prof = prof; else, par_est.prof = []; end
    tic
    [X_est, P_saved, R0] = run_UKF_rough(t, y, M1, v1, par_est, dt);
    fprintf('AUKF time %.1f s\n', toc);

    M_est = X_est(2*n+3, :);   v_est = X_est(2*n+4, :);
    sig_M = sqrt(squeeze(P_saved(2*n+3, 2*n+3, :)))';
    sig_v = sqrt(squeeze(P_saved(2*n+4, 2*n+4, :)))';
    fc_est = contact_force_rough(t, X_est(1:2*n+2, :), M_est, v_est, par_est);

    idx = t >= 10;   % steady-state window for summary statistics
    fprintf('M  : mean %.0f kg (true %d),  rel. error %.2f %%\n', ...
        mean(M_est(idx)), M1, 100*(mean(M_est(idx)) - M1)/M1);
    fprintf('v  : mean %.2f m/s (true %d), rel. error %.2f %%\n', ...
        mean(v_est(idx)), v1, 100*(mean(v_est(idx)) - v1)/v1);
    fprintf('f_c: RMS error %.1f %% of M g;  mean-force error %.2f %% (WIM-relevant)\n', ...
        100*rms(fc_est(idx) - fc_true(idx))/(M1*par.g), ...
        100*(mean(fc_est(idx)) - mean(fc_true(idx)))/(M1*par.g));

    results(c).class = classes{c};   results(c).prof = prof;
    results(c).t = t;  results(c).disp_true = disp_true;  results(c).y = y;
    results(c).X_true = X;  results(c).X_est = X_est;  results(c).R0 = R0;
    results(c).M_est = M_est;  results(c).v_est = v_est;
    results(c).sig_M = sig_M;  results(c).sig_v = sig_v;
    results(c).fc_true = fc_true;  results(c).fc_est = fc_est;
end
% --- every class again, now with the profile KNOWN to the filter ----------
%     (e.g. measured by a profilometer).  Two priors are compared:
%       known_paper : manuscript initialisation (0.7 v, P0_v = 10)
%       known_data  : data-driven initial speed with a prior std of 5 % of
%                     v0 (consistent with its accuracy), see initial_speed_estimate.m
%     Metrics: converged M, v errors and the low-frequency (< F_LP Hz)
%     contact-force agreement, since the high-frequency wheel-hop force is
%     not observable from a mid-span displacement sensor.
F_LP = 5;                                    % low-pass cut-off [Hz]
known_paper = struct([]); known_data = struct([]);
fprintf('\n=== profile known to the filter ===\n');
for c = 1:numel(classes)
    r = results(c); t = r.t; idx = t >= 10;
    par_k = par; par_k.prof = r.prof;
    v0 = initial_speed_estimate(t, r.y, par.L);
    for variant = 1:2
        if variant == 1, opts = struct();                                      tag = 'init 0.7v        ';
        else,            opts = struct('init', [0.7 v0/v1], 'P0_v', (0.05*v0)^2); tag = 'data v0, 5% prior'; end
        [X_k, P_k] = run_UKF_rough(t, r.y, M1, v1, par_k, dt, opts);
        k.M_est = X_k(2*n+3,:); k.v_est = X_k(2*n+4,:);
        k.fc_est = contact_force_rough(t, X_k(1:2*n+2,:), k.M_est, k.v_est, par_k);
        k.dM = 100*(mean(k.M_est(idx))-M1)/M1; k.dv = 100*(mean(k.v_est(idx))-v1)/v1;
        a = lowpass_bw(r.fc_true, F_LP, dt); b = lowpass_bw(k.fc_est, F_LP, dt);
        k.corr_lp = corr(b(idx)', a(idx)'); k.rms_lp = 100*rms(b(idx)-a(idx))/(M1*par.g);
        k.rms_full = 100*rms(k.fc_est(idx)-r.fc_true(idx))/(M1*par.g);
        fprintf('class %s, r known, %s : dM %+.2f %%  dv %+.3f %%  corr(<%d Hz) %.2f  RMS(<%d Hz) %.1f %%  RMS(full) %.1f %%\n', ...
            classes{c}, tag, k.dM, k.dv, F_LP, k.corr_lp, F_LP, k.rms_lp, k.rms_full);
        if variant == 1, if isempty(known_paper), known_paper = k; else, known_paper(end+1) = k; end
        else,            if isempty(known_data),  known_data  = k; else, known_data(end+1)  = k; end, end
    end
    % same low-frequency metrics for the roughness-unknown run
    a = lowpass_bw(r.fc_true, F_LP, dt); b = lowpass_bw(r.fc_est, F_LP, dt);
    results(c).corr_lp = corr(b(idx)', a(idx)'); results(c).rms_lp = 100*rms(b(idx)-a(idx))/(M1*par.g);
    fprintf('class %s, r unknown           : corr(<%d Hz) %.2f  RMS(<%d Hz) %.1f %%\n', classes{c}, F_LP, results(c).corr_lp, F_LP, results(c).rms_lp);
end
known = known_data(end);                      % backward compatibility (class C)

save('roughness_results.mat', 'results', 'known', 'known_paper', 'known_data', 'F_LP', 'par', 'M1', 'v1', 'dt', 'ROUGH_KNOWN', '-v7.3');

%% ------------------------------------------------------------- figures ---
plot_rough_figures
