% BMM_UKF_HIGHSPEED_ROUGH  Reviewer 2, C4: quarter-car on an ISO 8608 class C deck at
% higher speeds and masses; roughness NOT included in the filter model; data-driven
% initial speed with a 5 % prior.  Exports the first-crossing mid-span displacement
% (M = 5000 kg) and a table of converged errors for all (v, M) cases.
clc; clear; rng('default'); rng(1234);
par.n_modes = 4; par.m = 200; par.EI = 5.3e9; par.L = 60; par.g = 9.81; par.mu_w = 0.10;
f_v = 1.5; zeta_v = 0.20; dt = 0.004; tspan = [0 30]; n = par.n_modes;
phi_center = mode_shapesPP(par.L/2, n, par.L)';
prof = iso8608_profile('C', 0.01, 3, 300, 42);
speeds = [20 30 45 60 80]; masses = [5000 10000 15000];
amplitude_factor = 0.01; offset = 0.01; frequency = 5;
F = '../../PAPER-2026-BMM-FUSION/JSV_Submission/rebuttal/2026-09-11-round-1/figures/';
w1 = sqrt(par.EI/par.m)*(pi/par.L)^2; vcr = w1*par.L/pi;
xi = linspace(0, 1, 301)'; Dd = xi; dnames = {'xi'};
res = [];
fprintf('%6s %7s %6s %6s %8s %8s %8s %8s\n', 'v', 'M', 'v/vcr', 'DAF', 'v0 err%', 'dM %', 'dv %', 'fdyn %');
for M1 = masses
    m_s = (1-par.mu_w)*M1; par.k_v = m_s*(2*pi*f_v)^2; par.c_v = 2*zeta_v*m_s*(2*pi*f_v);
    for v1 = speeds
        par_true = par; par_true.prof = prof;
        [t, X] = ode4(@(tt,X) BMM_SS_rough(tt, X, M1, v1, par_true), tspan, zeros(2*n+2,1), dt); t = t(:)';
        d = phi_center*X(:,1:n)'; fc = contact_force_rough(t, X', M1, v1, par_true);
        rng(1234); Nt = numel(t); tau = linspace(0,2*pi,Nt);
        y = d + (offset + amplitude_factor*abs(d) + amplitude_factor*sin(frequency*tau)).*randn(1,Nt);
        v0 = initial_speed_estimate(t, y, par.L);
        par_est = par; par_est.prof = [];
        Xe = run_UKF_rough(t, y, M1, v1, par_est, dt, struct('init', [0.7 v0/v1], 'P0_v', (0.05*v0)^2));
        idx = t >= 10; Me = Xe(2*n+3,:); ve = Xe(2*n+4,:);
        d_static = M1*par.g*par.L^3/(48*par.EI);
        r = struct('v', v1, 'M', M1, 'vvcr', v1/vcr, 'DAF', max(d(idx))/d_static, 'v0err', 100*(v0-v1)/v1, ...
            'dM', 100*(mean(Me(idx))-M1)/M1, 'dv', 100*(mean(ve(idx))-v1)/v1, 'fdyn', 100*max(abs(fc(idx)-M1*par.g))/(M1*par.g), ...
            'Mhat', Me, 'vhat', ve, 't', t);
        if isempty(res), res = r; else, res(end+1) = r; end
        fprintf('%6d %7d %6.2f %6.2f %8.1f %8.2f %8.3f %8.1f\n', v1, M1, r.vvcr, r.DAF, r.v0err, r.dM, r.dv, r.fdyn);
        if M1 == 5000
            nc = round(par.L/v1/dt) + 1;
            Dd = [Dd, 1e3*interp1(linspace(0,1,nc)', d(1:nc)', xi)]; dnames{end+1} = sprintf('d%d', v1);
        end
    end
end
writetable(array2table(Dd, 'VariableNames', dnames), [F 'highspeed_disp.csv']);
save('highspeed_rough_results.mat', 'res', 'par', 'prof', 'speeds', 'masses', 'dt', '-v7.3');
disp('done')
