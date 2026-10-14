% Does knowing the profile change the parameter estimates?  All classes,
% paper initialisation vs data-driven v0 with a consistent (tighter) prior.
clc; clear; load('roughness_results.mat');
n = par.n_modes; Mg = M1*par.g;
lp = @(x, fc) lowpass_bw(x, fc, dt);
fprintf('%-8s %-28s %8s %8s %10s %12s\n', 'class', 'filter setting', 'dM [%]', 'dv [%]', 'corr<5Hz', 'RMS<5Hz [%Mg]');
for c = 1:numel(results)
    r = results(c); t = r.t; idx = t >= 10;
    v0 = initial_speed_estimate(t, r.y, par.L);
    settings = {
      'r unknown, init 0.7v',            [],     struct();
      'r known, init 0.7v',              r.prof, struct();
      'r unknown, data v0, P0v 0.1',     [],     struct('init', [0.7 v0/v1], 'P0_v', 0.1);
      'r known, data v0, P0v 0.1',       r.prof, struct('init', [0.7 v0/v1], 'P0_v', 0.1);
    };
    for k = 1:size(settings,1)
        pe = par; pe.prof = settings{k,2};
        Xe = run_UKF_rough(t, r.y, M1, v1, pe, dt, settings{k,3});
        fe = contact_force_rough(t, Xe(1:2*n+2,:), Xe(2*n+3,:), Xe(2*n+4,:), pe);
        a = lp(r.fc_true, 5); b = lp(fe, 5);
        fprintf('%-8s %-28s %+8.2f %+8.3f %10.2f %12.1f\n', r.class, settings{k,1}, ...
            100*(mean(Xe(2*n+3,idx))-M1)/M1, 100*(mean(Xe(2*n+4,idx))-v1)/v1, corr(b(idx)',a(idx)'), 100*rms(b(idx)-a(idx))/Mg);
    end
end
function y = lowpass_bw(x, fc, dt)
[b, a] = butter(4, fc/(0.5/dt)); y = filtfilt(b, a, x);
end
