% BMM_UKF_HIGHSPEED_BANK  R2-C4 cases (quarter-car, class C deck, roughness not in
% the filter) with the FILTER BANK initialisation instead of the data-driven v0:
% 0.7 M start, initial speeds 12:10:92 m/s (none equal to a true speed), selection
% by innovation log-likelihood over the first 10 s.
clc; clear; rng('default'); rng(1234);
diary('highspeed_bank_log.txt'); diary on
par.n_modes = 4; par.m = 200; par.EI = 5.3e9; par.L = 60; par.g = 9.81; par.mu_w = 0.10;
f_v = 1.5; zeta_v = 0.20; dt = 0.004; tspan = [0 30]; n = par.n_modes;
phi_center = mode_shapesPP(par.L/2, n, par.L)';
prof = iso8608_profile('C', 0.01, 3, 300, 42);
cases = [30 5000; 60 5000; 80 5000; 60 10000; 80 10000];
amplitude_factor = 0.01; offset = 0.01; frequency = 5; bank = 12:10:92;
fprintf('%6s %7s | %8s %8s | %6s\n','v','M','dM %','dv %','vsel');
res = struct([]);
for c = 1:size(cases,1)
    v1 = cases(c,1); M1 = cases(c,2);
    m_s = (1-par.mu_w)*M1; par.k_v = m_s*(2*pi*f_v)^2; par.c_v = 2*zeta_v*m_s*(2*pi*f_v);
    par_true = par; par_true.prof = prof;
    [t, X] = ode4(@(tt,X) BMM_SS_rough(tt, X, M1, v1, par_true), tspan, zeros(2*n+2,1), dt); t = t(:)';
    d = phi_center*X(:,1:n)'; rng(1234); Nt = numel(t); tau = linspace(0,2*pi,Nt);
    y = d + (offset + amplitude_factor*abs(d) + amplitude_factor*sin(frequency*tau)).*randn(1,Nt);
    par_est = par; par_est.prof = []; idx = t >= 10;
    LLs = zeros(size(bank)); dMs = LLs; dvs = LLs; Xs = cell(size(bank));
    for b = 1:numel(bank)
        [Xe,~,~,LL] = run_UKF_rough(t, y, M1, v1, par_est, dt, struct('v_init', bank(b))); Xs{b} = Xe;
        LLs(b) = sum(LL(t>0 & t<=10)); dMs(b) = 100*(mean(Xe(2*n+3,idx))-M1)/M1; dvs(b) = 100*(mean(Xe(2*n+4,idx))-v1)/v1;
    end
    [~,bs] = max(LLs);
    fprintf('%6d %7d | %8.2f %8.2f | %6d\n', v1, M1, dMs(bs), dvs(bs), bank(bs));
    fprintf('   bank dM: %s\n', mat2str(round(dMs,1)));
    res(c).v = v1; res(c).M = M1; res(c).dM = dMs(bs); res(c).dv = dvs(bs); res(c).sel = bank(bs); res(c).dMs = dMs; res(c).LLs = LLs; res(c).Mhat = Xs{bs}(2*n+3,:); res(c).vhat = Xs{bs}(2*n+4,:); res(c).t = t;
end
save('highspeed_bank_results.mat','res','bank'); fprintf('DONE\n'); diary off
