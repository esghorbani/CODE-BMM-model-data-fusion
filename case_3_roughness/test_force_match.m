% Compare reconstructed vs true contact force, classes A and C, roughness
% unknown / known to the filter: overall, per 1-s window, and by frequency band.
clc; clear; load('roughness_results.mat');
n = par.n_modes; Mg = M1*par.g; fs = 10; fn = 'Times New Roman';
sel = [1 numel(results)];                         % class A and class C
lp = @(x, fc) lowpass_bw(x, fc, dt);                     % low-pass, fc in Hz
figure('Position', [100 100 500 320]); k = 0;
for c = sel
    r = results(c); t = r.t; idx = t >= 10;
    par_k = par; par_k.prof = r.prof;
    if c == numel(results), fk_all = known.fc_est;
    else
        Xk = run_UKF_rough(t, r.y, M1, v1, par_k, dt);
        fk_all = contact_force_rough(t, Xk(1:2*n+2,:), Xk(2*n+3,:), Xk(2*n+4,:), par_k);
        fprintf('class %s, r known: dM %+.2f %%, dv %+.3f %%\n', r.class, 100*(mean(Xk(2*n+3,idx))-M1)/M1, 100*(mean(Xk(2*n+4,idx))-v1)/v1);
    end
    ft = r.fc_true; fu = r.fc_est; fk = fk_all;
    fprintf('\n=== class %s (t >= 10 s), std of true dynamic force %.1f %% of Mg ===\n', r.class, 100*std(ft(idx))/Mg);
    for tag = {'unknown', 'known'}
        if strcmp(tag{1}, 'unknown'), fe = fu; else, fe = fk; end
        % per-1-s-window correlation
        cw = []; for w0 = 10:29, ww = t >= w0 & t < w0+1; cw(end+1) = corr(fe(ww)', ft(ww)'); end
        fprintf('  r %-7s : RMS err %5.1f %% | corr overall %.2f, per-second corr min/median/max %.2f/%.2f/%.2f\n', ...
            tag{1}, 100*rms(fe(idx)-ft(idx))/Mg, corr(fe(idx)',ft(idx)'), min(cw), median(cw), max(cw));
        for fc = [2 5 10 20]
            a = lp(ft, fc); b = lp(fe, fc);
            fprintf('      low-passed at %2d Hz: RMS err %5.1f %% of Mg, corr %.2f (true dyn. std in band %.1f %%)\n', ...
                fc, 100*rms(b(idx)-a(idx))/Mg, corr(b(idx)',a(idx)'), 100*std(a(idx))/Mg);
        end
    end
    % zoomed panels
    w = t >= 13 & t <= 13.4;
    k = k+1; subplot(2,2,k); plot(t(w), 1e-3*ft(w), 'k--', 'LineWidth', 1); hold on; plot(t(w), 1e-3*fu(w), 'LineWidth', 1.2)
    title(sprintf('Class %s, r unknown to filter', r.class), 'FontName', fn, 'FontSize', fs, 'FontWeight', 'normal')
    ylabel('Contact force (kN)', 'FontName', fn, 'FontSize', fs); set(gca, 'FontName', fn, 'FontSize', fs)
    if k == 1, legend({'Ground Truth', 'UKF'}, 'FontSize', 7, 'Location', 'best'); end
    k = k+1; subplot(2,2,k); plot(t(w), 1e-3*ft(w), 'k--', 'LineWidth', 1); hold on; plot(t(w), 1e-3*fk(w), 'LineWidth', 1.2)
    title(sprintf('Class %s, r known to filter', r.class), 'FontName', fn, 'FontSize', fs, 'FontWeight', 'normal')
    ylabel('Contact force (kN)', 'FontName', fn, 'FontSize', fs); set(gca, 'FontName', fn, 'FontSize', fs)
    if k == 4, xlabel('Time (sec)', 'FontName', fn, 'FontSize', fs); subplot(2,2,3); xlabel('Time (sec)', 'FontName', fn, 'FontSize', fs); end
end
exportgraphics(gcf, 'Force_rough_zoom.pdf', 'ContentType', 'vector', 'BackgroundColor', 'none');

function y = lowpass_bw(x, fc, dt)
[b, a] = butter(4, fc/(0.5/dt));
y = filtfilt(b, a, x);
end
