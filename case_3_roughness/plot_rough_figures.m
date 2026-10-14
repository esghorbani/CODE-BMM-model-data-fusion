% PLOT_ROUGH_FIGURES  Figures of BMM_UKF_Rough_Main (standalone).
% Loads roughness_results.mat when the study variables are not in the workspace.
if ~exist('results', 'var') || ~exist('known_data', 'var')
    load('roughness_results.mat');
end
n = par.n_modes; classes = {results.class}; t = results(1).t; cR = numel(classes);
phi_center = mode_shapesPP(par.L/2, n, par.L)';
fs = 10; fn = 'Times New Roman';
cols = lines(numel(classes));

% roughness profiles
xg = linspace(0, par.L, 3001);
figure('Position', [100 100 500 170])
hold on
for c = 1:numel(classes)
    plot(xg, 1e3*eval_roughness(results(c).prof, xg), 'Color', cols(c,:), 'LineWidth', 1)
end
legend(strcat('ISO 8608 class ', classes), 'Location', 'best', 'NumColumns', 3)
xlabel('x (m)', 'FontName', fn, 'FontSize', fs)
ylabel('r(x) (mm)', 'FontName', fn, 'FontSize', fs)
set(gca, 'FontName', fn, 'FontSize', fs); xlim([0 par.L])
exportgraphics(gcf, 'Roughness_profiles.pdf', 'ContentType', 'vector', 'BackgroundColor', 'none');

% parameter estimates
figure('Position', [100 100 500 170])
subplot(1,2,1); hold on
for c = 1:numel(classes)
    plot(results(c).t, results(c).M_est, 'Color', cols(c,:), 'LineWidth', 1.2)
end
plot(t, M1*ones(size(t)), 'k--', 'LineWidth', 1.2)
legend([strcat('Class ', classes), {'Ground Truth'}], 'Location', 'best')
ylabel('Moving mass (kg)', 'FontName', fn, 'FontSize', fs)
xlabel('Time (sec)', 'FontName', fn, 'FontSize', fs)
ylim([0.5*M1 1.2*M1]); xlim([0 15]); set(gca, 'FontName', fn, 'FontSize', fs)
subplot(1,2,2); hold on
for c = 1:numel(classes)
    plot(results(c).t, results(c).v_est, 'Color', cols(c,:), 'LineWidth', 1.2)
end
plot(t, v1*ones(size(t)), 'k--', 'LineWidth', 1.2)
ylabel('Velocity (m/s)', 'FontName', fn, 'FontSize', fs)
xlabel('Time (sec)', 'FontName', fn, 'FontSize', fs)
ylim([0 1.4*v1]); xlim([0 15]); set(gca, 'FontName', fn, 'FontSize', fs)
exportgraphics(gcf, 'Parameters_rough.pdf', 'ContentType', 'vector', 'BackgroundColor', 'none');

% contact force, one passage after convergence
figure('Position', [100 100 750 170])
subplot(1,3,1)   % true force: roughest class drawn first so class A stays visible
plot(t, 1e-3*results(cR).fc_true, 'Color', cols(cR,:), 'LineWidth', 0.8); hold on
plot(t, 1e-3*results(1).fc_true, 'Color', cols(1,:), 'LineWidth', 1)
plot(t, 1e-3*M1*par.g*ones(size(t)), 'k--', 'LineWidth', 1)
legend({['Class ' classes{cR}], ['Class ' classes{1}], 'M g'}, 'Location', 'best')
ylabel('Contact force (kN)', 'FontName', fn, 'FontSize', fs)
xlabel('Time (sec)', 'FontName', fn, 'FontSize', fs)
title('(a) True force', 'FontName', fn, 'FontSize', fs, 'FontWeight', 'normal')
xlim([12 16]); set(gca, 'FontName', fn, 'FontSize', fs)
ft_lp = lowpass_bw(results(cR).fc_true, F_LP, dt);
subplot(1,3,2)   % filter assumes smooth deck, low-passed comparison
plot(t, 1e-3*ft_lp, '--', 'Color', [0.85 0.33 0.1], 'LineWidth', 1); hold on
plot(t, 1e-3*lowpass_bw(results(cR).fc_est, F_LP, dt), 'Color', [0 0.45 0.74], 'LineWidth', 1.2)
legend('Ground Truth', 'UKF, r unknown', 'Location', 'best')
ylabel(sprintf('Contact force < %d Hz (kN)', F_LP), 'FontName', fn, 'FontSize', fs)
xlabel('Time (sec)', 'FontName', fn, 'FontSize', fs)
title(['(b) Class ' classes{cR} ', r unknown to filter'], 'FontName', fn, 'FontSize', fs, 'FontWeight', 'normal')
xlim([12 16]); set(gca, 'FontName', fn, 'FontSize', fs)
subplot(1,3,3)   % filter knows the profile (data-driven v0, tight prior)
plot(t, 1e-3*ft_lp, '--', 'Color', [0.85 0.33 0.1], 'LineWidth', 1); hold on
plot(t, 1e-3*lowpass_bw(known_data(cR).fc_est, F_LP, dt), 'Color', [0 0.45 0.74], 'LineWidth', 1.2)
legend('Ground Truth', 'UKF, r known', 'Location', 'best')
ylabel(sprintf('Contact force < %d Hz (kN)', F_LP), 'FontName', fn, 'FontSize', fs)
xlabel('Time (sec)', 'FontName', fn, 'FontSize', fs)
title(['(c) Class ' classes{cR} ', r known to filter'], 'FontName', fn, 'FontSize', fs, 'FontWeight', 'normal')
xlim([12 16]); set(gca, 'FontName', fn, 'FontSize', fs)
exportgraphics(gcf, 'Force_rough.pdf', 'ContentType', 'vector', 'BackgroundColor', 'none');

% displacement, roughest class
figure('Position', [100 100 500 170])
subplot(1,2,1)
plot(t, results(cR).disp_true, 'LineWidth', 1); hold on
plot(t, results(cR).y, ':', 'LineWidth', 0.9)
legend('Ground Truth', 'Noisy Measurement', 'Location', 'best')
ylabel('Displacement (m)', 'FontName', fn, 'FontSize', fs)
xlabel('Time (s)', 'FontName', fn, 'FontSize', fs)
xlim([0 10]); set(gca, 'FontName', fn, 'FontSize', fs)
subplot(1,2,2)
plot(t, phi_center*results(cR).X_est(1:n,:), 'LineWidth', 1.2); hold on
plot(t, results(cR).disp_true, '--', 'LineWidth', 1.2)
legend('UKF estimated', 'Ground Truth', 'Location', 'best')
ylabel('Displacement (m)', 'FontName', fn, 'FontSize', fs)
xlabel('Time (s)', 'FontName', fn, 'FontSize', fs)
xlim([10 15]); set(gca, 'FontName', fn, 'FontSize', fs)
exportgraphics(gcf, 'Displacement_rough.pdf', 'ContentType', 'vector', 'BackgroundColor', 'none');

% zoom on the full-band force, classes A and C, r unknown / known
figure('Position', [100 100 500 320]); w = t >= 13 & t <= 13.4; k = 0; lab = 'abcd';
for c = [1 cR]
    for known_flag = [false true]
        k = k + 1; subplot(2,2,k)
        if known_flag, fe = known_data(c).fc_est; tag = 'r known to filter'; else, fe = results(c).fc_est; tag = 'r unknown to filter'; end
        plot(t(w), 1e-3*results(c).fc_true(w), 'k--', 'LineWidth', 1); hold on
        plot(t(w), 1e-3*fe(w), 'LineWidth', 1.2)
        title(sprintf('(%s) Class %s, %s', lab(k), classes{c}, tag), 'FontName', fn, 'FontSize', fs, 'FontWeight', 'normal')
        ylabel('Contact force (kN)', 'FontName', fn, 'FontSize', fs); set(gca, 'FontName', fn, 'FontSize', fs)
        if k == 1, legend({'Ground Truth', 'UKF'}, 'FontSize', 7, 'Location', 'best'); end
        if k > 2, xlabel('Time (sec)', 'FontName', fn, 'FontSize', fs); end
    end
end
exportgraphics(gcf, 'Force_rough_zoom.pdf', 'ContentType', 'vector', 'BackgroundColor', 'none');
