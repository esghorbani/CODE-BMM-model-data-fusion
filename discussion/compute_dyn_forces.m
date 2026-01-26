


clc
clear
close all

load('BMM_UKF_dataset_with_force_org.mat')

k = 76;  % <<< SELECT CASE INDEX

% -------------------------------
% Extract force data
% -------------------------------
tF = dataset(k).force_time;
Fd = dataset(k).force_dynamic_mean;
Fs = dataset(k).force_static_est;

% -------------------------------
% Extract state estimates
% -------------------------------
t    = dataset(k).time;
Xest = dataset(k).X_est;
P_saved=dataset(k).P_est;

M_true = dataset(k).mass;
V_true = dataset(k).velocity;

% infer n_modes safely
St_Len  = size(Xest,1);
n_modes = (St_Len - 2)/2;

M_est = Xest(2*n_modes+1, :);
V_est = Xest(2*n_modes+2, :);

% ===============================
% Extract covariance diagonals
% ===============================
Nt = length(t);

sigma_M = zeros(1,Nt);
sigma_V = zeros(1,Nt);

for i = 1:Nt-1
    sigma_M(i) = sqrt(P_saved(2*n_modes+1, 2*n_modes+1, i));
    sigma_V(i) = sqrt(P_saved(2*n_modes+2, 2*n_modes+2, i));
end


% -------------------------------
% Plot
% -------------------------------
figure('Position',[100 100 900 280])

% ===== Force =====
subplot(1,3,1)
plot(tF, Fd, 'b', 'LineWidth', 1.6); hold on
plot(tF, Fs, '--r', 'LineWidth', 1.6)
xlabel('Time (s)')
ylabel('Force (N)')
title('Dynamic vs static force')
legend('\mu_F(t)', 'M_{est} g', 'Location','best')
grid on
xlim([0 15])

% ===== Mass =====
subplot(1,3,2)
plot(t, M_est, 'b', 'LineWidth', 1.6); hold on
plot(t, M_true*ones(size(t)), '--r', 'LineWidth', 1.6)

tf=t(2:end);
t = t(:).';  
% ±2σ band
fill([t fliplr(t)], ...
     [M_est + 2*sigma_M, fliplr(M_est - 2*sigma_M)], ...
     'b', 'FaceAlpha',0.2, 'EdgeColor','none')

xlabel('Time (s)')
ylabel('Mass (kg)')
title('Estimated mass')
legend('Estimated','Ground truth','±2σ','Location','best')
grid on
xlim([0 15])


% ===== Velocity =====
subplot(1,3,3)
plot(t, V_est, 'b', 'LineWidth', 1.6); hold on
plot(t, V_true*ones(size(t)), '--r', 'LineWidth', 1.6)

% ±2σ band
fill([t fliplr(t)], ...
     [V_est + 2*sigma_V, fliplr(V_est - 2*sigma_V)], ...
     'b', 'FaceAlpha',0.2, 'EdgeColor','none')

xlabel('Time (s)')
ylabel('Velocity (m/s)')
title('Estimated velocity')
legend('Estimated','Ground truth','±2σ','Location','best')
grid on
xlim([0 15])


% -------------------------------
% Global title
% -------------------------------
sgtitle(sprintf('Case %d | M = %.0f kg | v = %.1f m/s', ...
        k, dataset(k).mass, dataset(k).velocity))

%%
nCases = length(dataset);

Fdyn_all = zeros(1,nCases);
Fstat_all = zeros(1,nCases);

for k = 1:nCases
    tF = dataset(k).force_time;
    Fd = dataset(k).force_dynamic_mean;
    Fs = dataset(k).force_static_est;

    idx = tF >= (tF(end) - 3);

    Fdyn_all(k)  = mean(Fd(idx));
    Fstat_all(k) = mean(Fs(idx));
end

figure('Position',[100 100 420 380])
scatter(Fstat_all, Fdyn_all, 40, 'filled'); hold on
plot([min(Fstat_all) max(Fstat_all)], ...
     [min(Fstat_all) max(Fstat_all)], 'k--', 'LineWidth',1.5)

xlabel('Average static force  M g  (N)')
ylabel('Average dynamic force  <F_{dyn}>  (N)')
title('Average dynamic force vs static force')
grid on
axis equal
%%
DAF = Fdyn_all ./ Fstat_all;

figure('Position',[100 100 420 280])
histogram(DAF, 20)
xlabel('Dynamic Amplification Factor  <F_{dyn}> / (Mg)')
ylabel('Count')
grid on

fprintf('DAF mean = %.4f\n', mean(DAF));
fprintf('DAF std  = %.4f\n', std(DAF));

%%
clc
clear
close all

load('BMM_UKF_dataset_with_force_org.mat')   % or BMM_UKF_dataset.mat

nCases = length(dataset);

M_est_all = zeros(1,nCases);
V_est_all = zeros(1,nCases);
M_true_all = zeros(1,nCases);
V_true_all = zeros(1,nCases);

% -------------------------------
% Extract steady-state estimates
% -------------------------------
for k = 1:nCases

    t    = dataset(k).time;
    Xest = dataset(k).X_est;

    % infer n_modes
    St_Len  = size(Xest,1);
    n_modes = (St_Len - 2)/2;

    M_est = Xest(2*n_modes+1, :);
    V_est = Xest(2*n_modes+2, :);

    % steady-state window (last 3 seconds)
    t_ss = 3;
    idx  = t >= (t(end) - t_ss);

    M_est_all(k) = mean(M_est(idx));
    V_est_all(k) = mean(V_est(idx));

    M_true_all(k) = dataset(k).mass;
    V_true_all(k) = dataset(k).velocity;
end

case_id = 1:nCases;

% -------------------------------
% Plot
% -------------------------------
figure('Position',[100 100 850 280])

% ===== Mass =====
subplot(1,2,1)
plot(case_id, M_est_all, 'bo', 'MarkerSize',5, 'LineWidth',1.2); hold on
plot(case_id, M_true_all, 'r--', 'LineWidth',1.6)
xlabel('Case number')
ylabel('Mass (kg)')
% title('Estimated mass vs ground truth (all cases)')
legend('Estimated','Ground truth','Location','best')
grid on

% ===== Velocity =====
subplot(1,2,2)
plot(case_id, V_est_all, 'bo', 'MarkerSize',5, 'LineWidth',1.2); hold on
plot(case_id, V_true_all, 'r--', 'LineWidth',1.6)
xlabel('Case number')
ylabel('Velocity (m/s)')
% title('Estimated velocity vs ground truth (all cases)')
legend('Estimated','Ground truth','Location','best')
grid on

% sgtitle('UKF parameter estimation accuracy across all cases')

% %%
% figure('Position',[100 100 900 280])
% 
% % ===== Force =====
% subplot(1,2,1)
% plot(tF, Fd, 'b', 'LineWidth', 1.6); hold on
% plot(tF, Fs, '--r', 'LineWidth', 1.6)
% xlabel('Time (s)')
% ylabel('Force (N)')
% % title('Dynamic vs static force')
% legend('dynamic force', 'static force', 'Location','best')
% grid on
% xlim([0 15])
% 
% subplot(1,2,2)
% scatter(Fstat_all, Fdyn_all, 40, 'filled'); hold on
% plot([min(Fstat_all) max(Fstat_all)], ...
%      [min(Fstat_all) max(Fstat_all)], 'k--', 'LineWidth',1.5)
% 
% xlabel('Average static force (N)')
% ylabel('Average dynamic force (N)')
% % title('Average dynamic force vs static force')
% grid on
% axis equal

%%
clc
clear
close all
load('BMM_UKF_dataset_with_force_org.mat')

% <<< SELECT TWO CASE INDICES >>>
k1 = 41;
k2 = 128;  % choose another case with different mass/velocity

% ===============================
% CASE 1
% ===============================
t1    = dataset(k1).time;
Xest1 = dataset(k1).X_est;
P_saved1 = dataset(k1).P_est;
M_true1 = dataset(k1).mass;
V_true1 = dataset(k1).velocity;

St_Len1  = size(Xest1,1);
n_modes1 = (St_Len1 - 2)/2;
M_est1 = Xest1(2*n_modes1+1, :);
V_est1 = Xest1(2*n_modes1+2, :);

% Align sizes safely
Nt1 = min([length(t1), size(Xest1,2), size(P_saved1,3)]);
t1 = t1(1:Nt1);
M_est1 = M_est1(1:Nt1);
V_est1 = V_est1(1:Nt1);

sigma_M1 = zeros(1,Nt1);
sigma_V1 = zeros(1,Nt1);
for i = 1:Nt1
    sigma_M1(i) = sqrt(P_saved1(2*n_modes1+1, 2*n_modes1+1, i));
    sigma_V1(i) = sqrt(P_saved1(2*n_modes1+2, 2*n_modes1+2, i));
end
t1 = t1(:).';

% ===============================
% CASE 2
% ===============================
t2    = dataset(k2).time;
Xest2 = dataset(k2).X_est;
P_saved2 = dataset(k2).P_est;
M_true2 = dataset(k2).mass;
V_true2 = dataset(k2).velocity;

St_Len2  = size(Xest2,1);
n_modes2 = (St_Len2 - 2)/2;
M_est2 = Xest2(2*n_modes2+1, :);
V_est2 = Xest2(2*n_modes2+2, :);

% Align sizes safely
Nt2 = min([length(t2), size(Xest2,2), size(P_saved2,3)]);
t2 = t2(1:Nt2);
M_est2 = M_est2(1:Nt2);
V_est2 = V_est2(1:Nt2);

sigma_M2 = zeros(1,Nt2);
sigma_V2 = zeros(1,Nt2);
for i = 1:Nt2
    sigma_M2(i) = sqrt(P_saved2(2*n_modes2+1, 2*n_modes2+1, i));
    sigma_V2(i) = sqrt(P_saved2(2*n_modes2+2, 2*n_modes2+2, i));
end
t2 = t2(:).';

% ===============================
% Define colors
% ===============================
color1 = [0 0.4470 0.7410];       % blue
color2 = [0.8500 0.3250 0.0980];  % orange

% ===============================
% Plot
% ===============================
figure('Position',[100 100 700 300])

% ===== Mass =====
subplot(1,2,1)
hold on

% Case 1
fill([t1 fliplr(t1)], [M_est1 + 2*sigma_M1, fliplr(M_est1 - 2*sigma_M1)], ...
    color1, 'FaceAlpha',0.15, 'EdgeColor','none')
p1 = plot(t1, M_est1, 'Color', color1, 'LineWidth', 1.6);
plot(t1, M_true1*ones(size(t1)), '--', 'Color', color1, 'LineWidth', 1.4)

% Case 2
fill([t2 fliplr(t2)], [M_est2 + 2*sigma_M2, fliplr(M_est2 - 2*sigma_M2)], ...
    color2, 'FaceAlpha',0.15, 'EdgeColor','none')
p2 = plot(t2, M_est2, 'Color', color2, 'LineWidth', 1.6);
plot(t2, M_true2*ones(size(t2)), '--', 'Color', color2, 'LineWidth', 1.4)

xlabel('Time (s)')
ylabel('Mass (kg)')
legend([p1 p2], {sprintf('Case %d (M=%.0f kg)', k1, M_true1), ...
                 sprintf('Case %d (M=%.0f kg)', k2, M_true2)}, ...
       'Location','best')
grid on
xlim([0 15])
box on

% ===== Velocity =====
subplot(1,2,2)
hold on

% Case 1
fill([t1 fliplr(t1)], [V_est1 + 2*sigma_V1, fliplr(V_est1 - 2*sigma_V1)], ...
    color1, 'FaceAlpha',0.15, 'EdgeColor','none')
p1 = plot(t1, V_est1, 'Color', color1, 'LineWidth', 1.6);
plot(t1, V_true1*ones(size(t1)), '--', 'Color', color1, 'LineWidth', 1.4)

% Case 2
fill([t2 fliplr(t2)], [V_est2 + 2*sigma_V2, fliplr(V_est2 - 2*sigma_V2)], ...
    color2, 'FaceAlpha',0.15, 'EdgeColor','none')
p2 = plot(t2, V_est2, 'Color', color2, 'LineWidth', 1.6);
plot(t2, V_true2*ones(size(t2)), '--', 'Color', color2, 'LineWidth', 1.4)

xlabel('Time (s)')
ylabel('Velocity (m/s)')
legend([p1 p2], {sprintf('Case %d (v=%.1f m/s)', k1, V_true1), ...
                 sprintf('Case %d (v=%.1f m/s)', k2, V_true2)}, ...
       'Location','best')
grid on
xlim([0 15])
box on
%%
% If you have innovations stored, or compute them:
% NIS = innovation' * inv(S) * innovation
% where S is the innovation covariance

% If innovations are saved in dataset:
% innovation = dataset(k).innovation;
% S = dataset(k).S;

% % Otherwise, approximate using estimation error (if ground truth available)
% M_true = dataset(k).mass;
% V_true = dataset(k).velocity;
% M_est = Xest(2*n_modes+1, 1:Nt);
% V_est = Xest(2*n_modes+2, 1:Nt);
% 
% % Normalized estimation error squared (NEES) - proxy for filter consistency
% NEES_M = ((M_est - M_true) ./ sigma_M).^2;
% NEES_V = ((V_est - V_true) ./ sigma_V).^2;
% 
% figure('Position', [100 100 800 400])
% 
% subplot(2,1,1)
% plot(t, NEES_M, 'b', 'LineWidth', 1.2)
% hold on
% yline(chi2inv(0.95, 1), '--r', '95% bound', 'LineWidth', 1.5)
% ylabel('NEES_M')
% title('Normalized Estimation Error Squared - Mass')
% grid on
% xlim([0 15])
% ylim([0 10])
% 
% subplot(2,1,2)
% plot(t, NEES_V, 'r', 'LineWidth', 1.2)
% hold on
% yline(chi2inv(0.95, 1), '--r', '95% bound', 'LineWidth', 1.5)
% ylabel('NEES_V')
% xlabel('Time (s)')
% title('Normalized Estimation Error Squared - Velocity')
% grid on
% xlim([0 15])
% ylim([0 10])
% 
% sgtitle('Filter Consistency Check (NEES)')
%%

