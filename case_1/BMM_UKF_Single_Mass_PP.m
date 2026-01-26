
clc
close all
clear all

% Set the random seed for reproducibility
rng('default');  % Resets the generator to the modern configuration
rng(1234);  % You can choose any seed number you like


global n_modes m EI L ;   % Global variables

n_modes=4;
M1 = 5000;                % First moving mass (assumed zero as in the original code)
m = 200;                % Beam unit mass
EI = 5.3*10^9;         % Young's modulus times moment of inertia
L = 60;                % Beam length
v1 = 15;               % Velocity of the first moving mass
tspan = [0 30];  % Time span for the simulation
X0 = 0*[ones(1,n_modes)  zeros(1,n_modes)];  % Initial state [x1; x2; ...; xdn]
% [t, X] = ode45(@(tt, X) beam_moving_mass_state_space(tt, X, n_modes,M1,M2,m,EI,L,v1,v2), tspan, X0);
dt=0.004;
% Solve the ODE using the custom ode4 solver with a fixed time step
[t, X] = ode4(@(tt, X) BMM_SS(tt, X, n_modes, M1,  m, EI, L, v1), tspan, X0, dt);


x_center = L / 2;
phi_center = mode_shapesPP(x_center, n_modes, L)';
% Calculate the total displacement at the center of the beam over time

% Plot the total displacement at the center of the beam over time
x = X(:,1:n_modes)';
total_displacement_center = (phi_center* x);
xd = X(:,n_modes+1:end)';

total_vel_center = (phi_center* xd);


FFT_signal=fft(total_displacement_center);
Fs=1/dt;
freq = (0:Fs/length(total_displacement_center):Fs/2);







% Define the noise level (relative to the signal amplitude)
num_steps = length(total_displacement_center);
time = linspace(0, 2 * pi, num_steps);  % Example time vector from 0 to 2*pi

% Define the parameters for time-varying noise level
amplitude_factor = 0.01;  % This factor controls the proportion of noise relative to the signal
offset = 0.01;           % Offset to ensure a minimum noise level
frequency = 5;            % Frequency of the sine wave

% Initialize the noisy displacement vector
total_displacement_center_noisy = zeros(size(total_displacement_center));

% Add time-varying noise proportional to the signal amplitude
for i = 1:num_steps
    % Calculate the noise level based on the signal amplitude and a sinusoidal variation
    
    noise_level(i) = (offset + amplitude_factor * abs(total_displacement_center(i)) + ...
                  amplitude_factor * sin(frequency * time(i)));
    
    % Add the noise to the signal
    total_displacement_center_noisy(i) = total_displacement_center(i) + noise_level(i)* randn() ;
end

% Add white noise to the total_displacement_center
% total_displacement_center_noisy = total_displacement_center + noise_level * randn(size(total_displacement_center));

% Create figure and set size

% Define color palettes
colors1 = [ ... % 'dark' palette
    0.2081, 0.1663, 0.5292; ...
    0.2617, 0.5725, 0.6494; ...
    0.3692, 0.7887, 0.3822; ...
    0.8686, 0.9517, 0.3667; ...
    0.9945, 0.7700, 0.2625; ...
    0.9626, 0.4795, 0.1954];
colors = [ ... % 'crest' palette
    0.212, 0.318, 0.368; ...
    0.343, 0.522, 0.553; ...
    0.518, 0.712, 0.686; ...
    0.674, 0.822, 0.767; ...
    0.839, 0.917, 0.878; ...
    0.957, 0.984, 0.941];
colors2 = [ ... % 'plasma' palette
    0.0504, 0.0298, 0.5271; ...
    0.2445, 0.0159, 0.6685; ...
    0.5365, 0.0894, 0.6397; ...
    0.7938, 0.3058, 0.4824; ...
    0.9574, 0.5714, 0.3155; ...
    0.9910, 0.8849, 0.0993];

figure('Position', [100, 100, 900, 300])  % [left, bottom, width, height]

% First subplot for displacement
subplot(2,2,1)
plot(t, total_displacement_center, 'Color', colors1(2, :), 'LineWidth', 1.5)
hold on
plot(t, total_displacement_center_noisy, 'Color', colors1(1, :), 'LineWidth', 1)
legend('Ground Truth', 'Noisy Measurement', 'Location', 'Best')
title('Total Displacement at the Center of the Beam', 'FontName', 'Times New Roman', 'FontSize', 10)
ylabel('Displacement (m)', 'FontName', 'Times New Roman', 'FontSize', 10)
xlabel('Time (s)', 'FontName', 'Times New Roman', 'FontSize', 10)
set(gca, 'FontName', 'Times New Roman', 'FontSize', 10)
xlim([0,15])
% Second subplot for PSD
subplot(2,2,2)
plot(freq(1:2000), log10(abs(FFT_signal(1:2000))), 'LineWidth', 1.5)
title('Total Displacement at the Center of the Beam', 'FontName', 'Times New Roman', 'FontSize', 10)
ylabel('PSD', 'FontName', 'Times New Roman', 'FontSize', 10)
xlabel('Frequency (rad/sec)', 'FontName', 'Times New Roman', 'FontSize', 10)
set(gca, 'FontName', 'Times New Roman', 'FontSize', 10)

subplot(2,2,3)
plot(freq(1:1000), log10(abs(FFT_signal(1:1000))), 'LineWidth', 1.5)
title('Total Displacement at the Center of the Beam', 'FontName', 'Times New Roman', 'FontSize', 10)
ylabel('PSD', 'FontName', 'Times New Roman', 'FontSize', 10)
xlabel('Frequency (rad/sec)', 'FontName', 'Times New Roman', 'FontSize', 10)
set(gca, 'FontName', 'Times New Roman', 'FontSize', 10)
xlim([1,5])

subplot(2,2,4)
plot(freq(1:1000), log10(abs(FFT_signal(1:1000))), 'LineWidth', 1.5)
title('Total Displacement at the Center of the Beam', 'FontName', 'Times New Roman', 'FontSize', 10)
ylabel('PSD', 'FontName', 'Times New Roman', 'FontSize', 10)
xlabel('Frequency (rad/sec)', 'FontName', 'Times New Roman', 'FontSize', 10)
set(gca, 'FontName', 'Times New Roman', 'FontSize', 10)
xlim([17,22])
% Adjust layout


figure

histogram(total_displacement_center_noisy )

hold on

histogram(total_displacement_center)

%%
% it is required to review the algprithm of UKF at first by the related
% paper to understand this code easily. it is tried to follow the nomation
% from the paper  "Application of the unscented Kalman filter for real-time
% nonlinear structural system identification" by A. Smyth. The example is
% used from the  thesis named "Identification of Civil Structural Parameters Using
% the Extended Kalman Filter".

clc
close all
% Number of steps
nSteps = floor((tspan(end) - tspan(1)) / dt);

% Initialize the solution
tt = linspace(tspan(1), tspan(end), nSteps+1)';
%     X = zeros(nSteps+1, length(X0));


% Using this code, i want to estimate the states, moving mass and
% speed.
% the state vector is [X, Xdot, M and V]

X_int = [zeros(1,n_modes)  zeros(1,n_modes) 0.7*M1 0.7*v1 ];

p=diag([ones(1,2*n_modes),10^4,10^1]);  %initial covariance

R0(1)=10^-1 ;            %covariance of measurement noise
Q=(10^-6.2 )*eye(2*n_modes+2); %covariance of process noise

St_Len=length(X_int);

% X_est = zeros(St_Len,nSteps+1);


X_est(:,1) = X_int';
gama=0;
r=1;
alfa=0.9;


% Time integration using 4th-order Runge-Kutta
for i = 1:nSteps


    [W,SP] = sigmapoint(X_est(:,i),p,gama,St_Len); %sigma points

    %     X_new = TF(@(tt, X)BMM_SS,X_temp,tt, dt);

    MF=@(X)phi_center* X(1:n_modes);

    % Calculation of trnasformed sigma points and its mean
    XS=0;
    YS=0;
    for g=1:2*St_Len+1
        x=SP(:,g);
        F(:,g)=TF(x,tt(i), dt);
        x1=F(:,g);
        H(:,g)=MF(x1);
        XS1(:,g)=F(:,g)*W(g);
        YS1(:,g)=H(:,g)*W(g);
        XS=XS+XS1(:,g);                   % average of transformed Sigma points
        YS=YS+YS1(:,g);                    % average of output of system considering sigma points
    end


    %adding adaptive Based on the following paper

    %Song, M., Astroza, R., Ebrahimian, H., Moaveni, B. and Papadimitriou, C., 2020.
    %Adaptive Kalman filters for nonlinear finite element model updating. Mechanical
    % Systems and Signal Processing, 143, p.106837.
    % we used forgeting factor method based on Innovation-based cov. matching
    % method

    %       d_temp=udd(t)-YS;

    %       R0(t+1)=alfa*R0(t)+(1-alfa)*(d_temp*d_temp'+YS);


    % Calculation of Covariance matrices
    % Initial covariance matrices
    PX=0;
    PY=0;
    PXY=0;

    for b=1:2*St_Len+1

        PX1=W(b)*(F(:,b)-XS)*(F(:,b)-XS)';
        PY1=W(b)*(H(:,b)-YS)*(H(:,b)-YS)';
        PXY1=W(b)*(F(:,b)-XS)*(H(:,b)-YS)';
        PX=PX+PX1;
        PY=PY+PY1;
        PXY=PXY+PXY1;
    end
    PY1=PY;
    PX=PX+Q;  %adding covariance of the process noise
    PY=PY+R0(i);  %%adding covariance of the measurement noise

    k=PXY/PY;    % Kalman Gain
    Res(i)=total_displacement_center_noisy(i)-YS;
    X_est(:,i+1)=XS+k*Res(i);  % Updated States
    p=(PX)-k*(PY)*k';         % updated states covariance matrix for new iteration

    P_saved(:,:,i)=p;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %adding adaptive Based on the following paper

    %Song, M., Astroza, R., Ebrahimian, H., Moaveni, B. and Papadimitriou, C., 2020.
    %Adaptive Kalman filters for nonlinear finite element model updating. Mechanical
    % Systems and Signal Processing, 143, p.106837.
    % we used forgeting factor method based on Residual-based cov. matching
    % method

    x_temp=X_est(:,i+1);
    %         F_temp=f(x_temp);
    H_temp=MF(x_temp);
    eps(i)= total_displacement_center(i)-H_temp;
    %
    R0(i+1)=alfa*R0(i)+(1-alfa)*(eps(i)*eps(i)'+PY1); %updating the measurement error cov.

    %        R0=mean(eps(t)*eps(t)')+YS;
    r=r+1;




end
%%
close all
clc

total_vel_center_estimated = (phi_center* X_est(n_modes+1:2*n_modes,:));

total_dis_center_estimated = (phi_center* X_est(1:n_modes,:));
figure('Position', [100, 100, 500, 170])  % [left, bottom, width, height]



% First subplot for displacement
subplot(1,2,1)
plot(t, total_displacement_center, 'Color', colors1(1, :), 'LineWidth', 1)
hold on
plot(t, total_displacement_center_noisy, 'Color', colors(3, :), 'LineWidth', 0.9,LineStyle=':')
legend('Ground Truth', 'Noisy Measurement', 'Location', 'Best')
% title('Total Displacement at the Center of the Beam', 'FontName', 'Times New Roman', 'FontSize', 10)
ylabel('Displacement (m)', 'FontName', 'Times New Roman', 'FontSize', 10)
xlabel('Time (s)', 'FontName', 'Times New Roman', 'FontSize', 10)
set(gca, 'FontName', 'Times New Roman', 'FontSize', 10)
xlim([0,10])
% Second subplot for PSD
subplot(1,2,2)
plot(freq(1:1000), log10(abs(FFT_signal(1:1000))), 'LineWidth', 1.5)
% title('Total Displacement at the Center of the Beam', 'FontName', 'Times New Roman', 'FontSize', 10)
ylabel('PSD', 'FontName', 'Times New Roman', 'FontSize', 10)
xlabel('Frequency (rad/sec)', 'FontName', 'Times New Roman', 'FontSize', 10)
set(gca, 'FontName', 'Times New Roman', 'FontSize', 10)
xlim([0,5])
exportgraphics(gcf, 'Displacement_dis.pdf', 'ContentType', 'vector', 'BackgroundColor', 'none');
% Adjust layout

var_mass = sqrt(P_saved(2*n_modes+1, 2*n_modes+1, :));  % Assuming P is 3D and this extracts the variances
var_mass = reshape(var_mass, [1, length(var_mass)]);
mass_ub = X_est(2*n_modes+1, 2:end) + 3 * var_mass;
mass_lb = X_est(2*n_modes+1, 2:end) - 3 * var_mass;

var_vel = sqrt(P_saved(2*n_modes+2, 2*n_modes+2, :));  % Assuming P is 3D and this extracts the variances
var_vel = reshape(var_vel, [1, length(var_mass)]);
vel_ub = X_est(2*n_modes+2, 2:end) + 3 * var_vel;
vel_lb = X_est(2*n_modes+2, 2:end) - 3 * var_vel;

%%

figure('Position', [100, 100, 500, 170])  % [left, bottom, width, height]

subplot(1,2,1)
plot(t(2:end), X_est(2*n_modes+1,2:end), 'LineWidth', 1.5)
hold on
Mass = M1 * ones(1, length(t(2:end)));
plot(t(2:end), Mass, '--', 'LineWidth', 1.5)
% Plot the shaded area for 2-sigma range
fill([t(2:end)' fliplr(t(2:end)')], [mass_ub fliplr(mass_lb)], 'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none')

legend('UKF estimated', 'Ground Truth', 'Location', 'Best')
ylabel('Moving mass (kg)', 'FontName', 'Times New Roman', 'FontSize', 10)
xlabel('Time (sec)', 'FontName', 'Times New Roman', 'FontSize', 10)
ylim([0.5*M1, 1.2*M1])  % Set y-axis limits as needed
xlim([0,15])
set(gca, 'FontName', 'Times New Roman', 'FontSize', 10)

subplot(1,2,2)
plot(t(2:end), X_est(2*n_modes+2,2:end), 'LineWidth', 1.5)
hold on
Vel = v1 * ones(1, length(tt(2:end)));
plot(t(2:end), Vel, '--', 'LineWidth', 1.5)
fill([t(2:end)' fliplr(t(2:end)')], [vel_ub fliplr(vel_lb)], 'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none')

% legend('UKF estimated', 'Ground Truth', 'Location', 'Best')
ylabel('Velocity (m/s)', 'FontName', 'Times New Roman', 'FontSize', 10)
xlabel('Time (sec)', 'FontName', 'Times New Roman', 'FontSize', 10)
set(gca, 'FontName', 'Times New Roman', 'FontSize', 10)
xlim([0,15])


% Save the figure
exportgraphics(gcf, 'Parameters_dis.pdf', 'ContentType', 'vector', 'BackgroundColor', 'none');

figure('Position', [100, 100, 500, 170])  % [left, bottom, width, height]

subplot(1,2,1)
plot(t, total_dis_center_estimated, 'LineWidth', 1.5)
hold on
plot(t, total_displacement_center, '--', 'LineWidth', 1.5)
legend('UKF estimated', 'Ground Truth', 'Location', 'Best')
ylabel(' Displacement (m)', 'FontName', 'Times New Roman', 'FontSize', 10)
xlabel('Time (sec)', 'FontName', 'Times New Roman', 'FontSize', 10)
% ylim([0.5*M1, 1.2*M1])  % Set y-axis limits as needed
xlim([10, 15])
set(gca, 'FontName', 'Times New Roman', 'FontSize', 10)



subplot(1,2,2)
plot(t, total_vel_center_estimated, 'LineWidth', 1.5)
hold on
plot(t, total_vel_center, '--', 'LineWidth', 1.5)
% legend('UKF estimated', 'Ground Truth', 'Location', 'Best')
ylabel(' Velocity (m/s)', 'FontName', 'Times New Roman', 'FontSize', 10)
xlabel('Time (sec)', 'FontName', 'Times New Roman', 'FontSize', 10)
set(gca, 'FontName', 'Times New Roman', 'FontSize', 10)
xlim([10, 15])

exportgraphics(gcf, 'states_dis.pdf', 'ContentType', 'vector', 'BackgroundColor', 'none');

FFT_R0=fft(R0);

figure('Position', [100, 100, 500, 170])  % [left, bottom, width, height]
subplot(1,2,1)
plot(t, 2*sqrt(R0), 'LineWidth', 1.5)
% min_line = (offset-amplitude_factor)^2 * ones(1, length(t));
% max_line= (offset+amplitude_factor)^2 * ones(1, length(t));
% ave_line=(offset)^2 * ones(1, length(t));
min_line = (min(noise_level)) * ones(1, length(t));
max_line= (max(noise_level)) * ones(1, length(t));
ave_line=(mean(noise_level)) * ones(1, length(t));
hold on
plot(t, min_line,'r','LineWidth', 1.5)
hold on
plot(t, max_line,'g','LineWidth', 1.5)
hold on
plot(t, ave_line,'k','LineWidth', 1.5)
ylabel('R', 'FontName', 'Times New Roman', 'FontSize', 10)
xlabel('Time (sec)', 'FontName', 'Times New Roman', 'FontSize', 10)
% Create the legend
lgd = legend({'R', 'Lb', 'Ub', 'Av'}, 'NumColumns', 2);

% Adjust font size
lgd.FontSize = 8;  % Smaller font size to reduce the box size

% Manually adjust the position and size of the legend box
lgd.Position = [0.257, 0.7, 0.15, 0.15];  % [left, bottom, width, height]

% Optionally, turn off the legend box outline for a cleaner look
lgd.Box = 'off';
xlim([5,20])
% ylim([-0.0001,0.0002121])

subplot(1,2,2)
plot(freq(1:3700), (abs(FFT_R0(1:3700))), 'LineWidth', 1.5)
% title('Power spectural density of Error', 'FontName', 'Times New Roman', 'FontSize', 10)
ylabel('PSD', 'FontName', 'Times New Roman', 'FontSize', 10)
xlabel('Frequency (rad/sec)', 'FontName', 'Times New Roman', 'FontSize', 10)
set(gca, 'FontName', 'Times New Roman', 'FontSize', 10)
exportgraphics(gcf, 'R.pdf', 'ContentType', 'vector', 'BackgroundColor', 'none');
%%
% Create the figure
figure('Position', [100, 100, 500, 170]);  % [left, bottom, width, height]

% Plot on the left
subplot(1,2,1)
plot(t(2:end), X_est(2*n_modes+1,2:end), 'LineWidth', 1.5)
hold on
Mass = M1 * ones(1, length(t(2:end)));
plot(t(2:end), Mass, '--', 'LineWidth', 1.5)
fill([t(2:end)' fliplr(t(2:end)')], [mass_ub fliplr(mass_lb)], 'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none')

legend('UKF estimated', 'Ground Truth', 'Location', 'Best')
ylabel('Moving mass (kg)', 'FontName', 'Times New Roman', 'FontSize', 10)
xlabel('Time (sec)', 'FontName', 'Times New Roman', 'FontSize', 10)
ylim([0.5*M1, 1.2*M1])
xlim([0,15])
set(gca, 'FontName', 'Times New Roman', 'FontSize', 10)

% Plot on the right
subplot(1,2,2)
plot(t(2:end), X_est(2*n_modes+2,2:end), 'LineWidth', 1.5)
hold on
Vel = v1 * ones(1, length(t(2:end)));
plot(t(2:end), Vel, '--', 'LineWidth', 1.5)
fill([t(2:end)' fliplr(t(2:end)')], [vel_ub fliplr(vel_lb)], 'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none')

ylabel('Velocity (m/s)', 'FontName', 'Times New Roman', 'FontSize', 10)
xlabel('Time (sec)', 'FontName', 'Times New Roman', 'FontSize', 10)
set(gca, 'FontName', 'Times New Roman', 'FontSize', 10)
xlim([0,15])

% Save the figure without borders using exportgraphics
exportgraphics(gcf, 'Parameters_dis.pdf', 'ContentType', 'vector', 'BackgroundColor', 'none');
%%


