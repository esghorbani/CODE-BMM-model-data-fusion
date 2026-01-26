clc
close all
clear all

% Set the random seed for reproducibility
rng(1234);  % You can choose any seed number you like


global n_modes m EI L ;   % Global variables

n_modes=4;
M1 = 5000;                % First moving mass (assumed zero as in the original code)
m = 200;                % Beam unit mass
EI = 5.3*10^9;         % Young's modulus times moment of inertia
L = 60;                % Beam length
v1 = 15;               % Velocity of the first moving mass
tspan = [0 4];  % Time span for the simulation
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







% Define the noise level (standard deviation of the noise)
num_steps = length(total_displacement_center);
time = linspace(0, 2 * pi, num_steps);  % Example time vector from 0 to 2*pi

% Define the time-varying noise level as a sinusoidal function of time
amplitude = 0.005;  % Amplitude of the sine wave
offset = 0.01;     % Offset for the noise level
frequency = 5;     % Frequency of the sine wave

noise_level = offset + amplitude * sin(frequency * time);

% Initialize the noisy displacement vector
total_displacement_center_noisy = zeros(size(total_displacement_center));

% Add time-varying white noise to the total_displacement_center
for i = 1:num_steps
    total_displacement_center_noisy(i) = total_displacement_center(i) + noise_level(i) * randn();
end

% Add white noise to the total_displacement_center
% total_displacement_center_noisy = total_displacement_center + noise_level * randn(size(total_displacement_center));

% Create figure and set size

% figure('Position', [100, 100, 900, 300])  % [left, bottom, width, height]
% 
% % First subplot for displacement
% subplot(2,2,1)
% plot(t, total_displacement_center, 'r', 'LineWidth', 1.5)
% hold on
% plot(t, total_displacement_center_noisy, 'b', 'LineWidth', 1)
% legend('Ground Truth', 'Noisy Measurement', 'Location', 'Best')
% title('Total Displacement at the Center of the Beam', 'FontName', 'Times New Roman', 'FontSize', 10)
% ylabel('Displacement (m)', 'FontName', 'Times New Roman', 'FontSize', 10)
% xlabel('Time (s)', 'FontName', 'Times New Roman', 'FontSize', 10)
% set(gca, 'FontName', 'Times New Roman', 'FontSize', 10)
% xlim([5,15])
% % Second subplot for PSD
% subplot(2,2,2)
% plot(freq(1:2000), log10(abs(FFT_signal(1:2000))), 'LineWidth', 1.5)
% title('Total Displacement at the Center of the Beam', 'FontName', 'Times New Roman', 'FontSize', 10)
% ylabel('PSD', 'FontName', 'Times New Roman', 'FontSize', 10)
% xlabel('Frequency (rad/sec)', 'FontName', 'Times New Roman', 'FontSize', 10)
% set(gca, 'FontName', 'Times New Roman', 'FontSize', 10)
% 
% subplot(2,2,3)
% plot(freq(1:1000), log10(abs(FFT_signal(1:1000))), 'LineWidth', 1.5)
% title('Total Displacement at the Center of the Beam', 'FontName', 'Times New Roman', 'FontSize', 10)
% ylabel('PSD', 'FontName', 'Times New Roman', 'FontSize', 10)
% xlabel('Frequency (rad/sec)', 'FontName', 'Times New Roman', 'FontSize', 10)
% set(gca, 'FontName', 'Times New Roman', 'FontSize', 10)
% xlim([1,5])
% 
% subplot(2,2,4)
% plot(freq(1:1000), log10(abs(FFT_signal(1:1000))), 'LineWidth', 1.5)
% title('Total Displacement at the Center of the Beam', 'FontName', 'Times New Roman', 'FontSize', 10)
% ylabel('PSD', 'FontName', 'Times New Roman', 'FontSize', 10)
% xlabel('Frequency (rad/sec)', 'FontName', 'Times New Roman', 'FontSize', 10)
% set(gca, 'FontName', 'Times New Roman', 'FontSize', 10)
% xlim([17,22])
% % Adjust layout



%%%%%%%%%%%%%%%%%
figure
[N,edges] = histcounts(total_displacement_center_noisy, 'Normalization','pdf');
edges = edges(2:end) - (edges(2)-edges(1))/2;
plot(edges, N);

hold on

[N,edges] = histcounts(total_displacement_center, 'Normalization','pdf');
edges = edges(2:end) - (edges(2)-edges(1))/2;
plot(edges, N);
