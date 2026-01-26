
clc
close all
clear all


% Example usage
n_modes = 5;  % Number of modes to consider
M1 = 400;                % First moving mass (assumed zero as in the original code)
m = 50;                % Beam unit mass
EI = 5*10^9;         % Young's modulus times moment of inertia
L = 60;                % Beam length
v1 = 5;               % Velocity of the first moving mass
tspan = [0 60];  % Time span for the simulation
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

FFT_signal=fft(total_displacement_center);
Fs=1/dt;
freq = (0:Fs/length(total_displacement_center):Fs/2);

figure
plot(freq(1:2000), log10(abs(FFT_signal(1:2000))));
title('Total Displacement at the Center of the Beam');
xlabel('Frequency (rad/sec)');
ylabel('PSD');


figure
plot(t, total_displacement_center);
title('Total Displacement at the Center of the Beam');
xlabel('Time (s)');
ylabel('Displacement (m)');









