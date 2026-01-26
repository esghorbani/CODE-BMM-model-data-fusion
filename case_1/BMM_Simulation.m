
clc
close all
clear all


% Example usage
n_modes = 8;  % Number of modes to consider
M1 = 250;                % First moving mass (assumed zero as in the original code)
M2 = 0;              % Second moving mass
m = 20;                % Beam unit mass
EI = 60000000;         % Young's modulus times moment of inertia
L = 30;                % Beam length
v1 = 5;               % Velocity of the first moving mass
v2 = 0;               % Velocity of the second moving mass    % Beam length
tspan = [0 60];  % Time span for the simulation
X0 = [ones(1,n_modes)  ones(1,n_modes)];  % Initial state [x1; x2; ...; xdn]
% [t, X] = ode45(@(tt, X) beam_moving_mass_state_space(tt, X, n_modes,M1,M2,m,EI,L,v1,v2), tspan, X0);
dt=0.001;
% Solve the ODE using the custom ode4 solver with a fixed time step
[t, X] = ode4(@(tt, X) beam_moving_mass_state_space(tt, X, n_modes, M1, M2, m, EI, L, v1, v2), tspan, X0, dt);


x_center = L / 2;
phi_center = mode_shapes(x_center, n_modes, L)';
% Calculate the total displacement at the center of the beam over time

% Plot the total displacement at the center of the beam over time
x = X(:,1:n_modes)';
total_displacement_center = (phi_center* x);

FFT_signal=fft(total_displacement_center);
Fs=1/dt;
freq = (0:Fs/length(total_displacement_center):Fs/2);

figure
plot(freq(1:6700), log10(abs(FFT_signal(1:6700))));
title('Total Displacement at the Center of the Beam');
xlabel('Frequency (rad/sec)');
ylabel('PSD');


figure
plot(t, total_displacement_center);
title('Total Displacement at the Center of the Beam');
xlabel('Time (s)');
ylabel('Displacement (m)');


function [t, X] = ode4(odefun, tspan, X0, dt)
    % Number of steps
    nSteps = floor((tspan(end) - tspan(1)) / dt);
    
    % Initialize the solution
    t = linspace(tspan(1), tspan(end), nSteps+1)';
    X = zeros(nSteps+1, length(X0));
    X(1,:) = X0';
    
    % Time integration using 4th-order Runge-Kutta
    for i = 1:nSteps
        tt = t(i);
        Xt = X(i,:)';
        
        k1 = odefun(tt, Xt);
        k2 = odefun(tt + dt/2, Xt + dt/2 * k1);
        k3 = odefun(tt + dt/2, Xt + dt/2 * k2);
        k4 = odefun(tt + dt, Xt + dt * k3);
        
        X(i+1,:) = Xt + (dt/6) * (k1 + 2*k2 + 2*k3 + k4);
    end
end



function [dX] = beam_moving_mass_state_space(tt, X, n,M1,M2,m,EI,L,v1,v2)
% Parameters


% Derived parameters
alpha1 = M1 / (m * L);
beta1 = pi * v1 / L;
alpha2 = M2 / (m * L);
beta2 = pi * v2 / L;

t1 = mod(tt, L / v1);  % Time modulation for the first moving mass
t2 = mod(tt, L / v2);  % Time modulation for the second moving mass

% Initialization of matrices
Mt = zeros(n, n);      % Mass matrix
Kt = zeros(n, n);      % Stiffness matrix
Ct = zeros(n, n);      % Damping matrix
ft = zeros(n, 1);      % Force vector

% Loop through each element to populate matrices and force vector
for i = 1:n
    for j = 1:n
        dd = (i == j);  % Kronecker delta

        % Mass matrix
        Mt(i, j) = dd + 2 * alpha1 * sin(i * beta1 * t1) * sin(j * beta1 * t1) ...
            + 2 * alpha2 * sin(i * beta2 * t2) * sin(j * beta2 * t2);

        % Damping matrix
        Ct(i, j) = 4 * alpha1 * beta1 * j * sin(i * beta1 * t1) * cos(j * beta1 * t1) ...
            + 4 * alpha2 * beta2 * j * sin(i * beta2 * t2) * cos(j * beta2 * t2);

        % Stiffness matrix
        Kt(i, j) = (EI / m) * (i^2) * (j^2) * (pi^4 / L^4) * dd ...
            - 2 * alpha1 * (j^2) * (beta1^2) * sin(i * beta1 * t1) * sin(j * beta1 * t1) ...
            - 2 * alpha2 * (j^2) * (beta2^2) * sin(i * beta2 * t2) * sin(j * beta2 * t2);
    end

    % Force vector
    ft(i) = alpha1 * sqrt(2 * L) * 9.81 * sin(i * beta1 * t1) ...
        + alpha2 * sqrt(2 * L) * 9.81 * sin(i * beta2 * t2);
end

% State vector
x = X(1:n);   % Displacement
xd = X(n+1:2*n);  % Velocity

% Calculating acceleration
xdd = Mt \ (ft - Ct * xd - Kt * x);

% State derivative
dX = [xd; xdd];
end


function phi = mode_shapes(x, n, L)
phi = zeros(n, 1);
for i = 1:n
    phi(i) = sin((i * pi * x) / L);
end
end

