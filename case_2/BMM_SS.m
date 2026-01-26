function [dX] = BMM_SS(tt, X, n,M1,m,EI,L,v1)
% Parameters


% Derived parameters
alpha1 = M1 / (m * L);
beta1 = pi * v1 / L;

t1 = mod(tt, L / v1);  % Time modulation for the first moving mass

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
        Mt(i, j) = dd + 2 * alpha1 * sin(i * beta1 * t1) * sin(j * beta1 * t1);

        % Damping matrix
        Ct(i, j) = 4 * alpha1 * beta1 * j * sin(i * beta1 * t1) * cos(j * beta1 * t1);

        % Stiffness matrix
        Kt(i, j) = (EI / m) * (i^2) * (j^2) * (pi^4 / L^4) * dd ...
            - 2 * alpha1 * (j^2) * (beta1^2) * sin(i * beta1 * t1) * sin(j * beta1 * t1);
    end

    % Force vector
    ft(i) = alpha1 * sqrt(2 * L) * 9.81 * sin(i * beta1 * t1);
end

% State vector
x = X(1:n);   % Displacement
xd = X(n+1:2*n);  % Velocity

% Calculating acceleration
xdd = Mt \ (ft - Ct * xd - Kt * x);

% State derivative
dX = [xd; xdd];
end