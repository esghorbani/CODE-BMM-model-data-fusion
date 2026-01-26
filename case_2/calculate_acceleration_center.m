function total_accel_center = calculate_acceleration_center(t, X, n_modes, M1, m, EI, L, v1, phi_center)
    % Calculate the total acceleration at the center of the beam over time
    % Inputs:
    % - t: Time vector
    % - X: State vector from ODE solver [displacement; velocity]
    % - n_modes: Number of modes
    % - M1: First moving mass
    % - m: Beam unit mass
    % - EI: Young's modulus times moment of inertia
    % - L: Beam length
    % - v1: Velocity of the first moving mass
    % - phi_center: Mode shapes at the center of the beam
    
    % Initialize matrices
    x = X(:,1:n_modes)';   % Displacement
    xd = X(:,n_modes+1:end)';  % Velocity
    xdd = zeros(size(xd));  % Initialize acceleration matrix
    
    % Derived parameters
    alpha1 = M1 / (m * L);
    beta1 = pi * v1 / L;

    % Loop through each time step
    for k = 1:length(t)
        % Compute the matrices at the current time step
        tt = t(k);
        x_k = x(:,k);
        xd_k = xd(:,k);

        t1 = mod(tt, L / v1);  % Time modulation for the first moving mass

        % Recalculate the matrices at time t1
        Mt = zeros(n_modes, n_modes);      % Mass matrix
        Kt = zeros(n_modes, n_modes);      % Stiffness matrix
        Ct = zeros(n_modes, n_modes);      % Damping matrix
        ft = zeros(n_modes, 1);            % Force vector

        for i = 1:n_modes
            for j = 1:n_modes
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

        % Compute acceleration at this time step
        xdd(:,k) = Mt \ (ft - Ct * xd_k - Kt * x_k);
    end

    % Total acceleration at the center of the beam
    total_accel_center = (phi_center * xdd);
end
