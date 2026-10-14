function [X_est, P_saved, R0] = run_UKF_displacement(t, y_disp, M1, v1, param)

global n_modes L

dt = param.dt;
Nt = length(t);
if ~isfield(param, 'R_covariance'), param.R_covariance = 'innovation'; end
if ~any(strcmpi(param.R_covariance, {'prior', 'posterior', 'innovation'}))
    error('R_covariance must be ''prior'', ''posterior'' or ''innovation''.');
end

X_int  = [zeros(1,n_modes), zeros(1,n_modes), 0.7*M1, 0.7*v1];
St_Len = length(X_int);

X_est   = zeros(St_Len, Nt);
P_saved = zeros(St_Len, St_Len, Nt);

X_est(:,1) = X_int';
P = diag([ones(1,2*n_modes), 1e4, 1e1]);
P_saved(:,:,1) = P;

if ~isfield(param, 'Q'),    param.Q    = 1e-6; end   % optional overrides (defaults unchanged)
if ~isfield(param, 'alfa'), param.alfa = 0.7;  end
Q = param.Q * eye(St_Len);
R0 = zeros(1,Nt);
R0(1) = 1e-1;

gama = 0;
alfa = param.alfa;

% Use global L to stay consistent with your dynamics
phi_center = mode_shapesPP(L/2, n_modes, L)';
MF = @(X) phi_center * X(1:n_modes);
measurement_row = [phi_center, zeros(1,n_modes+2)];

for k = 1:Nt-1

    [W,SP] = sigmapoint(X_est(:,k), P, gama, St_Len);

    XS = 0; YS = 0;
    F = zeros(St_Len, size(SP,2));
    H = zeros(1, size(SP,2));

    for j = 1:size(SP,2)
        F(:,j) = TF(SP(:,j), t(k), dt);
        H(j)   = MF(F(:,j));
        XS = XS + W(j)*F(:,j);
        YS = YS + W(j)*H(j);
    end

    PX = 0;
    for j = 1:size(SP,2)
        PX  = PX  + W(j)*(F(:,j)-XS)*(F(:,j)-XS)';
    end

    PX  = PX + Q;
    PX  = (PX + PX')/2;
    % The linear measurement sees the complete predicted covariance,
    % including the additive process noise introduced before correction.
    PXY = PX * measurement_row';
    PY1 = measurement_row * PX * measurement_row';
    PY  = PY1 + R0(k);

    K   = PXY / PY;
    Res = y_disp(k+1) - YS;

    X_est(:,k+1) = XS + K*Res;
    correction = eye(St_Len) - K*measurement_row;
    P = correction*PX*correction' + K*R0(k)*K';
    P = (P + P')/2;

    P_saved(:,:,k+1) = P;

    % Adaptive R (measurement-based)
    H_temp = MF(X_est(:,k+1));
    eps    = y_disp(k+1) - H_temp;

    % Default retains Song et al.'s prior-covariance approximation.
    % The exact posterior term is available for a separate sensitivity run.
    if strcmpi(param.R_covariance, 'innovation')
        % innovation form (Mohamed & Schwarz 1999): R = nu*nu' - H P- H', floored
        R0(k+1) = max(alfa*R0(k) + (1-alfa)*(Res^2 - PY1), 1e-8);
    else
        R_cov = PY1;
        if strcmpi(param.R_covariance, 'posterior')
            R_cov = measurement_row * P * measurement_row';
        end
        R0(k+1) = alfa*R0(k) + (1-alfa)*(eps^2 + R_cov);
    end
end

end
