function [X_est, P_saved, R0, LL] = run_UKF_rough(t, y_disp, M1, v1, par, dt, opts)
%RUN_UKF_ROUGH  Adaptive UKF for the beam / moving-oscillator model.
%
%   [X_est, P_saved, R0] = run_UKF_rough(t, y_disp, M1, v1, par, dt)
%
%   Port of case_1/run_UKF_displacement.m to the extended state
%       X = [q; qd; u_v; u_v'; M; v]      (2 n_modes + 4)
%   Measurement: mid-span beam displacement (same sensor as Case 1).
%   The filter is initialised at 70 % of the true M1 and v1, as in Case 1.
%   Residual-based forgetting-factor adaptation of R follows
%   Song et al. (2020), MSSP 143:106837.
%
%   opts (optional struct) : filter tuning, defaults = Case 1 script
%     .alfa   forgetting factor of the adaptive R      (default 0.9)
%     .Q      process-noise variance (scalar)          (default 10^-6.2)
%     .P0_veh initial variance of the u_v, u_v' states (default 1e-2)
%     .init   initial guess as fraction of truth [M v]  (default [0.7 0.7])
%     .P0_M   initial variance of M                      (default 1e4)
%     .P0_v   initial variance of v                      (default 1e1)
%     .R0_init initial measurement-noise variance         (default 1e-1)
%     .v_init  absolute initial speed [m/s] (overrides init(2)*v1)
%     .R_covariance 'innovation' (default: R = nu^2 - C P- C', Mohamed & Schwarz 1999,
%                    floored at 1e-8), 'posterior' (residual form R = eps^2 + C P+ C'),
%                    or 'prior' (Song et al. approximation)
%   LL : per-step innovation log-likelihood (for filter-bank selection)

n  = par.n_modes;
Nt = length(t);
if nargin < 7, opts = struct(); end
if ~isfield(opts, 'alfa'),   opts.alfa   = 0.9;      end
if ~isfield(opts, 'Q'),      opts.Q      = 10^-6.2;  end
if ~isfield(opts, 'P0_veh'), opts.P0_veh = 1e-2;     end
if ~isfield(opts, 'init'),   opts.init   = [0.7 0.7];end
if ~isfield(opts, 'P0_M'),   opts.P0_M   = 1e4;      end
if ~isfield(opts, 'P0_v'),   opts.P0_v   = 1e1;      end
if ~isfield(opts, 'R0_init'), opts.R0_init = 1e-1;   end
if ~isfield(opts, 'R_covariance'), opts.R_covariance = 'innovation'; end
R_covariance = validatestring(opts.R_covariance, {'prior', 'posterior', 'innovation'});
if isfield(opts, 'v_init') && ~isempty(opts.v_init), opts.init(2) = opts.v_init/v1; end

X_int  = [zeros(1,n), zeros(1,n), 0, 0, opts.init(1)*M1, opts.init(2)*v1];   % initial state: beam and vehicle at rest, M and v at 70 % of truth
St_Len = length(X_int);

X_est   = zeros(St_Len, Nt);
P_saved = zeros(St_Len, St_Len, Nt);

X_est(:,1) = X_int';
%            modal states     u_v    u_v'   M     v
P = diag([ones(1,2*n), opts.P0_veh, opts.P0_veh, opts.P0_M, opts.P0_v]);
P_saved(:,:,1) = P;

I  = eye(St_Len);
Q  = opts.Q * I;
R0 = zeros(1,Nt);
R0(1) = opts.R0_init;
LL = zeros(1,Nt);

gama = 0;   % sigma-point spread parameter (see sigmapoint.m)
alfa = opts.alfa;

phi_center = mode_shapesPP(par.L/2, n, par.L)';
MF = @(X) phi_center * X(1:n);          % measurement function
C = [phi_center, zeros(1,n+4)];         % linear measurement matrix

% ---- filter loop: predict with sigma points, update with the mid-span displacement ----
for k = 1:Nt-1

    [W,SP] = sigmapoint(X_est(:,k), P, gama, St_Len);   % sigma points and weights

    XS = 0; YS = 0;
    F = zeros(St_Len, size(SP,2));
    H = zeros(1, size(SP,2));

    % propagate each sigma point one RK4 step; XS, YS = predicted state and measurement
    for j = 1:size(SP,2)
        F(:,j) = TF_rough(SP(:,j), t(k), dt, par);
        H(j)   = MF(F(:,j));
        XS = XS + W(j)*F(:,j);
        YS = YS + W(j)*H(j);
    end

    % predicted state covariance
    PX = 0;
    for j = 1:size(SP,2)
        PX  = PX  + W(j)*(F(:,j)-XS)*(F(:,j)-XS)';
    end

    % Process noise belongs to the predicted state before observation.
    % For this linear measurement its contributions to PY and PXY are exact.
    PX  = PX + Q;
    PX  = (PX + PX')/2;
    PXY = PX*C';
    PY1 = C*PX*C';
    PY  = PY1 + R0(k);

    K   = PXY / PY;   % Kalman gain
    Res = y_disp(k+1) - YS;   % innovation

    LL(k+1) = -0.5*(log(2*pi*PY) + Res^2/PY);   % innovation log-likelihood (filter-bank score)
    X_est(:,k+1) = XS + K*Res;
    % Joseph form is equivalent to PX-K*PY*K' for this linear update and
    % better preserves covariance symmetry and positive semidefiniteness.
    A = I - K*C;
    P = A*PX*A' + K*R0(k)*K';
    P = (P + P')/2;

    P_saved(:,:,k+1) = P;

    % Adaptive R (residual-based, uses the noisy measurement only)
    H_temp = MF(X_est(:,k+1));
    eps    = y_disp(k+1) - H_temp;

    if strcmp(R_covariance, 'posterior')
        R_term = C*P*C';
        R0(k+1) = alfa*R0(k) + (1-alfa)*(eps^2 + R_term);
    elseif strcmp(R_covariance, 'innovation')
        % innovation form: R = nu*nu' - C P- C' (Mohamed & Schwarz 1999), floored
        R0(k+1) = max(alfa*R0(k) + (1-alfa)*(Res^2 - PY1), 1e-8);
    else
        R_term = PY1;  % Song's prior-covariance approximation
        R0(k+1) = alfa*R0(k) + (1-alfa)*(eps^2 + R_term);
    end
end

end
