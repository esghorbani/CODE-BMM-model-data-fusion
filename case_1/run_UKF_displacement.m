function [X_est, P_saved, R0] = run_UKF_displacement(t, y_disp, M1, v1, param)

global n_modes L

dt = param.dt;
Nt = length(t);

X_int  = [zeros(1,n_modes), zeros(1,n_modes), 0.7*M1, 0.7*v1];
St_Len = length(X_int);

X_est   = zeros(St_Len, Nt);
P_saved = zeros(St_Len, St_Len, Nt);

X_est(:,1) = X_int';
P = diag([ones(1,2*n_modes), 1e4, 1e1]);

Q = 1e-6 * eye(St_Len);
R0 = zeros(1,Nt);
R0(1) = 1e-1;

gama = 0;
alfa = 0.7;

% Use global L to stay consistent with your dynamics
phi_center = mode_shapesPP(L/2, n_modes, L)';
MF = @(X) phi_center * X(1:n_modes);

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

    PX = 0; PY = 0; PXY = 0;
    for j = 1:size(SP,2)
        PX  = PX  + W(j)*(F(:,j)-XS)*(F(:,j)-XS)';
        PY  = PY  + W(j)*(H(j)-YS)^2;
        PXY = PXY + W(j)*(F(:,j)-XS)*(H(j)-YS);
    end

    PY1 = PY;      % <-- IMPORTANT (same as your script)
    PX  = PX + Q;
    PY  = PY + R0(k);

    K   = PXY / PY;
    Res = y_disp(k) - YS;

    X_est(:,k+1) = XS + K*Res;
    P = PX - K*PY*K';

    P_saved(:,:,k) = P;

    % Adaptive R (measurement-based)
    H_temp = MF(X_est(:,k+1));
    eps    = y_disp(k) - H_temp;

    R0(k+1) = alfa*R0(k) + (1-alfa)*(eps^2 + PY1);  % <-- PY1 not PY
end

end
