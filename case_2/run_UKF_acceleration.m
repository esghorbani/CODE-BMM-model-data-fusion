function [X_est, P_saved] = run_UKF_acceleration(t, y_meas, M1, v1, param)

global n_modes

dt = param.dt;
Nt = length(t);

X_int = [zeros(1,n_modes), zeros(1,n_modes), 0.7*M1, 0.7*v1];
St_Len = length(X_int);

X_est = zeros(St_Len, Nt);
P_saved = zeros(St_Len, St_Len, Nt);

X_est(:,1) = X_int';
P = diag([1e-6*ones(1,2*n_modes), 1e4, 1e1]);   % retuned for the corrected filter
P_saved(:,:,1) = P;

Q = diag([1e-13*ones(1,2*n_modes), 1e-2, 1e-6]);   % modal, mass, speed
R0 = zeros(1,Nt);
R0(1) = 1e-1;

gama = 1;
alfa = 0.9;   % value used for the manuscript results

phi_center = mode_shapesPP(param.L/2, n_modes, param.L)';

for k = 1:Nt-1

    [W, SP] = sigmapoint(X_est(:,k), P, gama, St_Len);

    XS = 0; YS = 0;
    F = zeros(St_Len, size(SP,2));
    H = zeros(1, size(SP,2));

    for j = 1:size(SP,2)
        F(:,j) = TF(SP(:,j), t(k), dt);
        yj = BMM_UKF(t(k+1), F(:,j));
        H(j) = phi_center * yj(n_modes+1:2*n_modes);

        XS = XS + W(j)*F(:,j);
        YS = YS + W(j)*H(j);
    end

    PX = 0; PY = 0; PXY = 0;
    for j = 1:size(SP,2)
        PX  = PX  + W(j)*(F(:,j)-XS)*(F(:,j)-XS)';
        PY  = PY  + W(j)*(H(j)-YS)^2;
        PXY = PXY + W(j)*(F(:,j)-XS)*(H(j)-YS);
    end

    PX = PX + Q;
    PY = PY + R0(k);

    K = PXY / PY;
    res = y_meas(k+1) - YS;   % measurement at the predicted time t(k+1)

    X_est(:,k+1) = XS + K*res;
    P = PX - K*PY*K';

    P_saved(:,:,k+1) = P;

    % Adaptive R: innovation form (Mohamed & Schwarz 1999), floored
    R0(k+1) = max(alfa*R0(k) + (1-alfa)*(res^2 - (PY - R0(k))), 1e-8);
end

end
