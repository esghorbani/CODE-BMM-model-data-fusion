function dX = BMM_UKF_rough(tt, X, par)
%BMM_UKF_ROUGH  Estimator (process) model used inside the AUKF.
%
%   dX = BMM_UKF_rough(tt, X, par)
%
%   Augmented state  X = [q (n); qd (n); u_v; u_v'; M; v]   (2n+4 x 1)
%   M and v are treated as constant unknown parameters (dM = dv = 0),
%   exactly as in BMM_UKF.m of Case 1 / Case 2.
%
%   par.prof controls what the FILTER assumes about the deck:
%     par.prof = []    -> filter assumes a smooth deck; the true roughness
%                         acts as an unmodelled physical excitation
%                         (robustness test requested by Reviewer 1, C4)
%     par.prof = prof  -> roughness profile known to the filter
%                         (e.g. measured by a profilometer)

n   = par.n_modes;
q   = X(1:n);
qd  = X(n+1:2*n);
uv  = X(2*n+1);
uvd = X(2*n+2);
M   = X(2*n+3);
v   = X(2*n+4);

[qdd, uvdd] = vbi_rough_rhs(tt, q, qd, uv, uvd, M, v, par);

dX = [qd; qdd; uvd; uvdd; 0; 0];
end
