function fc = contact_force_rough(t, Xtraj, M, v, par)
%CONTACT_FORCE_ROUGH  Reconstruct the beam/vehicle contact force history.
%
%   fc = contact_force_rough(t, Xtraj, M, v, par)
%
%   t     : 1xNt time vector
%   Xtraj : (2n+2) x Nt trajectory of [q; qd; u_v; u_v']
%   M, v  : scalars (ground truth) or 1xNt vectors (AUKF estimates)
%   fc    : 1xNt contact force [N], positive downward on the beam
%
%   f_c = M g - m_s u_v'' - m_w y_w''   (see VBI_ROUGH_RHS)

n  = par.n_modes;
Nt = length(t);
if isscalar(M), M = M*ones(1,Nt); end
if isscalar(v), v = v*ones(1,Nt); end

fc = zeros(1,Nt);
for k = 1:Nt
    X = Xtraj(:,k);
    [~, ~, fc(k)] = vbi_rough_rhs(t(k), X(1:n), X(n+1:2*n), X(2*n+1), X(2*n+2), ...
                                  M(k), v(k), par);
end
end
