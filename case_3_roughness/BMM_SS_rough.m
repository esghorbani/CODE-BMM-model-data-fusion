function dX = BMM_SS_rough(tt, X, M, v, par)
%BMM_SS_ROUGH  Ground-truth dynamics: beam + moving oscillator on rough deck.
%
%   dX = BMM_SS_rough(tt, X, M, v, par)
%
%   State  X = [q (n); qd (n); u_v; u_v']   (2n+2 x 1)
%   M, v are the TRUE (known) mass and speed used to synthesise the data.
%   See VBI_ROUGH_RHS for the equations of motion.

n   = par.n_modes;
q   = X(1:n);
qd  = X(n+1:2*n);
uv  = X(2*n+1);
uvd = X(2*n+2);

[qdd, uvdd] = vbi_rough_rhs(tt, q, qd, uv, uvd, M, v, par);

dX = [qd; qdd; uvd; uvdd];
end
