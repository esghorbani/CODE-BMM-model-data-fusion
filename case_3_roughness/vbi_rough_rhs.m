function [qdd, uvdd, fc] = vbi_rough_rhs(tt, q, qd, uv, uvd, M, v, par)
%VBI_ROUGH_RHS  Beam / moving-oscillator interaction on a rough deck.
%
%   [qdd, uvdd, fc] = vbi_rough_rhs(tt, q, qd, uv, uvd, M, v, par)
%
%   Vehicle model (moving oscillator, Yang-Yau-Wu 2004):
%     * unsprung (wheel) mass  m_w = mu_w * M  in rigid contact with the deck
%     * sprung (body) mass     m_s = (1-mu_w) * M  on a spring k_v / dashpot c_v
%   Total vehicle mass M = m_w + m_s is the quantity estimated by the AUKF.
%
%   Special cases:
%     mu_w = 1                 -> pure moving mass of the paper (BMM_SS.m)
%     mu_w = 0                 -> pure SDOF sprung mass (reviewer's request)
%     par.prof = []            -> smooth deck (r = 0)
%
%   Kinematics of the contact point x_c = v t  (z positive downward):
%     y_w   = z(x_c,t) + r(x_c)                                  (wheel disp.)
%     y_w'  = z_t + v z_x + v r'(x_c)
%     y_w'' = z_tt + 2 v z_xt + v^2 z_xx + v^2 r''(x_c)
%
%   Sprung mass (u_v = dynamic body displacement about static equilibrium):
%     m_s u_v'' = -[ k_v (u_v - y_w) + c_v (u_v' - y_w') ] = -F_s
%
%   Contact force on the beam:
%     f_c = M g + F_s - m_w y_w''  =  M g - m_s u_v'' - m_w y_w''
%
%   Modal equations (phi_k = sqrt(2/L) sin(k pi x/L), modal mass = m):
%     T q'' + (Y + C) q' + Lam q = h
%     T_ij   = d_ij + (m_w/m)    phi_i phi_j
%     Y_ij   =        (2 m_w v/m) phi_i phi_j'
%     C_ij   = 2 zeta_i w_i d_ij          (structural modal damping; optional
%                                           field par.zeta, scalar or n x 1,
%                                           default 0 = undamped as in paper)
%     Lam_ij = w_i^2 d_ij + (m_w v^2/m) phi_i phi_j''
%     h_i    = (phi_i/m) [ M g + F_s - m_w v^2 r''(x_c) ]
%   With mu_w = 1 and r = 0 these are exactly the matrices in BMM_SS.m.
%
%   INPUTS  tt  time [s];  q,qd (n x 1) modal disp./vel.;  uv,uvd scalars;
%           M total mass [kg]; v speed [m/s];
%           par : struct with n_modes, m, EI, L, g, mu_w, k_v, c_v, prof
%   OUTPUTS qdd (n x 1), uvdd scalar, fc contact force [N] (down positive)

n = par.n_modes;  m = par.m;  L = par.L;  g = par.g;

% --- position of the contact point (mass re-enters the span periodically,
%     same convention as BMM_SS.m) ---------------------------------------
t1 = mod(tt, L / v);
xc = v * t1;

% --- mode shapes and their spatial derivatives at x_c --------------------
k     = (1:n).';
phi   = sqrt(2/L) * sin(k*pi*xc/L);          % phi_k(x_c)
dphi  = sqrt(2/L) * (k*pi/L) .* cos(k*pi*xc/L);
ddphi = -(k*pi/L).^2 .* phi;

% --- deck roughness at x_c -----------------------------------------------
[r, dr, ddr] = eval_roughness(par.prof, xc);

% --- vehicle masses -------------------------------------------------------
%   Default: fixed mass RATIO, m_w = mu_w*M (both parts scale with M).
%   If par.m_w_fixed is given, the wheel mass is a KNOWN vehicle property
%   [kg] and only the sprung mass m_s = M - m_w is unknown (M = m_w + m_s).
if isfield(par, 'm_w_fixed') && ~isempty(par.m_w_fixed)
    m_w = par.m_w_fixed;
    m_s = max(M - m_w, 0);
else
    m_w = par.mu_w * M;
    m_s = (1 - par.mu_w) * M;
end

% --- wheel kinematics (terms not involving q'') ---------------------------
yw  = phi.'*q  + r;
ywd = phi.'*qd + v*(dphi.'*q) + v*dr;

% --- suspension force and sprung-mass acceleration ------------------------
if m_s > 0
    Fs   = par.k_v*(uv - yw) + par.c_v*(uvd - ywd);
    uvdd = -Fs / m_s;
else                       % pure moving mass: no suspension DOF
    Fs   = 0;
    uvdd = 0;
end

% --- beam modal matrices --------------------------------------------------
omega2 = (par.EI/m) * (k*pi/L).^4;                 % undamped modal freq^2
T   = eye(n) + (m_w/m) * (phi*phi.');
Y   = (2*m_w*v/m) * (phi*dphi.');
Lam = diag(omega2) + (m_w*v^2/m) * (phi*ddphi.');
h   = (phi/m) * (M*g + Fs - m_w*v^2*ddr);

% structural (modal) damping, zero unless par.zeta is given
if isfield(par, 'zeta') && any(par.zeta ~= 0)
    zeta = par.zeta(:) .* ones(n,1);
    Cst  = diag(2*zeta.*sqrt(omega2));
else
    Cst  = 0;
end

qdd = T \ (h - (Y + Cst)*qd - Lam*q);

% --- contact force (needed for force reconstruction) ---------------------
ywdd = phi.'*qdd + 2*v*(dphi.'*qd) + v^2*(ddphi.'*q) + v^2*ddr;
fc   = M*g + Fs - m_w*ywdd;
end
