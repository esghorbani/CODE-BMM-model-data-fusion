function [r, dr, ddr] = eval_roughness(prof, x)
%EVAL_ROUGHNESS  Roughness elevation and its spatial derivatives at x.
%
%   [r, dr, ddr] = eval_roughness(prof, x)
%
%   prof : profile struct from ISO8608_PROFILE, or [] for a smooth surface
%   x    : position(s) along the beam [m], any shape
%
%   r   = sum_k A_k cos(2 pi n_k x + theta_k)             [m]
%   dr  = d r/dx  = -sum_k A_k (2 pi n_k)   sin(...)      [-]
%   ddr = d2r/dx2 = -sum_k A_k (2 pi n_k)^2 cos(...)      [1/m]
%
%   The analytic derivatives are needed because the moving-mass/oscillator
%   equations contain v*r'(vt) (velocity of the contact point) and
%   v^2*r''(vt) (acceleration of the contact point).

if isempty(prof)
    r = zeros(size(x)); dr = r; ddr = r;
    return
end

shp = size(x);
x   = x(:);                                  % column
kx  = 2*pi*prof.n;                           % 1xN wave numbers [rad/m]
arg = x*kx + prof.theta;                     % numel(x) x N
Acol = prof.A(:);

r   =  cos(arg)               * Acol;
dr  = -(sin(arg) .* kx)       * Acol;
ddr = -(cos(arg) .* (kx.^2))  * Acol;

r = reshape(r, shp); dr = reshape(dr, shp); ddr = reshape(ddr, shp);
end
