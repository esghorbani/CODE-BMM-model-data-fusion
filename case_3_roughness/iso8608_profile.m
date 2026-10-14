function prof = iso8608_profile(roadClass, n_min, n_max, N, seed)
%ISO8608_PROFILE  Stochastic deck-surface roughness profile per ISO 8608:2016.
%
%   prof = iso8608_profile(roadClass, n_min, n_max, N, seed)
%
%   The one-sided displacement PSD of the profile is the ISO 8608 power law
%
%       G_d(n) = G_d(n0) * (n / n0)^(-w),   n0 = 0.1 cycles/m,  w = 2,
%
%   where n is the spatial frequency [cycles/m] and G_d(n0) [m^3] sets the
%   roughness class (geometric means of ISO 8608, Table C.2):
%
%       A: 16e-6   B: 64e-6   C: 256e-6   D: 1024e-6   E: 4096e-6
%
%   A sample profile is synthesised by the sinusoidal superposition method
%   (Agostinacchio et al., 2014; Yang et al., 2004):
%
%       r(x) = sum_{k=1}^{N} A_k cos(2*pi*n_k*x + theta_k),
%       A_k  = sqrt( 2 * G_d(n_k) * dn ),   theta_k ~ U[0, 2*pi),
%
%   with n_k = n_min + (k - 1/2) dn,  dn = (n_max - n_min)/N.
%
%   INPUTS
%     roadClass : 'A'..'E' (char) or a numeric G_d(n0) in m^3
%     n_min     : lowest spatial frequency  [cycles/m]  (default 0.01)
%     n_max     : highest spatial frequency [cycles/m]  (default 3.00)
%     N         : number of harmonics                   (default 300)
%     seed      : rng seed for the random phases        (default 1)
%
%   OUTPUT  prof (struct) with fields
%     class, Gd0, n0, w, n (1xN), A (1xN), theta (1xN), dn, seed,
%     sigma_r : theoretical RMS of r(x) = sqrt(sum A_k^2 / 2)
%
%   Evaluate the profile and its spatial derivatives with EVAL_ROUGHNESS.
%
%   References
%     ISO 8608:2016, Mechanical vibration - Road surface profiles -
%       Reporting of measured data.
%     Agostinacchio, Ciampa, Olita (2014), Eur. Transp. Res. Rev. 6:267-275.
%     Yang, Yau, Wu (2004), Vehicle-Bridge Interaction Dynamics, World Sci.

if nargin < 2 || isempty(n_min), n_min = 0.01; end
if nargin < 3 || isempty(n_max), n_max = 3.00; end
if nargin < 4 || isempty(N),     N     = 300;  end
if nargin < 5 || isempty(seed),  seed  = 1;    end

% --- roughness coefficient G_d(n0) [m^3] -------------------------------
if ischar(roadClass) || isstring(roadClass)
    roadClass = upper(char(roadClass));
    switch roadClass
        case 'A', Gd0 = 16e-6;
        case 'B', Gd0 = 64e-6;
        case 'C', Gd0 = 256e-6;
        case 'D', Gd0 = 1024e-6;
        case 'E', Gd0 = 4096e-6;
        otherwise, error('iso8608_profile: unknown ISO 8608 class "%s".', roadClass);
    end
else
    Gd0 = roadClass;          % numeric G_d(n0) supplied directly
    roadClass = 'custom';
end

n0 = 0.1;                     % reference spatial frequency [cycles/m]
w  = 2;                       % waviness exponent

% --- discretise the spatial-frequency band ------------------------------
dn = (n_max - n_min) / N;
n  = n_min + (( 1:N ) - 0.5) * dn;         % band centres (1xN)
Gd = Gd0 * (n / n0).^(-w);                 % PSD at the band centres

% --- random phases and harmonic amplitudes ------------------------------
s = rng;                                   % preserve caller's rng state
rng(seed);
theta = 2*pi*rand(1, N);
rng(s);

A = sqrt(2 * Gd * dn);                     % amplitude of each cosine

prof = struct('class', roadClass, 'Gd0', Gd0, 'n0', n0, 'w', w, ...
              'n', n, 'A', A, 'theta', theta, 'dn', dn, 'seed', seed, ...
              'sigma_r', sqrt(sum(A.^2)/2));
end
