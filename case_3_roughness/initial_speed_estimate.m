function [v0, T0] = initial_speed_estimate(t, y, L, v_range)
%INITIAL_SPEED_ESTIMATE  Data-driven initial guess of the crossing speed.
%
%   [v0, T0] = initial_speed_estimate(t, y, L, v_range)
%
%   The mid-span response of a beam under a vehicle re-entering the span
%   periodically is periodic with the crossing time T = L/v (for a single
%   crossing, T is the time between entry and exit, equally visible in the
%   signal).  T is estimated in two steps from the noisy, detrended
%   measurement y:
%     1. coarse: the lowest-frequency peak of |FFT(y)| within the admissible
%        band [v_min, v_max]/L whose amplitude is at least 30 % of the
%        largest peak in the band (the fundamental of the crossing pattern;
%        higher harmonics and beam frequencies lie above it);
%     2. refine: the local maximum of the autocorrelation of y closest to
%        the coarse period.
%
%   v_range : [v_min v_max] admissible speeds (default [5 100] m/s)

if nargin < 4, v_range = [5 100]; end
dt = t(2) - t(1);
y  = y(:) - mean(y);
N  = numel(y);

% --- step 1: coarse period from the lowest strong spectral peak ----------
Nf = 2^nextpow2(8*N);                        % zero padding for resolution
Y  = abs(fft(y, Nf)); f = (0:Nf-1)'/(Nf*dt);
band = f >= v_range(1)/L & f <= v_range(2)/L;
fb = f(band); Yb = Y(band);
ipk = find(Yb(2:end-1) > Yb(1:end-2) & Yb(2:end-1) >= Yb(3:end)) + 1;
ipk = ipk(Yb(ipk) >= 0.3*max(Yb));
if isempty(ipk), [~, ipk] = max(Yb); end
T_coarse = 1/fb(ipk(1));

% --- step 2: refine on the autocorrelation --------------------------------
[c, lags] = xcorr(y, 'coeff');
c = c(lags >= 0); lags = lags(lags >= 0);
lag_c = round(T_coarse/dt);
win = max(2, round(0.7*lag_c)) : min(numel(c)-1, round(1.3*lag_c));
cw = c(win);
loc = find(cw(2:end-1) > cw(1:end-2) & cw(2:end-1) >= cw(3:end)) + 1;
if isempty(loc), [~, loc] = max(cw); end
[~, j] = min(abs(win(loc) - lag_c));         % local maximum nearest the coarse period
T0 = lags(win(loc(j))) * dt;
v0 = L / T0;
end
