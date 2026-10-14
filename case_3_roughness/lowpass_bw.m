function y = lowpass_bw(x, fc, dt)
%LOWPASS_BW  Zero-phase 4th-order Butterworth low-pass, cut-off fc [Hz].
[b, a] = butter(4, fc/(0.5/dt));
y = filtfilt(b, a, x);
end
