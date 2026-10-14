% VERIFY_ISO8608  Check the generated deck profiles against ISO 8608:2016.
%  1. Welch PSD of a long sample vs the target power law and the class limits.
%     NOTE: the sinusoidal-superposition profile is a sum of N discrete lines
%     spaced dn = 0.01 cycles/m (period 1/dn = 100 m > span).  The PSD window
%     must therefore be shorter than 1/dn so that the lines merge into the
%     continuous target spectrum (50-m windows -> 0.02 cycles/m resolution).
%  2. Band-integrated power vs the analytical integral of G_d(n).
%  3. Sample RMS vs analytical  sigma^2 = Gd(n0) n0^2 (1/n_min - 1/n_max).
%  4. Consistency of the class constants with the angular-frequency table.
%  5. IRI of each profile with the Golden Car (Sayers 1995) at 80 km/h.
clc; clear; close all;
L = 60; classes = {'A','B','C'}; Gd0 = [16 64 256]*1e-6; n0 = 0.1;
lims = [2 8 32 128 512 2048]*1e-6;          % Gd(n0) upper limits: A<32, B<128, C<512 ...
fs = 10; fn = 'Times New Roman'; cols = lines(3);
dx = 0.05; xl = (0:dx:6000-dx)'; win = round(50/dx);
fprintf('%-6s %9s %9s %9s | %8s %8s %8s | %6s %6s %6s | %s\n','class','sig_th','sig_60m','sig_6km','P/Gd@.05','P/Gd@0.1','P/Gd@0.5','B1','B2','B3','IRI [m/km]');
figure('Position', [100 100 750 200]); subplot(1,3,[2 3]); hold on
for c = 1:3
    prof = iso8608_profile(classes{c}, 0.01, 3, 300, 42);
    r60 = eval_roughness(prof, linspace(0, L, 3001));
    rl = eval_roughness(prof, xl); rl = rl - mean(rl);
    [P, n] = pwelch(rl, hann(win), [], [], 1/dx, 'onesided');
    band = n >= 0.02 & n <= 2.83;
    loglog(n(band), P(band), 'Color', cols(c,:), 'LineWidth', 1);
    Gd = @(nn) Gd0(c)*(nn/n0).^-2;
    ratio = @(nn) interp1(n, P, nn)/Gd(nn);
    % band-integrated power ratios (sample / analytical) in three bands
    bands = [0.02 0.1; 0.1 0.5; 0.5 2.83]; br = zeros(1,3);
    for b = 1:3
        ib = n >= bands(b,1) & n <= bands(b,2);
        br(b) = trapz(n(ib), P(ib)) / (Gd0(c)*n0^2*(1/bands(b,1) - 1/bands(b,2)));
    end
    sig_th = sqrt(Gd0(c)*n0^2*(1/0.01 - 1/3));
    IRI = golden_car_iri(xl, rl);
    fprintf('%-6s %7.2f mm %7.2f mm %7.2f mm | %8.2f %8.2f %8.2f | %6.2f %6.2f %6.2f | %5.1f\n', classes{c}, ...
        1e3*sig_th, 1e3*std(r60), 1e3*std(rl), ratio(0.05), ratio(0.1), ratio(0.5), br, IRI);
end
nn = logspace(log10(0.011), log10(2.83), 50);
for k = 1:numel(lims), loglog(nn, lims(k)*(nn/n0).^-2, 'k:', 'LineWidth', 0.6); end
for c = 1:3, loglog(nn, Gd0(c)*(nn/n0).^-2, '--', 'Color', cols(c,:), 'LineWidth', 1.2); end
text(0.012*ones(1,5), lims(1:5)*(0.012/n0)^-2*1.8, {'A','B','C','D','E'}, 'FontName', fn, 'FontSize', 8)
xlabel('Spatial frequency n (cycles/m)', 'FontName', fn, 'FontSize', fs); ylabel('G_d(n) (m^3)', 'FontName', fn, 'FontSize', fs)
legend({'sample PSD, class A', 'sample PSD, class B', 'sample PSD, class C'}, 'Location', 'northeast', 'FontSize', 7)
title('(b) PSD of the sampled profiles vs ISO 8608 (dotted: class limits)', 'FontName', fn, 'FontSize', fs, 'FontWeight', 'normal')
set(gca, 'FontName', fn, 'FontSize', fs, 'XScale', 'log', 'YScale', 'log'); xlim([0.01 3]); ylim([1e-9 1e-1]); grid on
subplot(1,3,1); hold on; xg = linspace(0, L, 3001);
for c = 1:3, prof = iso8608_profile(classes{c}, 0.01, 3, 300, 42); plot(xg, 1e3*eval_roughness(prof, xg), 'Color', cols(c,:), 'LineWidth', 1); end
xlabel('x (m)', 'FontName', fn, 'FontSize', fs); ylabel('r(x) (mm)', 'FontName', fn, 'FontSize', fs)
title('(a) Profiles on the span', 'FontName', fn, 'FontSize', fs, 'FontWeight', 'normal'); set(gca, 'FontName', fn, 'FontSize', fs); xlim([0 L])
legend(strcat('class ', classes), 'Location', 'southwest', 'FontSize', 7)
exportgraphics(gcf, 'ISO8608_PSD_check.pdf', 'ContentType', 'vector', 'BackgroundColor', 'none');
fprintf('\nGd(n0=0.1) implied by Gd(Omega0=1 rad/m) = 1e-6 (ISO 8608 class A, angular table): %.2e m^3 (used: 1.6e-5)\n', 2*pi*1e-6*(2*pi*0.1)^-2);

function IRI = golden_car_iri(x, y)
% IRI [m/km] from the Golden Car quarter-car (Sayers 1995): 80 km/h,
% k_t/m_s = 653, k_s/m_s = 63.3, c/m_s = 6.0, m_u/m_s = 0.15,
% 250-mm moving-average pre-filter, IRI = (1/L) int |zs' - zu'| dt.
v = 80/3.6; dx = x(2)-x(1); dt = dx/v;
y = movmean(y, max(1, round(0.25/dx)));
K1 = 653; K2 = 63.3; C = 6.0; mu = 0.15;
A = [0 1 0 0; -K2 -C K2 C; 0 0 0 1; K2/mu C/mu -(K1+K2)/mu -C/mu];
B = [0; 0; 0; K1/mu];
sys = ss(A, B, [0 1 0 -1], 0);
t = (0:numel(y)-1)'*dt;
slope = zeros(4,1); slope([2 4]) = (y(round(11/dx)) - y(1))/(11/v);   % standard initialisation
rs = lsim(sys, y, t, slope);
IRI = trapz(t, abs(rs)) / (x(end)/1000);
end
