% PLOT_GENERALITY_ERROR_MAP  Converged mass and speed errors over the (mass, speed) grid of the
% 147-case study, without and with model error, as filled contours (generality_error_map.pdf).
close all
C = readtable('generality_cases_clean.csv'); E = readtable('generality_cases_error.csv');   % per-case summary tables: clean and model-error runs
Ms = unique(C.M_true); vs = unique(C.v_true); [MM, VV] = meshgrid(Ms, vs);   % mass-speed grid
grid_of = @(T, col) reshape(arrayfun(@(m, v) T.(col)(T.M_true==m & T.v_true==v), MM, VV), size(MM));   % column of T arranged on the grid
Z = {grid_of(C,'dM'), grid_of(E,'dM'); grid_of(C,'dv'), grid_of(E,'dv')};   % (a) dM clean, (b) dM error; (c) dv clean, (d) dv error
ttl = {'Mass error (%), no model error', 'Mass error (%), model error'; 'Speed error (%), no model error', 'Speed error (%), model error'};
lv = {-12:1:0, -0.5:0.05:0.5}; fn = 'Times New Roman'; fs = 14.4;   % contour levels for mass (%) and speed (%)
figure('Position',[100 100 900 560])
for r = 1:2
  for c = 1:2
    subplot(2,2,(r-1)*2+c)
    % filled contours with the case grid as dots
    contourf(MM/1000, VV, Z{r,c}, lv{r}, 'LineColor', [0.3 0.3 0.3], 'LineWidth', 0.3); hold on
    plot(MM(:)/1000, VV(:), 'k.', 'MarkerSize', 4)
    if r == 1 && c == 1   % cases A (4500 kg, 28 m/s) and B (5800 kg, 16 m/s) of the convergence figure, panel (a) only
      plot([4.5 5.8], [28 16], 'o', 'MarkerSize', 11, 'LineWidth', 2.2, 'MarkerEdgeColor', [0.85 0 0])
      plot([4.5 5.8], [28 16], '+', 'MarkerSize', 7, 'LineWidth', 1.4, 'Color', [0.85 0 0])
      text(4.5+0.09, 28, 'A', 'FontName', fn, 'FontSize', fs+1, 'FontWeight', 'bold', 'Color', [0.85 0 0], 'VerticalAlignment', 'middle', 'BackgroundColor', 'w', 'Margin', 1)
      text(5.8-0.09, 16, 'B', 'FontName', fn, 'FontSize', fs+1, 'FontWeight', 'bold', 'Color', [0.85 0 0], 'VerticalAlignment', 'middle', 'HorizontalAlignment', 'right', 'BackgroundColor', 'w', 'Margin', 1)
    end
    caxis([lv{r}(1) lv{r}(end)]); colormap(gca, flipud(parula)); cb = colorbar; cb.Label.String = '%';
    xlabel('Moving mass (t)'); ylabel('Speed (m/s)'); title(ttl{r,c}, 'FontWeight', 'normal')
    set(gca, 'FontName', fn, 'FontSize', fs); xlim([4 6]); ylim([13 31])
    lab = {'(a)','(b)';'(c)','(d)'}; text(-0.13, -0.2, lab{r,c}, 'Units', 'normalized', 'FontName', fn, 'FontSize', fs)
  end
end
exportgraphics(gcf, 'generality_error_map.pdf', 'ContentType', 'vector', 'BackgroundColor', 'none');
fprintf('mass error, no model error: mean %.2f, worst %.2f at M=%d v=%g\n', mean(C.dM), min(C.dM), C.M_true(C.dM==min(C.dM)), C.v_true(C.dM==min(C.dM)));
fprintf('mass error, model error   : mean %.2f, worst %.2f at M=%d v=%g\n', mean(E.dM), min(E.dM), E.M_true(E.dM==min(E.dM)), E.v_true(E.dM==min(E.dM)));
for v = vs', fprintf('v=%2d: mean dM clean %6.2f  error %6.2f\n', v, mean(C.dM(C.v_true==v)), mean(E.dM(E.v_true==v))); end
disp('map written')
