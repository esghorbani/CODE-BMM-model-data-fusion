% EXPORT_SPEED_FIGURE  CSV data for the pgfplots version of the speed-study figure (letter, R2-C4).
load('reviewer2_results.mat');
F = '../../PAPER-2026-BMM-FUSION/JSV_Submission/rebuttal/2026-09-11-round-1/figures/';
xi = linspace(0, 1, 301)';                        % normalised position vt/L during the first crossing
D = xi; names = {'xi'};
for k = 1:numel(C)
    n = round(par.L/C(k).v/dt) + 1;               % samples in the first crossing
    D = [D, 1e3*interp1(linspace(0,1,n)', C(k).disp(1:n)', xi)];
    names{end+1} = sprintf('d%d', C(k).v);
end
writetable(array2table(D, 'VariableNames', names), [F 'speed_disp.csv']);
T = table([C.v]', [C.dM2]', [C.dv2]', [C.DAF]', [C.fdyn_max]', 'VariableNames', {'v','dM','dv','DAF','fdyn'});
writetable(T, [F 'speed_summary.csv']); disp(T); disp('csv ok')
