% EXPORT_TUNING_FIGURES  CSV data for the pgfplots Q- and alpha-sensitivity figures (letter, R2-C3).
% Data-driven initialisation runs (init == 2) of BMM_UKF_Reviewer2_Studies.
load('reviewer2_results.mat');
F = '../../PAPER-2026-BMM-FUSION/JSV_Submission/rebuttal/2026-09-11-round-1/figures/';
t = D.t; k = find(t <= 15); k = k(1:5:end);
Ad = A([A.init] == 2); Bd = B([B.init] == 2);
M = t(k)'; names = {'t'};
for i = 1:numel(Ad), M = [M, Ad(i).M(k)']; names{end+1} = sprintf('M_Q%d', i); end
writetable(array2table(M, 'VariableNames', names), [F 'q_sensitivity.csv']);
fprintf('Q values: %s\n', mat2str([Ad.Q], 3));
M = t(k)'; names = {'t'};
for i = 1:numel(Bd)
    if any(abs(Bd(i).alfa - [0.5 0.7 0.9]) < 1e-9)
        tag = strrep(sprintf('%.1f', Bd(i).alfa), '.', '');
        M = [M, Bd(i).M(k)', Bd(i).v(k)']; names{end+1} = ['M_a' tag]; names{end+1} = ['v_a' tag];
    end
end
writetable(array2table(M, 'VariableNames', names), [F 'alpha_sensitivity.csv']);
disp(strjoin(names, ' ')); disp('csv ok')
% alpha runs with the ORIGINAL initialisation (0.7*M, 0.7*v, prior std 3.2 m/s), init == 1
Bp = B([B.init] == 1); M = t(k)'; names = {'t'};
for i = 1:numel(Bp)
    if any(abs(Bp(i).alfa - [0.5 0.7 0.9]) < 1e-9)
        tag = strrep(sprintf('%.1f', Bp(i).alfa), '.', '');
        M = [M, Bp(i).M(k)', Bp(i).v(k)']; names{end+1} = ['M_a' tag]; names{end+1} = ['v_a' tag];
        fprintf('alpha %.1f (init 0.7): dM %+.2f %%, dv %+.3f %%\n', Bp(i).alfa, Bp(i).dM, Bp(i).dv);
    end
end
writetable(array2table(M, 'VariableNames', names), [F 'alpha_sensitivity_init07.csv']); disp('csv2 ok')
