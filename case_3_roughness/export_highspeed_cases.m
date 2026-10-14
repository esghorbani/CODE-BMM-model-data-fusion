% EXPORT_HIGHSPEED_CASES  Mass-estimate histories (Mhat/M) of the tabulated high-speed cases.
load('highspeed_rough_results.mat');
F = '../../PAPER-2026-BMM-FUSION/JSV_Submission/rebuttal/2026-09-11-round-1/figures/';
cases = [30 5000; 60 5000; 80 5000; 60 10000; 80 10000];
t = res(1).t; k = find(t <= 15); k = k(1:5:end); Mx = t(k)'; names = {'t'};
for i = 1:size(cases,1)
    r = res([res.v] == cases(i,1) & [res.M] == cases(i,2));
    Mx = [Mx, (r.Mhat(k)/r.M)']; names{end+1} = sprintf('r%d', i);
    fprintf('case %d: v=%d M=%d DAF %.2f v0err %.1f dM %.2f dv %.3f\n', i, r.v, r.M, r.DAF, r.v0err, r.dM, r.dv);
end
writetable(array2table(Mx, 'VariableNames', names), [F 'highspeed_cases.csv']); disp('csv ok')
