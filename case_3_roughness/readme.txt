Case 3 - Robustness analysis (Section 6) and generality study (Section 7)
=========================================================================

Model
-----
Ground truth: Euler-Bernoulli beam + moving oscillator on a rough deck.
  * sprung mass m_s on spring k_v / dashpot c_v over an unsprung wheel mass
    m_w = mu_w*M in rigid contact with the deck; estimated parameters are the
    total mass M = m_s + m_w and the speed v
  * mu_w = 1 reproduces the rigid moving mass of Section 2 (case_1/BMM_SS.m)
  * optional modal damping par.zeta (scalar or per mode)
  * deck profile r(x): ISO 8608:2016 PSD G_d(n) = G_d(n0) (n/n0)^-2,
    n0 = 0.1 cycles/m, sum of N = 300 cosines with uniform random phases,
    band 0.01-3 cycles/m, classes A/B/C (G_d(n0) = 16/64/256e-6 m^3)
Filter: adaptive UKF of Case 1 on the state [q; qd; u_v; u_v'; M; v]. In all
  studies of the paper the filter assumes a rigid mass, an undamped beam and a
  smooth deck (par.prof = []), so the suspension, damping and roughness of the
  data act as unmodelled physics.

Files
-----
iso8608_profile.m        ISO 8608 profile generator (harmonic amplitudes/phases)
eval_roughness.m         r(x), r'(x), r''(x) evaluated analytically
vbi_rough_rhs.m          shared right-hand side: modal matrices, sprung-mass
                         equation and contact force
BMM_SS_rough.m           truth dynamics   (state [q; qd; u_v; u_v'])
BMM_UKF_rough.m          filter dynamics  (state [q; qd; u_v; u_v'; M; v])
TF_rough.m               RK4 step of the filter model
run_UKF_rough.m          adaptive UKF loop (port of case_1/run_UKF_displacement.m);
                         opts struct: alfa, Q, v_init, P0 settings; returns the
                         innovation log-likelihood used by the speed bank
contact_force_rough.m    contact-force reconstruction from a state trajectory
ode4.m, sigmapoint.m, mode_shapesPP.m   copied unchanged from case_1
verify_iso8608.m         checks the generated profiles against ISO 8608
                         (Welch PSD, band power, RMS, class limits, IRI)
test_reduction_to_case1.m  checks that the extended model reduces exactly to
                         case_1/BMM_SS.m for mu_w = 1 on a smooth deck
                         (run with ../case_1 on the MATLAB path)

Section 6, Robustness analysis
------------------------------
BMM_UKF_Robustness_Combined.m    quarter-car truth (mu_w = 0.1, f_v = 1.5 Hz,
    zeta_v = 0.2), roughness {none, A, B, C} x beam damping {0, 2, 5} %,
    original filter model. Also computes the model-form error (filter model
    driven with the true M, v minus the true displacement).
    -> robustness_combined_results.mat      (Fig. robustness_params, Fig. robustness_errors)
BMM_UKF_Robustness_Divergence.m  roughness beyond class C (G_d(n0) = 256e-6 .. 8192e-6 m^3)
BMM_UKF_Robustness_Divergence2.m           ... and 16384e-6 .. 65536e-6 m^3
    -> robustness_divergence*_results.mat   (Fig. robustness_errors, breakdown)
BMM_UKF_Robustness_Dynamics.m    class C data with beam stiffness EI, beam mass m
    and beam damping different from the filter values
    -> robustness_dynamics_results.mat      (Fig. robustness_dynamics)
Each script runs in about 1-5 minutes and writes a *_log.txt diary with the
converged errors quoted in the text. The figures of Section 6 are drawn with
pgfplots from CSV exports of these .mat files (time histories of M and v,
model-error and R variances, converged errors); the plotting sources are not
part of this repository.

Section 7, Generality study (147 cases) and dynamic force
---------------------------------------------------------
BMM_UKF_Generality_Forces.m (part, nparts, mode)
    M = 4000:100:6000 kg x v = 13:3:31 m/s, Case-1 noise, initial mass 0.7 M,
    speed bank 10:8:34 m/s selected by the innovation log-likelihood over the
    first 10 s. mode = 'clean': truth = rigid mass, smooth, undamped;
    mode = 'error': quarter-car + class C roughness + 5 % beam damping with
    the filter unchanged. The mass list is split over nparts parallel MATLAB
    processes; the paper used 4 parts, e.g.
        for p = 1:4, BMM_UKF_Generality_Forces(p, 4, 'clean'); end
        for p = 1:4, BMM_UKF_Generality_Forces(p, 4, 'error'); end
    -> generality_forces[_error]_part<p>_results.mat
plot_generality_forces.m (mode)   Fig. dyn_force_vs_statics (mode 'clean'):
    reconstructed force history of the 5000 kg / 28 m/s case and the mean
    dynamic force vs the static weight of all cases. Also writes the per-case
    summary table generality_cases[_error]_from_forces.csv.
plot_two_cases_both.m             Fig. two_cases_both: convergence of cases
    A (4500 kg, 28 m/s) and B (5800 kg, 16 m/s) without and with model error,
    both started from the speed the bank selected in the model-error run.
plot_generality_error_map.m       Fig. generality_error_map: converged mass and
    speed errors over the grid, without and with model error. Reads the
    per-case summary tables generality_cases_clean.csv and
    generality_cases_error.csv (columns case, M_true, M_est, v_true, v_est,
    start, dM, dv), i.e. the quantities written by plot_generality_forces.m.
