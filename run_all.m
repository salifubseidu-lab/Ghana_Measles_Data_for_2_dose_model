%RUN_ALL  Reproduce every numerical result, table and figure of the manuscript
%   "Modeling and Optimizing Control of Measles in West Africa".
%   Output: figures/ (PDF and EPS), results/ (CSV/MAT) and the console log.
%   Set QUICK = true before running for a fast check with reduced settings.
%
%   Approximate run times (single core): Sections 3-6 minutes; Section 9 about 20 min
%   (PRCC, basin); Section 10 about 30-60 min; Section 8 several hours (DE, 10,000-step
%   MCMC, profiles, sensitivity refits).  Octave is several times slower.
%   Section numbers follow the manuscript; folder names keep the original numbering.
if ~exist('QUICK', 'var'), QUICK = false; end
root = fileparts(mfilename('fullpath')); cd(root); setup_paths();
diary(fullfile(root, 'results', 'run_all_log.txt'));
banner = @(s) fprintf('\n==================== %s ====================\n', s);
banner('Section 3: thresholds, elasticities');            s03_thresholds_elasticities; s03_fig_rc;
banner('Section 5: classification of equilibria');        if QUICK, NS = 300; else, NS = 2000; end
MODE = 'measles'; s05_verify_classification; MODE = 'wide'; s05_verify_classification; s05_fig_h;
banner('Section 6: bifurcation analysis');                s06_hopf_table; s06_fig_rho_bif; s06_regime_map;
if QUICK, NK = 61; else, NK = 241; end;                   s06_bif2p_codim2;
s06_fig_ts; s06_sweep_hysteresis;
if QUICK, TEND = 2000; else, TEND = 8000; end;            s06_hopf_amplitude_sim;
banner('Section 9: numerical simulations');              s10_baseline_ts; s10_vaccination_maps; s10_equilibrium_burden; s10_honeymoon;
if QUICK, NA = 6; NB = 8; NS = 200; else, NA = 14; NB = 36; NS = 2000; end
s10_basin; s10_prcc;
banner('Section 10: optimal control');                    s11_static_controls;
if QUICK, MAXIT = 50; else, MAXIT = 3000; end;           s11_optimal_control; s11_cost_effectiveness;
banner('Section 8: data');                                s08_data_summary;
banner('Section 8: parameter estimation');
if QUICK, NSTEPS = 400; DEITER = 60; SEEDS = 1; NMEVAL = 300; else, NSTEPS = 10000; DEITER = 600; SEEDS = [1 2 3]; NMEVAL = 5000; end
s09_fit_mcmc;
if QUICK, NDRAW = 300; NCF = 50; else, NDRAW = 1500; NCF = 300; end
s09_posterior_analysis; s09_solver_identity_checks;
if ~QUICK, s09_profiles; s09_sensitivity; end
if QUICK, NSTEPS_N = 600; else, NSTEPS_N = 12000; end;    s09_niamey;
s09_figures;
diary off;
