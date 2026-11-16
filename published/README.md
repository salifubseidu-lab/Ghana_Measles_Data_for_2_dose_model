# Published estimation outputs

`posterior_draws.csv`
: 8,140 posterior draws from the affine-invariant ensemble sampler (44 walkers, 10,000 iterations, first 5,000 discarded, thinned by half the shortest autocorrelation time). One row per draw, parameters on their natural scale: `R0`, `A_s` (seasonal amplitude), `t_s` (seasonal phase, years), `iota` (imported infectious-equivalents), `rho` (reporting probability), `s0` (initial effective susceptible fraction), `omega1` (first-dose efficacy), `phi_NB` (negative binomial size), `lambda_ref` (reference force of infection, per year), `c_npi` (contact reduction under COVID-19 measures), `h_SIA` (campaign reach factor).

`posterior_summary.csv`
: Posterior median and 95% credible interval (2.5% and 97.5% quantiles) of each parameter, plus `theta = 1 - omega1` and the seasonal peak month `1 + 12 t_s`. These are the values in Table 4 of the paper. They were computed from the full retained chain, so quantiles recomputed from the thinned draws differ slightly (Monte Carlo error).

`python_reference_estimation.json`
: Posterior mode, log-likelihood (-721.259) and log-posterior (-759.127) from the original Python run, with a few intermediate quantities. Use it to check that the MATLAB code reproduces the reference fit.
