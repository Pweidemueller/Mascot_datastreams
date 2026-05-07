# MASCOT Datastreams

## Max-interval ODE subdivision

`MascotLogPflag` ships an opt-in alternative outer loop, enabled by setting
`useMaxInterval="true" maxInterval="<value>" implementation="java"` on the
distribution element. Instead of stepping through the dynamics' `mascotshifts`
grid, each tree-event interval is split into `ceil(L / maxInterval)` equal
sub-intervals; per-cell rates are precomputed once per likelihood call so the
inner loop is a simple array index. On the SARS-CoV-2 999-seq dataset at the
user's published `mascotshifts` step (0.005), an equilibrated MCMC tree, and
`maxInterval = 0.02` (the tree-event-floor sweet spot), this gives ~5%
wall-time speedup at bitwise-matching `logP`. Design rationale,
testing strategy, and known follow-ups live in `MAX_INTERVAL_PLAN.md` and
`Thingstoreturnto.md`.

## Benchmarks

End-to-end reproducibility scripts under `scripts/`:

- `bash scripts/tier1.sh` — correctness sweep on a small fixed-tree XML.
  Confirms the new path converges to the same `logP` as the old grid as
  `maxInterval → 0`. Outputs in `sandbox/tier1*`.
- `bash scripts/tier2.sh` — timing sweep across two grid baselines plus a
  multi-parameter bias-stability check. Reports per-call wall time,
  doEuler-call count, and matched-accuracy speedup. Outputs in
  `sandbox/tier2*`.
- `ant intervals -Dintervals.xml=... -Dintervals.mode={old,new} -Dintervals.maxInterval=...`
  — dump the `(length, t_start, t_end, t_start_type, t_end_type)` interval
  list MascotLogPflag would iterate; useful for inspecting where each method
  spends ODE work. See `scripts/plot_intervals.py` for the histogram view.
- `ant repeat-only -Drepeat.xml=... -Drepeat.n=N` — single-XML determinism
  check; prints `logP` and `doEuler` count for each of N consecutive calls.
  Catches issues like seed leakage and class-loading non-determinism.

The harness sources are `test/mascotdatastreams/benchmark/{Tier1Runner,
Tier2Runner, IntervalDumper, RepeatabilityCheck}.java`. They live under
`test/` so they don't ship in the BEAST package jar.

## TODO

- **`GammaPoisson`: remove defensive clamp on `p`** (`GammaPoisson.java:89, :136`)
  `p = r/(r+μ)` is mathematically in (0,1) for positive r and μ, but floating-point underflow/overflow can push it to exactly 0 or 1 when μ >> r or μ << r (extreme mean or dispersion proposals during MCMC). The clamp masks this instead of handling it correctly. Fix: reparameterise `GammaPoissonDistributionImpl` to accept (r, μ) directly and compute log-probabilities in log space — `logP = log(r) − log(r+μ)` and `log1mP = log(μ) − log(r+μ)` — which are always finite for positive inputs and make the clamp unnecessary.

