# Max-interval ODE subdivision — implementation & testing plan

## Goal

Replace the user-defined fixed grid (`gridRateShifts`) as the *ODE-step subdivider* with a single parameter `maxInterval`: the maximum length of any single Euler/RK4 integration step inside `MascotLogPflag.calculateLogP()`. Tree events (coalescent / sampling) remain natural breakpoints; between events, an interval longer than `maxInterval` is split into `ceil(L / maxInterval)` equal pieces.

Rates are still queried at the **midpoint** of each subinterval (matches the convention used today by `StructuredSkylinePrevalence.getCoalescentRate(int i)`), so behaviour converges to the existing implementation as resolution increases.

The change is opt-in. The current code path is preserved verbatim and remains the default.

---

## Scope of code changes

Three files. No changes to the abstract `mascot.dynamics.Dynamics` interface or to the Mascot JNI / native paths.

### 1. `src/mascotdatastreams/distribution/StructuredSkylinePrevalence.java`

Add **two** public methods that take time directly. They duplicate the rate-construction logic from the existing index-based methods, but skip the `i → midpoint(i)` step and use the supplied `t`.

```java
public double[] getCoalescentRateAtTime(double t) {
    double[] coal = new double[getDimension()];
    for (int j = 0; j < coal.length; j++) {
        if (parametricFunction.get(j).isTime) {
            coal[j] = Math.min(maxRateInput.get(),
                               1.0 / parametricFunction.get(j).getNeTime(t));
        } else {
            int intervalNr = findIntervalIndex(t);
            coal[j] = Math.min(maxRateInput.get(),
                               1.0 / parametricFunction.get(j).getNeInterval(intervalNr));
        }
    }
    return coal;
}

public double[] getBackwardsMigrationAtTime(double t) {
    int dim = dimensionInput.get();
    double[] m = new double[dim * dim];
    double[] NeInt = new double[dim];

    for (int j = 0; j < dim; j++) {
        if (parametricFunction.get(j).isTime) {
            if (parametricFunction.get(j) instanceof SplinePrevalenceToNe) {
                NeInt[j] = ((SplinePrevalenceToNe) parametricFunction.get(j))
                           .getPrevalenceTime(t);
            } else {
                NeInt[j] = parametricFunction.get(j).getNeTime(t);
            }
        } else {
            NeInt[j] = parametricFunction.get(j).getNeInterval(findIntervalIndex(t));
        }
    }
    // identical matrix construction loop as in getBackwardsMigration(int i)
    for (int a = 0; a < dim; a++) {
        for (int b = 0; b < dim; b++) {
            if (a != b) {
                if (indicatorInput.get().getArrayValue(dirs[b][a]) > 0.5)
                    m[a * dim + b] = Math.min(maxRateInput.get(),
                        NeInt[b] * migration.getArrayValue(dirs[b][a]) / NeInt[a]);
                else
                    m[a * dim + b] = 0.0;
            }
        }
    }
    return m;
}

private int findIntervalIndex(double t) {
    // binary search over rateShiftsInput on absolute backward time
    // returns intervalNr clamped to [firstlargerzero, dim-2] to match the
    // convention used in getCoalescentRate(int i)
}
```

Refactor existing `getCoalescentRate(int i)` and `getBackwardsMigration(int i)` to delegate:

```java
@Override
public double[] getCoalescentRate(int i) {
    int intervalNr = mapIndexToIntervalNr(i);
    return getCoalescentRateAtTime(rateShiftsInput.get().getIntervalMidpoint(intervalNr));
}
```

This guarantees zero behavioural drift on the old path — same time, same result.

### 2. `src/mascotdatastreams/distribution/MascotLogPflag.java`

**Two new inputs:**

```java
public Input<Boolean> useMaxIntervalInput = new Input<>("useMaxInterval",
    "If true, ignore dynamics grid for ODE subdivision; subdivide each tree-event interval into ceil(L / maxInterval) equal pieces and evaluate rates at each subinterval midpoint.",
    false);

public Input<Double> maxIntervalInput = new Input<>("maxInterval",
    "Maximum length of any single ODE step when useMaxInterval=true.",
    Double.POSITIVE_INFINITY);
```

**Branch in `calculateLogP()`:**

```java
public double calculateLogP() {
    if (useMaxIntervalInput.get()) {
        return calculateLogP_maxInterval();
    }
    // existing body unchanged
    ...
}
```

**New method `calculateLogP_maxInterval()`** (forked from the existing loop, ~150 lines):

```java
double currentTime = 0.0;                 // backward from present
int treeInterval = 0;
double nextTreeEvent = treeIntervals.getInterval(treeInterval);
nrLineages = activeLineages.size();
linProbsLength = nrLineages * states;

if (!computeLikelihood) { first++; return 0; }

double maxInt = maxIntervalInput.get();
StructuredSkylinePrevalence d = (StructuredSkylinePrevalence) dynamics;

do {
    double L = nextTreeEvent;
    if (L > 0) {
        int n = (L > maxInt && maxInt > 0)
                ? (int) Math.ceil(L / maxInt) : 1;
        double subL = L / n;
        for (int k = 0; k < n; k++) {
            double tMid = currentTime + (k + 0.5) * subL;
            coalescentRates = d.getCoalescentRateAtTime(tMid);
            double[] migrationRates = d.getBackwardsMigrationAtTime(tMid);
            int[] indicators = d.getIndicators(0); // structurally constant
            logP += doEulerAtTime(subL, coalescentRates, migrationRates, indicators);
            if (logP == Double.NEGATIVE_INFINITY) return logP;
        }
        currentTime += L;
    }

    // handle the tree event at currentTime
    IntervalType type = treeIntervals.getIntervalType(treeInterval);
    if (type == IntervalType.COALESCENT) {
        nrLineages--;
        // ratesInterval is no longer meaningful here; pass 0 or refactor coalesce()
        logP += coalesceAtTime(treeInterval, currentTime);
    } else if (type == IntervalType.SAMPLE) {
        nrLineages++;
        sampleAtTime(treeInterval, currentTime);
    }
    treeInterval++;
    try { nextTreeEvent = treeIntervals.getInterval(treeInterval); }
    catch (Exception e) { break; }
} while (nextTreeEvent <= Double.POSITIVE_INFINITY);

first++;
return logP;
```

**Helper methods to add:**
- `doEulerAtTime(double dt, double[] coalRates, double[] migRates, int[] indicators)` — wraps the Euler/RK4 integrator without referring to a `ratesInterval` index.
- `coalesceAtTime(int treeInterval, double t)`, `sampleAtTime(int treeInterval, double t)` — versions of `coalesce`/`sample` that don't depend on `ratesInterval` or `nextRateShift`. Most likely the existing methods only need a small refactor (any reference to `ratesInterval` removed or replaced by the absolute time).

**Cache & implementation guards (initAndValidate or top of the new method):**

```java
if (useMaxIntervalInput.get()) {
    if (cacheInput.get())
        throw new IllegalArgumentException("useMaxInterval=true is incompatible with useCache=true (initial version)");
    if (implementationInput.get() != MascotImplementation.java)
        throw new IllegalArgumentException("useMaxInterval=true requires implementation=\"java\"");
    if (!(dynamics instanceof StructuredSkylinePrevalence))
        throw new IllegalArgumentException("useMaxInterval=true requires StructuredSkylinePrevalence dynamics");
}
```

Cache disablement keeps the first version safe; the `coalRatesInterval` / `nextRateShifts` arrays are simply not exercised on the new path. We can revisit re-enabling caching once correctness is verified.

### 3. `src/mascotdatastreams/dynamics/Spline.java` and `SplinePrevalenceToNe.java`

**No changes required.** Their public time-based API (`getPrevalence(t)`, `getTransmissionRate(t)`, `getNeTime(t)`) is already what the new path needs.

---

## Backward compatibility

- Default `useMaxInterval=false` ⇒ identical code path, identical results.
- Existing XMLs run unchanged. The two example XMLs do **not** need to be modified to test the new path; instead, copies that add `useMaxInterval="true" maxInterval="..."` to the `<distribution id="Mascot...">` element are created in `examples/benchmarks/` (see Tier-2 below).
- The refactor of `getCoalescentRate(int i)` / `getBackwardsMigration(int i)` to delegate to the new `*AtTime` methods is a pure refactor — same midpoint, same maxRate clipping. A quick unit test should confirm bitwise equality on a fixed input before touching the loop.

---

## Testing strategy

Three tiers, escalating in cost. Stop at any tier if results don't justify continuing.

### Tier 1 — Correctness, single-evaluation (minutes, on laptop)

**Dataset:** `examples/6_2_simulation_datastreams.xml` (125 seqs, 2 demes).

**Procedure:**

1. Fix all parameters at their initial values; disable MCMC (use BEAST's `-validate` or just compute `logP` once).
2. Compute `logP_old` on the original XML — call this our baseline at the user's chosen grid resolution.
3. Compute `logP_old_fine` on a copy of the XML with `gridRateShifts` densified by ~10× — this is the "converged old" reference; if `|logP_old - logP_old_fine|` is small, the user's grid was already adequate.
4. Compute `logP_new(maxInterval)` on a copy of the XML with `useMaxInterval="true"` for `maxInterval ∈ {1e-2, 1e-3, 1e-4, 1e-5}` (tree height ≈ 0.15, current grid spacing ≈ 2e-4, so this brackets the user's grid both ways).
5. Plot `logP_new(maxInterval)` vs `maxInterval`. Expect monotone convergence to `logP_old_fine` as `maxInterval → 0`.

**Acceptance:**
- At `maxInterval ≈ tree_height / 1000` (~1.5e-4): `|logP_new - logP_old_fine| < 0.1` (well below per-MCMC-step Hastings-ratio noise).
- Convergence curve is monotone — no oscillation.

**Add as a JUnit test:** `test/mascotdatastreams/distribution/MaxIntervalConvergenceTest.java`. Loads a fixed minimal XML, computes `logP` on both paths at fixed parameters, asserts agreement.

### Tier 2 — Per-call timing benchmark (tens of minutes, on laptop)

**Datasets:** both XMLs.
- Small: `6_2_simulation_datastreams.xml` (125 seqs, 2 demes) — sanity check, confirms the new path isn't *slower* on small problems.
- Large: `SARSCoV2_Epsilon_BayArea_results_1000seq_datastreams_noclip.xml` (999 seqs, 4 demes) — primary speedup target.

**Procedure:**

1. Build a small Java harness `test/mascotdatastreams/benchmark/MascotLogPBenchmark.java` that:
   - Loads an XML, initialises `MascotLogPflag`.
   - Warm-up: 10 calls to `calculateLogP()` (JIT settles).
   - Timed: 100 calls. Record mean ± stdev wall time.
   - Also records the **number of `doEuler` invocations per `calculateLogP()` call** — this is the cleanest mechanism-level metric and avoids JIT noise.
2. Run on each XML at:
   - **Old path:** the XML's existing grid, plus a 2× densified and 0.5× sparsified grid.
   - **New path:** `maxInterval ∈ {1e-3, 1e-4, 1e-5}` (small XML) and `maxInterval ∈ {tree_height/100, /500, /2000}` (large XML).
3. Cross-tabulate: for each `maxInterval`, what's the *equivalent* old-grid resolution that gives matching accuracy? Compare wall time at matched accuracy, not at matched parameter.

**Output:** a small `benchmarks/results.csv` with columns `dataset, mode, parameter, n_doEuler_per_call, mean_ms, std_ms, logP`.

**Decision:** speedup at matched accuracy must be ≥ 1.5× on the large XML to justify Tier 3. If it's < 1.2× the idea isn't worth pursuing further.

**Total runtime estimate:** 100 calls on the 999-seq tree at 4 demes is roughly 5–30 minutes per configuration depending on grid density. Eight configurations ⇒ 1–4 hours total. No remote machine needed.

### Tier 3 — Short MCMC equivalence (hours, only if Tier 2 passes)

**Dataset:** `SARSCoV2_Epsilon_BayArea_results_1000seq_datastreams_noclip.xml`.

**Procedure:**

1. Reduce `chainLength` from 1e7 to 1e6 (or whatever fits in ~2–4 h wall time on the available machine), keep everything else identical.
2. Run **two** chains in parallel, identical seeds:
   - Run A: original XML.
   - Run B: copy with `useMaxInterval="true" maxInterval=<value chosen from Tier 2>`.
3. Compare:
   - **ESS/hour** for `posterior`, `Mascot.t:SimDataset` likelihood, all migration rates, all `NeScaler.*`, `caseCounts.dispersion`. Use Tracer or `coda::effectiveSize()`.
   - **Posterior summaries** (mean, median, 95% HPD) for the same parameters. Should overlap within MCMC noise.
   - **Trace-plot visual check** for any obvious bias.

**Acceptance:**
- ESS/hour for the new path is at least as good as the old path (ideally 1.5×+, matching Tier-2 timing).
- All 95% HPDs overlap; mean estimates within ~1 SE of each other.

**If accepted:** the new path is correct *and* faster end-to-end on a real workload — the change is ready to consider as default-on (separate decision; not part of this plan).

---

## Files to create

| File | Purpose |
|---|---|
| `src/mascotdatastreams/distribution/StructuredSkylinePrevalence.java` | **Modify**: add `getCoalescentRateAtTime(double)`, `getBackwardsMigrationAtTime(double)`, `findIntervalIndex(double)`. Refactor existing index-based methods to delegate. |
| `src/mascotdatastreams/distribution/MascotLogPflag.java` | **Modify**: add `useMaxInterval` and `maxInterval` inputs, branch at top of `calculateLogP()`, add `calculateLogP_maxInterval()`, `doEulerAtTime(...)`, `coalesceAtTime(...)`, `sampleAtTime(...)`. Add validation guards. |
| `test/mascotdatastreams/distribution/MaxIntervalConvergenceTest.java` | **Create**: Tier-1 JUnit test. Asserts new path converges to old-fine baseline as `maxInterval → 0` on the 6_2 simulation. |
| `test/mascotdatastreams/benchmark/MascotLogPBenchmark.java` | **Create**: Tier-2 timing harness. Loads XML, runs warm-up + timed calls, prints CSV row. Run from `ant` or a small shell script. |
| `examples/benchmarks/6_2_maxInterval_*.xml` | **Create**: copies of `6_2_simulation_datastreams.xml` with `useMaxInterval="true"` and varied `maxInterval`. One per Tier-1 / Tier-2 config. |
| `examples/benchmarks/SARSCoV2_maxInterval_*.xml` | **Create**: same for the large XML, plus a `chainLength=1e6` variant for Tier 3. |
| `examples/benchmarks/results.csv` | **Output** of Tier 2; not committed until results stable. |

## Files NOT to change

- `src/mascotdatastreams/dynamics/Spline.java`
- `src/mascotdatastreams/dynamics/SplinePrevalenceToNe.java`
- `mascot.dynamics.Dynamics` (abstract class in the upstream Mascot jar)
- Any of the JNI / native (`MascotImplementation.allnative`, `indicators`) paths — they remain untouched and the new flag rejects them at validation time.

---

## Open questions / decisions deferred

1. **Midpoint vs. left-endpoint rate evaluation.** Plan defaults to midpoint to match the existing convention. Worth comparing once on Tier 1 to confirm midpoint really is more accurate per ODE step than left-endpoint — if not, left-endpoint simplifies bookkeeping slightly.
2. **Re-enabling cache on the new path.** Disabled in v1. Re-enabling requires re-keying `coalLogP[]` / `coalLinProbs[]` by tree node rather than by `ratesInterval`. Worth doing only if Tier 3 shows ESS/hour close to baseline despite cache being off.
3. **`maxInterval` default value.** `Double.POSITIVE_INFINITY` is safe (means "never subdivide between tree events") but useless; in practice we'd want a sensible default like `tree_height / 1000`. Defer until Tier-2 results suggest a robust rule of thumb.
4. **Whether `getInterval(int)` should still return real-time durations on the new path.** The new path doesn't call it, so no decision needed. Old path is unchanged.
