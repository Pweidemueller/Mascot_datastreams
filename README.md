# MASCOT Datastreams


Mascot Datastreams is an extension to [BEAST2](http://beast2.org)’s MASCOT package that lets you incorporate time-stamped case count data as an additional likelihood term tied to per-deme prevalence I(t) (via Mascot’s `Skygrowth`). It provides:

- A prevalence-based case-count likelihood that links counts to I(t) from Mascot `Skygrowth` (`CaseCountLikelihood`).
- A Gamma–Poisson (Negative Binomial) distribution parameterized by mean and dispersion (`GammaPoisson`).
- A prevalence-to-Ne mapping for downstream coalescent models (`PrevalenceToNeSkygrowth`).
- A spline-based prevalence-to-Ne mapping using cubic spline interpolation (`PrevalenceToNeSpline`).
- Updated example XMLs demonstrating both workflows under `examples/`.

This README summarizes the modeling assumptions, parameterization, and how to use the package from Java and from BEAST XMLs.

> Note: The prevalence-based workflow is now the default. `CaseCountLikelihood` expects a per-deme `prevalence` (`mascot.parameterdynamics.Skygrowth`) and uses I(t) as the distribution mean (optionally scaled). The older Ne-based input and legacy helper classes have been removed.


## What this package assumes

- Case counts are observed at specific times (years before the most recent sample) and associated to demes consistent with your MASCOT analysis.
- For each observation i at time $t_i$ for deme $d_i$, the expected count is linked to the corresponding prevalence: $E[count_i] = s · I_{d_i}(t_i)$, where $s$ is an optional positive scaling factor (default 1.0).
- Observations are conditionally independent given I(t) and the chosen count distribution and its overdispersion parameter(s).
- Times are on the same scale as your MASCOT tree and dynamics (years before most recent sample).
- Trait/deme indexing must be consistent across your analysis when aggregating multiple single-deme likelihoods.

Limitations to keep in mind
- The link $E[count] = s · I(t)$ is pragmatic; ensure the scaling $s$ and dispersion $\alpha$ provide adequate fit via posterior predictive checks.


## Statistical model and parameterization

Let $X$ be the observed case count for a given (deme, time) pair and let $\mu = s · I_{deme}(time)$.

- Distribution: Gamma–Poisson (Negative Binomial) with mean $\mu$ and dispersion $\alpha > 0$.
  - Shape $r = 1/\alpha$, success probability $p = r / (r + \mu)$.
  - $E[X] = \mu$.
  - $Var[X] = \mu + \alpha \mu^2$.
  - PMF: $P(X = k) = Γ(r + k) / (Γ(k+1) Γ(r)) · p^r · (1 − p)^k$.
  - CDF: $I(p; r, k + 1)$, the regularized incomplete beta.

The overall log-likelihood is the sum over observations of the log PMF at the realized counts with mean set to I(t) (optionally scaled) at each observation’s time and deme.


## How prevalence I(t) is provided

`CaseCountLikelihood` expects a single-deme Mascot `Skygrowth` component via the `prevalence` input. Its `getNeTime(t)` is interpreted as I(t) in forward time (values are exp(log I(t))). You instantiate one likelihood per deme and sum their log-probabilities.

At evaluation time, for each observation the likelihood sets the distribution mean μ to I_deme(t) multiplied by the optional positive `scaling` factor.

### Skygrowth expectations

- **Rate-shift handling**: Provide `rateShifts` and a log-scale `RealParameter` for control points (dimension = shifts + 1).
- **Time indexing**: `getNeTime(t)` expects forward times in years before the most recent sample.
- **Interpolation**: Piecewise exponential interpolation between log-prevalence control points, per Mascot’s `Skygrowth` semantics.
- **Per-deme composition**: Use one `Skygrowth` per deme and one `CaseCountLikelihood` per deme; sum the log-likelihoods.

## Mapping prevalence to Ne(t)

For downstream phylodynamic components that require Ne(t), you have two options:

### PrevalenceToNeSkygrowth (piecewise-exponential)

Use `mascotdatastreams.dynamics.PrevalenceToNeSkygrowth` for piecewise-exponential interpolation:

- **Inputs**: `prevalence` (Skygrowth, interpreted as log-prevalence), `uninfectiousRate` γ, optional `coalescentScale` c (default 2).
- **Mapping**: With forward-time slope of log I(t) denoted g_fwd, transmission_rate(t) = g_fwd + γ, and
  Ne(t) = I(t) / (c · transmission_rate(t)). Numerical clamps ensure robust behavior after the last rate shift and at extremes.

### PrevalenceToNeSpline (cubic spline)

Use `mascotdatastreams.dynamics.PrevalenceToNeSpline` for smooth cubic spline interpolation:

- **Inputs**: `logInfected` (RealParameter with log-prevalence values), `rateShifts` (RealParameter with time points), `uninfectiousRate` γ, optional `coalescentScale` c (default 2).
- **Mapping**: Uses cubic spline interpolation between rate shift points, then computes Ne(t) = I(t) / (c · transmission_rate(t)) where transmission_rate(t) = dlogI/dt + γ.
- **Advantages**: Smooth interpolation, analytical derivatives, direct log-prevalence storage, precomputed grid points for efficient and consistent lookup of both prevalence and derivatives.
- **Note**: TO DO - The numerical clamping should be revisited as it might interfere with inference/convergence. Ideally, unreasonable Ne values should be rejected by the sampler rather than clamped.


## Provided classes (API)

- `mascotdatastreams.distribution.GammaPoisson`
  - Inputs: `mean` (RealParameter), `dispersion` (RealParameter; α > 0).
  - Parameterization implements r = 1/α, p = r/(r+μ). Provides PMF, log-PMF, and CDF via Apache Commons Math.

- `mascotdatastreams.distribution.CaseCountLikelihood`
  - Inputs: `prevalence` (`mascot.parameterdynamics.Skygrowth`), `caseCounts` (RealParameter), `caseTimes` (RealParameter), `distribution` (`GammaPoisson`), optional `scaling` (RealParameter; s > 0), optional `uninfectiousRate` (unused placeholder).
  - For each observation: μ ← s · I(t) from `prevalence`, then contribute log PMF at the observed count.

- `mascotdatastreams.dynamics.PrevalenceToNeSkygrowth`
  - Inputs: `prevalence` (`mascot.parameterdynamics.Skygrowth`), `uninfectiousRate` γ (RealParameter), optional `coalescentScale` c (RealParameter; default 2).
  - Output: an `NeDynamics` compatible component where `getNeTime(t)` returns Ne(t) derived from I(t).

- `mascotdatastreams.dynamics.PrevalenceToNeSpline`
  - Inputs: `logInfected` (RealParameter with log-prevalence values), `rateShifts` (RealParameter with time points), `uninfectiousRate` γ (RealParameter), optional `coalescentScale` c (RealParameter; default 2).
  - Output: an `NeDynamics` compatible component using cubic spline interpolation for smooth prevalence trajectories with precomputed grid points for efficient and consistent lookup of both prevalence and derivatives.


## Usage from Java (authoritative example)

See `test/mascotdatastreams/distribution/CaseCountLikelihoodPrevalenceTest.java` for a complete, up-to-date example that:

- Builds a minimal Mascot setup and per-deme `Skygrowth` prevalence.
- Provides `caseCounts` and `caseTimes` directly as `RealParameter`s (no wrapper class).
- Instantiates a single `GammaPoisson` (mean is overwritten per observation by the likelihood) and optional `scaling`.
- Constructs per-deme `CaseCountLikelihood`s and sums their log probabilities.


## Conceptual BEAST XML wiring (guidance)

The example XMLs under `examples/` reflect the current API. Conceptually, a per-deme case-count likelihood looks like:

1. Define Mascot demes and trait mapping (as in standard MASCOT setups).
2. Define prevalence per deme using `mascot.parameterdynamics.Skygrowth` with `mascot.dynamics.RateShifts`.
3. Add the case-count likelihood with per-deme observations and count distribution, for example:

```xml
<distribution id="CaseCountLikelihood.Deme0" spec="mascotdatastreams.distribution.CaseCountLikelihood">
    <prevalence idref="PrevDeme0"/>
    <caseCounts idref="caseCounts.Deme0"/>
    <caseTimes idref="caseTimes.Deme0"/>
    <distribution idref="CaseCountDist"/>
    <!-- Optional scaling of the mean: mu = s * I(t) -->
    <!-- <scaling spec="beast.base.inference.parameter.RealParameter" value="0.1"/> -->
</distribution>
```

`GammaPoisson` in the same XML:

```xml
<GammaPoisson id="CaseCountDist" spec="mascotdatastreams.distribution.GammaPoisson">
    <!-- The mean is overwritten per observation from I(t); provide any positive placeholder. -->
    <mean spec="beast.base.inference.parameter.RealParameter" value="1.0"/>
    <dispersion spec="beast.base.inference.parameter.RealParameter" value="0.5"/>
</GammaPoisson>
```

Prevalence-to-Ne mapping (see `examples/prevalence_to_ne_demo.xml`):

```xml
<PrevalenceToNeSkygrowth id="NeFromPrevDeme1" spec="mascotdatastreams.dynamics.PrevalenceToNeSkygrowth">
    <prevalence idref="PrevDeme1"/>
    <uninfectiousRate spec="beast.base.inference.parameter.RealParameter" value="0.2"/>
    <coalescentScale spec="beast.base.inference.parameter.RealParameter" value="2.0"/>
</PrevalenceToNeSkygrowth>
```

Important
- The exact XML attributes (`id`, `idref`, `spec`) may vary depending on your BEAST/MASCOT versions and package registration. The semantic mapping (input names) is authoritative: `prevalence`, `caseCounts`, `caseTimes`, `distribution`, and optional `scaling` for `CaseCountLikelihood`; `mean` and `dispersion` for `GammaPoisson`; `prevalence`, `uninfectiousRate`, and optional `coalescentScale` for `PrevalenceToNeSkygrowth`.


## Compatibility notes

- The example XML files under `examples/` are up to date with the current API.
- This package depends on MASCOT classes such as `RateShifts` and `Skygrowth`. Ensure the MASCOT plugin is available on the BEAST 2 classpath.


## Priors and inference tips

- **Dispersion (α) prior**
  - Use a weakly informative prior to allow overdispersion beyond Poisson, e.g., LogNormal on α (median ~0.5, wide sd), Exponential(rate ≈ 1), or Half-Normal(σ ≈ 1) truncated to α > 0.
  - Constrain α to a reasonable domain (e.g., 1e-6 to 10) to avoid pathological proposals; in BEAST specify lower/upper bounds on the `RealParameter`.
  - Current implementation treats `dispersion` effectively as a scalar; if you pass a vector, only the first value is used. Start with a single shared α across demes.

- **Prevalence and scaling**
  - If counts are systematically larger/smaller than I(t), adjust `scaling` (s) and/or use larger α. Validate via prior/posterior predictive checks.

- **Practical tuning**
  - Initialize α in a moderate range (e.g., 0.3–2.0).
  - Ensure observation times are on the same scale (years before most recent sample) as your Mascot dynamics.
  - Use standard MASCOT priors/smoothing on log-prevalence (Skygrowth). The case-count likelihood adds information; check sensitivity so it does not dominate unduly.

## Not a Knot Spline
**OVERVIEW:**
This class provides cubic spline interpolation with not-a-knot boundary conditions,
which ensures smooth interpolation between data points. It also precomputes
values at grid points for efficient lookup.

**CUBIC SPLINE MATHEMATICS:**
A cubic spline consists of cubic polynomials S_i(t) on each interval [x[i], x[i+1]]:
S_i(t) = a[i] + b[i]*(t-x[i]) + c[i]*(t-x[i])^2 + d[i]*(t-x[i])^3

**COEFFICIENT INTERPRETATION:**
- a[i] = y[i]: Function value at knot x[i] (given data)
- b[i]: Coefficient of linear term (t-x[i])
- c[i]: Coefficient of quadratic term (t-x[i])^2 
- d[i]: Coefficient of cubic term (t-x[i])^3

The actual derivatives at knots are computed from these coefficients:
- S'(x[i]) = b[i] (first derivative at left endpoint of interval i)
- S''(x[i]) = 2*c[i] (second derivative at left endpoint of interval i)
- S'''(x[i]) = 6*d[i] (third derivative at left endpoint of interval i)

**NOT-A-KNOT CONDITION:**
The not-a-knot condition is a boundary condition that eliminates the need
for additional constraints at the endpoints. It requires:
- At x[1]: S'''_0(x[1]) = S'''_1(x[1]) (third derivative continuity)
- At x[n-1]: S'''_{n-2}(x[n-1]) = S'''_{n-1}(x[n-1]) (third derivative continuity)

This creates a tridiagonal system of equations that is solved using
LU decomposition to determine the second derivatives at all knots.

## Building and installing

- This repository contains an Ant `build.xml`. Build and add the resulting JAR to your BEAST 2 installation, alongside the MASCOT plugin. Alternatively, integrate into your development environment so that BEAST sees the package on the classpath.


## License
TBD


