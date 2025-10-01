# MASCOT Datastreams


Mascot Datastreams is an extension to [BEAST2](http://beast2.org)’s MASCOT package that lets you incorporate time-stamped case count data as an additional likelihood term tied to the effective population size Ne(t) per deme. It provides:

- A light-weight data container for case counts with times and deme indices (`CaseCountData`).
- A likelihood that connects case counts to Ne(t) via a parametric count distribution (`CaseCountLikelihood`).
- A recommended Gamma–Poisson (Negative Binomial) distribution parameterized by mean and dispersion (`GammaPoisson`).
- A legacy/compatibility Negative Binomial wrapper (`NegativeBinomialDistribution`).

This README summarizes the modeling assumptions, parameterization, and how to use the package from Java and conceptually from BEAST XMLs.

> Note: A prevalence-based workflow is now available and preferred for case counts. `CaseCountLikelihood` expects prevalence inputs (`prevalence`) and treats the mean as I(t) = exp(log I(t)). The legacy Ne-based input is deprecated.


## What this package assumes

- Case counts are observed at specific times (in years, relative to the most recent sample) and assigned to demes (traits) consistent with your MASCOT analysis.
- For each observation i with time t_i and trait/deme d_i, the expected count is linked to the corresponding effective population size: E[count_i] = Ne_d_i(t_i).
- Observations are conditionally independent given Ne(t) and the chosen count distribution and its overdispersion parameter(s).
- Times are provided on the same absolute time scale as your MASCOT tree and Ne dynamics (years before the most recent sample).
- Trait indices are 0-based integers that refer to the demes in your MASCOT analysis.

Limitations to keep in mind
- Using Ne directly as the mean for case counts can underfit very large observed counts if plausible Ne is constrained by the phylogeny. In such settings a different link (e.g., via prevalence) may be preferable in future versions.


## Statistical model and parameterization

Let X be the observed case count for a given (deme, time) pair and let μ = Ne_deme(time).

- Default distribution: Gamma–Poisson (Negative Binomial) with mean μ and dispersion α > 0.
  - Shape r = 1/α, success probability p = r / (r + μ).
  - E[X] = μ.
  - Var[X] = μ + α μ².
  - PMF: P(X = k) = Γ(r + k) / (Γ(k+1) Γ(r)) · p^r · (1 − p)^k.
  - CDF: I(p; r, k + 1), the regularized incomplete beta.
- Legacy option: `NegativeBinomialDistribution` (Pascal) parameterized by (mean, α) but internally converting to integer r ≈ round(1/α). Prefer `GammaPoisson` because it supports non-integer r.

The overall log-likelihood is the sum over observations of the log PMF at the realized counts with mean set to Ne at the observation’s time and deme.


## How Ne(t) is used

`CaseCountLikelihood` expects an `NeDynamicsList` from MASCOT that provides one Ne dynamic per deme and exposes `getNeTime(t)` to query Ne at arbitrary times. In practice you will construct this using MASCOT components, e.g. a Skyline/Skygrowth parameterization with `RateShifts` and then assemble them into an `InitializedNeDynamicsList`.

At evaluation time, for each observation the likelihood sets the distribution mean μ to the current Ne for that deme and time.

### Current Ne dynamics expectations

- **Rate-shift handling**: `NeDynamicsList` stores one `NeDynamics` per deme; each `Skygrowth` instance binds a log-scale `RealParameter` and matching `RateShifts`. The parameter vector is dimension `rateShifts + 1`, describing log Ne values at breakpoints.
- **Time indexing**: `getNeTime(t)` expects *forward* times measured in years before the most recent sample (identical scale as tree heights). Internally, `Skygrowth.getNeTime(t)` locates the interval containing `t` by scanning `RateShifts` in ascending order.
- **Interpolation**: `Skygrowth` converts the log-scale control points into per-interval growth rates (piecewise exponential interpolation). Within interval `i`, it evaluates `Ne(t) = exp(logNe_i - growth_i · (t - rateShift_{i-1}))`, matching MASCOT’s convention of backward-time smoothing.
- **Per-deme access**: callers select the relevant deme with `NeDynamicsList.get(demeIndex)`; demes must be ordered consistently with trait indices. `InitializedNeDynamicsList` relaxes the input requirement so BEAUti can instantiate empty lists before values are bound.


## Provided classes (API)

- `mascotdatastreams.distribution.CaseCountData`
  - Inputs: `caseCounts` (RealParameter), `observationTimes` (RealParameter), `traitIndices` (RealParameter; converted to 0-based ints).
  - All three arrays must have identical length. Non-negative counts are required.
  - Times are in years relative to the most recent sample.

- `mascotdatastreams.distribution.GammaPoisson`
  - Inputs: `mean` (RealParameter), `dispersion` (RealParameter; α > 0).
  - Parameterization implements r = 1/α, p = r/(r+μ). Provides PMF, log-PMF, and CDF via Apache Commons Math special functions.

- `mascotdatastreams.distribution.NegativeBinomialDistribution` (legacy/compatibility)
  - Inputs: `mean` (RealParameter), `dispersion` (RealParameter; α > 0).
  - Internally casts r ≈ round(1/α) to satisfy the Pascal distribution’s integer requirement.

- `mascotdatastreams.distribution.CaseCountLikelihood`
  - Inputs: `prevalence` (`mascotdatastreams.dynamics.PrevalenceDynamicsList`), `caseCounts` (`CaseCountData`), `distribution` (`ParametricDistribution`, typically `GammaPoisson`).
  - For each observation: μ ← I_deme(t) = exp(log I(t)), then contributes log PMF at the observed count.

- `mascotdatastreams.dynamics.PrevalenceDynamicsList`
  - Container for per-deme log-prevalence dynamics.

- `mascotdatastreams.dynamics.PrevalenceSkygrowth`
  - Log-prevalence skyline with `mascot.dynamics.RateShifts`, mirroring `mascot.parameterdynamics.Skygrowth` semantics.


## Usage from Java (authoritative example)

See `test/mascotdatastreams/distribution/CaseCountLikelihoodTest.java` for a complete, up-to-date example that:

- Builds a small tree, demes (traits), and MASCOT Ne dynamics using `RateShifts` and `Skygrowth`.
- Creates `CaseCountData` for two demes across multiple time points.
- Instantiates `GammaPoisson` with a dispersion parameter (the mean is overwritten per observation).
- Constructs `CaseCountLikelihood` and evaluates its log probability.


## Conceptual BEAST XML wiring (guidance)

While the provided example XMLs in this repository are outdated, the conceptual wiring in a BEAST 2 XML looks like this (class names shown as `spec` for clarity; adjust to your BEAST 2 version and package registration):

1. Define MASCOT demes and trait mapping (as in standard MASCOT setups).
2. Define Ne dynamics per deme (e.g., `mascot.parameterdynamics.Skygrowth` with `mascot.dynamics.RateShifts`) and assemble them into an `mascot.util.InitializedNeDynamicsList`.
3. Add the case-count likelihood with its data and count distribution, for example:

```xml
<distribution id="CaseCountLikelihood" spec="mascotdatastreams.distribution.CaseCountLikelihood">
    <NeDynamics idref="NeList"/>
    <caseCounts id="CaseData" spec="mascotdatastreams.distribution.CaseCountData">
        <caseCounts spec="beast.base.inference.parameter.RealParameter" value="5 7 6 8 10 3 4 5 6 8"/>
        <observationTimes spec="beast.base.inference.parameter.RealParameter" value="0.0 0.1 0.11 0.24 0.3  0.0 0.15 0.16 0.25 0.37"/>
        <traitIndices spec="beast.base.inference.parameter.RealParameter" value="0 0 0 0 0  1 1 1 1 1"/>
    </caseCounts>
    <distribution id="CaseCountDist" spec="mascotdatastreams.distribution.GammaPoisson">
        <!-- The mean is overwritten per observation from Ne_deme(t); provide any positive placeholder. -->
        <mean spec="beast.base.inference.parameter.RealParameter" value="1.0"/>
        <dispersion spec="beast.base.inference.parameter.RealParameter" value="0.5"/>
    </distribution>
</distribution>
```

Prevalence-based wiring (preferred):

```xml
<distribution id="CaseCountLikelihood" spec="mascotdatastreams.distribution.CaseCountLikelihood">
    <prevalence idref="PrevList"/>
    <caseCounts id="CaseData" spec="mascotdatastreams.distribution.CaseCountData">
        <caseCounts spec="beast.base.inference.parameter.RealParameter" value="5 7 6 8 10 3 4 5 6 8"/>
        <observationTimes spec="beast.base.inference.parameter.RealParameter" value="0.0 0.1 0.11 0.24 0.3  0.0 0.15 0.16 0.25 0.37"/>
        <traitIndices spec="beast.base.inference.parameter.RealParameter" value="0 0 0 0 0  1 1 1 1 1"/>
    </caseCounts>
    <distribution id="CaseCountDist" spec="mascotdatastreams.distribution.GammaPoisson">
        <!-- The mean is overwritten per observation from I_deme(t); provide any positive placeholder. -->
        <mean spec="beast.base.inference.parameter.RealParameter" value="1.0"/>
        <dispersion spec="beast.base.inference.parameter.RealParameter" value="0.5"/>
    </distribution>
</distribution>
```

Important
- The exact XML attributes (`id`, `idref`, `spec`) may vary depending on your BEAST/MASCOT versions and package registration. The semantic mapping (input names) is authoritative: `NeDynamics`, `caseCounts`, and `distribution` for `CaseCountLikelihood`; `caseCounts`, `observationTimes`, and `traitIndices` for `CaseCountData`; `mean` and `dispersion` for `GammaPoisson`.


## Compatibility notes

- The example XML files under `examples/` correspond to an older API and should not be used as-is. Prefer the Java test as a reference for current semantics and adapt your XML accordingly.
- This package depends on MASCOT classes such as `RateShifts`, `Skygrowth`, and an `NeDynamicsList` implementation. Ensure the MASCOT plugin is available on the BEAST 2 classpath.


## Priors and inference tips

- **Dispersion (α) prior**
  - Use a weakly informative prior to allow overdispersion beyond Poisson: for example, LogNormal on α (median ~0.5, wide sd), Exponential(rate ≈ 1), or Half-Normal(σ ≈ 1) truncated to α > 0.
  - Constrain α to a reasonable domain (e.g., 1e-6 to 10) to avoid pathological proposals; in BEAST specify lower/upper bounds on the `RealParameter`.
  - Current implementation treats `dispersion` effectively as a scalar; if you pass a vector, only the first value is used. Start with a single shared α across demes.

- **Interpreting Ne vs. prevalence**
  - Ne is a coalescent effective population size, not prevalence. Using μ = Ne(t) as the expected count is a pragmatic link; for very large observed counts this can be inadequate if phylogenetically plausible Ne is small.
  - If you frequently observe very large counts, prefer generous overdispersion (larger α) and perform prior/posterior predictive checks. A future model may link counts to prevalence or include an explicit scaling, which is not currently in this package.

- **Practical tuning**
  - Initialize α in a moderate range (e.g., 0.3–2.0) to avoid near-Poisson behavior that can hinder mixing.
  - Ensure observation times are on the same scale (years before most recent sample) as your MASCOT Ne dynamics, and verify trait indices are 0-based and match deme ordering.
  - Use standard MASCOT priors/smoothing on logNe (e.g., Skygrowth smoothing). The case-count likelihood adds information; check that it does not dominate unduly via sensitivity analyses.


## Building and installing

- This repository contains an Ant `build.xml`. Build and add the resulting JAR to your BEAST 2 installation, alongside the MASCOT plugin. Alternatively, integrate into your development environment so that BEAST sees the package on the classpath.


## License
TBD

