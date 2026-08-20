# MASCOT-DS Tutorial

This tutorial walks through the example analysis in [`examples/2deme_fixedtree.xml`](examples/2deme_fixedtree.xml) to explain how a MASCOT-DS XML is put together and how to adapt it to your own data.

We are still working on BEAUti support. Instead, you start from a working example XML and edit it directly in a text editor. This tutorial explains what each block in the example does, so you know exactly what to change for your own dataset.

## 1. What the example analysis is
`2deme_fixedtree.xml` analyses a simulated outbreak with two demes (`I0`/`I1`, referred to below as Deme1/Deme2) on a **fixed genealogy** (161 tips, supplied as a Newick string, not sampled during MCMC). For each deme, three surveillance data streams are fit jointly with the tree (aka geneology) to infer prevalence trajectories and between deme migration rates:

- case counts (`caseCounts.Deme1`/`Deme2`)
- seroprevalence surveys (`seroTestedCounts.Deme1`/`Deme2`, `seroWithAntibodiesCounts.Deme1`/`Deme2`)
- wastewater viral concentrations (`wastewaterConcentration.Deme1`/`Deme2`)

All three are linked to a per-deme **prevalence spline**, which also determines each deme's effective population size in the structured-coalescent tree likelihood. Migration between demes is estimated as a forward-in-time rate (`migrationRatesSkyline`, dimension 2 for a 2-deme model).

Because the tree is fixed, this example isolates how well the data streams constrain the shared prevalence trajectory and the migration rate. Users can easily add tree inference on provided sequences by using separate functions as part of the wider BEAST2 ecosystem (not part of MASCOT-DS). For example, users could use BEAUTi to set up an XML for tree inference based on their pathogen genetic sequences (e.g. https://github.com/nicfel/targetedbeast/) and then add the MASCOT-DS functionality straight into that produced XML.

## 2. Anatomy of the XML

### 2.1 Sequence data and tree

```xml
<data id="SimDataset" spec="Alignment" name="alignment">
  <sequence id="seq_leaf_0" spec="Sequence" taxon="leaf_0" totalcount="4" value="????"/>
  ...
</data>
```
The alignment here is a placeholder — every sequence is `????` (missing data). MASCOT-DS uses the tree for its structured-coalescent likelihood, not the sequence data directly, so when the tree is fixed the alignment only needs to define the taxon set. If you *are* estimating the tree from sequence data, replace this block with your real alignment and add a substitution/site model + tree likelihood as in a normal BEAST2 analysis.

The tree itself is fixed via a `TreeParser`:

```xml
<init spec="beast.base.evolution.tree.TreeParser" id="NewickTree.t:SimDataset"
      adjustTipHeights="false" initial="@Tree.t:SimDataset" taxa="@SimDataset"
      IsLabelledNewick="true" newick="(...);"/>
```

To use your own fixed tree, replace the `newick` string (tip labels must match your taxon names) or remove this `<init>` block entirely and add tree operators if you want to co-estimate the tree.

Each tip's deme assignment is given by a `typeTraitSet`:

```xml
<typeTrait id="typeTraitSet.t:SimDataset" spec="mascot.util.InitializedTraitSet"
           traitname="type" value="leaf_0=I0,leaf_1=I0,...">
  <taxa id="TaxonSet.1" spec="TaxonSet" alignment="@SimDataset"/>
</typeTrait>
```

`value` is a comma-separated `taxon=deme` list. Deme names (`I0`, `I1`, ...) are arbitrary labels — just be consistent with them elsewhere (e.g. logger names, `dimension` counts).

### 2.2 Epidemiological data streams (per deme)

For each deme, three `RealParameter` groups hold the raw surveillance data as **fixed** (`estimate="false"`) state nodes — a *counts* vector and a matching *times* vector (times are backwards in time relative to the most recent sample in the tree, same units as the tree, usually years):

```xml
<parameter id="caseCounts.Deme1:SimDataset" ... estimate="false">5 13 4 6 ...</parameter>
<parameter id="caseTimes.Deme1:SimDataset" ... estimate="false">0.0010 0.0038 ...</parameter>

<parameter id="seroTestedCounts.Deme1:SimDataset" ... estimate="false">51 693 831</parameter>
<parameter id="seroWithAntibodiesCounts.Deme1:SimDataset" ... estimate="false">43 133 0</parameter>
<parameter id="seroTestedTimes.Deme1:SimDataset" ... estimate="false">0.0339 0.1435 0.2147</parameter>

<parameter id="wastewaterConcentration.Deme1:SimDataset" ... estimate="false">0.162 0.050 ...</parameter>
<parameter id="wastewaterConcentrationTimes.Deme1:SimDataset" ... estimate="false">0.0010 0.0038 ...</parameter>
```

To use your own data: replace these vectors with your observed values and observation times (same length, one time per observation). You can omit any of the three data streams for a deme simply by deleting its `<parameter>` state nodes *and* the corresponding likelihood/prior/operator/logger blocks described below — MASCOT-DS does not require all three streams to be present, however excluding certain data streams might alter the interpreation of e.g. the inferred magnitude of the prevalence trajectories.

`PopSize.Deme1`/`Deme2` is the (fixed) census population size of each deme, used to convert seroprevalence/wastewater signal to absolute prevalence.

### 2.3 The prevalence spline and Ne dynamics

The core dynamics object is the `spline`, one per deme:

```xml
<spline id="splinePrev.Deme1.t:SimDataset" spec="mascotdatastreams.dynamics.Spline" clipTransRate="false">
  <logInfected idref="SkylinePrev.Deme1.t:SimDataset"/>
  <rateShifts id="SkygrowthRateShifts" spec="mascot.dynamics.RateShifts">0.0000 0.0268 ... 0.2678</rateShifts>
  <gridRateShifts id="SplineGridRateShifts" spec="mascot.dynamics.RateShifts">0.0000 0.0003 ... 0.2678</gridRateShifts>
  <uninfectiousRate idref="uninfectiousRate.t:SimDataset"/>
</spline>
```

- `logInfected` (`SkylinePrev.Deme1`) is the estimated state node: log-prevalence at each of the coarse `rateShifts` knots (dimension 11 here, one value per `SkygrowthRateShifts` entry).
- `rateShifts` are the (coarse) knot times where `logInfected` is estimated; `gridRateShifts` is the (finer) grid the spline is evaluated on internally — the same grid is shared across demes.
- `uninfectiousRate` is the recovery/becoming uninfectious rate (make sure to provide this in the same units as the other parameters, most likely in years).

This spline feeds `SplinePrevalenceToNe`, which converts the prevalence trajectory into an effective population size for the structured coalescent, given migration rates into the deme:

```xml
<neDynamics id="NeDynamics.Deme1.t:SimDataset" spec="mascotdatastreams.dynamics.SplinePrevalenceToNe">
  <spline idref="splinePrev.Deme1.t:SimDataset"/>
</neDynamics>
```

Both demes' `NeDynamics` are collected under `StructuredSkylinePrevalence`, which is the `dynamics` fed to the tree likelihood (`MascotLogPflag`, MASCOT-DS's drop-in replacement for MASCOT's standard `Mascot` distribution):

```xml
<distribution id="Mascot.t:SimDataset" spec="mascotdatastreams.distribution.MascotLogPflag"
              tree="@Tree.t:SimDataset" compute_likelihood="true">
  <dynamics id="StructuredSkyline.t:SimDataset" spec="mascotdatastreams.distribution.StructuredSkylinePrevalence"
            dimension="2" forwardsMigration="@migrationRatesSkyline.t:SimDataset" maxRate="1000">
    <NeDynamics id="NeDynamicsList.t:SimDataset" spec="mascot.util.InitializedNeDynamicsList">
      <neDynamics idref="NeDynamics.Deme1.t:SimDataset"/>
      <neDynamics idref="NeDynamics.Deme2.t:SimDataset"/>
    </NeDynamics>
    <rateShifts idref="SplineGridRateShifts"/>
    <indicators id="indicatorsSkyline.t:SimDataset" spec="parameter.BooleanParameter" dimension="2" estimate="false">true</indicators>
    <typeTrait idref="typeTraitSet.t:SimDataset"/>
  </dynamics>
  <structuredTreeIntervals id="StructuredTreeIntervals.t:SimDataset" spec="mascot.distribution.StructuredTreeIntervals" tree="@Tree.t:SimDataset"/>
</distribution>
```

`dimension="2"` must match the number of demes. For an *n*-deme model, add one `neDynamics` block per extra deme, extend `migrationRatesSkyline` to `dimension="n*(n-1)"`.

### 2.4 Data stream likelihoods

Each likelihood plugs the deme's spline into a stream-specific submodel:

```xml
<distribution id="caseCountLikelihood.Deme1.t:SimDataset" spec="mascotdatastreams.distribution.CaseCountLikelihood">
  <prevalenceSpline idref="splinePrev.Deme1.t:SimDataset"/>
  <caseCounts idref="caseCounts.Deme1:SimDataset"/>
  <caseTimes idref="caseTimes.Deme1:SimDataset"/>
  <scaling idref="caseCounts.scaling.Deme1:SimDataset"/>
  <distribution id="gammaPoisson.Deme1.t:SimDataset" spec="mascotdatastreams.distribution.GammaPoisson"
                dispersion="@caseCounts.dispersion:SimDataset"/>
</distribution>
```

`caseCounts.scaling.Deme1` is an estimated ascertainment-rate-like scaling factor (how many reported cases per unit prevalence); `caseCounts.dispersion` is the shared overdispersion parameter of the Gamma-Poisson observation model.

```xml
<distribution id="seroprevalenceLikelihood.Deme1.t:SimDataset" spec="mascotdatastreams.distribution.SeroprevalenceLikelihood">
  <prevalenceSpline idref="splinePrev.Deme1.t:SimDataset"/>
  <seroPeopleTested idref="seroTestedCounts.Deme1:SimDataset"/>
  <seroPeopleSeropositive idref="seroWithAntibodiesCounts.Deme1:SimDataset"/>
  <seroTimes idref="seroTestedTimes.Deme1:SimDataset"/>
  <populationSize idref="PopSize.Deme1:SimDataset"/>
  <scaling idref="seroprevalence.scaling.Deme1:SimDataset"/>
  <distribution id="binomial.seroprevalenceLikelihood.Deme1.t:Spec" spec="mascotdatastreams.distribution.Binomial"/>
</distribution>
```

```xml
<distribution id="wastewaterLikelihood.Deme1.t:SimDataset" spec="mascotdatastreams.distribution.WastewaterLikelihood">
  <prevalenceSpline idref="splinePrev.Deme1.t:SimDataset"/>
  <concentrations idref="wastewaterConcentration.Deme1:SimDataset"/>
  <concentrationTimes idref="wastewaterConcentrationTimes.Deme1:SimDataset"/>
  <populationSize idref="PopSize.Deme1:SimDataset"/>
  <scaling idref="wastewater.scaling.Deme1:SimDataset"/>
  <distribution id="logNormal.wastewaterLikelihood.Deme1.t:SimDataset" spec="mascotdatastreams.distribution.LogNormal"
                sd="@wastewater.sigma:SimDataset"/>
</distribution>
```

To drop a data stream for a deme (e.g. no wastewater surveillance), delete its `<distribution>` block, its state-node parameters, its prior and operators, and its logger entries. Nothing else needs to change — the tree likelihood only depends on the spline, not on which observation likelihoods are attached to it.

### 2.5 Regularising priors

One extra "prior" enforces internal consistency rather than expressing genuine prior belief:

```xml
<distribution id="regularizeTransmissionRate.Deme1.t:SimDataset" spec="mascotdatastreams.util.TransmissionSmallerThan">
  <spline idref="splinePrev.Deme1.t:SimDataset"/>
</distribution>
```

These reject states where the implied transmission rate becomes non-positive at any grid point of the prevalence trajectory. One per deme.

Standard `LogNormal` priors are placed on the scaling/dispersion/migration parameters, e.g.:

```xml
<prior id="migrationRatesSkylinePrior.t:SimDataset" name="distribution" x="@migrationRatesSkyline.t:SimDataset">
  <LogNormal id="LogNormal.3" name="distr">
    <parameter spec="parameter.RealParameter" estimate="false" name="M">0.5</parameter>
    <parameter spec="parameter.RealParameter" estimate="false" name="S">0.5</parameter>
  </LogNormal>
</prior>
```

Adjust `M`/`S` to reflect your own prior beliefs about migration rate, ascertainment scaling, wastewater scaling, dispersion, etc.

### 2.6 Operators

The AVMN (adaptable variance multivariate normal) operator jointly proposes the "global" scalar parameters (migration rates, scalings, dispersion), while dedicated operators handle the spline vectors and an up/down operator co-scales a deme's log-prevalence against its own case-count and wastewater scalings (since prevalence and ascertainment scaling are only identifiable up to their product):

```xml
<operator id="UpDownPrevScaling.Deme1:SimDataset" spec="mascotdatastreams.operators.UpDownLogScaleOperator"
          scaleFactor="0.2" weight="1.5" optimise="true">
  <upLogParameter idref="SkylinePrev.Deme1.t:SimDataset"/>
  <downRealParameter idref="caseCounts.scaling.Deme1:SimDataset"/>
  <downRealParameter idref="wastewater.scaling.Deme1:SimDataset"/>
</operator>
```

If a scaling parameter doesn't exist for your setup (e.g. you dropped wastewater), remove its `<downRealParameter>` entry from this operator.

Additional operators for specific parameters might be necessary for appropriate mixing of the MCMC chain. This might need some trial and error. The [example](examples/2deme_fixedtree.xml) is a good starting point though. 

If the tree is not fixed, you'll additionally need the usual tree operators (`SubtreeSlide`, `Exchange`, `WilsonBalding`, etc.) and a tree-height/root-time scaler — none are present here since `Tree.t:SimDataset` has no operators acting on it in this example.

### 2.7 Loggers

- `tracelog` — the main `.log` file: posterior/prior, per-deme likelihoods, scaling/dispersion parameters, the tree-likelihood component, and the log-prevalence skyline vectors.
- `screenlog` — posterior/prior progress printed to the console.
- `treelog` — the (untyped) sampled trees.
- `typedTreelogger` (`mascotdatastreams.logger.StructuredTreeLogger`) — trees annotated with deme state along branches, for visualising migration history.
- `typedEventsTreelogger` — trees annotated with individual migration events (`mascot.logger.MigrationCountLogger` via a `MappedMascot`).
- `NeDynamicsLogger.Deme1`/`Deme2` — per-deme effective-population-size trajectory over the grid.
- `cumulativeIncidenceLogger.Deme1`/`Deme2` (`mascotdatastreams.logger.CumulativeIncidenceLogger`) — cumulative incidence implied by the fitted prevalence spline, comparable against seroprevalence data.

## 3. Running the analysis

Install the package (see the main [README](README.md) for build instructions), then run:

```bash
beast examples/2deme_fixedtree.xml
```

Outputs (`.log`, `.trees` files) are written next to the XML, named from its `$(filebase)` (`2deme_fixedtree`). Inspect the `.log` file in [Tracer](https://github.com/beast-dev/tracer) to check ESS and convergence of the posterior, per-deme scaling/dispersion parameters, and the log-prevalence skyline; use the `NeDynamics.*.log` and `cumulativeIncidence.*.log` files to plot fitted trajectories against your raw case-count/seroprevalence/wastewater data as a posterior-predictive check.

## 4. Adapting the example to your own data

A minimal checklist for turning this example into your own analysis:

1. Replace the alignment (or keep a placeholder if using a fixed tree) and the `newick` fixed-tree string with your taxa, or remove the fixed-tree `<init>` and add tree operators if the tree should be estimated.
2. Update `typeTraitSet` with your taxa's true deme assignments.
3. For each deme, replace the `caseCounts`/`caseTimes`, `seroTestedCounts`/`seroWithAntibodiesCounts`/`seroTestedTimes`, and `wastewaterConcentration`/`wastewaterConcentrationTimes` vectors with your own observations (dropping any stream you don't have, per §2.2/§2.4).
4. Set `PopSize.DemeX` to each deme's actual population size.
5. Adjust `SkygrowthRateShifts`/`SplineGridRateShifts` to span your study period at appropriate resolution (coarse knots vs. fine evaluation grid), and resize `SkylinePrev.DemeX` (`dimension`) to match the number of coarse knots.
6. Review/adjust priors (`M`/`S` of each `LogNormal`) to reflect realistic ranges for your migration rates, ascertainment scaling, wastewater scaling, and dispersion.
7. If you add a third (or more) deme, extend `dimension` on `StructuredSkylinePrevalence` and `migrationRatesSkyline`, add a `neDynamics` block per deme, and add matching data-stream/likelihood/logger blocks.
8. Sanity-check `chainLength`, `storeEvery`, and `logEvery` for your dataset size and run time budget.

For background on the structured-coalescent MASCOT model itself (deme structure, migration rates, effective population sizes), see the [MASCOT repository](https://github.com/CompEvol/Mascot) and the [taming-the-BEAST MASCOT tutorial](https://github.com/taming-the-BEAST/Mascot-Tutorial) — MASCOT-DS reuses MASCOT's coalescent machinery and adds the datastream-linked prevalence dynamics described above.
