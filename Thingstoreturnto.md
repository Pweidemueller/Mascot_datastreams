# Things to return to

Running list of issues, workarounds, and follow-ups discovered along the way.
Each entry: what, where, why it's deferred, what to do later.

---

## Mascot upstream `Euler2ndOrder.initWithIndicators` does not allocate `sumDotStates`

**Where:** `mascot.ode.Euler2ndOrder.initWithIndicators(double[], int[], double[], int)` (upstream Mascot library, `deps/Mascot/Mascot.v3.0.7.src.jar`).

**What:** The `init(double[], double[], int)` overload allocates three temporary arrays:
```java
sumStates = new double[states];
tCR = new double[states];
sumDotStates = new double[states];
```
The sibling `initWithIndicators(...)` allocates only the first two — `sumDotStates` is left null. `Euler2ndOrder2.computeSecondDerivate` then calls `clearArray(sumDotStates, states)` and NPEs.

**How we currently fix it:** `src/mascotdatastreams/ode/Euler2ndOrder.java` is a vendored fork of the upstream class with the missing allocation added. The 2/3/4-state specialisations (`Euler2ndOrder2/3/4.java`) are also vendored so they inherit from our fixed parent. `MascotLogPflag.initAndValidate()` swaps in the vendored classes when `useMaxInterval=true`, and `doEulerAtTime` now uses `initWithIndicators(...)` to keep the indicator-based sparsity optimisation.

**To do later:** File a one-line PR upstream to add `sumDotStates = new double[states];` to `initWithIndicators`. Once a fixed Mascot release ships and is referenced via `deps/Mascot`, delete `src/mascotdatastreams/ode/Euler2ndOrder*.java` and revert the swap in `MascotLogPflag.initAndValidate` to use `mascot.ode.Euler2ndOrderN`.

---
## Why is `java` instead of `allnative` necessary for implementation? 

```
public Input<MascotImplementation> implementationInput = new Input<>("implementation", "implementation, one of " + MascotImplementation.values().toString(),
        MascotImplementation.allnative, MascotImplementation.values());
```

## RandomTree-initialised XMLs need a pinned RNG seed for reproducible benchmarks

**Where:** `examples/SARSCoV2_Epsilon_BayArea_results_1000seq_datastreams_noclip.xml` line 1076 — `<init id="RandomTree.t:SimDataset" spec="RandomTree" ...>` (no seed specified). The small `examples/6_2_simulation_datastreams.xml` is fine because it uses a fixed Newick string via `TreeParser`.

**What:** When `XMLParser.parseFile()` instantiates a `RandomTree`, it uses BEAST's process-wide `Randomizer` whose seed defaults to the system clock. Two JVM invocations on the same XML produced different `logP` and different `doEuler` counts (e.g. -1714.5 / 1396 vs -1713.0 / 1442).

**How we sidestep:** `Tier1Runner`, `Tier2Runner`, and `RepeatabilityCheck` all call `Randomizer.setSeed(42)` before each `parser.parseFile(...)`. This pins the random tree across the whole sweep so configs are comparable, and pins it across separate JVM invocations so re-runs are byte-identical.

**To do later:**
- If we ever do Tier 3 (real MCMC equivalence runs) on the SARS XML, the operators move the tree, so the random init only matters for the starting point. Either accept that or re-roll the operator seed alongside the init seed for a fully reproducible chain.
- Consider exposing the seed via a CLI flag in the harnesses (`-Dtier2.seed=...`) — already done for Tier2Runner; mirror in Tier1Runner if we expand its sweep.

---

## Class-loading consumes Randomizer values on the first XML parse in a JVM

**Where:** Inside BEAST's `XMLParser.parseFile(...)` — concretely whichever lazily-loaded class first uses `Randomizer` (operators, RandomTree, etc.). Affects any XML that drives a `RandomTree` `<init>`.

**What:** With `Randomizer.setSeed(42)` called immediately before each `parser.parseFile(...)`, the **first** XML in a fresh JVM consumed N random values during class loading and then ran `RandomTree.initAndValidate()` with the leftover RNG state. The **second** XML (classes already loaded) consumed 0 random values during class loading and got the full fresh seeded sequence — different tree, different `logP`, different `n_doEuler`.

We caught this comparing the SARS rerun (new first → old second) to the original main sweep (max=5e-3 processed mid-sweep): bitwise-identical XML gave `logP=-1699.45 / 1618 doEuler` in the rerun but `-1721.89 / 1432 doEuler` in the main sweep. Confirmed by reordering the rerun pair: whichever XML went first received the "first parse" tree.

**How we currently fix it:** `Tier2Runner.main` does an explicit class-loading primer parse on `args[1]` before entering the timing loop. After the primer, every real parse with `Randomizer.setSeed(42)` sees the same post-class-loading RNG state. Verified: both `(new, old)` and `(old, new)` orderings now produce identical `logP` and `n_doEuler` for each XML.

**Caveat:** `Tier1Runner` has the same theoretical exposure but the small-XML examples use `TreeParser` (fixed Newick), so it never bit us there. Worth adding the same primer if Tier 1 ever expands to RandomTree XMLs.

**To do later:**
- Mirror the primer in `Tier1Runner` defensively.
- Better long-term fix: replace `<init spec="RandomTree" ...>` in the SARS XML with a `<init spec="TreeParser" newick="..."/>` block that bakes in the seed-42 starting tree. Eliminates the random init dependency entirely and makes the whole pipeline reproducible without harness gymnastics. Feasible — just means dumping the seed-42 Newick once and editing the XML.

---

## SARS XML at user's published settings appears under-resolved

**Where:** `examples/SARSCoV2_Epsilon_BayArea_results_1000seq_datastreams_noclip.xml`, the `<rateShifts id="mascotshifts" ...>` element (line 1211, step 0.005, range 0..1.5).

**What:** In Tier 2 we found that:
- The OLD path's outer iteration grid is `mascotshifts` (step 0.005, ~300 points), separate from `SplineGridRateShifts` (the spline eval grid, ~1000 pts).
- At the user's grid setting (mascotshifts step 0.005), the OLD path returns `logP ≈ -1730`.
- The NEW path at `maxInterval=1e-4` (~10 000 sub-events) returns `logP ≈ -1750` — a **20 logP unit gap**.

That's a much bigger discrepancy than we saw on the small (6_2) XML, where the user's grid was already well-tuned. It strongly suggests that the user's published `mascotshifts` grid for the 1000-seq SARS analysis is too coarse and the existing posterior has a systematic bias of order tens of logP. Worth flagging in the manuscript and rerunning the published analyses with a finer grid (or with the new max-interval path).

**To do later:**
1. Confirm by running OLD with `mascotshifts` densified 5×, 10×, 20×; check whether logP plateaus.
2. If confirmed, audit which posteriors in `manuscript/` were generated with the under-resolved grid.
3. Decide: re-run those analyses with finer grid, or with `useMaxInterval=true maxInterval=1e-4` once the new path is validated end-to-end (Tier 3).

---
