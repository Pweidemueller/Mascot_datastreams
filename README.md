# MASCOT Datastreams

**MASCOT-DS** (MASCOT Datastreams) is a [BEAST2](https://www.beast2.org/) package that extends [MASCOT](https://github.com/CompEvol/Mascot) to jointly infer structured phylodynamic histories from genomic data (pathogen sequences) *together with* epidemiological data streams such as case counts, seroprevalence surveys, and wastewater measurements.

Like MASCOT, MASCOT-DS models how a pathogen moves between demes (e.g. locations or subpopulations) along a phylogeny, using the structured coalescent approximation to relate migration rates and effective population sizes to tree topology and branch lengths. MASCOT-DS in contrast is built around infering prevalence and forward migration rates, improving interpretability of transmission dynamics estimates. MASCOT-DS parameterises each deme through a "prevalence spline" — a trajectory of underlying infections over time. It then calculates the effective population size explicitly from the prevalence and derived transmission rate to evaluate the structured coalescent likelihood.  This same prevalence trajectory is plugged into likelihood terms for other data streams:

- **Case counts** — a negative binomial distribution serving as a discrete, non-negative likelihood function (`CaseCountLikelihood` + `GammaPoisson`) linking reported counts to the prevalence spline through a scaling factor and overdispersion parameter.
- **Seroprevalence** — a binomial distribution serving as the seroprevalence likelihood (`SeroprevalenceLikelihood` + `Binomial`) linking survey testing results to cumulative incidence implied by the spline.
- **Wastewater concentrations** — a log-normal distribution serving as a continuous, strictly positive likelihood function (`WastewaterLikelihood` + `LogNormal`) linking normalised wastewater pathogen DNA concentrations to the prevalence spline through a scaling factor and standard deviation to model the noisiness of the data.

Because all these likelihoods share the same underlying prevalence dynamics as the tree likelihood, MASCOT-DS lets the genealogy and the epidemiological data streams jointly inform a single, consistent epidemic trajectory per deme, rather than analysing each data source separately.

MASCOT-DS does not yet have BEAUti support — analyses are set up by editing an example XML file directly (see [tutorial.md](tutorial.md)).

## Installation
We are still working on a Beauti integration. For now, you can use MASCOT-DS by simply placing the provided zipped jar file (under Releases) in the package folder of your BEAST2 installation (for MACOS this is often: `/Users/<username>/Library/Application Support/BEAST/2.7/`) and unzip it there (make sure it's in a dedicated folder called `MASCOTDS`). Make sure to also install MASCOT either via the instructions [here](https://github.com/CompEvol/Mascot) or via Beauti.

You can also build MASCOT-DS from source using ant:
```
ant build
```
## How to run
We are still working on a Beauti integration. For now, you can use the xml file in the `examples` folder as a reference on how to add different data types, associated priors and operators. Running the `2deme_fixedtree.xml` with beast after you've installed the Mascot-DS package the analysis should progress with approximately 40min per 1M samples (tested in MacBook Air M3 using one thread). 1 million samples should be sufficient to investigate convergence in tracer.



## Citation
MASCOT-DS improves transmission dynamics inference by integrating multiple epidemiological data streams with phylodynamic inference
Paula H. Weidemüller, Luis R. Esquivel Gomez, Isabel Rodriguez-Barraquer, Nicola F. Müller
medRxiv 2026.08.21.26361056; doi: https://doi.org/10.64898/2026.08.21.26361056

## License
The java source code is licensed under the [GNU General Public License v3.0](LICENSE)
