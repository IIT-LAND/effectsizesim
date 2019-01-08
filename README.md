# effectsizesim

Code to reproduce simulations in Lombardo, Lai, & Baron-Cohen (in press). Big
data approaches to decomposing heterogeneity across the autism spectrum.
```Molecular Psychiatry```. doi:10.1038/s41380-018-0321-0.
https://www.nature.com/articles/s41380-018-0321-0

Simulations on effect size variability and inflation for simple case-control
comparisons. Simulations will generate population-level data, and then over
repeated simulated experiments, will randomly sample from those populations.
Samples will be of different sizes. Simulations will also run across many
population effect sizes. Finally, there is a simulation for calculating bias in
sample prevalence of nested subtypes within one population. This is meant to
illustrate how specific strata within the population can be enriched in studies
that have small sample size.

Code implements the simulations in ```python``` or ```matlab```.

```effect_size_simulation.m``` or ```effect_size_simulation.py``` will do the
main simulation to reproduce Figure 3A-E.

```effect_size_inflation_sim.m``` reproduces Figure 3F.

```heterogeneity_simulation.m``` or ```heterogeneity_simulation.py``` will do
the simulation to reproduce Figure 4.
