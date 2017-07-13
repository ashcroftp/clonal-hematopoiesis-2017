# Clonal-hematopoiesis-2017

Supplementary software and data files for manuscript "Clonal Dominance and Transplantation Dynamics in Hematopoietic Stem Cell Compartments", by Peter Ashcroft, Markus G. Manz, & Sebastian Bonhoeffer (2017).




## Stochastic simulations
All simulations can be performed using the [`gillespie-simulation.cpp`](https://github.com/ashcroftp/clonal-hematopoiesis-2017/blob/master/gillespie-simulation.cpp) program, which is written in c++.
It uses the `<random>` header, which requires c++11 or later.
Compilation may require:
```
g++ -std=c++11 gillespie-simulation.cpp
```

Execution of the compiled object `a.out` follows:
```
./a.out paramIndex simindex nRuns
```
where `paramIndex` is used to sweep over parameter values, `simIndex` is used for output file labeling, and `nRuns` is the number of realisations of the stochastic process to run.
This program is designed to be run on a cluster.
Say we want to run 10,000 realisations of a process for a given parameter set (`paramIndex=0`, say).
We can spread this as 1,000 realisations run on 10 nodes (`nRuns=1000` and `simIndex` takes values 0-9). 

The code can be manipulated in many ways:

* Changing parameters to study clonal hematopoiesis in different organisms. The variable `paramindex` is an input to the executable, and can be used to perform parameter sweeps.

* Sampling can be achieved within the actual Gillespie algorithm. Currently the code is set up to record the time when different levels of chimerism are reached.

* The algorithm can be run until a specified number of cells is reached (useful for reconstitution of preconditioned hosts, or investigating different levels of clonality) or it can be run until one of the cell types reaches extinction.

* You can output the data from sampling, or the end state of the system.




## Mathematica notebook
The notebook [`clonal-hem-notebook.nb`](https://github.com/ashcroftp/clonal-hematopoiesis-2017/blob/master/clonal-hem-notebook.nb) is originally written in Mathematica 10.3.0.0.
It contains multiple sections with analytical steps, data analysis, and graph production, as detailed below:

1. Evaluate steady-state parameters, residence times of cells in each anatomical compartment, and fluxes of cells between compartments.

2. Analytical calculations of fixation probability and mean fixation times. This is handled separately for neutral and advantageous clones.

3. Analysis of clonal dominance in mice.

4. Analysis of clonal dominance in humans.

5. What is the probability that preconditioned mice are reconstituted by a dose of donor cells.

6. How successful is transplantation into a non-preconditioned mouse.