# smurfs: **S**ampling **M**ultivariate **U**niform **R**ecords with **F**ixed **S**um

`smurfs` offers R implementations of algorithms for the generation of unbiased task utilization vectors. A task utilization vector defines the utilization of a set of $n$ tasks with a fixed total utilization. The total utilization, as well as the range possible utilization values for each task, may vary case by case.

From a sampling perspective, an unbiased task utilization vector amounts to a uniform sample over the $n$-dimensional space $(X_1, ..., X_N)$ where:

  - $x_i + ... + x_n = U$
  
  - $x_i \in [a_i, b_i], \forall i$
  
Without any intervention, the $x_i$ will be bounded by $[0,U]$. However, a proper general-purpose solution will allow each value to be bounded uniquely.

This has a lot of applications, even outside of the task utilization space! So, this package aims to offer R users access to the convenient algorithms that do such sampling. It features the suite of unbiased task utilization generators discussed in "Generating Utilization Vectors for the Systematic
Evaluation of Schedulability Tests" (Griffin, Bate and Davis, 2020). This includes four previously-defined case-specific solutions and the authors' proposed general-purpose solution, called the Dirichlet-Rescale (DRS) algorithm.

Of the methods in this package, DRS stands out as the number-one option due to its flexibility; in fact, it may be used as a general-purpose replacement for them. The other methods are included for primarily for completeness, but may be more efficient in their specific use cases.

## Installation

You can install the development version of smurfs like so:

``` r
# devtools::install_github("aaron-nlt/smurfs")
```

