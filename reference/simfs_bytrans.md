# Reformat simulated multi-state data with one row per simulated transition

Reformat simulated multi-state data with one row per simulated
transition

## Usage

``` r
simfs_bytrans(simfs)
```

## Arguments

- simfs:

  Output from
  [`sim.fmsm`](http://chjackson.github.io/flexsurv-dev/reference/sim.fmsm.md)
  representing simulated histories from a multi-state model.

## Value

Data frame with four columns giving transition start state, transition
end state, transition name and the time taken by the transition.
