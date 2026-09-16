# Convert Simulation Results to Array Format

Converts a list of simulation results (data.frames) to a 3-dimensional
array format for easier analysis and manipulation. This is a utility
function for post-processing simulation results from
[`run_parallel_simulation`](https://t-momozaki.github.io/DRBayes/reference/run_parallel_simulation.md).

## Usage

``` r
convert_to_array(simulation_results)
```

## Arguments

- simulation_results:

  List of data.frames, typically the output from
  [`run_parallel_simulation`](https://t-momozaki.github.io/DRBayes/reference/run_parallel_simulation.md).

## Value

A 3-dimensional array with dimensions:

- Dimension 1:

  Observations (rows of each data.frame)

- Dimension 2:

  Variables (columns of each data.frame)

- Dimension 3:

  Simulations (list elements)

The array includes appropriate dimension names for easy indexing.

## Details

This function is useful when you need to:

- Perform array-based operations across simulations

- Calculate statistics across the simulation dimension

- Export results to other software that expects array format

- Use vectorized operations for faster computation

**Memory warning:** The resulting array can be very large for many
simulations or high-dimensional datasets. Consider the memory
requirements before conversion.

## See also

[`run_parallel_simulation`](https://t-momozaki.github.io/DRBayes/reference/run_parallel_simulation.md),
[`generate_dataset`](https://t-momozaki.github.io/DRBayes/reference/generate_dataset.md)

## Examples

``` r
# Generate simulation results with only primary confounders
results <- run_parallel_simulation(num_simulations = 5, nn = 100, pp = 0,
                                   seed = 1)

# Convert to array format
result_array <- convert_to_array(results)

# Dimensions: observations, variables, simulations
dim(result_array)  # 100 x 6 x 5 for pp = 0
#> [1] 100   6   5
dimnames(result_array)[[2]]
#> [1] "AA"  "YY"  "WW1" "WW2" "WW3" "WW4"

# Access a specific simulation
sim_1 <- result_array[, , 1]  # First simulation as a matrix

# The unadjusted difference in means for each simulation. The unit level
# treatment effect is 110, so the gap is the confounding.
ate_by_sim <- apply(result_array, 3, function(sim) {
  treated <- sim[sim[, "AA"] == 1, "YY"]
  control <- sim[sim[, "AA"] == 0, "YY"]
  mean(treated) - mean(control)
})
round(ate_by_sim, 1)
#> sim_1 sim_2 sim_3 sim_4 sim_5 
#> 126.8 129.9 125.8 117.6 119.7 
```
