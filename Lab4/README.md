# Lab 4: Euclidean Distance and Parallel Computing

## Overview

In this lab, we calculate the same distance
matrix using nested loops, vectorization, and `foreach`, then compare their
running times.

The data contain 3,109 counties. For counties i and j, calculate:

```text
D[i, j] = sqrt((latitude[i] - latitude[j])^2 +
               (longitude[i] - longitude[j])^2)
```

The output is an N by N matrix. These are Euclidean distances in coordinate
degrees.

## Learning objectives

1. Calculate pairwise distances with nested loops.
2. Replace the inner loop with vectorized arithmetic.
3. Use `foreach` for sequential and parallel calculations.
4. Compare runtime, speedup, and parallel efficiency.

## Project structure

```text
Lab4/
|-- data-raw/
|   |-- counties.csv
|-- R/
|   |-- lab4.r
|   |-- run_lab4.r
|-- results/
|   |-- timing.csv
|-- tests/
|   |-- lab4_test.r
|-- README.md
```

## Requirements

Install the packages once:

```r
install.packages(c("foreach", "doParallel", "microbenchmark"))
```

## In-class exercise

```r
counties <- read.csv(file.path("data-raw", "counties.csv"))
mat <- as.matrix(counties[seq_len(n_counties), c("latitude", "longitude")])
```

### Task 1: Nested loops

Write `distance_loop(mat)`. Preallocate an N by N matrix, use one loop over
county i and another over county j, and fill each entry with its Euclidean
distance.

### Task 2: Vectorization

Write `distance_vectorized(mat)`. Keep the loop over i, but calculate distances
from county i to all counties in one vector expression. Store that vector in
`distances[i, ]`.

Check that the two methods return the same distances before comparing time.

### Task 3: Foreach

Write `distance_foreach(mat, parallel = FALSE)`. Each iteration calculates one
row of the distance matrix. Combine the returned rows with `.combine = "rbind"`.

Use `%do%` for sequential execution. Then register a two-worker cluster with
`doParallel`, switch to `%dopar%`, and calculate the same matrix in parallel.
Stop the cluster when finished.

Reference implementations are in `R/lab4.r`. The cluster setup and timing
comparison are in `R/run_lab4.r`.

### Compare the methods

Time all four versions with `microbenchmark(..., times = 5)`:

- nested loops;
- the vectorized loop;
- sequential `foreach`;
- parallel `foreach` with two workers.

Use the median elapsed time. The supplied runner saves the comparison to
`results/timing.csv`. It starts the cluster before timing, so parallel timings
include calculation, communication, and combining results, but exclude cluster
startup and shutdown.

Calculate:

```text
Vectorization gain = nested-loop time / vectorized time
Parallel speedup   = sequential foreach time / parallel foreach time
Parallel efficiency = parallel speedup / number of workers
```

Discuss together:

1. How much faster is vectorization than the nested loop?
2. Does vectorization change the O(N^2) amount of work?
3. Is sequential `foreach` faster than the vectorized `for` loop?
4. Does parallel `foreach` beat both sequential `foreach` and vectorization?
   How do communication and combining results affect the comparison?


## How to run

Set `Lab4` as the working directory. From the course repository root:

```r
setwd("Lab4")
```

Run the correctness checks:

```r
source(file.path("tests", "lab4_test.r"))
```

Run the timing comparison:

```r
source(file.path("R", "run_lab4.r"))
```

## Author

Xiangyu Song

## Date

2026-09-15
