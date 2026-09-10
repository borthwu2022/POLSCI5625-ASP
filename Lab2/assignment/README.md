# Lab 2 Assignment: Thresholded Pairwise Distances

## Problem

Suppose `x` and `y` contain one-dimensional policy positions. For every pair
`(x[i], y[j])`, calculate the Euclidean distance

```r
abs(x[i] - y[j])
```

If a distance is strictly less than `0.5`, record it as `0`. Preserve every
distance greater than or equal to `0.5`.

Use the following reproducible vectors:

```r
set.seed(123)
x <- runif(3000)
y <- runif(5000)
threshold <- 0.5
```

The resulting matrix has 3,000 rows, 5,000 columns, and 15,000,000 entries.

## Requirements

- R 4.0 or later.
- The recommended package `Matrix`.

Install `Matrix` once if it is not available:

```r
install.packages("Matrix")
```

## Project structure

```text
assignment/
|-- data-raw/
|   |-- lab2_vectors.rds
|-- R/
|   |-- assignment_helpers.r
|   |-- gendata.r
|   |-- lab2_assignment.r
|   |-- lab2_assignment_solution.r
|   |-- run_assignment.r
|-- results/
|   |-- large_benchmark.csv
|-- tests/
|   |-- lab2_assignment_test.r
|-- README.md
```



## Task 3: Nested loops

Write `thresholded_distance_loop(x, y, threshold)`.

- Preallocate a dense `length(x)` by `length(y)` matrix.
- Use two nested loops.
- Calculate `abs(x[i] - y[j])` once for each pair.
- Store `0` when the distance is strictly less than `threshold`.
- Otherwise, store the original distance.
- Do not grow the result inside the loops.

## Task 4: Vectorized dense matrix

Write `thresholded_distance_vectorized(x, y, threshold)` without nested loops.
The primary version may use `outer()`.

The instructor solution also provides
`thresholded_distance_vectorized_handcoded()` immediately below the primary
version. This alternative does not use `outer()`. Its strategy is:

1. use `matrix()` to repeat `x` across `length(y)` columns;
2. use `sweep(..., MARGIN = 2, STATS = y, FUN = "-")` to subtract `y[j]`
   from column `j`;
3. take the absolute value; and
4. replace distances strictly below the threshold with zero.

The function must return a regular dense R matrix.

## Task 5: Vectorized sparse matrix

Write `thresholded_distance_sparse(x, y, threshold)`.

- Reuse the vectorized calculation from Task 4.
- Convert the thresholded matrix with
  `Matrix::Matrix(distances, sparse = TRUE)`.
- Return a sparse matrix from the `Matrix` package.

The vectors are generated on `[0, 1]`, so the threshold creates many zeros.
This is important: sparse storage is useful only when enough entries are zero.

## Correctness before benchmarking

Check the implementations on a small example before measuring speed:

```r
x_test <- c(0, 1)
y_test <- c(0, 0.5, 2)

expected <- matrix(
  c(0, 0.5, 2,
    1, 0.5, 1),
  nrow = 2,
  byrow = TRUE
)
```

Verify that:

1. Tasks 3 and 4 equal `expected`;
2. `as.matrix()` of Task 5 equals `expected`;
3. the output dimensions are correct;
4. values exactly equal to `0.5` are preserved; and
5. the dense and sparse versions contain the same number of nonzero entries.

Use a numerical tolerance when appropriate. Run the supplied tests with the
assignment folder as the working directory:

```r
source(file.path("tests", "lab2_assignment_test.r"))
```

## Performance comparison

After the correctness checks pass, benchmark all three functions using the
same large vectors and computer. Report at least:

```text
method | elapsed seconds | stored result bytes | nonzero entries
```

Use `system.time()` for elapsed time and `object.size()` for the stored size.
For the sparse result, use `Matrix::nnzero()` to count nonzero entries.

Calculate the time gain from Task 3 to Task 4:

```r
time_gain <- time_loop / time_vectorized
time_reduction_percent <- 100 * (1 - time_vectorized / time_loop)
```

Calculate the stored-result memory reduction from Task 4 to Task 5:

```r
memory_reduction_percent <-
  100 * (1 - as.numeric(object.size(d_sparse)) /
           as.numeric(object.size(d_dense)))
```


## Deliverables

Submit:

- the completed `R/lab2_assignment.r` file;
- a timing and memory table;
- the correctness checks; and
- short answers to the six questions.


## Files

- `R/assignment_helpers.r` contains shared input checks and timing support.
- `R/gendata.r` creates the reproducible assignment vectors.
- `R/lab2_assignment.r` is the student starter file.
- `R/lab2_assignment_solution.r` is the instructor solution, including the
  additional hand-coded vectorized alternative.
- `tests/lab2_assignment_test.r` tests the Tasks 3-5 implementations.
- `R/run_assignment.r` runs the full assignment benchmark.
- `data-raw/` and `results/` contain only assignment inputs and outputs.

## Running the assignment materials

Start with `Lab2/assignment` as the working directory, then run:

```r
source(file.path("R", "gendata.r"))
source(file.path("tests", "lab2_assignment_test.r"))
source(file.path("R", "run_assignment.r"))
```


