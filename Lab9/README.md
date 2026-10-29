# Lab 9: Getting Started with Rcpp

## Overview

Write small C++ functions and call them from R with `Rcpp::cppFunction()`.
Practice scalar arithmetic, numeric vectors, zero-based indexing, and logical
comparisons using vote totals and voting eligibility. Then put the three
exercise functions in separate C++ files and check them from R with `sourceCpp()`.

This lab follows Lecture 9: compilation and `cppFunction()` (slides 3-11),
scalar types and operators (slides 12-15), vectors and indexing (slides 16-18),
and an optional loop and `sourceCpp()` example (slides 19-24).

## Learning objectives

1. Compile a C++ function and call it from R using `cppFunction()` and `sourceCpp()`.
2. Distinguish scalar `double`, vector `NumericVector`, and logical `bool`.
3. Declare local variables and use `return` and semicolons.
4. Translate between R's one-based and C++'s zero-based indexing.
5. Write and check `voteDifference()`, `secondValue()`, and `isEligible()`.
6. Explain why compiled code does not guarantee a faster analysis.

## Project structure

```text
Lab9/
|-- R/
|   |-- lab9_exercises.r       # student examples and exercise prompts
|   |-- lab9.r                 # reference implementations with cppFunction()
|-- src/
|   |-- voteDifference.cpp    # Exercise 1
|   |-- secondValue.cpp       # Exercise 2
|   |-- isEligible.cpp        # Exercise 3
|   |-- sumSquares.cpp        # optional lecture example
|-- tests/
|   |-- lab9_test.r           # compile and check the three exercise files
|-- README.md
```

No external data are needed. Small inputs are defined directly in the scripts.
As in earlier labs, paths are relative to the lab folder.

## Requirements and setup

You need R, the **Rcpp** package, and a working C++ compiler. Install Rcpp once:

```r
install.packages("Rcpp")
```

A compiler is a separate requirement: Windows uses Rtools matched to its R
version; macOS uses Xcode Command Line Tools (`xcode-select --install` in a
terminal); Linux needs the C++ toolchain and development files for R.
If compilation fails, read the first compiler error and check that toolchain.
The optional timing exercise also uses `microbenchmark`.

From the course repository root in R:

```r
setwd("Lab9")
library(Rcpp)
```

`cppFunction()` takes C++ source code inside an R string, creates a wrapper,
compiles and loads it, and makes the function available in the current R
session. Rcpp does not translate arbitrary R code into C++.

For the core exercises, use finite, nonmissing numeric inputs. Scalar `double`
arguments take one number; they do not automatically apply a function to a
whole vector. Vote totals are nonnegative and `totalVotes` must be positive.
A scalar C++ `bool` has only `true` and `false`; it cannot represent R's `NA`.
These short teaching functions do not implement a general missing-data policy.

## In-class exercise

Open `R/lab9_exercises.r` and work through it from top to bottom. It contains
all five supplied examples and commented templates for the three exercises.
Each example function has a short R implementation example directly below
its definition, with input values, a function call, and the expected output.
Fill each template and uncomment its `cppFunction()` block. Write your
functions before consulting `R/lab9.r`, which contains the complete solutions.

### Task 1: Scalars, arithmetic, and return values

Start with two numbers:

```r
cppFunction('
  double addNumbers(double x, double y) {
    return x + y;
  }
')

addNumbers(2, 3)       # 5
addNumbers(2.5, 1.25)  # 3.75
```

The first `double` declares the return type. Each argument also has a declared
type. Curly braces enclose the function body, `return` supplies the result,
and a semicolon ends a statement. C++ does not automatically return the last
expression as R does.

The next example introduces a local variable:

```r
cppFunction('
  double votePercentage(double votes, double totalVotes) {

    double proportion = votes / totalVotes;

    return proportion * 100.0;
  }
')

votePercentage(450, 600)  # 75
```

Because the arguments are doubles, this uses decimal division. In contrast,
C++ integer division `5 / 2` gives `2`; `5.0 / 2.0` gives `2.5`.

**Exercise 1.** Write `voteDifference(double votesA, double votesB)`, returning
a `double`. Define the difference as **the first total minus the second**,
`votesA - votesB`. This is a signed vote margin, measured in votes.

| R call | Expected result |
|---|---:|
| `voteDifference(620, 580)` | 40 |
| `voteDifference(580, 620)` | -40 |
| `voteDifference(600, 600)` | 0 |

Discuss: What does a negative margin mean? How does a difference in vote counts
differ from a percentage or a difference in percentage points?

### Task 2: Numeric vectors and indexing

Rcpp provides vector classes that connect to R objects:

```r
cppFunction('
  NumericVector addTen(NumericVector x) {
    return x + 10.0;
  }
')

addTen(c(1, 2, 3))  # 11 12 13
```

`NumericVector` represents a vector of doubles. Here Rcpp's vector arithmetic
adds 10 to every element and returns a new result vector.

Access one element with square brackets:

```r
cppFunction('
  double firstValue(NumericVector x) {
    return x[0];
  }
')

firstValue(c(12, 18, 21, 30))  # 12
```

| Position | R expression | C++ expression | Value |
|---|---|---|---:|
| First | `x[1]` | `x[0]` | 12 |
| Second | `x[2]` | `x[1]` | 18 |
| Third | `x[3]` | `x[2]` | 21 |
| Fourth | `x[4]` | `x[3]` | 30 |

**Exercise 2.** Write `secondValue(NumericVector x)`, returning a `double`.
Return the second element using its C++ index.

```r
secondValue(c(12, 18, 21, 30))  # 18
secondValue(c(7, 9))           # 9
```

Before reading the element, add this length check:

```cpp
if (x.size() < 2) {
  stop("x must contain at least two values.");
}
```

Rcpp's `[]` access does not generally check bounds. Valid indexes for a vector
of length `n` run from `0` to `n - 1`. Do not call the unguarded demonstration
of `firstValue()` on an empty vector. The reference implementation adds a
similar length check to `firstValue()` so both accessors fail with an R error
on inputs that are too short.

Discuss: Why is the second C++ index `1`? Why does `secondValue()` return a
`double` even though its input is a `NumericVector`?

### Task 3: Logical results

A comparison can return a Boolean value:

```r
cppFunction('
  bool isGreater(double x, double y) {
    return x > y;
  }
')

isGreater(5, 3)  # TRUE
isGreater(3, 3)  # FALSE
```

C++ uses `true` and `false`; Rcpp converts the returned `bool` into R's
`TRUE` or `FALSE`. Use `==` for equality and `=` for assignment.

**Exercise 3.** Write `isEligible(double age)`, returning a `bool`. Return
`TRUE` when the person's age is **at least 18**, including age 18.
This is a simplified age rule for the programming exercise.

| R call | Expected result |
|---|---|
| `isEligible(16)` | `FALSE` |
| `isEligible(17.9)` | `FALSE` |
| `isEligible(18)` | `TRUE` |
| `isEligible(21)` | `TRUE` |

Discuss: Which comparison includes the boundary? Why does testing age 18 help
find a mistake that testing only ages 16 and 21 would miss?

## Compile the exercise functions from C++ files

The two teaching scripts use `cppFunction()`: C++ definitions live inside R
strings. After completing the exercises there, put each C++ function in its
matching file under `src/`. The supplied `.cpp` files contain reference
solutions; try writing your functions before consulting them.

| Exercise | C++ file | Function exported to R |
|---|---|---|
| 1 | `src/voteDifference.cpp` | `voteDifference()` |
| 2 | `src/secondValue.cpp` | `secondValue()` |
| 3 | `src/isEligible.cpp` | `isEligible()` |

Copy the C++ function definition itself, without the surrounding
`cppFunction('...')` call. Each separate C++ file includes:

```cpp
#include <Rcpp.h>
using namespace Rcpp;
```

The header makes Rcpp definitions available; the namespace directive lets you
write `NumericVector` instead of `Rcpp::NumericVector`. Place
`// [[Rcpp::export]]` immediately before the function so `sourceCpp()` creates
an R-callable wrapper. `cppFunction()` supplies that setup automatically.

With `Lab9` as the R working directory, compile a file and call its function:

```r
Rcpp::sourceCpp(file.path("src", "voteDifference.cpp"))
voteDifference(620, 580)  # 40

Rcpp::sourceCpp(file.path("src", "secondValue.cpp"))
secondValue(c(12, 18, 21, 30))  # 18

Rcpp::sourceCpp(file.path("src", "isEligible.cpp"))
isEligible(18)  # TRUE
```

## What the test script does

`tests/lab9_test.r` is a correctness check for the three `.cpp` exercise files:

1. Compile each file with `Rcpp::sourceCpp()` in a fresh environment.
2. Check that all three functions were exported and are callable from R.
3. Call the functions on inputs with known answers.
4. Print `PASS` messages when checks succeed; stop with an error if a check fails.

For example, `voteDifference(620, 580)` must return 40,
`secondValue(c(12, 18, 21, 30))` must return 18, and `isEligible(18)` must return
`TRUE`. Additional checks cover reversed vote totals, ties, short vectors,
and ages below 18. `stopifnot()` enforces these expectations. A fresh
environment prevents an existing R function from hiding a missing C++ export.

The test script reads the exercise implementations from `src/`; it does not
source the teaching scripts. `sumSquares.cpp` is a separate optional exercise.
No runner or result tables are needed for this workflow.

## How to run the tests

From the course repository root in R:

```r
setwd("Lab9")
source(file.path("tests", "lab9_test.r"))
```

Or, with `Lab9` as the shell working directory:

```sh
Rscript tests/lab9_test.r
```

Compilation requires temporary-file access. The script does not install
packages, download data, or save analysis outputs. The `PASS` messages and
final success message appear in the console. To try the teaching examples
independently, run `source(file.path("R", "lab9_exercises.r"))` or
`source(file.path("R", "lab9.r"))`.

## Optional extension: A loop in a C++ source file

Lecture 9 combines a numeric vector, zero-based indexing, and a running total.
Open `src/sumSquares.cpp` and find the input type, output type, accumulator,
loop condition, and return statement. Then run:

```r
Rcpp::sourceCpp(file.path("src", "sumSquares.cpp"))
sumSquares(c(2, 3))  # 13
sum(c(2, 3)^2)      # 13
```

`#include <Rcpp.h>` loads the definitions, `using namespace Rcpp;` makes Rcpp
names available, and `// [[Rcpp::export]]` marks the function for export to R.
For the earlier one-function strings, `cppFunction()` supplies this setup.
The loop uses `i < x.size()` to visit indexes `0` through `n - 1`.
An empty vector produces a sum of zero because the loop never runs.

After compiling, optionally compare with efficient vectorized R:

```r
# Install microbenchmark once if needed:
# install.packages("microbenchmark")
x <- seq(-1, 1, length.out = 100000)
stopifnot(isTRUE(all.equal(sumSquares(x), sum(x * x))))

microbenchmark::microbenchmark(
  R = sum(x * x),
  Rcpp = sumSquares(x),
  times = 30L
)
```

Compile **before** timing repeated function calls. Compilation is a separate
startup cost that still matters for one-off work. Do not assume Rcpp will win:
many R functions already use compiled code, and calling C++ has overhead.
Small scalar examples teach syntax and are not evidence of a speed gain.

Lecture 9's overall-speedup example is also a useful limit: if a function takes
80% of runtime and becomes 10 times faster, the entire analysis speeds up by
`1 / (0.20 + 0.80 / 10)`, about 3.57 times. Even eliminating that function's
runtime entirely would give at most a 5-fold speedup.

## Reproducibility and checks

All supplied examples are deterministic; no random seed or generated dataset
is needed. The test script compiles the three exercise files and checks signed
margins, ties, zero-based access, short-vector errors, logical return values,
and the exact age-18 boundary. It compares numeric answers with known values
using a floating-point tolerance. The optional sum-of-squares function can be
compiled and evaluated separately with the commands above.

## Author

Xiangyu Song

## Date

2026-10-27
