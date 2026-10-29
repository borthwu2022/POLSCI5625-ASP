######################################################
#                                                    #
# Lab 9: Getting started with Rcpp                   #
#                                                    #
# R/lab9_exercises.r                                 #
#                                                    #
# In-class examples and exercises for:               #
#   Task 1. scalar arithmetic and vote difference    #
#   Task 2. numeric vectors and element access       #
#   Task 3. comparisons and age eligibility          #
#                                                    #
######################################################

# Work from Lab9. See README.md for compiler setup.

if (!requireNamespace("Rcpp", quietly = TRUE)) {
  stop("Install Rcpp with install.packages('Rcpp') before starting Lab 9.",
       call. = FALSE)
}
library(Rcpp)


# Task 1. Scalars and arithmetic

cppFunction('
  double addNumbers(double x, double y) {
    return x + y;
  }
')

# Example: combine vote totals from two precincts.
precinct_one <- 120
precinct_two <- 80
addNumbers(precinct_one, precinct_two)  # 200

cppFunction('
  double votePercentage(double votes, double totalVotes) {

    double proportion = votes / totalVotes;

    return proportion * 100.0;
  }
')

# Example: calculate a candidate's percentage of all votes.
candidate_votes <- 450
total_votes <- 600
votePercentage(candidate_votes, total_votes)  # 75
# Use a positive total_votes value.

# Exercise 1. Return the signed margin: first total minus second total.
# Complete the function, then uncomment the entire cppFunction() block.
# cppFunction('
#   double voteDifference(double votesA, double votesB) {
#     // TODO: Return votesA minus votesB.
#   }
# ')

# Check your answer after compiling it:
# voteDifference(620, 580)  # 40
# voteDifference(580, 620)  # -40
# voteDifference(600, 600)  # 0


# Task 2. Numeric vectors and zero-based indexing

cppFunction('
  NumericVector addTen(NumericVector x) {
    return x + 10.0;
  }
')

# Example: add ten to every value in a numeric vector.
values <- c(1, 2, 3)
addTen(values)  # 11 12 13

cppFunction('
  double firstValue(NumericVector x) {
    return x[0];
  }
')

# Example: retrieve the first recorded vote total.
vote_totals <- c(12, 18, 21, 30)
firstValue(vote_totals)  # 12
# This supplied demonstration requires a nonempty vector.
# Do not call it with numeric(0): [] does not generally check bounds.

# Exercise 2. Return the second value of a numeric vector.
# Keep the length guard; fill in the return statement and uncomment the block.
# cppFunction('
#   double secondValue(NumericVector x) {
#     if (x.size() < 2) {
#       stop("x must contain at least two values.");
#     }
#     // TODO: Return the second value using its C++ index.
#   }
# ')

# secondValue(c(12, 18, 21, 30))  # 18
# secondValue(c(7, 9))           # 9
# secondValue(c(7))              # should raise an informative R error


# Task 3. Comparisons and Boolean return values

cppFunction('
  bool isGreater(double x, double y) {
    return x > y;
  }
')

# Example: check whether candidate A received more votes.
candidate_a <- 620
candidate_b <- 580
isGreater(candidate_a, candidate_b)  # TRUE
isGreater(580, 580)                 # FALSE

# Exercise 3. Return TRUE when age is at least 18, including age 18.
# cppFunction('
#   bool isEligible(double age) {
#     // TODO: Return the comparison for the age rule.
#   }
# ')

# isEligible(16)    # FALSE
# isEligible(17.9)  # FALSE
# isEligible(18)    # TRUE
# isEligible(21)    # TRUE

# After completing the exercises, put the three C++ function definitions in
# their matching src/*.cpp files, then run tests/lab9_test.r.
# See README.md for the sourceCpp() workflow and commands.
