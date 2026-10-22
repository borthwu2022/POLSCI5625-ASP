# Lab 8: Shingles, MinHash, and Locality-Sensitive Hashing

## Overview

Find similar BBC News articles using the Week 8 pipeline: represent documents
as shingle sets, summarize their contents with MinHash signatures, and use
locality-sensitive hashing (LSH) to select candidate pairs. Compare **2-, 4-,
and 10-character shingles**. Then use **4-character shingles** for MinHash and
LSH, and calculate the average bucket size for several band sizes.

This lab follows Week 8's random-permutation MinHash (slides 3-7), candidate
buckets and banding (slides 9-14), and the 360-row example (slides 15-16).

## Learning objectives

1. Create sets of overlapping character shingles and compare their sizes.
2. Calculate exact Jaccard similarity between shingle sets.
3. Write a MinHash function using a common random ordering for all documents.
4. Write an LSH banding function and remove repeated candidate pairs.
5. Calculate per-band and pooled average occupied-bucket sizes.
6. Explain how rows per band affect candidate counts and similarity screening.

## Project structure

```text
Lab8/
|-- data-raw/
|   |-- bbc-news-data.csv
|-- R/
|   |-- lab8.r
|   |-- run_lab8.r
|-- tests/
|   |-- lab8_test.r
|-- results/
|   |-- documents.csv
|   |-- category_counts.csv
|   |-- document_shingle_counts.csv
|   |-- shingle_summary.csv
|   |-- shingle_dictionary_4.csv
|   |-- minhash_signature_4.csv
|   |-- minhash_accuracy.csv
|   |-- pair_comparison.csv
|   |-- bucket_sizes_by_band.csv
|   |-- bucket_size_distribution.csv
|   |-- lsh_summary.csv
|   |-- candidate_probabilities.csv
|   |-- lsh_candidates_4.csv
|   |-- candidate_verification_4.csv
|   |-- shingles_minhash_lsh.pdf
|-- README.md
```

## Requirements

R only; no additional packages. Paths are relative to `Lab8`.

## Data and objective

The supplied BBC file contains **2,225 articles**: 510 business, 386
entertainment, 417 politics, 511 sport, and 401 tech. Its columns are
`category`, `filename`, `title`, and `content`. The original attachment is
copied unchanged into `data-raw/`.

## In-class exercise

Write your functions before consulting the reference implementations. As in
the previous labs, `R/lab8.r` contains functions, `R/run_lab8.r` runs the
analysis and saves outputs, and `tests/lab8_test.r` checks the calculations.

### Setup: Read and inspect the data

Although the filename ends in `.csv`, the separator is a **tab**. Read it with:

```r
news <- read.delim(file.path("data-raw", "bbc-news-data.csv"),
                   sep = "\t", quote = "", comment.char = "",
                   colClasses = "character", fileEncoding = "UTF-8")
news$doc_id <- paste(news$category, news$filename, sep = "/")
```

Using `quote = ""` preserves quotation marks within article text. Filenames
repeat across categories, so `category/filename` identifies a document.

Combine `title` and `content` with a space. Convert ASCII uppercase letters to
lowercase, replace each run of characters outside `a-z` and `0-9` with one
space, and trim leading/trailing spaces. Retain internal spaces, numbers, and
common words. These rules also replace accented characters with spaces.
Category labels are used only to describe results. Retain all records,
including repeated articles.

### Task 1: Create 2-, 4-, and 10-shingle representations

Write `normalize_text(text)` and `make_shingles(text, k)`.
A **k-shingle is k consecutive characters**, including spaces. Move one
character at a time and retain each distinct shingle once per document.
For example:

```r
make_shingles("banana", 2)   # "ba", "an", "na"
make_shingles("banana", 4)   # "bana", "anan", "nana"
make_shingles("a b", 2)      # "a ", " b"
```

For each $k \in \{2,4,10\}$, build one set per document. Report the number of
unique shingles across the corpus, the total document-shingle incidences,
and the minimum, median, mean, and maximum number of unique shingles per
article. Use sets rather than shingle frequencies.

Write `encode_shingles(shingle_sets)` to create one shared dictionary and map
each document to its integer shingle IDs. These ID lists represent the
nonzero entries of a characteristic matrix whose rows are shingles and
columns are documents. A dense matrix is unnecessary.

Write `jaccard(a, b)`:

```math
J(A,B)=\frac{|A\cap B|}{|A\cup B|}.
```

Using the same 1,000 random document pairs (seed 123), compare exact Jaccard
similarities for all three values of $k$. Save the representation summary
and pair comparisons.

### Task 2: Write a MinHash function for 4-shingles

Use the 4-shingle dictionary and set **m = 360** and **seed = 123**.
Write `minhash_row(shingle_ids, ordering)` and
`minhash_signature(shingle_ids, n_shingles, m, seed)`.

For each signature row:

1. Generate a random permutation of all shingle row IDs.
2. Use that **same permutation for every document**.
3. Find the first shingle in the ordering that each document contains.
4. Store that shingle's **original row ID** in the signature.

A convenient implementation assigns each ID its rank in the permutation,
then finds the present ID with the smallest rank. The returned signature has
360 rows and 2,225 document columns. Sample a fresh ordering for each row.
Generate one ordering at a time.

Check the lecture's example:

```text
A = {1, 2, 3, 4}       B = {1, 2, 3, 5}
Ordering               h(A)   h(B)
3, 1, 5, 2, 4           3      3
4, 2, 1, 5, 3           4      2
5, 3, 4, 1, 2           3      5
2, 5, 1, 3, 4           2      2
1, 4, 3, 5, 2           1      1
```

Three of five rows match, giving an estimate of 0.6. The minimum original ID
would be 1 for both sets in every row; random ordering is what makes MinHash
informative.

For the BBC evaluation pairs, calculate

```math
\widehat J(A,B)=\frac{1}{360}\sum_{h=1}^{360}
\mathbf{1}\{\operatorname{sig}_h(A)=\operatorname{sig}_h(B)\}.
```

Compare this estimate with exact 4-shingle Jaccard. On the 1,000 random pairs,
report mean error, mean absolute error, RMSE, and maximum absolute error.
Save the full signature and its dictionary so each signature entry can be
traced back to a shingle.

### Task 3: Write LSH banding for the 4-shingle signature

Write `lsh_bands(signature, rows_per_band)`.
Let **r be rows per band**, **b be the number of bands**, and **m = br = 360**.
Start with the lecture's example, **r = 4 and b = 90**.

For each band, use the tuple of its $r$ signature entries as a document's
bucket key. Keep band identities separate. Two documents become a candidate
pair when all entries match in **at least one common band**. Matching only
some rows of a band is insufficient.

Build the bucket membership lists, generate unordered pairs within buckets
containing at least two documents, and retain each pair only once across
bands. Equal tuples are used directly as keys, so the implementation adds
no extra collisions from mapping keys into a smaller hash table.

Reuse the **same 360-row signature** with:

| Rows per band $r$ | Number of bands $b=360/r$ |
|---:|---:|
| 2 | 180 |
| 4 | 90 |
| 10 | 36 |
| 20 | 18 |
| 30 | 12 |
| 60 | 6 |

For each setting, report the unique candidate count and the fraction of all
$\binom{2225}{2}=2,474,200$ possible document pairs selected. The sum of
within-bucket pair counts can exceed the unique candidate count because a
pair may collide in multiple bands.

### Task 4: Calculate average bucket size given band size

For band $t$, let $B_t$ be its number of occupied buckets and $n_{tj}$ the
number of documents in occupied bucket $j$. Every document enters exactly
one bucket in each band, so $\sum_{j=1}^{B_t} n_{tj}=N=2225$.

Calculate the average occupied-bucket size **within each band**:

```math
\overline n_t=\frac{\sum_{j=1}^{B_t}n_{tj}}{B_t}=\frac{N}{B_t}.
```

Then calculate the **pooled mean across all occupied buckets**, the primary
average reported in `lsh_summary.csv`:

```math
\overline n_{\mathrm{pooled}}
=\frac{\sum_{t=1}^{b}\sum_{j=1}^{B_t}n_{tj}}{\sum_{t=1}^{b}B_t}
=\frac{bN}{\sum_{t=1}^{b}B_t}.
```

Include singleton buckets of size 1. We use exact tuple keys, so there is no
fixed table capacity or count of unused bucket slots in this calculation.

Also report the unweighted average of the band means,
$\frac{1}{b}\sum_t N/B_t$. It can differ from the pooled mean because the
bands can have different numbers of occupied buckets.

For example, suppose four documents produce bucket sizes `(3, 1)` in band 1
and `(2, 1, 1)` in band 2. The band means are 2 and 4/3. The pooled mean is
`8/5 = 1.6`; the average of band means is `5/3 = 1.6667`.

Save one row per band with its occupied-bucket count, mean, largest bucket,
singleton count, and within-bucket pair occurrences. Also save a frequency
table of bucket sizes. Verify that the bucket sizes sum to $N$ in every band,
and to $bN$ within each band-size setting.

### Task 5: Verify candidates and explain the tradeoff

For **r = 4**, calculate exact Jaccard on the original 4-shingle sets for
**every unique candidate pair**. Save document IDs, exact similarities, and
MinHash estimates, sorted by exact similarity. Report the number and fraction
of candidates with exact Jaccard at least **0.8**. Inspect several article
pairs using `documents.csv` and the original data.

Using the lecture's independent MinHash rows, calculate and plot:

```math
P(\text{candidate}\mid s)=1-(1-s^r)^b,\qquad b=360/r.
```

Calculate this probability at $s=0.25,0.50,0.75,0.90$ for every tested band
size. For $r=4,b=90$, the probabilities at 0.25 and 0.75 are approximately
0.297 and 1, as in the lecture.

## How to run

From the course repository root in R:

```r
setwd("Lab8")
source(file.path("tests", "lab8_test.r"))
source(file.path("R", "run_lab8.r"))
```

Or, with `Lab8` as the shell working directory:

```sh
Rscript tests/lab8_test.r
Rscript R/run_lab8.r
```

The tests use small hand-checkable examples plus checks on the source data.
The runner performs the full-corpus calculations, saves 14 CSVs and a two-page
PDF, and checks bucket totals and candidate deduplication. No data-generation
script is needed because the supplied BBC file is the input.

## Reference results

Run `R/run_lab8.r` to regenerate the saved calculations. The reference tables
below report the accompanying output files.

Reference calculations use all 2,225 articles, 360 MinHash rows, and R 4.5.2.

| Shingle length | Corpus vocabulary | Mean unique shingles per article | Mean Jaccard, random pairs |
|---:|---:|---:|---:|
| 2 | 974 | 281.46 | 0.626154 |
| 4 | 54,594 | 1,321.87 | 0.128251 |
| 10 | 2,254,375 | 2,068.49 | 0.002820 |

Short 2-shingles give a mean random-pair Jaccard of 0.6262. The same
pairs have mean similarity 0.1283 with 4-shingles and 0.00282 with
10-shingles, illustrating how short character sequences are widely shared.

For 4-shingle MinHash, mean absolute error is **0.013949** and
RMSE is **0.017567** on the 1,000 random pairs. The observed mean
error is -0.007576; maximum absolute error is 0.068998.

| Rows per band | Bands | Pooled mean bucket size | Mean of band means | Unique candidates | All pairs selected |
|---:|---:|---:|---:|---:|---:|
| 2 | 180 | 3.258270 | 3.531970 | 2,175,151 | 87.9133% |
| 4 | 90 | 1.192340 | 1.197243 | 47,717 | 1.9286% |
| 10 | 36 | 1.074115 | 1.074117 | 170 | 0.0069% |
| 20 | 18 | 1.071629 | 1.071631 | 164 | 0.0066% |
| 30 | 12 | 1.069883 | 1.069884 | 160 | 0.0065% |
| 60 | 6 | 1.066549 | 1.066550 | 151 | 0.0061% |

Each average includes singleton buckets. With more rows per band, the
buckets are more selective and fewer pairs reach exact comparison.

At **r = 4**, **164 of 47,717 candidates** have exact 4-shingle
Jaccard at least 0.8, giving candidate precision **0.3437%**. LSH selects
**1.9286%** of all possible pairs for exact verification.
The saved candidate table lets you inspect which articles share text.


## Discuss the results

1. How do vocabulary size, average set size, and random-pair Jaccard change
   between 2-, 4-, and 10-shingles? Why do short shingles occur in many articles?
2. Why must every document use the same ordering within a MinHash row? Why
   must each signature row use a fresh random ordering?
3. Why does the fraction of equal signature entries estimate Jaccard? Does
   agreement in all 360 rows prove the original sets are equal?
4. How do average bucket size and candidate count change as rows per band
   increase? Account for the corresponding change in the number of bands.
5. Why do pooled bucket means and averages of band means differ? Why does
   excluding singleton buckets answer a different question?
6. Why must candidate pairs be deduplicated before counting comparisons?
7. How do repeated articles affect large buckets? Would LSH still generate
   few candidates if many documents had identical contents?
8. At a threshold of 0.8, what fraction of the r = 4 candidates pass exact
   verification? Candidate verification measures precision; measuring recall
   also requires finding highly similar pairs that LSH did not select.
9. Why can candidate generation still produce quadratic work when a bucket
   is large, even though bucket assignment avoids comparing every pair?

## Reproducibility and checks

The runner uses all 2,225 records in their original order. Both the 1,000
sampled evaluation pairs and the 360 independent random shingle orderings use
seed 123. The shared shingle dictionary is sorted using radix order. Seeded
sampling uses R's Mersenne-Twister generator and Rejection sampler.

Tests reproduce the lecture's five MinHash signatures and six LSH candidates,
verify the MinHash identity by enumerating all 120 permutations of five IDs,
and compare banding with an independent pairwise definition on small
examples. They also check singleton buckets, separate band identities,
repeated candidate pairs, unequal band bucket counts, and source-data totals.

## Author

Xiangyu Song

## Date

2026-10-20
