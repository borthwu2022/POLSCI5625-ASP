######################################################
#                                                    #
# Lab 8: Shingles, MinHash, and LSH                  #
#                                                    #
# tests/lab8_test.r                                  #
#                                                    #
# Check shingles and the hand-coded MinHash updates. #
# Test LSH candidates and average bucket sizes.      #
#                                                    #
######################################################

source(file.path("R", "lab8.r"))

# Character windows, spaces, and set membership.
stopifnot(
  identical(normalize_text("  BBC: NEWS!! 2026  "), "bbc news 2026"),
  identical(make_shingles("banana", 2L), c("ba", "an", "na")),
  identical(make_shingles("banana", 4L), c("bana", "anan", "nana")),
  identical(make_shingles("a b", 2L), c("a ", " b")),
  length(make_shingles("short", 10L)) == 0L,
  jaccard(c(1, 2, 2), c(2, 3)) == 1 / 3,
  jaccard(integer(), integer()) == 1
)
encoded <- encode_shingles(list(c("ba", "an"), c("an", "na")))
stopifnot(identical(encoded$vocabulary, c("an", "ba", "na")),
          identical(encoded$ids, list(c(1L, 2L), c(1L, 3L))))
cat("PASS: character shingles, deduplication, dictionary, and exact Jaccard.\n")

# Week 8, slides 4-5. Store original IDs of first-present shingles.
sets <- list(A = 1:4, B = c(1L, 2L, 3L, 5L))
orders <- rbind(c(3L, 1L, 5L, 2L, 4L), c(4L, 2L, 1L, 5L, 3L),
                c(5L, 3L, 4L, 1L, 2L), c(2L, 5L, 1L, 3L, 4L),
                c(1L, 4L, 3L, 5L, 2L))
sig <- t(vapply(seq_len(nrow(orders)), function(h) {
  minhash_row(sets, orders[h, ])
}, integer(2)))
stopifnot(all(sig[, 1L] == c(3L, 4L, 3L, 2L, 1L)),
          all(sig[, 2L] == c(3L, 2L, 5L, 2L, 1L)),
          signature_similarity(sig, 1L, 2L) == 0.6,
          jaccard(sets$A, sets$B) == 0.6)

# Enumerate all 5! permutations: equality probability is exactly 3/5.
permutations <- function(x) {
  if (length(x) == 1L) return(matrix(x, nrow = 1L))
  do.call(rbind, lapply(seq_along(x), function(i) cbind(x[i], permutations(x[-i]))))
}
all_orders <- permutations(1:5)
agree <- apply(all_orders, 1L, function(ordering) {
  # Independent direct scan, with no rank vector.
  first_a <- ordering[ordering %in% sets$A][1L]
  first_b <- ordering[ordering %in% sets$B][1L]
  actual <- minhash_row(sets, ordering)
  stopifnot(all(actual == c(first_a, first_b)))
  actual[1L] == actual[2L]
})
stopifnot(sum(agree) == 72L, length(agree) == 120L)
cat("PASS: slide signatures and exact MinHash identity over all 120 orderings.\n")

one <- minhash_signature(sets, 5L, m = 120L, seed = 123L)
two <- minhash_signature(sets, 5L, m = 120L, seed = 123L)
duplicate <- minhash_signature(list(1:4, 1:4, 5L), 5L, m = 120L)
stopifnot(identical(one, two), all(duplicate[, 1L] == duplicate[, 2L]),
          all(duplicate[, 1L] != duplicate[, 3L]))
cat("PASS: seeded reproducibility, identical sets, and disjoint sets.\n")

# Week 8, slide 11: records 2-5 form all six candidates; record 1 has none.
slide_sig <- cbind(c(30, 8, 6, 2, 102, 16, 8, 112, 1, 76),
                   matrix(rep(c(30, 8, 1, 2, 102, 16, 8, 161, 1, 27), 4L), ncol = 4L))
slide_lsh <- lsh_bands(slide_sig, rows_per_band = 5L)
expected <- t(combn(2:5, 2L))
stopifnot(nrow(slide_lsh$candidates) == 6L,
          setequal(pair_codes(slide_lsh$candidates$i, slide_lsh$candidates$j, 5L),
                   pair_codes(expected[, 1L], expected[, 2L], 5L)),
          slide_lsh$summary$pair_occurrences == 12,
          slide_lsh$summary$mean_bucket_size == 2.5)

# Independent pairwise criterion on a small fixture for every divisor of m.
fixture <- cbind(c(1, 2, 3, 4), c(1, 2, 8, 9), c(5, 6, 3, 4), c(1, 2, 3, 4))
for (r in c(1L, 2L, 4L)) {
  fit <- lsh_bands(fixture, r)
  all_pairs <- t(combn(seq_len(ncol(fixture)), 2L))
  selected <- apply(all_pairs, 1L, function(ij) {
    matches <- fixture[, ij[1L]] == fixture[, ij[2L]]
    any(vapply(split(matches, rep(seq_len(nrow(fixture) / r), each = r)),
               all, logical(1)))
  })
  stopifnot(setequal(pair_codes(fit$candidates$i, fit$candidates$j, 4L),
                     pair_codes(all_pairs[selected, 1L], all_pairs[selected, 2L], 4L)),
            !anyDuplicated(pair_codes(fit$candidates$i, fit$candidates$j, 4L)),
            sum(fit$distribution$bucket_size * fit$distribution$bucket_count) ==
              ncol(fixture) * nrow(fixture) / r)
}
# The same tuple appearing in DIFFERENT bands does not create a candidate.
separate_bands <- lsh_bands(cbind(c(1, 2, 3, 4), c(3, 4, 1, 2)), 2L)
stopifnot(nrow(separate_bands$candidates) == 0L,
          separate_bands$summary$mean_bucket_size == 1)
cat("PASS: slide candidates, pair deduplication, band identity, and bucket accounting.\n")

# Means differ when bands contain different numbers of occupied buckets.
# Band 1: (3, 1), mean 2; band 2: (2, 1, 1), mean 4/3.
# Pooled mean = 8/5, while mean of the two band means = 5/3.
unequal <- lsh_bands(rbind(c(1, 1, 1, 2), c(1, 1, 2, 3)), 1L)
stopifnot(unequal$summary$occupied_buckets == 5L,
          unequal$summary$mean_bucket_size == 8 / 5,
          abs(unequal$summary$mean_of_band_means - 5 / 3) < 1e-12,
          candidate_probability(0, 4L) == 0,
          candidate_probability(1, 4L) == 1,
          abs(candidate_probability(0.25, 4L) - (1 - (1 - 0.25^4)^90)) < 1e-12,
          abs(candidate_probability(0.75, 4L) - (1 - (1 - 0.75^4)^90)) < 1e-12)
cat("PASS: both definitions of mean bucket size and the lecture probability formula.\n")

news <- read_bbc()
stopifnot(nrow(news) == 2225L, length(unique(news$category)) == 5L,
          !anyDuplicated(news$doc_id),
          all(table(news$category) == c(business = 510L, entertainment = 386L,
                                      politics = 417L, sport = 511L, tech = 401L)))
text <- normalize_text(paste(news$title, news$content))
stopifnot(all(nchar(text) >= 10L))
pairs <- sample_document_pairs(nrow(news), size = 1000L, seed = 123L)
stopifnot(nrow(pairs) == 1000L, all(pairs$i < pairs$j),
          !anyDuplicated(pair_codes(pairs$i, pairs$j, nrow(news))))
cat("PASS: all 2,225 BBC records, category counts, document IDs, and sampled pairs.\n")
cat("All Lab 8 tests passed.\n")
