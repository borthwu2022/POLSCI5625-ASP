######################################################
#                                                    #
# Lab 8: Shingles, MinHash, and LSH                  #
#                                                    #
# R/run_lab8.r                                       #
#                                                    #
# Read the BBC data and compare shingle sets.        #
# Calculate MinHash signatures and LSH candidates.   #
# Save bucket-size summaries and comparison plots.   #
#                                                    #
######################################################

source(file.path("R", "lab8.r"))
dir.create("results", showWarnings = FALSE)
save_csv <- function(x, name) {
  write.csv(x, file.path("results", name), row.names = FALSE, na = "NA")
}

news <- read_bbc()
text <- normalize_text(paste(news$title, news$content))
n <- nrow(news)
cat("Read", n, "BBC articles. Creating character shingle sets.\n")
metadata <- news[, c("doc_id", "category", "filename", "title")]
metadata$normalized_characters <- nchar(text)
save_csv(metadata, "documents.csv")
save_csv(as.data.frame(table(category = news$category)), "category_counts.csv")

# Use the same 1,000 uniformly sampled pairs for every shingle size.
pairs <- sample_document_pairs(n, size = 1000L, seed = 123L)
pair_comparison <- pairs
pair_comparison$doc_id_i <- news$doc_id[pairs$i]
pair_comparison$doc_id_j <- news$doc_id[pairs$j]
pair_comparison$same_category <- news$category[pairs$i] == news$category[pairs$j]
representation_summary <- vector("list", 3L)
shingle_sizes <- data.frame(doc_id = news$doc_id)

# Task 1. Build each representation in turn. Keep integer sets only for k = 4.
for (index in seq_along(c(2L, 4L, 10L))) {
  k <- c(2L, 4L, 10L)[index]
  sets <- lapply(text, make_shingles, k = k)
  names(sets) <- news$doc_id
  sizes <- lengths(sets)
  vocabulary_size <- length(unique(unlist(sets, use.names = FALSE)))
  exact <- vapply(seq_len(nrow(pairs)), function(p) {
    jaccard(sets[[pairs$i[p]]], sets[[pairs$j[p]]])
  }, numeric(1))
  shingle_sizes[[paste0("shingles_", k)]] <- sizes
  pair_comparison[[paste0("jaccard_", k)]] <- exact
  representation_summary[[index]] <- data.frame(
    k = k, documents = n, vocabulary_size = vocabulary_size,
    total_incidences = sum(sizes), min_shingles = min(sizes),
    median_shingles = median(sizes), mean_shingles = mean(sizes),
    max_shingles = max(sizes),
    mean_random_jaccard = mean(exact),
    median_random_jaccard = median(exact)
  )
  if (k == 4L) encoded4 <- encode_shingles(sets)
  rm(sets)
  cat("  Finished k =", k, "; vocabulary =", vocabulary_size, "\n")
}
representation_summary <- do.call(rbind, representation_summary)
save_csv(representation_summary, "shingle_summary.csv")
save_csv(shingle_sizes, "document_shingle_counts.csv")
save_csv(data.frame(shingle_id = seq_along(encoded4$vocabulary),
                    shingle = encoded4$vocabulary), "shingle_dictionary_4.csv")

# Task 2. A 360 x n signature: one independent ordering per row, shared by all docs.
cat("Calculating the 4-shingle MinHash signature (m = 360, seed = 123).\n")
signature4 <- minhash_signature(encoded4$ids, length(encoded4$vocabulary),
                                m = 360L, seed = 123L, verbose = TRUE)
save_csv(data.frame(hash_row = seq_len(nrow(signature4)), signature4,
                    check.names = FALSE), "minhash_signature_4.csv")
pair_comparison$minhash_4 <- vapply(seq_len(nrow(pairs)), function(p) {
  signature_similarity(signature4, pairs$i[p], pairs$j[p])
}, numeric(1))
pair_comparison$error_4 <- pair_comparison$minhash_4 - pair_comparison$jaccard_4
random_errors <- pair_comparison$error_4
minhash_accuracy <- data.frame(
  k = 4L, signature_rows = 360L, seed = 123L, random_pairs = length(random_errors),
  mean_error = mean(random_errors), mean_absolute_error = mean(abs(random_errors)),
  rmse = sqrt(mean(random_errors^2)), max_absolute_error = max(abs(random_errors))
)
save_csv(minhash_accuracy, "minhash_accuracy.csv")

# Tasks 3-4. Reuse ONE signature and vary rows per band, keeping m = 360 fixed.
# r = 4 reproduces the lecture's 90-band example.
band_sizes <- c(2L, 4L, 10L, 20L, 30L, 60L)
summaries <- band_statistics <- bucket_distributions <- vector("list", length(band_sizes))
for (index in seq_along(band_sizes)) {
  r <- band_sizes[index]
  cat("LSH with", r, "rows per band and", 360L / r, "bands.\n")
  lsh <- lsh_bands(signature4, rows_per_band = r)
  summaries[[index]] <- lsh$summary
  band_statistics[[index]] <- lsh$per_band
  bucket_distributions[[index]] <- lsh$distribution
  candidate_codes <- pair_codes(lsh$candidates$i, lsh$candidates$j, n)
  pair_comparison[[paste0("candidate_r", r)]] <-
    pair_codes(pairs$i, pairs$j, n) %in% candidate_codes
  if (r == 4L) candidates4 <- lsh$candidates
  cat("  Mean occupied bucket size:", format(lsh$summary$mean_bucket_size, digits = 7),
      "; unique candidate pairs:", nrow(lsh$candidates), "\n")
  rm(lsh, candidate_codes)
}
lsh_summary <- do.call(rbind, summaries)
per_band <- do.call(rbind, band_statistics)
bucket_distribution <- do.call(rbind, bucket_distributions)
save_csv(lsh_summary, "lsh_summary.csv")
save_csv(per_band, "bucket_sizes_by_band.csv")
save_csv(bucket_distribution, "bucket_size_distribution.csv")
save_csv(pair_comparison, "pair_comparison.csv")

# The lecture's probability calculation, plus other band sizes for comparison.
probabilities <- expand.grid(similarity = c(0.25, 0.5, 0.75, 0.9),
                             rows_per_band = band_sizes)
probabilities$bands <- 360 / probabilities$rows_per_band
probabilities$candidate_probability <- mapply(candidate_probability,
                                              probabilities$similarity,
                                              probabilities$rows_per_band)
save_csv(probabilities, "candidate_probabilities.csv")

# Task 5. Verify ALL default-r=4 candidates with their original 4-shingle sets.
cat("Calculating exact 4-shingle Jaccard for", nrow(candidates4), "r = 4 candidates.\n")
verified <- score_pairs(candidates4, encoded4$ids, signature4, verbose = TRUE)
verified$doc_id_i <- news$doc_id[verified$i]
verified$doc_id_j <- news$doc_id[verified$j]
verified$same_category <- news$category[verified$i] == news$category[verified$j]
verified <- verified[order(-verified$jaccard, verified$i, verified$j), ]
rownames(verified) <- NULL
save_csv(verified, "lsh_candidates_4.csv")
threshold <- 0.8
verification <- data.frame(
  rows_per_band = 4L, bands = 90L, threshold = threshold,
  candidate_pairs = nrow(verified),
  pairs_at_or_above_threshold = sum(verified$jaccard >= threshold),
  candidate_precision = if (nrow(verified)) mean(verified$jaccard >= threshold) else NA_real_,
  fraction_of_all_pairs_compared = nrow(verified) / choose(n, 2)
)
save_csv(verification, "candidate_verification_4.csv")

# Independent accounting checks on the full corpus.
stopifnot(
  all(dim(signature4) == c(360L, n)),
  !anyNA(signature4),
  all(lsh_summary$assignments == n * lsh_summary$bands),
  all(lsh_summary$candidate_pairs <= choose(n, 2)),
  all(lsh_summary$pair_occurrences >= lsh_summary$candidate_pairs),
  all(abs(lsh_summary$mean_bucket_size -
          lsh_summary$assignments / lsh_summary$occupied_buckets) < 1e-12),
  all(vapply(split(bucket_distribution, bucket_distribution$rows_per_band),
             function(x) sum(x$bucket_size * x$bucket_count) == n * x$bands[1L],
             logical(1))),
  !anyDuplicated(pair_codes(verified$i, verified$j, n))
)

# One two-page PDFs
save_figures <- function() {
  pdf(file.path("results", "shingles_minhash_lsh.pdf"), width = 11, height = 8.5)
  on.exit(dev.off(), add = TRUE)
  colors <- c("#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#555555")
  par(mfrow = c(2, 2), mar = c(4.2, 5.2, 3.3, 1), oma = c(1, 0, 2, 0),
      mgp = c(3.3, 0.8, 0), las = 1)
  boxplot(pair_comparison[, c("jaccard_2", "jaccard_4", "jaccard_10")],
          names = c("k = 2", "k = 4", "k = 10"), col = colors[1:3],
          ylab = "Exact Jaccard similarity", main = "Same 1,000 random article pairs",
          ylim = c(0, 1), outline = FALSE)
  plot(pair_comparison$jaccard_4, pair_comparison$minhash_4,
       pch = 16, cex = 0.6, col = adjustcolor(colors[1], 0.35),
       xlim = c(0, 1), ylim = c(0, 1), xlab = "Exact 4-shingle Jaccard",
       ylab = "Fraction of matching MinHash rows", main = "360-row MinHash estimate")
  abline(0, 1, col = colors[2], lwd = 2)
  barplot(representation_summary$mean_shingles, names.arg = c("k = 2", "k = 4", "k = 10"),
          col = colors[1:3], ylab = "Mean unique shingles per article",
          main = "Size of each set representation")
  hist(random_errors, breaks = 25, col = colors[1], border = "white",
       xlab = "MinHash estimate - exact Jaccard", main = "Error on 1,000 random pairs")
  abline(v = 0, col = colors[2], lwd = 2)
  mtext("BBC News: character shingles and MinHash", outer = TRUE, font = 2, cex = 1.3)

  mtext("1 / 2", side = 1, outer = TRUE, line = 0, cex = 0.7)

  # Reset the panel layout for page 2.
  par(mfrow = c(2, 2), mar = c(4.2, 5.2, 3.3, 1), oma = c(1, 0, 2, 0),
      mgp = c(3.3, 0.8, 0), las = 1)
  s <- seq(0, 1, length.out = 501L)
  plot(s, candidate_probability(s, band_sizes[1L]), type = "l", lwd = 2,
       col = colors[1], ylim = c(0, 1), yaxt = "n", xlab = "True Jaccard similarity",
       ylab = "Probability of becoming a candidate", main = "Banding curves: m = 360")
  for (index in seq_along(band_sizes)[-1L]) {
    lines(s, candidate_probability(s, band_sizes[index]), col = colors[index], lwd = 2)
  }
  axis(2, at = seq(0, 1, by = 0.2), labels = sprintf("%.1f", seq(0, 1, by = 0.2)))
  box()
  legend("bottomright", legend = paste0("r = ", band_sizes, ", b = ", 360 / band_sizes),
         col = colors, lwd = 2, cex = 0.75, bty = "n")
  plot(lsh_summary$rows_per_band, lsh_summary$mean_bucket_size, type = "b", pch = 16,
       log = "x", col = colors[1], xaxt = "n", xlab = "Rows per band (log scale)",
       ylab = "Mean occupied-bucket size", main = "Pooled mean: b * n / occupied buckets")
  axis(1, at = band_sizes, labels = band_sizes)
  plot(lsh_summary$rows_per_band, lsh_summary$candidate_fraction, type = "b", pch = 16,
       log = "x", col = colors[2], xaxt = "n", ylim = c(0, 1),
       xlab = "Rows per band (log scale)", ylab = "Fraction of all unordered pairs",
       main = "Unique candidate pairs / choose(n, 2)")
  axis(1, at = band_sizes, labels = band_sizes)
  candidate_hist <- hist(verified$jaccard, breaks = seq(0, 1, by = 0.05), plot = FALSE)
  candidate_hist$counts <- candidate_hist$counts / 1000
  plot(candidate_hist, freq = TRUE, col = colors[3], border = "white",
       xlab = "Exact Jaccard among r = 4 candidates", ylab = "Candidate pairs (thousands)",
       main = "Verify candidates with original sets")
  abline(v = threshold, col = colors[2], lwd = 2, lty = 2)
  mtext("BBC News: LSH and bucket size", outer = TRUE, font = 2, cex = 1.3)
  mtext("2 / 2", side = 1, outer = TRUE, line = 0, cex = 0.7)
}
save_figures()

cat("\nShingle representations:\n")
print(representation_summary, row.names = FALSE, digits = 6)
cat("\nMinHash accuracy on random pairs:\n")
print(minhash_accuracy, row.names = FALSE, digits = 6)
cat("\nLSH and mean occupied-bucket sizes:\n")
print(lsh_summary[, c("rows_per_band", "bands", "mean_bucket_size",
                      "mean_of_band_means", "candidate_pairs", "candidate_fraction")],
      row.names = FALSE, digits = 7)
cat("\nCandidate verification at exact Jaccard >= 0.8 (r = 4):\n")
print(verification, row.names = FALSE, digits = 6)
cat("\nAll accounting checks passed. Saved CSVs and a two-page PDF in results/.\n")
