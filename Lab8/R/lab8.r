######################################################
#                                                    #
# Lab 8: Shingles, MinHash, and LSH                  #
#                                                    #
# R/lab8.r                                           #
#                                                    #
# In-class reference implementations for:            #
#   Task 1. shingles and exact Jaccard               #
#   Task 2. random-permutation MinHash               #
#   Tasks 3-4. LSH and average bucket sizes          #
#   Task 5. candidate verification                   #
#                                                    #
######################################################

# Setup. Use the title and content; category is metadata, not an input feature.
read_bbc <- function(path = file.path("data-raw", "bbc-news-data.csv")) {
  news <- read.delim(path, sep = "\t", quote = "", comment.char = "",
                     colClasses = "character", fileEncoding = "UTF-8",
                     check.names = FALSE)
  stopifnot(all(c("category", "filename", "title", "content") %in% names(news)))
  news$doc_id <- paste(news$category, news$filename, sep = "/")
  stopifnot(!anyDuplicated(news$doc_id), !anyNA(news))
  news
}


# Task 1. A shingle is k consecutive CHARACTERS, including internal spaces.

# This makes shingle comparisons insensitive to capitalization, punctuation,
# and extra whitespace while preserving numbers and single spaces between words.
normalize_text <- function(text) {
  # Explicit ASCII case conversion makes preprocessing independent of locale.
  text <- chartr("ABCDEFGHIJKLMNOPQRSTUVWXYZ", "abcdefghijklmnopqrstuvwxyz", text)
  text <- gsub("[^a-z0-9]+", " ", text, perl = TRUE) # replaces punctuation, extra spaces, tabs, and newlines with a single space
  trimws(text)
}


make_shingles <- function(text, k = 4L) {
  stopifnot(length(text) == 1L, !is.na(text),
            length(k) == 1L, is.finite(k), k >= 1L, k == floor(k))
  if (nchar(text, type = "chars") < k) return(character())
  starts <- seq_len(nchar(text, type = "chars") - k + 1L)
  unique(substring(text, starts, starts + k - 1L))
}

# Store the nonzero row IDs of each document, rather than a dense 0/1 matrix.
# Rows use one common dictionary for the whole corpus, sorted in ASCII order.
encode_shingles <- function(shingle_sets) {
  vocabulary <- sort(unique(unlist(shingle_sets, use.names = FALSE)), method = "radix")
  ids <- lapply(shingle_sets, function(x) sort(match(x, vocabulary)))
  names(ids) <- names(shingle_sets)
  list(vocabulary = vocabulary, ids = ids)
}

jaccard <- function(a, b) {
  a <- unique(a)
  b <- unique(b)
  intersection <- sum(a %in% b)
  union <- length(a) + length(b) - intersection
  if (union == 0L) return(1)
  intersection / union
}


# Task 2. Exactly follow Week 8: select the ORIGINAL shingle ID that comes first
# in a common random ordering. The signature entry is not the minimum ID.
minhash_row <- function(shingle_ids, ordering) {
  stopifnot(length(ordering) > 0L, all(lengths(shingle_ids) > 0L))
  ranks <- integer(length(ordering))
  ranks[ordering] <- seq_along(ordering)
  vapply(shingle_ids, function(ids) ids[which.min(ranks[ids])], integer(1))
}

minhash_signature <- function(shingle_ids, n_shingles, m = 360L,
                              seed = 123L, verbose = FALSE) {
  stopifnot(length(shingle_ids) > 0L, all(lengths(shingle_ids) > 0L),
            length(m) == 1L, is.finite(m), m >= 1L, m == floor(m),
            n_shingles >= 1L,
            all(vapply(shingle_ids, function(x) {
              all(x >= 1L & x <= n_shingles & x == floor(x))
            }, logical(1))))
  set.seed(seed)
  signature <- matrix(0L, nrow = m, ncol = length(shingle_ids),
                      dimnames = list(paste0("h", seq_len(m)), names(shingle_ids)))
  for (h in seq_len(m)) {
    ordering <- sample.int(n_shingles)
    signature[h, ] <- minhash_row(shingle_ids, ordering)
    if (verbose && h %% 60L == 0L) {
      cat("  Completed", h, "of", m, "MinHash rows.\n")
    }
  }
  signature
}

signature_similarity <- function(signature, i, j) {
  mean(signature[, i] == signature[, j])
}


# Task 3. Put equal tuples into the same bucket WITHIN each band.
# Numeric pair codes allow us to remove repeated unordered pairs across bands.
# For i < j, code = i + n * (j - 1). No document-pair similarity is needed here.
pair_codes <- function(i, j, n) {
  pmin(i, j) + as.double(n) * (pmax(i, j) - 1L)
}

decode_pairs <- function(codes, n) {
  data.frame(i = as.integer((codes - 1) %% n + 1),
             j = as.integer((codes - 1) %/% n + 1))
}

lsh_bands <- function(signature, rows_per_band = 4L) {
  m <- nrow(signature)
  n <- ncol(signature)
  r <- rows_per_band
  stopifnot(is.matrix(signature), !anyNA(signature), m > 0L, n > 0L,
            length(r) == 1L, is.finite(r), r >= 1L, r == floor(r), m %% r == 0L)
  b <- m %/% r
  buckets <- vector("list", b)
  per_band <- distributions <- band_pairs <- vector("list", b)

  for (band in seq_len(b)) {
    rows <- seq.int((band - 1L) * r + 1L, band * r)
    # Delimiters make tuples unambiguous. Each band has its own key space.
    keys <- apply(signature[rows, , drop = FALSE], 2L, paste, collapse = ":")
    members <- split(seq_len(n), keys)
    sizes <- lengths(members)
    buckets[[band]] <- members

    # Task 4. Include occupied singleton buckets when calculating the mean.
    per_band[[band]] <- data.frame(
      rows_per_band = r, bands = b, band = band, documents = n,
      occupied_buckets = length(sizes),
      mean_bucket_size = mean(sizes),
      collision_buckets = sum(sizes >= 2L),
      singleton_buckets = sum(sizes == 1L),
      largest_bucket = max(sizes),
      pair_occurrences = sum(sizes * (as.double(sizes) - 1) / 2)
    )
    counts <- table(sizes)
    distributions[[band]] <- data.frame(
      rows_per_band = r, bands = b, band = band,
      bucket_size = as.integer(names(counts)), bucket_count = as.integer(counts)
    )
    collision_buckets <- members[sizes >= 2L]
    band_pairs[[band]] <- unlist(lapply(collision_buckets, function(ids) {
      pairs <- combn(ids, 2L)
      pair_codes(pairs[1L, ], pairs[2L, ], n)
    }), use.names = FALSE)
  }

  per_band <- do.call(rbind, per_band)
  distribution <- do.call(rbind, distributions)
  codes <- sort(unique(as.double(unlist(band_pairs, use.names = FALSE))))
  summary <- data.frame(
    signature_rows = m, rows_per_band = r, bands = b, documents = n,
    occupied_buckets = sum(per_band$occupied_buckets),
    assignments = b * n,
    mean_bucket_size = b * n / sum(per_band$occupied_buckets),
    mean_of_band_means = mean(per_band$mean_bucket_size),
    largest_bucket = max(per_band$largest_bucket),
    singleton_fraction = sum(per_band$singleton_buckets) / sum(per_band$occupied_buckets),
    pair_occurrences = sum(per_band$pair_occurrences),
    candidate_pairs = length(codes), total_pairs = choose(n, 2),
    candidate_fraction = if (n > 1L) length(codes) / choose(n, 2) else 0
  )
  list(buckets = buckets, per_band = per_band, distribution = distribution,
       candidates = decode_pairs(codes, n), summary = summary)
}

candidate_probability <- function(s, rows_per_band, signature_rows = 360L) {
  stopifnot(all(s >= 0 & s <= 1), rows_per_band >= 1L,
            signature_rows %% rows_per_band == 0L)
  b <- signature_rows / rows_per_band
  # Stable evaluation of 1 - (1 - s^r)^b, including the endpoints s = 0 and 1.
  -expm1(b * log1p(-s^rows_per_band))
}


# Task 5. Calculate exact Jaccard only after LSH selects candidates.
score_pairs <- function(pairs, shingle_ids, signature = NULL, verbose = FALSE) {
  scores <- estimates <- numeric(nrow(pairs))
  for (p in seq_len(nrow(pairs))) {
    i <- pairs$i[p]
    j <- pairs$j[p]
    scores[p] <- jaccard(shingle_ids[[i]], shingle_ids[[j]])
    if (!is.null(signature)) estimates[p] <- signature_similarity(signature, i, j)
    if (verbose && p %% 25000L == 0L) cat("  Scored", p, "candidate pairs.\n")
  }
  pairs$jaccard <- scores
  if (!is.null(signature)) pairs$minhash <- estimates
  pairs
}

# Sample distinct unordered pairs without constructing all choose(n, 2) pairs.
sample_document_pairs <- function(n, size = 1000L, seed = 123L) {
  stopifnot(n >= 2L, size >= 1L, size <= choose(n, 2))
  set.seed(seed)
  codes <- numeric()
  while (length(codes) < size) {
    i <- sample.int(n, 2L * size, replace = TRUE)
    j <- sample.int(n, 2L * size, replace = TRUE)
    codes <- unique(c(codes, pair_codes(i[i != j], j[i != j], n)))
  }
  decode_pairs(codes[seq_len(size)], n)
 }
