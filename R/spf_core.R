# Transparent review benchmark utilities. Coordinates and all windows are in um.
# The cluster sampler preserves component-wise cluster sizes at placement, not
# exact shapes, pair distances, or the connected-component sizes after placement.
suppressPackageStartupMessages({library(RANN); library(igraph); library(mgcv)})

spf_distance <- function(receiver, candidate, source) {
  ans <- rep(Inf, nrow(receiver))
  for (cc in intersect(unique(receiver$component), unique(candidate$component[source]))) {
    ii <- which(receiver$component == cc)
    jj <- source[candidate$component[source] == cc]
    ans[ii] <- RANN::nn2(as.matrix(candidate[jj, c("x", "y")]),
                         as.matrix(receiver[ii, c("x", "y")]), k = 1)$nn.dists[, 1]
  }
  ans
}

spf_contrast <- function(distance, near = 100, far = c(400, 800), minimum = 5L) {
  ni <- which(is.finite(distance) & distance <= near)
  fi <- which(is.finite(distance) & distance > far[1] & distance <= far[2])
  if (length(ni) < minimum || length(fi) < minimum) return(rep(NA_real_, length(distance)))
  v <- numeric(length(distance)); v[ni] <- 1 / length(ni); v[fi] <- -1 / length(fi)
  v
}

spf_cluster_spec <- function(source, threshold = 100) {
  groups <- list()
  for (cc in unique(source$component)) {
    ii <- which(source$component == cc)
    dd <- as.matrix(dist(source[ii, c("x", "y")]))
    ed <- which(dd <= threshold & dd > 0, arr.ind = TRUE)
    g <- igraph::make_empty_graph(length(ii), directed = FALSE)
    if (nrow(ed)) g <- igraph::add_edges(g, as.vector(t(ed)))
    for (ids in split(ii, igraph::components(g)$membership)) {
      center <- colMeans(source[ids, c("x", "y")])
      anchor <- ids[which.min(rowSums(sweep(as.matrix(source[ids, c("x", "y")]), 2, center)^2))]
      groups[[length(groups) + 1L]] <- data.frame(component = cc, size = length(ids),
                                                 niche = source$niche[anchor])
    }
  }
  do.call(rbind, groups)
}

spf_sampler_cache <- function(candidate, spec) {
  # Nearest-neighbour order is cached once. Fall back to exact distances if all
  # cached neighbours are already used, so component/size matching is retained.
  k <- min(nrow(candidate), max(spec$size) + sum(spec$size) + 10L)
  cache <- list()
  for (cc in unique(candidate$component)) {
    ii <- which(candidate$component == cc)
    nn <- RANN::nn2(as.matrix(candidate[ii, c("x", "y")]), k = min(k, length(ii)))$nn.idx
    cache[[as.character(cc)]] <- list(ii = ii, nn = matrix(ii[nn], nrow = length(ii)))
  }
  cache
}

spf_sample_sources <- function(candidate, spec, method = "cluster_component", cache = NULL) {
  total <- sum(spec$size)
  if (method == "label_unrestricted") return(sample.int(nrow(candidate), total))
  if (method == "count_component") return(unlist(lapply(unique(spec$component), function(cc) {
    pool <- which(candidate$component == cc)
    pool[sample.int(length(pool), sum(spec$size[spec$component == cc]))]
  }), use.names = FALSE))
  if (is.null(cache)) cache <- spf_sampler_cache(candidate, spec)
  chosen <- integer()
  for (j in sample.int(nrow(spec))) {
    cc <- spec$component[j]; size <- spec$size[j]
    pool <- which(candidate$component == cc & !seq_len(nrow(candidate)) %in% chosen)
    anchors <- pool
    if (method == "cluster_anchor_niche") anchors <- pool[candidate$niche[pool] == spec$niche[j]]
    if (!length(anchors) || length(pool) < size) stop("Insufficient exchangeability pool")
    seed <- anchors[sample.int(length(anchors), 1)]
    z <- cache[[as.character(cc)]]
    order <- z$nn[match(seed, z$ii), ]
    order <- order[!order %in% chosen]
    if (length(order) < size) {
      d <- rowSums(sweep(as.matrix(candidate[pool, c("x", "y")]), 2,
                         as.numeric(candidate[seed, c("x", "y")]))^2)
      order <- pool[order(d)]
    }
    chosen <- c(chosen, head(order, size))
  }
  chosen
}

spf_adjust <- function(score, receiver, model = "spatial", return_fitted = FALSE) {
  z <- receiver; z$score <- score; z$identity <- factor(z$identity)
  if (!"log_lib" %in% names(z)) z$log_lib <- 0
  lib_ok <- length(unique(z$log_lib)) >= 6L
  if (model == "none") fitted <- rep(mean(score), length(score))
  if (model == "linear") {
    f <- if (lib_ok) score ~ identity + log_lib + x + y else score ~ identity + x + y
    fitted <- fitted(stats::lm(f, data = z))
  }
  if (model == "spatial") {
    f <- if (lib_ok) score ~ identity + s(log_lib, k = 5) + s(x, y, k = 30) else
      score ~ identity + s(x, y, k = 30)
    fitted <- fitted(mgcv::gam(f, data = z, method = "REML"))
  }
  if (return_fitted) fitted else score - fitted
}

spf_pvalue <- function(observed, null) {
  null <- null[is.finite(null)]
  if (length(null) < 20L || !is.finite(observed)) return(c(p = NA, effect = NA, valid = length(null)))
  # A common centre treats the observed and simulated statistics symmetrically.
  # This fixes null-only centring; source-placement exchangeability is still required.
  mu <- mean(c(observed, null))
  c(p = (1 + sum(abs(null - mu) >= abs(observed - mu))) / (length(null) + 1),
    effect = observed - mean(null), valid = length(null))
}

spf_null_bank <- function(receiver, candidate, spec, B = 199L,
                          near = 100, far = c(400, 800)) {
  cache <- spf_sampler_cache(candidate, spec)
  point_key <- function(z) paste(z$component,sprintf('%.17g',z$x),sprintf('%.17g',z$y),sep=':')
  receiver_key <- point_key(receiver); candidate_key <- point_key(candidate)
  methods <- c("label_unrestricted", "count_component", "cluster_component", "cluster_anchor_niche")
  setNames(lapply(methods, function(method) vapply(seq_len(B), function(b) {
    src <- spf_sample_sources(candidate, spec, method, cache)
    d <- spf_distance(receiver, candidate, src)
    # A selected pseudo-source is excluded from its own receiver statistic.
    d[receiver_key %in% candidate_key[src]] <- Inf
    spf_contrast(d, near, far)
  }, numeric(nrow(receiver)))), methods)
}

spf_wilson <- function(k, n, level = .95) {
  if (!n) return(c(lower = NA_real_, upper = NA_real_))
  z <- qnorm((1 + level) / 2); p <- k / n; den <- 1 + z^2 / n
  mid <- (p + z^2 / (2 * n)) / den
  rad <- z * sqrt(p * (1-p) / n + z^2 / (4*n^2)) / den
  c(lower = max(0, mid-rad), upper = min(1, mid+rad))
}
