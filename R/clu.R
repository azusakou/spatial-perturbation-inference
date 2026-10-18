# Full-chip strict Clu field analysis from frozen gene-score/identity inputs.
# The background is refit once per observed chip and then frozen for all relocations.
suppressPackageStartupMessages({library(data.table); library(RANN); library(igraph); library(mgcv)})
component_grid <- function(x, y, bin = 500) {
  bx <- floor(x / bin); by <- floor(y / bin); key <- paste(bx, by, sep = ":")
  occupied <- unique(data.table(key, bx, by))
  lookup <- setNames(seq_len(nrow(occupied)), occupied$key)
  edges <- rbindlist(lapply(seq_len(nrow(occupied)), function(i) {
    nb <- paste(occupied$bx[[i]] + c(-1,-1,-1,0,0,1,1,1),
                occupied$by[[i]] + c(-1,0,1,-1,1,-1,0,1), sep = ":")
    j <- unname(lookup[nb]); j <- j[!is.na(j)]
    if (!length(j)) return(NULL)
    data.table(from = i, to = j)
  }))
  g <- make_empty_graph(n = nrow(occupied), directed = FALSE)
  if (nrow(edges)) g <- add_edges(g, as.vector(t(as.matrix(edges))))
  components(g)$membership[lookup[key]]
}

nearest_component <- function(coords, source, component, direction = FALSE) {
  d <- rep(Inf, nrow(coords)); angle <- rep(NA_real_, nrow(coords))
  for (cc in intersect(unique(component), unique(component[source]))) {
    query <- which(component == cc)
    src <- source[component[source] == cc]
    nn <- nn2(coords[src, , drop = FALSE], coords[query, , drop = FALSE], k = 1)
    d[query] <- nn$nn.dists[, 1]
    if (direction) {
      si <- src[nn$nn.idx[, 1]]
      angle[query] <- atan2(coords[query, 2] - coords[si, 2],
                            coords[query, 1] - coords[si, 1])
    }
  }
  list(distance = d, angle = angle)
}

source_cluster_spec <- function(coords, source, component, threshold_um = 100) {
  out <- list()
  for (cc in unique(component[source])) {
    s <- source[component[source] == cc]
    if (length(s) == 1L) {
      sizes <- 1L
    } else {
      dd <- as.matrix(dist(coords[s, , drop = FALSE]))
      edge <- which(dd <= threshold_um & dd > 0, arr.ind = TRUE)
      g <- make_empty_graph(n = length(s), directed = FALSE)
      if (nrow(edge)) g <- add_edges(g, as.vector(t(edge)))
      sizes <- as.integer(table(components(g)$membership))
    }
    out[[as.character(cc)]] <- sizes
  }
  out
}

sample_cluster_preserving_sources <- function(coords, neg, component, spec) {
  chosen <- integer()
  for (cc_name in names(spec)) {
    cc <- as.integer(cc_name)
    pool0 <- neg[component[neg] == cc]
    used <- integer()
    for (size in spec[[cc_name]][sample.int(length(spec[[cc_name]]))]) {
      pool <- setdiff(pool0, used)
      if (length(pool) < size) return(integer())
      seed <- pool[sample.int(length(pool), 1L)]
      nn <- nn2(coords[pool, , drop = FALSE], coords[seed, , drop = FALSE],
                k = min(size, length(pool)))$nn.idx[1, ]
      pick <- pool[nn]
      chosen <- c(chosen, pick); used <- c(used, pick)
    }
  }
  chosen
}

read_sparse_genes <- function(path, genes_wanted, cell, gene) {
  genes <- intersect(genes_wanted, gene$geneName)
  parts <- lapply(seq_along(genes), function(i) {
    g <- gene[geneName == genes[[i]]][1]
    z <- as.data.table(h5read(path, "/cellBin/geneExp",
      start = as.numeric(g$offset) + 1, count = as.numeric(g$cellCount)))
    j <- match(z$cellID, cell$id); keep <- !is.na(j) & z$count > 0
    data.table(i = i, j = j[keep], x = as.numeric(z$count[keep]))
  })
  trip <- rbindlist(parts)
  sparseMatrix(i = trip$i, j = trip$j, x = trip$x,
               dims = c(length(genes), nrow(cell)), dimnames = list(genes, NULL))
}

anisotropy <- function(score, distance, angle) {
  keep <- is.finite(score) & is.finite(distance) & distance <= 800 & is.finite(angle)
  if (sum(keep) < 100) return(c(amplitude = NA_real_, preferred_angle = NA_real_))
  d_bin <- findInterval(distance[keep], c(0, 50, 100, 200, 400, 800), all.inside = TRUE)
  y <- score[keep]
  y <- y - ave(y, d_bin, FUN = mean)
  X <- cbind(cos(angle[keep]), sin(angle[keep]))
  beta <- tryCatch(qr.solve(crossprod(X), crossprod(X, y)), error = function(e) c(NA, NA))
  c(amplitude = sqrt(sum(beta^2)), preferred_angle = atan2(beta[[2]], beta[[1]]))
}

effect_by_identity <- function(score, distance, identity, ids, scale=1) {
  vapply(c("All", ids), function(id) {
    keep <- if (id == "All") rep(TRUE, length(score)) else identity == id
    near <- keep & distance <= 100*scale
    far <- keep & distance > 400*scale & distance <= 800*scale
    if (sum(near) < 30 || sum(far) < 30) return(NA_real_)
    mean(score[near], na.rm = TRUE) - mean(score[far], na.rm = TRUE)
  }, numeric(1))
}


# A finite reference fraction is unavailable without a finite observation and
# at least one finite reference statistic. No new minimum-bank threshold is imposed.
clu_reference_fraction <- function(observed, reference, absolute=TRUE) {
  reference <- reference[is.finite(reference)]
  if(length(observed)!=1L || !is.finite(observed) || !length(reference))return(NA_real_)
  if(absolute) {
    centre <- mean(c(observed,reference))
    return((1+sum(abs(reference-centre)>=abs(observed-centre)))/(length(reference)+1))
  }
  (1+sum(reference>=observed))/(length(reference)+1)
}

strict_clu_chip <- function(z, n_perm=99L, source_bank=NULL) {
  chip <- as.character(unique(z$chip_id)); stopifnot(length(chip)==1L)
  coords <- as.matrix(z[, .(x_um,y_um)])
  cell <- z; source <- which(z$guide_call=="sgrna-clu"); neg <- which(z$guide_call=="Neg")
  raw_score <- z$raw_clu_score; identity <- as.character(z$identity)
  lib <- z$library_size
  bg <- data.table(score = raw_score[neg], identity = factor(identity[neg]),
                   log_lib = log1p(lib[neg]), x_um = coords[neg, 1], y_um = coords[neg, 2])
  bg_fit <- bam(score ~ identity + s(log_lib, k = 5) + s(x_um, y_um, k = 30),
                data = bg, method = "fREML", discrete = TRUE, nthreads = 2)
  adjusted <- residuals(bg_fit)

  near <- nearest_component(coords, source, cell$component, direction = TRUE)
  ids_eval <- sort(unique(identity[neg]))
  observed_effect <- effect_by_identity(adjusted, near$distance[neg], identity[neg], ids_eval)
  scales <- c(.75,1,1.25)
  window_observed <- sapply(scales,function(s)effect_by_identity(adjusted,near$distance[neg],identity[neg],ids_eval,s))
  observed_aniso <- anisotropy(adjusted, near$distance[neg], near$angle[neg])
  bins <- cut(near$distance[neg], c(-Inf, 50, 100, 200, 400, 800, Inf),
              labels = c("0-50", "50-100", "100-200", "200-400", "400-800", ">800"))
  observed_profile <- data.table(distance_bin = bins, score = adjusted)[,
    .(observed = mean(score), receiver_cells = .N), by = distance_bin]

  spec <- source_cluster_spec(coords, source, cell$component)
  null_effect <- matrix(NA_real_, nrow = length(observed_effect), ncol = n_perm,
                        dimnames = list(names(observed_effect), NULL))
  null_aniso <- rep(NA_real_, n_perm)
  null_profile <- matrix(NA_real_, nrow = 6, ncol = n_perm,
                         dimnames = list(levels(bins), NULL))
  mask_audit <- vector("list", n_perm)
  window_null <- array(NA_real_,c(length(observed_effect),n_perm,length(scales)))
  if(!is.null(source_bank))stopifnot(length(source_bank)==n_perm)
  set.seed(24000 + sum(utf8ToInt(chip)))
  for (b in seq_len(n_perm)) {
    pseudo <- if(is.null(source_bank))sample_cluster_preserving_sources(coords, neg, cell$component, spec)else source_bank[[b]]
    stopifnot(all(pseudo%in%neg),!anyDuplicated(pseudo))
    if (length(pseudo) != length(source)) next
    pd <- nearest_component(coords, pseudo, cell$component, direction = TRUE)
    # Observed guide sources cannot be receivers. Apply the same exclusion to
    # the pseudo-source configuration, while retaining the frozen GAM fit.
    use <- !neg %in% pseudo
    stopifnot(!any(neg[use] %in% pseudo))
    null_effect[, b] <- effect_by_identity(adjusted[use], pd$distance[neg][use], identity[neg][use], ids_eval)
    for(w in seq_along(scales))window_null[,b,w] <- if(scales[w]==1)null_effect[,b]else
      effect_by_identity(adjusted[use],pd$distance[neg][use],identity[neg][use],ids_eval,scales[w])
    null_aniso[[b]] <- anisotropy(adjusted[use], pd$distance[neg][use], pd$angle[neg][use])[["amplitude"]]
    pb <- cut(pd$distance[neg][use], c(-Inf, 50, 100, 200, 400, 800, Inf), labels = levels(bins))
    tmp <- data.table(bin = pb, score = adjusted[use])[, mean(score), by = bin]
    null_profile[match(as.character(tmp$bin), rownames(null_profile)), b] <- tmp$V1
    mask_audit[[b]] <- data.table(chip_id=chip,draw=b,base_receivers=length(neg),
      pseudo_source_spots=length(pseudo),excluded_pseudo_sources=sum(!use),effective_receivers=sum(use),source_receiver_overlap_after=0L,
      near_receivers=sum(pd$distance[neg][use]<=100),far_receivers=sum(pd$distance[neg][use]>400&pd$distance[neg][use]<=800),
      pseudo_indices=paste(sort(pseudo),collapse=","))
  }

  effect_out <- data.table(receiver_identity = names(observed_effect),
                           observed_effect = observed_effect,
                           null_mean = rowMeans(null_effect, na.rm = TRUE),
                           null_sd = apply(null_effect, 1, sd, na.rm = TRUE))
  effect_out[, valid_nulls:=rowSums(is.finite(null_effect))]
  effect_out[, null_adjusted_effect := observed_effect - null_mean]
  effect_out[, permutation_p := vapply(seq_len(.N), function(i) {
    z <- null_effect[i, is.finite(null_effect[i, ])]
    clu_reference_fraction(observed_effect[[i]],z)
  }, numeric(1))]
  effect_out[, `:=`(chip_id = chip, source_cells = length(source))]

  profile_out <- merge(observed_profile,
    data.table(distance_bin = rownames(null_profile),
               null_mean = rowMeans(null_profile, na.rm = TRUE),
               null_sd = apply(null_profile, 1, sd, na.rm = TRUE)), by = "distance_bin")
  profile_out[, excess := observed - null_mean]
  finite_excess <- profile_out$excess[is.finite(profile_out$excess)]
  peak <- if (length(finite_excess)) max(finite_excess) else NA_real_
  profile_out[, above_quarter_peak := is.finite(excess) & is.finite(peak) &
                peak > sqrt(.Machine$double.eps) & excess >= 0.25 * peak]
  radius_right_censored <- any(profile_out$above_quarter_peak &
                                 as.character(profile_out$distance_bin) == ">800")
  radius <- if (radius_right_censored) {
    NA_real_
  } else if (any(profile_out$above_quarter_peak)) {
    upper <- c(`0-50` = 50, `50-100` = 100, `100-200` = 200,
               `200-400` = 400, `400-800` = 800)
    max(upper[as.character(profile_out$distance_bin[profile_out$above_quarter_peak])])
  } else NA_real_
  profile_out[, `:=`(chip_id = chip, estimated_radius_um = radius,
                     radius_right_censored = radius_right_censored)]

  geometry <- data.table(chip_id = chip, source_cells = length(source),
    source_clusters = sum(lengths(spec)), largest_source_cluster = max(unlist(spec)),
    anisotropy_amplitude = observed_aniso[["amplitude"]],
    preferred_angle_radians = observed_aniso[["preferred_angle"]],
    anisotropy_null_mean = mean(null_aniso, na.rm = TRUE),
    anisotropy_p = clu_reference_fraction(observed_aniso[["amplitude"]],null_aniso,absolute=FALSE),
    estimated_radius_um = radius,
    radius_right_censored = radius_right_censored,
    direct_source_score = mean(raw_score[source]))

  reference <- data.table(receiver_identity=rep(rownames(null_effect),times=n_perm),
    draw=rep(seq_len(n_perm),each=nrow(null_effect)),statistic=as.vector(null_effect),chip_id=chip)
  profile_reference <- data.table(distance_bin=rep(rownames(null_profile),times=n_perm),
    draw=rep(seq_len(n_perm),each=nrow(null_profile)),statistic=as.vector(null_profile),chip_id=chip)
  window_effect <- rbindlist(lapply(seq_along(scales),function(w) {
   rbindlist(lapply(seq_along(observed_effect),function(i) {
    v<-window_null[i,,w];v<-v[is.finite(v)];obs<-window_observed[i,w];mu<-mean(c(obs,v))
    keep<-if(names(observed_effect)[i]=="All")rep(TRUE,length(neg))else identity[neg]==names(observed_effect)[i]
    nn<-sum(keep&near$distance[neg]<=100*scales[w]);fn<-sum(keep&near$distance[neg]>400*scales[w]&near$distance[neg]<=800*scales[w])
    data.table(chip_id=chip,receiver_identity=names(observed_effect)[i],scale=scales[w],
      near_um=100*scales[w],far_start_um=400*scales[w],far_end_um=800*scales[w],
      observed_near_receivers=nn,observed_far_receivers=fn,
      estimable=is.finite(obs)&&length(v)>0L,
      reason=if(nn<30L||fn<30L)"fewer than 30 observed near/far receivers"else if(!length(v))"no finite reference"else"conditional sensitivity only",
      observed_effect=obs,null_mean=if(length(v))mean(v)else NA_real_,null_sd=if(length(v)>1L)sd(v)else NA_real_,
      null_adjusted_effect=if(length(v))obs-mean(v)else NA_real_,valid_nulls=length(v),
      permutation_p=clu_reference_fraction(obs,v))
   }))
  }))
  list(effect=effect_out, profile=profile_out, geometry=geometry,reference=reference,
       profile_reference=profile_reference,receiver_mask_audit=rbindlist(mask_audit),window_effect=window_effect)
}

run_strict_clu <- function(cells, out, B=99L) {
  stopifnot(all(c("chip_id","component","library_size","raw_clu_score") %in% names(cells)))
  result <- lapply(unique(cells$chip_id), function(chip) strict_clu_chip(cells[chip_id==chip], B))
  names(result) <- unique(cells$chip_id)
  for(kind in c("effect","profile","geometry","reference","profile_reference","receiver_mask_audit","window_effect")) {
    z <- rbindlist(lapply(result, `[[`, kind), fill=TRUE)
    if(kind=="effect") z[, fdr_within_chip:=p.adjust(permutation_p,"BH",n=.N), by=chip_id]
    if(kind=="window_effect") {
      z[,fdr_within_chip_window:=p.adjust(permutation_p,"BH",n=.N),by=.(chip_id,scale)]
      z[,fdr_all_chip_windows:=p.adjust(permutation_p,"BH",n=nrow(z))]
    }
    fwrite(z,file.path(out,paste0("clu_",kind,".tsv")),sep="\t")
  }
  invisible(result)
}
