# Perturb-DBiT within-block section stress tests; no animal-level joint inference.
# ATG4A Hallmark membership is frozen as MSigDB 2026.1.Hs (CC BY 4.0).
run_dbit_atg4a <- function(input,out,B=500L) {
 suppressPackageStartupMessages({library(data.table);library(mgcv);library(RANN);library(igraph)})
 gene_sets <- setNames(vector("list",length(input$programme_names)),input$programme_names)
adjust_programmes <- function(dt) {
  receivers <- !dt$any_guide
  for (nm in names(gene_sets)) {
    kxy <- min(60L, max(10L, floor(sum(receivers) / 40)))
    fit <- gam(reformulate(c(sprintf("s(x,y,k=%d)", kxy), "s(log1p(lib_proxy),k=5)"), nm),
               data = dt[receivers], method = "REML")
    dt[[paste0(nm, "_adjusted")]] <- dt[[nm]] - as.numeric(predict(fit, newdata = dt))
  }
  dt
}

cluster_sizes <- function(dt) {
  src <- dt[atg4a == TRUE, .(x = as.integer(x / 50), y = as.integer(y / 50))]
  if (!nrow(src)) return(integer())
  key <- paste(src$x, src$y, sep = "x")
  edges <- rbindlist(lapply(seq_len(nrow(src)), function(i) {
    nb <- paste(src$x[[i]] + c(-1,-1,-1,0,0,1,1,1),
                src$y[[i]] + c(-1,0,1,-1,1,-1,0,1), sep = "x")
    j <- match(nb, key, nomatch = 0L)
    data.table(from = i, to = j[j > 0L])
  }))
  g <- graph_from_data_frame(edges[from != to], directed = FALSE,
                             vertices = data.frame(name = seq_len(nrow(src))))
  sort(as.integer(table(components(g)$membership)), decreasing = TRUE)
}

grow_pseudo_sources <- function(dt, sizes) {
  sample_one <- function(x) x[[sample.int(length(x), 1L)]]
  grid <- dt[, .(gx = as.integer(x / 50), gy = as.integer(y / 50))]
  key <- paste(grid$gx, grid$gy, sep = "x")
  used <- logical(nrow(grid))
  chosen_all <- integer()
  for (sz in sizes) {
    success <- FALSE
    for (attempt in seq_len(100L)) {
      seed <- sample_one(which(!used))
      chosen <- seed
      frontier <- integer()
      while (length(chosen) < sz) {
        cur <- sample_one(chosen)
        nb <- paste(grid$gx[[cur]] + c(-1,-1,-1,0,0,1,1,1),
                    grid$gy[[cur]] + c(-1,0,1,-1,1,-1,0,1), sep = "x")
        cand <- match(nb, key, nomatch = 0L)
        cand <- cand[cand > 0L]
        cand <- cand[!used[cand] & !cand %in% chosen]
        frontier <- unique(c(frontier, cand))
        if (!length(frontier)) break
        add <- sample_one(frontier)
        chosen <- c(chosen, add)
        frontier <- setdiff(frontier, add)
      }
      if (length(chosen) == sz) { success <- TRUE; break }
    }
    if (!success) stop("Could not place pseudo-source cluster")
    used[chosen] <- TRUE
    chosen_all <- c(chosen_all, chosen)
  }
  chosen_all
}

field_stat <- function(dt, score, source_idx) {
  base <- which(!dt$any_guide);recv <- setdiff(base,source_idx)
  stopifnot(!any(recv %in% source_idx))
  d <- nn2(as.matrix(dt[source_idx, .(x, y)]), as.matrix(dt[recv, .(x, y)]), k = 1)$nn.dists[, 1]
  v <- dt[[score]][recv]
  if (!any(d <= 100) || !any(d > 400 & d <= 800))
    return(c(effect = NA_real_, near_n = sum(d <= 100), far_n = sum(d > 400 & d <= 800),
      base_receivers=length(base),excluded_pseudo_sources=length(base)-length(recv),effective_receivers=length(recv)))
  c(effect = mean(v[d <= 100], na.rm = TRUE) -
             mean(v[d > 400 & d <= 800], na.rm = TRUE),
    near_n = sum(d <= 100), far_n = sum(d > 400 & d <= 800),
    base_receivers=length(base),excluded_pseudo_sources=length(base)-length(recv),effective_receivers=length(recv))
}

set.seed(340319)
s6 <- list(dt=copy(input$s6));s8 <- list(dt=copy(input$s8))
s6$dt <- adjust_programmes(s6$dt); s8$dt <- adjust_programmes(s8$dt)
coverage <- input$coverage
fwrite(coverage, file.path(out,"hallmark_gene_coverage.tsv"), sep = "\t")

# Programme selection uses only direct source contrast in Bru.LC.6.
training <- rbindlist(lapply(names(gene_sets), function(nm) {
  score <- paste0(nm, "_adjusted")
  src <- which(s6$dt$atg4a)
  recv <- which(!s6$dt$any_guide)
  d <- nn2(as.matrix(s6$dt[src, .(x, y)]), as.matrix(s6$dt[recv, .(x, y)]), k = 1)$nn.dists[, 1]
  data.table(programme = nm,
             direct_effect = mean(s6$dt[[score]][src]) -
               mean(s6$dt[[score]][recv][d > 400 & d <= 800]))
}))
selected <- training[which.max(abs(direct_effect))]
sign_frozen <- ifelse(selected$direct_effect >= 0, 1, -1)
score_name <- paste0(selected$programme, "_adjusted")
s6$dt[, frozen_score := sign_frozen * get(score_name)]
s8$dt[, frozen_score := sign_frozen * get(score_name)]
fwrite(training, file.path(out,"atg4a_training_programme_selection.tsv"), sep = "\t")

null_audit <- list()
evaluate <- function(obj, role) {
  dt <- obj$dt
  src <- which(dt$atg4a)
  sizes <- cluster_sizes(dt)
  observed <- field_stat(dt, "frozen_score", src)
  draws <- rbindlist(lapply(seq_len(B), function(b) {
    pseudo <- grow_pseudo_sources(dt, sizes)
    cbind(data.table(sample=unique(dt$sample),draw=b,pseudo_source_spots=length(pseudo),pseudo_indices=paste(sort(pseudo),collapse=","),
      source_receiver_overlap_after=0L),as.data.table(as.list(field_stat(dt,"frozen_score",pseudo))))
  }))
  null_audit[[as.character(unique(dt$sample))]] <<- draws
  null <- draws$effect[is.finite(draws$effect)]
  direct <- mean(dt$frozen_score[src]) - mean(dt$frozen_score[!dt$any_guide])
  data.table(sample = unique(dt$sample), role, selected_programme = selected$programme,
             frozen_sign = sign_frozen, source_spots = length(src),
             source_clusters = length(sizes), largest_cluster = max(sizes),
             direct_effect = direct, receiver_field_effect = observed[["effect"]],
             near_receivers = observed[["near_n"]], far_receivers = observed[["far_n"]],
             null_mean = mean(null), null_sd = sd(null),
             empirical_p_two_sided = if (length(null)) (1 + sum(abs(null - mean(c(observed[["effect"]], null))) >=
               abs(observed[["effect"]] - mean(c(observed[["effect"]], null))))) / (length(null) + 1) else NA_real_,
             null_permutations = length(null))
}
results <- rbind(evaluate(s6, "programme-selection section"),
                 evaluate(s8, "held-out section"))
fwrite(results, file.path(out,"atg4a_cross_section_field.tsv"), sep = "\t")
fwrite(rbindlist(null_audit),file.path(out,"atg4a_receiver_mask_audit.tsv.gz"),sep="\t",compress="gzip")

profile_one <- function(obj) {
  dt <- obj$dt; src <- which(dt$atg4a); recv <- which(!dt$any_guide)
  d <- nn2(as.matrix(dt[src, .(x, y)]), as.matrix(dt[recv, .(x, y)]), k = 1)$nn.dists[, 1]
  data.table(sample = unique(dt$sample), distance = d, score = dt$frozen_score[recv])[, .(
    mean = mean(score), se = sd(score) / sqrt(.N), n = .N),
    by = .(sample, distance_bin = cut(distance, c(-Inf,100,200,400,800,Inf),
      labels = c("0-100","100-200","200-400","400-800",">800")))]
}
profile <- rbind(profile_one(s6), profile_one(s8))
fwrite(profile, file.path(out,"perturb_dbit_atg4a_profiles.tsv"), sep = "\t")


invisible(results)
}

run_dbit_tgfbr2 <- function(input,out,B=1000L) {
 suppressPackageStartupMessages({library(data.table);library(mgcv);library(RANN);library(igraph)})
 make_section <- function(sample) {
  obj <- input[[if(sample=="mTSG.LV.2")"t2" else "t3"]]
  dt <- copy(obj$dt);recv <- !dt$any_guide
  kxy <- min(40L,max(10L,floor(sum(recv)/40)))
  fit <- gam(score~s(x,y,k=kxy)+s(log1p(lib_proxy),k=5),data=dt[recv],method="REML")
  dt[,adjusted:=score-as.numeric(predict(fit,newdata=dt))];obj$dt<-dt;obj
 }
cluster_sizes <- function(dt) {
  src <- dt[tgfbr2 == TRUE, .(gx = as.integer(x / 50), gy = as.integer(y / 50))]
  key <- paste(src$gx, src$gy, sep = "x")
  edges <- rbindlist(lapply(seq_len(nrow(src)), function(i) {
    nb <- paste(src$gx[[i]] + c(-1,-1,-1,0,0,1,1,1),
                src$gy[[i]] + c(-1,0,1,-1,1,-1,0,1), sep = "x")
    j <- match(nb, key, nomatch = 0L)
    data.table(from = i, to = j[j > 0L])
  }))
  g <- graph_from_data_frame(edges[from != to], directed = FALSE,
                             vertices = data.frame(name = seq_len(nrow(src))))
  sort(as.integer(table(components(g)$membership)), decreasing = TRUE)
}

grow <- function(dt, sizes) {
  sample_one <- function(x) x[[sample.int(length(x), 1L)]]
  grid <- dt[, .(gx = as.integer(x / 50), gy = as.integer(y / 50))]
  key <- paste(grid$gx, grid$gy, sep = "x"); used <- logical(nrow(grid)); out <- integer()
  for (sz in sizes) {
    for (attempt in seq_len(100L)) {
      chosen <- sample_one(which(!used)); frontier <- integer()
      while (length(chosen) < sz) {
        cur <- sample_one(chosen)
        nb <- paste(grid$gx[[cur]] + c(-1,-1,-1,0,0,1,1,1),
                    grid$gy[[cur]] + c(-1,0,1,-1,1,-1,0,1), sep = "x")
        cand <- match(nb, key, nomatch = 0L)
        cand <- cand[cand > 0L]
        cand <- cand[!used[cand] & !cand %in% chosen]
        frontier <- unique(c(frontier, cand))
        if (!length(frontier)) break
        add <- sample_one(frontier); chosen <- c(chosen, add); frontier <- setdiff(frontier, add)
      }
      if (length(chosen) == sz) break
    }
    if (length(chosen) != sz) stop("pseudo-source placement failed")
    used[chosen] <- TRUE; out <- c(out, chosen)
  }
  out
}

effect <- function(dt, src) {
  base <- which(!dt$any_guide);recv <- setdiff(base,src)
  stopifnot(!any(recv %in% src))
  d <- nn2(as.matrix(dt[src, .(x, y)]), as.matrix(dt[recv, .(x, y)]), k = 1)$nn.dists[, 1]
  v <- dt$adjusted[recv]
  e <- if (!any(d <= 100) || !any(d > 400 & d <= 800)) NA_real_ else mean(v[d <= 100]) - mean(v[d > 400 & d <= 800])
  c(effect=e,near_receivers=sum(d<=100),far_receivers=sum(d>400&d<=800),base_receivers=length(base),
    excluded_pseudo_sources=length(base)-length(recv),effective_receivers=length(recv))
}

set.seed(350319)
objs <- lapply(c("mTSG.LV.2", "mTSG.LV.3"), make_section)
section_results <- list(); nulls <- list(); profiles <- list();null_audit<-list()
for (i in seq_along(objs)) {
  dt <- objs[[i]]$dt; src <- which(dt$tgfbr2); sizes <- cluster_sizes(dt)
  obs <- effect(dt, src)[["effect"]]
  draws<-rbindlist(lapply(seq_len(B),function(b) {
   pseudo<-grow(dt,sizes)
   cbind(data.table(sample=unique(dt$sample),draw=b,pseudo_source_spots=length(pseudo),pseudo_indices=paste(sort(pseudo),collapse=","),
     source_receiver_overlap_after=0L),as.data.table(as.list(effect(dt,pseudo))))
  }))
  null_audit[[i]]<-draws;nul<-draws$effect[is.finite(draws$effect)]
  section_results[[i]] <- data.table(sample = unique(dt$sample), source_spots = length(src),
    source_clusters = length(sizes), largest_cluster = max(sizes), target_genes = length(objs[[i]]$observed_targets),
    effect = obs, null_mean = mean(nul), null_sd = sd(nul),
    p_two_sided = (1 + sum(abs(nul - mean(c(obs, nul))) >= abs(obs - mean(c(obs, nul))))) / (length(nul) + 1),
    valid_nulls = length(nul))
  nulls[[i]] <- nul
  recv <- which(!dt$any_guide)
  d <- nn2(as.matrix(dt[src, .(x, y)]), as.matrix(dt[recv, .(x, y)]), k = 1)$nn.dists[,1]
  profiles[[i]] <- data.table(sample = unique(dt$sample), distance = d, score = dt$adjusted[recv])[,
    .(mean = mean(score), se = sd(score)/sqrt(.N), n=.N), by=.(sample,
      distance_bin=cut(distance,c(-Inf,100,200,400,800,Inf),labels=c("0-100","100-200","200-400","400-800",">800")))]
}
section_results <- rbindlist(section_results); profile <- rbindlist(profiles)

section_results[,scope:="same published block; section-conditional source-relocation reference only"]
fwrite(section_results,file.path(out,"tgfbr2_section_fields.tsv"),sep="\t")
fwrite(profile,file.path(out,"perturb_dbit_tgfbr2_profiles.tsv"),sep="\t")
fwrite(rbindlist(lapply(seq_along(nulls),function(i)data.table(sample=section_results$sample[i],draw=seq_along(nulls[[i]]),statistic=nulls[[i]]))),file.path(out,"tgfbr2_section_nulls.tsv"),sep="\t")
fwrite(rbindlist(null_audit),file.path(out,"tgfbr2_receiver_mask_audit.tsv.gz"),sep="\t",compress="gzip")
invisible(section_results)
}

run_dbit <- function(input,out) {
 dir.create(out,recursive=TRUE,showWarnings=FALSE)
 a<-run_dbit_atg4a(input$atg4a,out)
 b<-run_dbit_tgfbr2(input$tgfbr2,out)
 invisible(list(atg4a=a,tgfbr2=b))
}
