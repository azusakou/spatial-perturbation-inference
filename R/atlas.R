# Ordinal shape atlas. B1..B5 do not equate physical distances across assays.
run_atlas <- function(d,out) {
suppressPackageStartupMessages(library(cluster))
common_bins <- paste0("B", 1:5)
curves <- list()

# Spatial Perturb-Seq: combine the two sub-100-um bins using receiver counts.
sps <- copy(d$sps_profiles)
sps[, common_bin := fcase(distance_bin %in% c("0-50", "50-100"), "0-100",
                          default = as.character(distance_bin))]
sps <- sps[, .(value = weighted.mean(mean_score, receiver_cells), units = uniqueN(chip_id)),
           by = .(perturbation, common_bin)]
curves[["sps"]] <- sps[, .(dataset = "Spatial Perturb-Seq",
                            field_id = paste0("SPS:", perturbation),
                            condition = perturbation, programme = perturbation,
                            distance_bin = common_bin, value, units)]

pm <- copy(d$perturb_map_profiles)
pm <- pm[, .(value = mean(mean_score), units = uniqueN(lesion)),
         by = .(lesion_class, program, distance_bin)]
curves[["pm"]] <- pm[, .(dataset = "Perturb-map",
                          field_id = paste("PM", lesion_class, program, sep = ":"),
                          condition = lesion_class, programme = program,
                          distance_bin = as.character(distance_bin), value, units)]

mi <- copy(d$mi_profiles)
mi <- mi[, .(value = mean(mean_score), units = uniqueN(sample_id)),
         by = .(primary_zone, program, distance_bin)]
curves[["mi"]] <- mi[, .(dataset = "Myocardial infarction",
                          field_id = paste("MI", primary_zone, program, sep = ":"),
                          condition = as.character(primary_zone), programme = program,
                          distance_bin = as.character(distance_bin), value, units)]

skin <- copy(d$skin_profiles)
skin <- skin[, .(value = mean(mean_score), units = uniqueN(sample_id)),
             by = .(time, program, distance_bin)]
curves[["skin"]] <- skin[, .(dataset = "Skin injury",
                              field_id = paste("Skin", time, program, sep = ":"),
                              condition = as.character(time), programme = program,
                              distance_bin = as.character(distance_bin), value, units)]

hcc <- copy(d$hcc_internal_profiles)
hcc <- hcc[, .(value = mean(mean_score), units = uniqueN(patient_id)),
           by = .(Response, program, distance_bin)]
hcc[, distance_bin := c(`0-1` = "0-100", `1-2` = "100-200", `2-4` = "200-400",
                         `4-8` = "400-800", `>8` = ">800")[as.character(distance_bin)]]
curves[["hcc"]] <- hcc[, .(dataset = "Paired HCC",
                            field_id = paste("HCC", Response, program, sep = ":"),
                            condition = Response, programme = program,
                            distance_bin, value, units)]

ext <- copy(d$hcc_external_profiles)
ext <- ext[, .(value = mean(mean_score), units = uniqueN(sample_id)),
           by = .(Response, program, distance_bin)]
ext[, distance_bin := c(`0-1` = "0-100", `1-2` = "100-200", `2-4` = "200-400",
                         `4-8` = "400-800", `>8` = ">800")[as.character(distance_bin)]]
curves[["ext"]] <- ext[, .(dataset = "External HCC",
                            field_id = paste("ExternalHCC", Response, program, sep = ":"),
                            condition = Response, programme = program,
                            distance_bin, value, units)]

long <- rbindlist(curves, fill = TRUE)
# Cross-platform comparison is ordinal: retain each native interval/unit and
# compare only nearest-to-farthest bin order, without inventing physical HCC um.
native_labels <- c("0-100", "100-200", "200-400", "400-800", ">800")
grid_labels <- c("0-1", "1-2", "2-4", "4-8", ">8")
is_hcc <- long$dataset %chin% c("Paired HCC", "External HCC")
long[, original_distance_bin := as.character(distance_bin)]
long[is_hcc, original_distance_bin := grid_labels[match(as.character(distance_bin), native_labels)]]
long[, original_distance_unit := fifelse(dataset %chin% c("Paired HCC", "External HCC"),
  "isotropic Visium hex distance / median nearest-neighbour hex spacing", "micrometres")]
long[, distance_bin := paste0("B", match(as.character(distance_bin), native_labels))]
long[, common_bin_interpretation := "ordinal order only; B1 nearest and B5 farthest; unequal native widths"]
long <- long[distance_bin %in% common_bins & is.finite(value)]
long[, distance_bin := factor(distance_bin, levels = common_bins)]
wide <- dcast(long, dataset + field_id + condition + programme + units ~ distance_bin,
              value.var = "value")
wide <- wide[complete.cases(wide[, ..common_bins])]
Y <- as.matrix(wide[, ..common_bins])
Y_center <- t(scale(t(Y)))
Y_center[!is.finite(Y_center)] <- 0
mid <- seq_len(5)
features <- cbind(Y_center,
  local_contrast = rowMeans(Y[, 1:2, drop = FALSE]) - rowMeans(Y[, 4:5, drop = FALSE]),
  monotonicity = apply(Y, 1, function(v) cor(mid, v, method = "spearman")),
  interface_peak = apply(Y, 1, function(v) max(v[2:4]) - max(v[c(1,5)])))
features[!is.finite(features)] <- 0
Z <- scale(features); Z[!is.finite(Z)] <- 0

set.seed(280319)
ks <- 2:min(6, nrow(Z) - 1)
sil <- rbindlist(lapply(ks, function(k) {
  fit <- kmeans(Z, centers = k, nstart = 100)
  sizes <- table(fit$cluster)
  data.table(k, mean_silhouette = mean(silhouette(fit$cluster, dist(Z))[, 3]),
             min_cluster_size = min(sizes))
}))
eligible <- sil[min_cluster_size >= 3]
best_k <- if (nrow(eligible)) eligible[which.max(mean_silhouette), k] else sil[which.max(mean_silhouette), k]
fit <- kmeans(Z, centers = best_k, nstart = 500)
wide[, cluster := fit$cluster]

centroid <- data.table(cluster = seq_len(best_k), fit$centers)
curve_cols <- common_bins
name_cluster <- function(k) {
  y <- colMeans(Y_center[fit$cluster == k, , drop = FALSE])
  delta <- mean(y[1:2]) - mean(y[4:5])
  peak <- which.max(y)
  if (abs(delta) < 0.35) return("flat/heterogeneous")
  if (delta < 0) return("local-suppressed")
  if (peak %in% 2:4 && y[[1]] < max(y[2:4]) - 0.25) return("interface-peaked")
  if (y[[5]] > 0) return("long-range sustained")
  if (y[[1]] - y[[2]] > 0.7) return("local-contained")
  "radial-decaying"
}
labels <- setNames(vapply(seq_len(best_k), name_cluster, character(1)), seq_len(best_k))
wide[, archetype := labels[as.character(cluster)]]

# Distance margin is an item-level stability diagnostic; it is not treated as
# a bootstrap confidence interval.
d_to_centroid <- sapply(seq_len(best_k), function(k) rowSums((Z - fit$centers[k, ])^2))
sorted_d <- t(apply(d_to_centroid, 1, sort))
wide[, assignment_margin := sorted_d[, 2] - sorted_d[, 1]]
wide[, evidence_status := fcase(
  grepl("^SPS:Clu$", field_id), "not supported by strict null",
  dataset == "Perturb-map" & condition == "Tgfbr2" & programme == "TGF-beta/stromal", "conditional lesion profile",
  dataset == "Perturb-map" & condition == "Tgfbr2" & programme == "T-cell", "directional trend",
  dataset %in% c("Myocardial infarction", "Skin injury", "Paired HCC", "External HCC"), "proxy/exploratory",
  default = "screening"
)]

atlas <- merge(long, wide[, .(field_id, cluster, archetype, assignment_margin, evidence_status)],
               by = "field_id")
fwrite(wide, file.path(out,"spf_atlas_field_embeddings.tsv"), sep = "\t")
fwrite(atlas, file.path(out,"spf_atlas_distance_curves.tsv"), sep = "\t")
fwrite(sil, file.path(out,"spf_archetype_model_selection.tsv"), sep = "\t")


invisible(wide)
}
