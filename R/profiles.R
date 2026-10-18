# Descriptive ordinal profile holdouts. Entire Perturb-map animals are excluded.
run_profiles <- function(d,out) {
bins <- paste0("B", 1:5)
native_bins <- c("0-100", "100-200", "200-400", "400-800", ">800")
mid <- seq_len(5)

energy_distance <- function(a, b) {
  pa <- abs(a - tail(a, 1)); pb <- abs(b - tail(b, 1))
  if (sum(pa) < 1e-12) pa <- rep(1, length(pa))
  if (sum(pb) < 1e-12) pb <- rep(1, length(pb))
  pa <- pa / sum(pa); pb <- pb / sum(pb)
  D <- abs(outer(mid, mid, "-"))
  sqrt(max(0, 2 * sum(outer(pa, pb) * D) -
             sum(outer(pa, pa) * D) - sum(outer(pb, pb) * D)))
}

radius <- function(y) {
  z <- abs(y - tail(y, 1)); peak <- max(z)
  if (!is.finite(peak) || peak < 1e-10) return(NA_real_)
  hit <- which(z <= 0.25 * peak & seq_along(z) > which.max(z))
  if (length(hit)) mid[min(hit)] else max(mid)
}

metrics <- function(obs, pred) {
  obs <- as.numeric(obs); pred <- as.numeric(pred)
  ao <- abs(obs - tail(obs, 1)); ap <- abs(pred - tail(pred, 1))
  to <- 0.25 * max(ao); tp <- 0.25 * max(ap)
  active_o <- ao >= to & max(ao) > 1e-10
  active_p <- ap >= tp & max(ap) > 1e-10
  data.table(
    response_field_reconstruction_error = sqrt(mean((obs - pred)^2)),
    ordinal_energy_distance = energy_distance(obs, pred),
    active_bin_agreement = mean(active_o == active_p),
    quarter_peak_bin_error = abs(radius(obs) - radius(pred)),
    radial_rank_concordance = suppressWarnings(cor(obs, pred, method = "spearman")),
    near_bin_score_error = abs(mean(obs[1:2]) - mean(pred[1:2]))
  )
}

make_wide <- function(dt, ids) {
  x <- dcast(dt, paste(paste(ids, collapse = " + "), "~ distance_bin"),
             value.var = "value", fun.aggregate = mean)
  x[complete.cases(x[, ..bins])]
}

evaluate <- function(w, scheme, hold_col, strata, scope) {
  ans <- vector("list", nrow(w))
  for (i in seq_len(nrow(w))) {
    test <- w[i]
    train <- w[get(hold_col) != test[[hold_col]]]
    matched <- train
    for (s in strata) matched <- matched[get(s) == test[[s]]]
    fallback_used <- !nrow(matched)
    if (fallback_used) matched <- train
    pred <- colMeans(as.matrix(matched[, ..bins]), na.rm = TRUE)
    obs <- as.numeric(test[, ..bins])
    ans[[i]] <- cbind(data.table(validation = scheme, scope = scope,
                                  held_out = as.character(test[[hold_col]]),
                                  field_id = as.character(test[["field_id"]]),
                                  training_fields = nrow(matched),
                                  fallback_used = fallback_used,
                                  held_out_animal = if ("mouse_id" %in% names(test)) as.character(test[["mouse_id"]]) else NA_character_,
                                  training_animals = if ("mouse_id" %in% names(matched)) uniqueN(matched$mouse_id) else NA_integer_,
                                  training_animal_ids = if ("mouse_id" %in% names(matched)) paste(sort(unique(matched$mouse_id)), collapse=";") else NA_character_), metrics(obs, pred))
  }
  rbindlist(ans, fill = TRUE)
}

# Leave-chip and leave-perturbation validation in Spatial Perturb-Seq.
sps <- copy(d$sps_profiles)
sps[, distance_bin := as.character(distance_bin)]
sps[, distance_bin := fcase(distance_bin %in% c("0-50", "50-100"), "0-100",
                            default = distance_bin)]
sps <- sps[, .(value = weighted.mean(mean_score, receiver_cells)),
           by = .(chip_id, perturbation, distance_bin)]
sps[, distance_bin := paste0("B", match(distance_bin, native_bins))]
sps[, field_id := paste(chip_id, perturbation, sep = ":")]
spsw <- make_wide(sps, c("field_id", "chip_id", "perturbation"))
v_chip <- evaluate(spsw, "leave-chip", "chip_id", "perturbation",
                   "Spatial Perturb-Seq; matched perturbation")
v_gene <- evaluate(spsw, "leave-gene/perturbation", "perturbation", "chip_id",
                   "Spatial Perturb-Seq; matched chip")

# Whole-animal and receiver-programme holdouts in Perturb-map.
pm <- copy(d$perturb_map_profiles)
pm[, `:=`(distance_bin = as.character(distance_bin), value = mean_score,
          field_id = paste(section_id, lesion, program, sep = ":"))]
pm[, distance_bin := paste0("B", match(distance_bin, native_bins))]
pm[, mouse_id := c(KP_1="mouse4", KP_2="mouse4", KP_3="mouse5", KP_4="mouse7")[as.character(section_id)]]
stopifnot(!anyNA(pm$mouse_id))
pmw <- make_wide(pm, c("field_id", "section_id", "lesion", "lesion_class", "program", "mouse_id"))
v_lesion <- evaluate(pmw, "leave-Perturb-map-animal", "mouse_id", c("lesion_class", "program"),
                     "Perturb-map; whole-animal exclusion, class/programme match when available")
v_niche <- evaluate(pmw, "leave-programme", "program", "lesion_class",
                    "Perturb-map; receiver-programme exclusion, not spatial-niche exclusion")

# Cross-organ/platform/dataset generalization on the common-bin atlas. These
# are intentionally hard zero-shot tests: the predictor never sees the held-out
# dataset and uses the training-atlas mean because programme labels do not align.
atlas <- copy(d$atlas)
# Group the actual measurement technology, rather than splitting studies that
# all use Visium. MI and skin are Visium in their source study/GEO descriptions.
platform_map <- c("Spatial Perturb-Seq"="Stereo-seq", "Perturb-map"="Visium",
  "Myocardial infarction"="Visium", "Skin injury"="Visium",
  "Paired HCC"="Visium", "External HCC"="Visium")
atlas[, platform := unname(platform_map[as.character(dataset)])]
stopifnot(!anyNA(atlas$platform), uniqueN(atlas$platform)==2L)
fwrite(atlas[, .(platform=unique(platform), annotation_curves=.N), by=dataset],
       file.path(out,"spf_platform_mapping.tsv"), sep="\t")
atlas[, organ := fcase(dataset == "Spatial Perturb-Seq", "brain",
                       dataset == "Perturb-map", "lung tumour",
                       dataset == "Myocardial infarction", "heart",
                       dataset == "Skin injury", "skin", default = "HCC patient tissues (liver and node)")]
v_dataset <- evaluate(atlas, "leave-dataset", "dataset", character(),
                      "Common-bin atlas; zero-shot dataset")
v_platform <- evaluate(atlas, "leave-platform", "platform", character(),
                       "Common-bin atlas; zero-shot platform")
v_organ <- evaluate(atlas, "leave-tissue-context", "organ", character(),
                    "Common-bin atlas; zero-shot tissue context")

val <- rbindlist(list(v_chip, v_gene, v_lesion, v_niche,
                      v_dataset, v_platform, v_organ), fill = TRUE)
unit_audit <- val[, .(test_annotation_curves=.N, held_out_groups=uniqueN(held_out),
  fallback_curves=sum(fallback_used, na.rm=TRUE),
  training_animal_count_min=if(all(is.na(training_animals))) NA_integer_ else min(training_animals,na.rm=TRUE),
  training_animal_count_max=if(all(is.na(training_animals))) NA_integer_ else max(training_animals,na.rm=TRUE)), by=validation]
animal_rows <- val[validation=="leave-Perturb-map-animal"]
stopifnot(nrow(animal_rows)>0L, all(vapply(seq_len(nrow(animal_rows)),function(i)
  !animal_rows$held_out_animal[i] %in% strsplit(animal_rows$training_animal_ids[i],";",fixed=TRUE)[[1]], logical(1))))
fwrite(unit_audit, file.path(out,"spf_reconstruction_unit_audit.tsv"), sep="\t")
fwrite(val, file.path(out,"spf_leaveout_validation_metrics.tsv"), sep = "\t")
summary <- val[, lapply(.SD, function(x) mean(x, na.rm = TRUE)),
               by = .(validation, scope),
               .SDcols = c("response_field_reconstruction_error", "ordinal_energy_distance",
                           "active_bin_agreement", "quarter_peak_bin_error",
                           "radial_rank_concordance", "near_bin_score_error")]
summary[, fields := val[, .N, by = validation]$N[match(validation, val[, unique(validation)])]]
fwrite(summary, file.path(out,"spf_leaveout_validation_summary.tsv"), sep = "\t")

# Receiver-identity holdout is scalar because identity-resolved radial curves
# were not estimable at sufficient counts; retain NAs rather than invent them.
id <- copy(d$clu_strict_expected)[receiver_identity != "All"]
id[, prediction := vapply(seq_len(.N), function(i) {
  train <- null_adjusted_effect[receiver_identity != receiver_identity[i]]
  if (any(is.finite(train))) mean(train, na.rm = TRUE) else NA_real_
}, numeric(1)), by = chip_id]
id[, training_count := vapply(seq_len(.N), function(i)
  sum(receiver_identity != receiver_identity[i] & is.finite(null_adjusted_effect)),
  integer(1)), by = chip_id]
idval <- id[, .(validation = "leave-cell-type", scope = "SPS Clu cluster-null scalar effect",
                held_out = receiver_identity, field_id = paste(chip_id, receiver_identity, sep = ":"),
                training_fields = training_count,
                response_field_reconstruction_error = abs(null_adjusted_effect - prediction),
                ordinal_energy_distance = NA_real_, active_bin_agreement = NA_real_,
                quarter_peak_bin_error = NA_real_, radial_rank_concordance = NA_real_,
                near_bin_score_error = abs(null_adjusted_effect - prediction))]
fwrite(idval, file.path(out,"spf_leave_cell_type_validation.tsv"), sep = "\t")

coverage <- data.table(
  requested_validation = c("leave-gene", "leave-Perturb-map-animal", "leave-cell-type", "leave-programme",
                           "leave-tissue-context", "leave-platform", "leave-dataset"),
  implementation = c("SPS perturbation holdout", "Perturb-map whole-animal exclusion; three reconstructed mouse groups",
                     "SPS broad receiver identity; scalar endpoint", "Perturb-map programme holdout",
                     "common-bin tissue context; HCC pools liver and node", "two actual technologies: Stereo-seq versus Visium", "common-bin zero-shot"),
  limitation = c("genes are guide-defined programmes", "annotation curves remain correlated within each of three held-out mouse groups; no new intervention experiment",
                 "radial metrics unavailable", "receiver-programme holdout only; spatial niche and lineage holdouts unavailable",
                 "descriptive profile model", "two platform groups; correlated curves are not independent technology replications", "descriptive profile model")
)
fwrite(coverage, file.path(out,"spf_validation_coverage.tsv"), sep = "\t")


invisible(summary)
}
