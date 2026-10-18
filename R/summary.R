# Fast endpoint-level statistics. Frozen scores are inputs, not raw count data.
wilson <- function(k,n) {
 z <- qnorm(.975); p <- k/n; den <- 1+z*z/n
 mid <- (p+z*z/(2*n))/den; rad <- z*sqrt(p*(1-p)/n+z*z/(4*n*n))/den
 c(max(0,mid-rad),min(1,mid+rad))
}
run_summary <- function(d,out) {
 z <- copy(d$synthetic_replicates)
 syn <- z[, {
  ok<-is.finite(p); n<-sum(ok); k<-sum(p[ok]<=.05); ci<-wilson(k,n)
  list(n_replicates=n,rejections=k,rejection_rate=k/n,mc_lower=ci[1],mc_upper=ci[2],
       mean_effect=mean(null_adjusted_effect[ok]),
       mean_injected_near_far=mean(injected_near_far[ok]),
       min_valid_permutations=min(valid_permutations))
 },by=.(scenario,amplitude,adjustment,null_model)]
 fwrite(syn,file.path(out,'synthetic_summary.tsv'),sep='\t')
 z <- copy(d$spatial_screening)
 screen <- z[is.finite(directional_z),.(chips=.N,
   Z=sum((near_far_effect-null_mean)/null_sd)/sqrt(.N),
   mean_null_adjusted_effect=mean(null_adjusted_effect)),by=perturbation]
 screen[,`:=`(p=2*pnorm(-abs(Z)),scope='chip-block screening approximation; not a four-mouse population test')]
 screen[,fdr:=p.adjust(p,'BH')]
 fwrite(screen,file.path(out,'spatial_screening.tsv'),sep='\t')
 for(kind in c('mi','skin')) {
  z <- copy(if(kind=='mi')d$mi_embedding else d$skin_effects)
  if(kind=='skin')z[,is_control:=time=='Unwounded']
  zz <- z[, {
   stopifnot(any(is_control),any(!is_control))
   obs <- mean(near_far_effect[!is_control])-mean(near_far_effect[is_control])
   nul <- combn(.N,sum(is_control),function(ii)mean(near_far_effect[-ii])-mean(near_far_effect[ii]))
   list(effect=obs,exact_p_two_sided=mean(abs(nul)>=abs(obs)-1e-12),exact_assignments=length(nul))
  },by=program]
  fwrite(zz,file.path(out,paste0(kind,'_exact.tsv')),sep='\t')
 }
 # This scalar holdout diagnosis does not re-estimate programme weights.
 fwrite(as.data.table(d$heldout_metrics),file.path(out,'reference_weight_holdout_frozen_metrics.tsv'),sep='\t')
 fwrite(as.data.table(d$native_expected),file.path(out,'native_frozen_summary.tsv'),sep='\t')
 invisible(list(synthetic=syn,screening=screen))
}

verify_results <- function(d,out) {
 checks<-list()
 compare<-function(label,file,expected,keys,columns,tolerance=1e-8) {
  actual<-fread(file.path(out,file));expected<-as.data.table(expected)
  k<-function(x)do.call(paste,c(lapply(x[,..keys],as.character),sep='|'))
  m<-match(k(expected),k(actual));stopifnot(!anyNA(m),nrow(actual)==nrow(expected))
  for(nm in columns) {
   x<-as.numeric(actual[[nm]][m]);y<-as.numeric(expected[[nm]])
   same_na<-identical(is.na(x),is.na(y));v<-abs(x-y);err<-if(any(is.finite(v)))max(v[is.finite(v)])else 0
   checks[[length(checks)+1L]]<<-data.table(check=paste(label,nm),max_absolute_difference=err,
     tolerance=tolerance,pass=same_na&&is.finite(err)&&err<=tolerance)
  }
 }
 compare('strict Clu','clu_effect.tsv',d$clu_strict_expected,c('chip_id','receiver_identity'),
  c('observed_effect','null_mean','null_sd','null_adjusted_effect','permutation_p','fdr_within_chip'),1e-10)
 compare('synthetic summary','synthetic_summary.tsv',d$synthetic_expected,
  c('scenario','amplitude','adjustment','null_model'),c('n_replicates','rejections','rejection_rate','mc_lower','mc_upper'),1e-12)
 compare('internal clinical','hcc_nested_loo_predictions.tsv',d$hcc_internal_predictions_expected,
  c('patient_id','model'),c('outcome','probability','selected_lambda'),1e-8)
 compare('external clinical','hcc_external_locked_predictions.tsv',d$hcc_external_predictions_expected,
  c('sample_id','model'),c('outcome','probability','lambda'),1e-8)
 compare('ordinal profile','spf_leaveout_validation_metrics.tsv',d$profiles_expected,
  c('validation','field_id'),c('response_field_reconstruction_error','ordinal_energy_distance',
    'active_bin_agreement','quarter_peak_bin_error','radial_rank_concordance','near_bin_score_error'),1e-10)
 compare('reference weight holdout','reference/heldout_Clu_metrics.tsv',d$heldout_metrics,
  c('heldout_mouse','programme'),c('n_programme_genes','auc','mean_score_difference'),1e-10)
 compare('DBiT ATG4A','dbit/atg4a_cross_section_field.tsv',d$dbit_atg4a,'sample',
  c('direct_effect','receiver_field_effect','null_mean','null_sd','empirical_p_two_sided'),1e-10)
 compare('DBiT Tgfbr2','dbit/tgfbr2_section_fields.tsv',d$dbit_tgfbr2,'sample',
  c('effect','null_mean','null_sd','p_two_sided'),1e-10)
 z<-rbindlist(checks);fwrite(z,file.path(out,'verification.tsv'),sep='\t')
 if(!all(z$pass))stop('Reproduction mismatch; see verification.tsv')
 message('PASS: ',nrow(z),' numerical comparisons against archived/current audit snapshots.')
 invisible(z)
}
