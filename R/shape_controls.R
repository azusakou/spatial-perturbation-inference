# Conditional Perturb-map shape/location diagnostics; no animal-level P value.
# Rigidly translate/rotate the original physical point footprint on the native grid.
# All mapped source spots must be in tissue and outside every named lesion.
run_shape_controls <- function(spot_path,outdir) {
 suppressPackageStartupMessages({library(data.table);library(RANN)})
 dir.create(outdir,recursive=TRUE,showWarnings=FALSE)
 put<-function(x,n)fwrite(x,file.path(outdir,n),sep="\t",na="NA")
 z0<-fread(spot_path,na.strings=c("","NA"));stopifnot(nrow(z0)==6363L,
   all(c("array_row","array_col")%in%names(z0)),all(z0$eligible_receiver==is.na(z0$lesion_class)))
 defs<-list(primary_stromal=c("Serpine1","Ctgf","Col1a1","Col1a2","Col3a1","Fn1","Tnc","Timp1","Tgfb1","Spp1"),
   t_cell=c("Cd3d","Cd3e","Cd3g","Trac","Cd8a","Cd8b1","Lck","Lat"),
   tgfb_target=c("Serpine1","Ctgf","Smad7","Acta2","Tagln","Timp1","Col1a1","Col1a2","Col3a1","Fn1","Tnc"))
 # All thresholds are fixed before examining programme outcomes. Nested levels
 # are all retained; none is selected from the resulting programme contrast.
 levels<-data.table(level=c("strict","moderate","broad"),edge_tolerance_um=c(100,200,400),
   normal_fraction_tolerance=c(.10,.20,.30),source_normal_tolerance=c(.10,.20,.30),
   receiver_count_max_ratio=c(2,3,4),minimum_candidates=20L)
 put(levels,"matching_design.tsv")
 put(data.table(parameter=c("orientation_degrees","maximum_grid_centre_error","near_um","far_interval_um","minimum_stratum_near_and_far","shape_mask","receiver_rule","minimum_matched_placements","observed_KP_size_ratio_calipers","observed_KP_distance"),
   value=c("0;60;120;180;240;300","0.25 median nearest-neighbour spacing","200","(400,800]","5 each",
     "published annotated capture positions; not a histological boundary","published non-lesion receivers minus transported source positions",
     "20","2;3;4","sqrt[(log(size ratio)/log2)^2+(edge difference/200)^2+(normal fraction difference/0.2)^2]")),"analysis_design.tsv")
 vals<-list();cands<-list();allcands<-list();allvals<-list();meta<-list();geom<-list();mem<-list();actual<-list();checks<-list()
 for(sec in c("KP_1","KP_2","KP_4")) {
  z<-copy(z0[section_id==sec]);xy<-as.matrix(z[,.(x_um,y_um)]);grid<-as.matrix(z[,.(array_row,array_col)])
  keys<-paste(grid[,1],grid[,2],sep=":");stopifnot(!anyDuplicated(keys))
  ix<-which(z$eligible_receiver);base_norm<-z$kmeans=="normal"
  X<-cbind(1,grid);aff<-qr.solve(X,xy);resid<-xy-X%*%aff
  spacing<-median(nn2(xy,xy,k=2)$nn.dists[,2])
  step<-matrix(c(-1,-1,-1,1,0,-2,0,2,1,-1,1,1),ncol=2,byrow=TRUE)
  edge<-Reduce(`|`,lapply(seq_len(nrow(step)),function(i)
    !(paste(grid[,1]+step[i,1],grid[,2]+step[i,2],sep=":")%in%keys)))
  stopifnot(any(edge));edge_dist<-nn2(xy[edge,,drop=FALSE],xy,k=1)$nn.dists[,1]
  geom[[sec]]<-data.table(section_id=sec,mouse=unique(z$mouse),spots=nrow(z),eligible_receivers=length(ix),
    boundary_spots=sum(edge),median_neighbour_spacing_um=spacing,
    lattice_fit_max_coordinate_residual_um=max(abs(resid)),lattice_fit_max_vector_residual_um=max(sqrt(rowSums(resid^2))))
  for(pr in names(defs)) {
   genes<-intersect(defs[[pr]],sub("^gene__","",grep("^gene__",names(z),value=TRUE)))
   genes<-genes[vapply(genes,function(g)all(!is.na(z[[paste0("gene__",g)]])),logical(1))]
   gc<-paste0("gene__",genes);v<-log1p(1e4*rowMeans(as.matrix(z[,..gc]))/pmax(z$whole_library_size,1))
   ref<-z$score_reference_receiver;s<-mad(v[ref],constant=1.4826);if(!is.finite(s)||s==0)s<-sd(v[ref])
   z[,(pr):=(v-median(v[ref]))/s]
   mem[[length(mem)+1L]]<-data.table(section_id=sec,programme=pr,requested=paste(defs[[pr]],collapse=";"),used=paste(genes,collapse=";"))
  }
  # Broad kmeans labels are frozen published normal/tumour strata, not cell types.
  qnorm<-mean(base_norm[ix]);eval_shape<-function(srcxy,source_indices) {
   r<-setdiff(ix,source_indices);d<-nn2(srcxy,xy[r,,drop=FALSE],k=1)$nn.dists[,1]
   nr<-r[d<=200];fr<-r[d>400&d<=800]
   feature<-c(edge_distance_um=mean(edge_dist[source_indices]),source_normal_fraction=mean(base_norm[source_indices]),
     near_normal_fraction=if(length(nr))mean(base_norm[nr]) else NA_real_,
     near_receivers=length(nr),far_receivers=length(fr),
     near_normal_receivers=sum(base_norm[nr]),near_tumour_receivers=sum(!base_norm[nr]),
     far_normal_receivers=sum(base_norm[fr]),far_tumour_receivers=sum(!base_norm[fr]),
     standardization_normal_weight=qnorm)
   out<-lapply(names(defs),function(pr) {
    a<-if(length(nr)&&length(fr))mean(z[[pr]][nr])-mean(z[[pr]][fr]) else NA_real_
    ss<-vapply(c(TRUE,FALSE),function(g) {
      n<-nr[base_norm[nr]==g];f<-fr[base_norm[fr]==g]
      if(length(n)<5L||length(f)<5L)NA_real_ else mean(z[[pr]][n])-mean(z[[pr]][f])
    },numeric(1))
    data.table(programme=pr,endpoint=c("raw","normal_tumour_standardized"),
      effect=c(a,if(all(is.finite(ss)))sum(c(qnorm,1-qnorm)*ss) else NA_real_))
   })
   list(features=feature,effects=rbindlist(out))
  }
  lesions<-unique(z[lesion_class%chin%c("Tgfbr2","KP control"),.(lesion=phenotypes,lesion_class)])
  for(j in seq_len(nrow(lesions))) {
   src<-which(z$phenotypes==lesions$lesion[j]);stopifnot(length(src)>=5L)
   observed<-eval_shape(xy[src,,drop=FALSE],src)
   base<-data.table(section_id=sec,mouse=unique(z$mouse),lesion=lesions$lesion[j],lesion_class=lesions$lesion_class[j],source_spots=length(src))
   actual[[length(actual)+1L]]<-cbind(base,as.data.table(as.list(observed$features)),observed$effects)
   # Anchoring one member of the shape at each tissue grid spot enumerates every
   # possible integer translation once for each fixed orientation, without sampling.
   cc<-list();vv<-list();k<-0L;mask_count<-0L;nonlesion_count<-0L
   translation_mask_count<-0L;translation_nonlesion_count<-0L
   max_shape_error<-0;max_center_error<-0
   # Six fixed 60-degree orientations preserve every physical pair distance.
   # Native hex-grid addresses define mask membership. Rotated physical centres
   # must remain within one-quarter spacing of those addresses; the maximum
   # discrepancy is reported rather than called exact capture-centre alignment.
   relative<-sweep(grid[src,,drop=FALSE],2,grid[src[1],],"-")
   axial<-cbind((relative[,2]-relative[,1])/2,relative[,1])
   physical<-sweep(xy[src,,drop=FALSE],2,xy[src[1],],"-")
   for(rot in 0:5)for(a in seq_len(nrow(z))) {
    ax<-axial
    if(rot>0L)for(rr in seq_len(rot))ax<-cbind(-ax[,2],ax[,1]+ax[,2])
    rotated<-cbind(ax[,2],2*ax[,1]+ax[,2])
    mapped_grid<-sweep(rotated,2,grid[a,],"+")
    shift<-grid[a,]-grid[src[1],];mapped<-match(paste(mapped_grid[,1],mapped_grid[,2],sep=":"),keys)
    if(anyNA(mapped))next
    mask_count<-mask_count+1L
    if(rot==0L)translation_mask_count<-translation_mask_count+1L
    if(any(!z$eligible_receiver[mapped]))next
    nonlesion_count<-nonlesion_count+1L
    if(rot==0L)translation_nonlesion_count<-translation_nonlesion_count+1L
    theta<-rot*pi/3;rotation<-matrix(c(cos(theta),sin(theta),-sin(theta),cos(theta)),nrow=2)
    translated<-sweep(physical%*%t(rotation),2,xy[a,],"+")
    # Inverse rotation recovers the original physical pair vectors exactly.
    shape_error<-max(abs(sweep(translated,2,translated[1,],"-")%*%rotation-physical))
    center_error<-max(sqrt(rowSums((translated-xy[mapped,,drop=FALSE])^2)))
    if(center_error>.25*spacing)next
    stopifnot(length(unique(mapped))==length(src),shape_error<1e-8,
      max(abs(grid[mapped,,drop=FALSE]-mapped_grid))==0)
    max_shape_error<-max(max_shape_error,shape_error);max_center_error<-max(max_center_error,center_error)
    ev<-eval_shape(translated,mapped);k<-k+1L
    cc[[k]]<-cbind(base,data.table(candidate_id=k,rotation_degrees=rot*60,anchor_row=grid[a,1],anchor_col=grid[a,2],
      shift_row=shift[1],shift_col=shift[2],max_grid_center_mismatch_um=center_error),as.data.table(as.list(ev$features)))
    vv[[k]]<-cbind(base,data.table(candidate_id=k),ev$effects)
   }
   cand<-rbindlist(cc);evs<-rbindlist(vv)
   if(nrow(cand)){allcands[[length(allcands)+1L]]<-copy(cand);allvals[[length(allvals)+1L]]<-copy(evs)}
   source_key<-keys[src];adj<-lapply(seq_len(length(src)),function(a)
     match(paste(grid[src[a],1]+step[,1],grid[src[a],2]+step[,2],sep=":"),source_key))
   remaining<-seq_along(src);component_size<-integer()
   while(length(remaining)) {
    reached<-remaining[1];front<-reached
    while(length(front)) {nxt<-setdiff(unique(unlist(adj[front])),reached);nxt<-nxt[!is.na(nxt)];reached<-union(reached,nxt);front<-nxt}
    component_size<-c(component_size,length(reached));remaining<-setdiff(remaining,reached)
   }
   meta[[length(meta)+1L]]<-cbind(base,data.table(mask_valid_placements=mask_count,
     nonlesion_grid_placements=nonlesion_count,physical_alignment_placements=k,
     translation_only_mask_placements=translation_mask_count,translation_only_nonlesion_placements=translation_nonlesion_count,
     source_components=length(component_size),largest_component_spots=max(component_size),
     source_max_x_span_um=diff(range(xy[src,1])),source_max_y_span_um=diff(range(xy[src,2])),
     original_edge_um=observed$features["edge_distance_um"],
     original_source_normal_fraction=observed$features["source_normal_fraction"],
     original_near_normal_fraction=observed$features["near_normal_fraction"],
     original_near_receivers=observed$features["near_receivers"],original_far_receivers=observed$features["far_receivers"]))
   checks[[length(checks)+1L]]<-cbind(base,data.table(tested_placements=nonlesion_count,
     maximum_pair_vector_error_um=max_shape_error,maximum_grid_center_mismatch_um=max_center_error,
     original_sources_excluded=TRUE,all_shapes_in_original_section=TRUE,
     independent_animal_count=length(unique(z0[section_id%chin%c("KP_1","KP_2","KP_4"),mouse]))))
   if(nrow(cand))cand[,`:=`(edge_difference_um=edge_distance_um-observed$features["edge_distance_um"],
     near_normal_difference=near_normal_fraction-observed$features["near_normal_fraction"],
     source_normal_difference=source_normal_fraction-observed$features["source_normal_fraction"],
     near_count_ratio=pmax(near_receivers/observed$features["near_receivers"],observed$features["near_receivers"]/near_receivers),
     far_count_ratio=pmax(far_receivers/observed$features["far_receivers"],observed$features["far_receivers"]/far_receivers))]
   for(l in seq_len(nrow(levels))) {
    good<-if(nrow(cand))cand[abs(edge_difference_um)<=levels$edge_tolerance_um[l]&
      abs(near_normal_difference)<=levels$normal_fraction_tolerance[l]&
      abs(source_normal_difference)<=levels$source_normal_tolerance[l]&
      near_count_ratio<=levels$receiver_count_max_ratio[l]&far_count_ratio<=levels$receiver_count_max_ratio[l]]else data.table()
    if(nrow(good))cands[[length(cands)+1L]]<-cbind(good,data.table(level=levels$level[l]))
    for(p in seq_len(nrow(observed$effects))) {
     pr<-observed$effects$programme[p];ep<-observed$effects$endpoint[p];obs<-observed$effects$effect[p]
     v<-if(nrow(evs)&&nrow(good))evs[programme==pr&endpoint==ep&candidate_id%in%good$candidate_id,effect]else numeric()
     v<-v[is.finite(v)]
     enough<-length(v)>=levels$minimum_candidates[l]&&is.finite(obs)
     centre<-if(length(v)&&is.finite(obs))mean(c(obs,v))else NA_real_
     vals[[length(vals)+1L]]<-cbind(base,data.table(level=levels$level[l],programme=pr,endpoint=ep,
       actual_effect=obs,accepted_candidates=nrow(good),finite_candidates=length(v),
       candidate_mean=if(length(v))mean(v)else NA_real_,candidate_sd=if(length(v)>1L)sd(v)else NA_real_,
       candidate_q025=if(length(v))unname(quantile(v,.025))else NA_real_,
       candidate_q975=if(length(v))unname(quantile(v,.975))else NA_real_,
       conditional_excess=if(enough)obs-mean(v)else NA_real_,
       upper_reference_fraction=if(enough)mean(v>=obs)else NA_real_,
       lower_reference_fraction=if(enough)mean(v<=obs)else NA_real_,
       pooled_absolute_reference_fraction=if(enough)(1+sum(abs(v-centre)>=abs(obs-centre)-1e-12))/(length(v)+1)else NA_real_,
       legacy_null_centred_absolute_fraction=if(enough)(1+sum(abs(v-mean(v))>=abs(obs-mean(v))-1e-12))/(length(v)+1)else NA_real_,
       estimable=enough,reason=if(!is.finite(obs))"actual coarse-stratum endpoint unestimable"else if(length(v)<20L)"fewer than 20 finite matched translations"else"conditional spatial diagnostic only",
       mean_edge_difference_um=if(nrow(good))mean(good$edge_difference_um)else NA_real_,
       mean_near_normal_difference=if(nrow(good))mean(good$near_normal_difference)else NA_real_,
       mean_source_normal_difference=if(nrow(good))mean(good$source_normal_difference)else NA_real_))
    }
   }
  }
 }
 v<-rbindlist(vals);mt<-rbindlist(meta);act<-rbindlist(actual);checked<-rbindlist(checks)
 put(v,"lesion_matched_diagnostics.tsv");put(mt,"candidate_support.tsv");put(act,"actual_endpoints.tsv")
 put(rbindlist(cands),"accepted_placement_geometry.tsv.gz");put(rbindlist(geom),"native_grid_geometry.tsv")
 put(rbindlist(allcands),"all_nonlesion_placement_geometry.tsv.gz");put(rbindlist(allvals),"all_nonlesion_placement_endpoints.tsv.gz")
 put(rbindlist(mem),"programme_membership.tsv");put(checked,"translation_integrity.tsv")
 # Complete-case section contrasts only. Missing a lesion is not silently a
 # different target; the actual and matched result use precisely the same set.
 section<-v[,{
   ok<-all(estimable)&all(c("Tgfbr2","KP control")%in%lesion_class)
   .(lesions=.N,estimable_lesions=sum(estimable),complete=ok,
     mutant_minus_KP_actual=if(ok)mean(actual_effect[lesion_class=="Tgfbr2"])-mean(actual_effect[lesion_class=="KP control"])else NA_real_,
     mutant_minus_KP_matched_excess=if(ok)mean(conditional_excess[lesion_class=="Tgfbr2"])-mean(conditional_excess[lesion_class=="KP control"])else NA_real_)
 },by=.(level,programme,endpoint,section_id,mouse)]
 # Missing sections are not discarded from the mouse4 two-section estimand.
 animal<-section[,.(informative_sections=.N,complete=all(complete),
   actual=if(all(complete))mean(mutant_minus_KP_actual)else NA_real_,
   matched_excess=if(all(complete))mean(mutant_minus_KP_matched_excess)else NA_real_),
   by=.(level,programme,endpoint,mouse)]
 summary<-animal[,.(animals=.N,complete_animals=sum(complete),
   equal_animal_actual=if(all(complete))mean(actual)else NA_real_,
   equal_animal_matched_excess=if(all(complete))mean(matched_excess)else NA_real_),by=.(level,programme,endpoint)]
 put(section,"section_contrasts.tsv");put(animal,"animal_contrasts.tsv");put(summary,"equal_animal_summary.tsv")
 # This composition sensitivity uses all original lesions. It is independent of
 # whether a transported shape had sufficient placement support.
 actual_section<-act[,.(effect=mean(effect[lesion_class=="Tgfbr2"])-mean(effect[lesion_class=="KP control"])),
   by=.(section_id,mouse,programme,endpoint)]
 actual_animal<-actual_section[,.(effect=mean(effect)),by=.(mouse,programme,endpoint)]
 actual_summary<-actual_animal[,.(equal_animal_effect=mean(effect),mouse4_effect=effect[mouse=="mouse4"],mouse7_effect=effect[mouse=="mouse7"]),by=.(programme,endpoint)]
 put(actual_section,"actual_composition_section_contrasts.tsv");put(actual_animal,"actual_composition_animal_contrasts.tsv")
 put(actual_summary,"actual_composition_equal_animal.tsv")
 # A complementary observed-KP geometry match does not claim identical shape.
 # Match within section without replacement; maximise eligible pair count then
 # minimise the fixed geometry distance, before looking at programme outcomes.
 gg<-unique(act[,.(section_id,mouse,lesion,lesion_class,source_spots,edge_distance_um,near_normal_fraction)])
 pairs<-list();supports<-list()
 for(l in seq_len(nrow(levels)))for(sec in c("KP_1","KP_2","KP_4")) {
  mut<-gg[section_id==sec&lesion_class=="Tgfbr2"][order(lesion)]
  kp<-gg[section_id==sec&lesion_class=="KP control"][order(lesion)]
  distances<-outer(seq_len(nrow(mut)),seq_len(nrow(kp)),Vectorize(function(m,k) {
   ds<-abs(log(mut$source_spots[m]/kp$source_spots[k]));de<-abs(mut$edge_distance_um[m]-kp$edge_distance_um[k]);
   dn<-abs(mut$near_normal_fraction[m]-kp$near_normal_fraction[k])
   if(ds>log(levels$receiver_count_max_ratio[l])||de>levels$edge_tolerance_um[l]||dn>levels$normal_fraction_tolerance[l])Inf
   else sqrt((ds/log(2))^2+(de/200)^2+(dn/.2)^2)
  }))
  options<-do.call(expand.grid,c(rep(list(0:nrow(kp)),nrow(mut)),list(KEEP.OUT.ATTRS=FALSE)))
  legal<-vapply(seq_len(nrow(options)),function(i) {
   a<-as.integer(options[i,]);s<-which(a>0L)
   !anyDuplicated(a[s])&&all(is.finite(distances[cbind(s,a[s])]))
  },logical(1));options<-options[legal,,drop=FALSE]
  matched<-rowSums(options>0);cost<-vapply(seq_len(nrow(options)),function(i) {
    a<-as.integer(options[i,]);s<-which(a>0L);sum(distances[cbind(s,a[s])])
  },numeric(1));best<-as.integer(options[order(-matched,cost)[1],])
  supports[[length(supports)+1L]]<-data.table(level=levels$level[l],section_id=sec,mouse=mut$mouse[1],
    mutants=nrow(mut),matched_mutants=sum(best>0L),complete=all(best>0L))
  for(m in seq_len(nrow(mut))) {
   k<-best[m];pairs[[length(pairs)+1L]]<-data.table(level=levels$level[l],section_id=sec,mouse=mut$mouse[m],
    mutant=mut$lesion[m],control=if(k>0L)kp$lesion[k]else NA_character_,matched=k>0L,
    geometry_distance=if(k>0L)distances[m,k]else NA_real_,
    source_size_ratio=if(k>0L)mut$source_spots[m]/kp$source_spots[k]else NA_real_,
    edge_difference_um=if(k>0L)mut$edge_distance_um[m]-kp$edge_distance_um[k]else NA_real_,
    near_normal_difference=if(k>0L)mut$near_normal_fraction[m]-kp$near_normal_fraction[k]else NA_real_)
  }
 }
 pairtab<-rbindlist(pairs);pairend<-rbindlist(lapply(seq_len(nrow(pairtab)),function(i) {
  p<-pairtab[i];m<-act[lesion==p$mutant];k<-act[lesion==p$control]
  cbind(p,m[,.(programme,endpoint)],data.table(effect=if(p$matched)m$effect-k$effect[match(paste(m$programme,m$endpoint),paste(k$programme,k$endpoint))]else NA_real_))
 }))
 real_section<-pairend[,.(complete=all(matched),matched_mutants=sum(matched),mutants=.N,
   effect=if(all(matched))mean(effect)else NA_real_),by=.(level,section_id,mouse,programme,endpoint)]
 real_animal<-real_section[,.(complete=all(complete),sections=.N,effect=if(all(complete))mean(effect)else NA_real_),
   by=.(level,mouse,programme,endpoint)]
 real_summary<-real_animal[,.(complete_animals=sum(complete),equal_animal_effect=if(all(complete))mean(effect)else NA_real_,
   mouse4_effect=effect[mouse=="mouse4"],mouse7_effect=effect[mouse=="mouse7"]),by=.(level,programme,endpoint)]
 put(pairtab,"observed_KP_geometry_pairs.tsv");put(rbindlist(supports),"observed_KP_geometry_support.tsv")
 put(pairend,"observed_KP_geometry_pair_endpoints.tsv");put(real_section,"observed_KP_geometry_section_contrasts.tsv")
 put(real_animal,"observed_KP_geometry_animal_contrasts.tsv");put(real_summary,"observed_KP_geometry_equal_animal.tsv")
 invisible(list(lesion=v,section=section,animal=animal,summary=summary))
}
