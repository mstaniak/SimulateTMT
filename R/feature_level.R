#' @export
simulate_unique_peptides = function(protein_profiles, num_unique,
                                    noise_sd = 0.1, y_axis_noise_rate = 1.5,
                                    feature_effects_range = c(-1, 1)) {
  psms_unique = data.table::rbindlist(
    lapply(split(protein_profiles, protein_profiles$Protein, drop = FALSE),
           function(protein_profile) {
             data.table::rbindlist(lapply(seq_len(num_unique), function(feature_id) {
               x = protein_profile
               x$feature_eff = runif(1,
                                     feature_effects_range[1],
                                     feature_effects_range[2])
               x$PSM = paste0("P", feature_id)
               x
             }))
           })
  )
  psms_unique[, Log2Intensity := Abundance + feature_eff +
                rnorm(1, sd = y_axis_noise_rate) + rnorm(.N, sd = noise_sd),
              by = c("Protein", "PSM")]
  psms_unique[, PSM := paste(PSM, as.numeric(Protein), sep = "")]
  psms_unique[, feature_eff := NULL]
  psms_unique[, MixCond := NULL]
  psms_unique[, ProtCond := NULL]
  psms_unique[, IsSignificant := NULL]
  psms_unique[, Abundance := NULL]
  data.table::setnames(psms_unique, "Protein", "ProteinName")
  psms_unique
}

#' @export
simulate_shared_profiles = function(profiles, combinations, num_shared,
                                    weights, noise_sd = 0.1,
                                    y_axis_noise_rate = 1.5,
                                    feature_effects_range = c(-1, 1)) {
  proteins = as.character(unique(profiles$Protein))
  channels = unique(profiles$Channel)
  shared_profiles = vector("list", nrow(combinations))
  for (i in seq_len(nrow(combinations))) {
    combination = unlist(combinations[i, ], TRUE, FALSE)
    weights_set = weights[[i]]
    num_shared_i = num_shared[i]

    selected_profiles = profiles[Protein %in% proteins[combination]]
    proteins_dt = unique(selected_profiles[, .(ProteinName = Protein,
                                               ProteinPresent = 1)])
    selected_dt = data.table::dcast(selected_profiles, Channel ~ Protein,
                                    value.var = "Abundance") # TODO: multi-mixture
    channels = unique(selected_dt$Channel)
    proteins_matrix = as.matrix(selected_dt[, -1, with = FALSE])
    shared_profile = proteins_matrix %*% weights_set
    shared_profile = data.table::data.table(Protein = 0,
                                            Channel = channels,
                                            Abundance = shared_profile[, 1])
    shared_annot = unique(profiles[, .(Run, Mixture, TechRepMixture,
                                       Channel, BioReplicate, Condition)])
    shared_profile = merge(shared_annot, shared_profile, by = "Channel")
    
    shared_profiles_i = suppressWarnings(simulate_unique_peptides(shared_profile, 10, 0.2, 2.1))
    shared_profiles_i[, PSM := as.numeric(as.factor(PSM))]
    protein_ids = paste(proteins[combination], sep = "_", collapse = "_")
    shared_profiles_i[, PSM := paste("PS", PSM, protein_ids, sep = "_")]
    shared_profiles_i = data.table::rbindlist(lapply(proteins[combination],
                                                    function(prot) {
                                                      x = shared_profiles_i
                                                      x$ProteinName = prot
                                                      x
                                                    }))  
    shared_profiles[[i]] = shared_profiles_i
  }

  data.table::rbindlist(shared_profiles)
}

#' @export
combine_features = function(unique_features, shared_features) {
  all_features = rbind(unique_features, shared_features)
  all_features[, Intensity := 2 ^ Log2Intensity]
  all_features[, PSM := gsub("_", "", PSM)]
  all_features[, PeptideSequence := PSM]
  all_features[, PSM := paste(PSM, 1, sep = "_")]
  all_features[,
               list(ProteinName, PeptideSequence, Charge = 1, PSM, Run, Mixture,
                    TechRepMixture, Channel, BioReplicate, Condition, Intensity)]
}
