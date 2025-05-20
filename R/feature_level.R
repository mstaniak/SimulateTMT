#' @export
simulate_unique_peptides = function(protein_profiles, num_unique, noise_sd,
                                    y_axis_range) {
  by_protein_profiles = split(protein_profiles, protein_profiles[["Protein"]],
                              drop = FALSE)
  if (length(num_unique) == 1) {
    num_unique = rep(num_unique, length(by_protein_profiles))
  }
  psms_unique = data.table::rbindlist(
    lapply(seq_along(num_unique),
           function(i) {
             protein_profile = by_protein_profiles[[i]]
             data.table::rbindlist(lapply(seq_len(num_unique[i]),
                                          function(feature_id) {
                                            x = protein_profile
                                            x$PSM = paste0("P", feature_id)
                                            x
                                          }))
           })
  )
  psms_unique[, Log2Intensity := Abundance + runif(1, y_axis_range[1],
                                                   y_axis_range[2]) + rnorm(.N, sd = noise_sd),
              by = c("Protein", "PSM")]
  psms_unique[, PSM := paste(PSM, as.numeric(Protein), sep = "")]
  psms_unique[, Abundance := NULL]
  data.table::setnames(psms_unique, "Protein", "ProteinName")
  psms_unique
}

#' @export
simulate_shared_profiles = function(profiles, combinations, num_shared,
                                    weights, noise_sd,
                                    y_axis_range) {
  proteins = sort(as.character(unique(profiles$Protein)))
  channels = unique(profiles$Channel)
  shared_profiles = vector("list", nrow(combinations))
  for (i in seq_len(nrow(combinations))) {
    combination = as.logical(unlist(combinations[i, ], TRUE, FALSE))
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
    shared_profile = data.table::data.table(Protein = "0",
                                            Channel = channels,
                                            Abundance = shared_profile[, 1])
    shared_annot = unique(profiles[, .(Run, Mixture, TechRepMixture,
                                       Channel, BioReplicate, Condition)])
    shared_profile = merge(shared_annot, shared_profile, by = "Channel")

    shared_profiles_i = suppressWarnings(simulate_unique_peptides(shared_profile,
                                                                  num_shared_i,
                                                                  noise_sd,
                                                                  y_axis_range))
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
  if (nrow(unique_features) > 0 | nrow(shared_features) > 0) {
    charge = 1
  } else {
    charge = numeric(0)
  }
  all_features = rbind(unique_features, shared_features)
  all_features[, Intensity := 2 ^ Log2Intensity]
  all_features[, PSM := gsub("_", "", PSM)]
  all_features[, PeptideSequence := PSM]
  all_features[, PSM := paste(PSM, 1, sep = "_")]
  all_features[,
               list(ProteinName, PeptideSequence, Charge = charge, PSM, Run, Mixture,
                    TechRepMixture, Channel, BioReplicate, Condition, Log2Intensity, Intensity)]
}
