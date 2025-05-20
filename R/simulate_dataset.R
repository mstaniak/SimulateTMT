#' @export
simulate_dataset = function(protein_level_data, protein_groups, combinations,
                            num_unique, num_shared, weights, channel_sd,
                            noise_sd, y_axis_range, num_cores) {
  protein_level_data = protein_level_data[, .(Protein, Mixture, BioReplicate,
                                              Condition, TechRepMixture,
                                              Channel, Run, Abundance)]
  all_grouped_proteins = unlist(protein_groups, TRUE, FALSE)
  all_profiles_shared = parallel::mclapply(seq_along(protein_groups), function(i) {
    profiles_source = protein_level_data[Protein %in% protein_groups[[i]]]
    data.table::rbindlist(
      lapply(split(profiles_source, profiles_source[, Run]), function(profiles) {
        unique_profiles = simulate_unique_peptides(profiles, num_unique,
                                                   noise_sd, y_axis_range)
        shared_profiles = simulate_shared_profiles(profiles, combinations, num_shared,
                                                   weights, noise_sd, y_axis_range)
        all_features = combine_features(unique_profiles[, .(ProteinName, PSM, Run, Mixture, TechRepMixture,
                                                            BioReplicate, Condition, Channel, Log2Intensity)],
                                        shared_profiles[, .(ProteinName, PSM, Run, Mixture, TechRepMixture,
                                                            BioReplicate, Condition, Channel, Log2Intensity)])
        all_features
      }
      ))
  }, mc.cores = num_cores)
  all_profiles_shared = data.table::rbindlist(all_profiles_shared,
                                              use.names = TRUE, fill = TRUE)
  all_profiles_no_shared = combine_features(simulate_unique_peptides(protein_level_data[!(Protein %in% all_grouped_proteins)],
                                                                     num_unique, noise_sd, y_axis_range),
                                            data.table())
  all_features = rbind(all_profiles_shared, all_profiles_no_shared,
                       use.names = TRUE, fill = TRUE)
  all_features[, PSM := paste(PSM, Run, sep = "_")]
  all_features[, ChannelVariability := add_variability(Channel, channel_sd),
               by = "ProteinName"] # ProteinGroup
  all_features[, Log2Intensity := Log2Intensity + ChannelVariability]
  all_features[, Intensity := 2 ^ Log2Intensity]
  all_features[, ChannelVariability := NULL]
  all_features[, IsUnique := uniqueN(ProteinName) == 1, by = "PSM"]
  all_features
}
