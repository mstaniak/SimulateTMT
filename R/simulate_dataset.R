#' @export
simulate_dataset = function(protein_level_data, protein_groups, combinations,
                            num_unique, num_shared, weights, channel_sd,
                            noise_sd, y_axis_noise_sd, num_cores) {
  unique_proteins = unique(protein_level_data$Protein)
  run_channels = unique(protein_level_data$RunChannel)
  channel_variation = rnorm(length(run_channels), 0, channel_sd)
  names(channel_variation) = run_channels

  all_grouped_proteins = unlist(protein_groups, TRUE, FALSE)
  all_profiles_shared = parallel::mclapply(seq_along(protein_groups), function(i) {
    profiles = protein_level_data[Protein %in% protein_groups[[i]]]
    unique_profiles = simulate_unique_peptides(profiles, num_unique,
                                               channel_variation,
                                               noise_sd, y_axis_noise_sd)
    shared_profiles = simulate_shared_profiles(profiles, combs, num_shared,
                                               weights, channel_variation,
                                               noise_sd, y_axis_noise_sd)
    all_features = combine_features(unique_profiles, shared_profiles)
    all_features
  }, mc.cores = num_cores)
  all_profiles_shared = data.table::rbindlist(all_profiles_shared,
                                              use.names = TRUE, fill = TRUE)
  all_profiles_no_shared = combine_features(simulate_unique_peptides(protein_level_data[!(Protein %in% all_grouped_proteins)],
                                                                     num_unique, channel_variation, noise_sd, y_axis_noise_sd),
                                            data.table())
  all_features = rbind(all_profiles_shared, all_profiles_no_shared,
                       use.names = TRUE, fill = TRUE)
  all_features[, IsUnique := uniqueN(ProteinName) == 1, by = "PSM"]
  all_features
}
