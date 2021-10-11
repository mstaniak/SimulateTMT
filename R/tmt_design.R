#' @import data.table
#' @export
create_TMT_design = function(num_proteins, num_significant,
                             num_mixtures, num_conditions,
                             num_subjects, experiment_type = "time") {
  proteins = data.table::data.table(Protein = seq_len(num_proteins),
                                    Join = TRUE)
  proteins[, Protein := as.character(Protein)]
  max_length = max(stringr::str_length(proteins$Protein))
  proteins[, Protein := stringr::str_pad(Protein, width = 4, pad = "0", side = "left")]
  proteins[, Protein := paste("Prot", Protein, sep = "_")]
  num_channels = num_mixtures * num_subjects

  experiment_design = data.table::data.table(
    Mixture = rep(seq_len(num_mixtures),
                  each = num_channels),
    BioReplicate = rep(rep(seq_len(num_subjects),
                           times = num_conditions),
                       times = num_mixtures),
    Condition = rep(rep(seq_len(num_conditions),
                        each = num_subjects),
                    times = num_mixtures),
    Join = TRUE)
  if (experiment_type == "groupComparison") {
    experiment_design[, BioReplicate := paste(Condition, BioReplicate, sep = "_")]
  }
  experiment_design[, BioReplicate := paste(Mixture, BioReplicate, sep = "_")]
  simulated_data = merge(proteins, experiment_design,
                         allow.cartesian = TRUE, by = "Join")
  simulated_data[, Join := NULL]
  simulated_data[, TechRepMixture := 1]
  simulated_data[, Channel := 1:.N, by = c("Protein", "Mixture")]
  simulated_data[, Run := paste(Mixture, TechRepMixture, sep = "_")]
  simulated_data = simulated_data[, lapply(.SD, as.character)]
  simulated_data[, MixCond := paste(Mixture, Condition, sep = "_")]
  simulated_data[, ProtCond := paste(Protein, Condition, sep = "_")]
  simulated_data[, RunChannel := paste(Run, Channel, sep = "_")]
  significant_proteins = unique(simulated_data$Protein)[seq_len(num_significant)]
  simulated_data[, IsSignificant := Protein %in% significant_proteins]
  simulated_data[, NestedBioRep := paste(Mixture, Condition, BioReplicate, sep = "_")]
  simulated_data
}


#' @export
simulate_log_abundances = function(tmt_design, baseline, log2FC, sd_mix,
                                   sd_cond_mix, sd_sub, sd_error) {
  num_mixtures = data.table::uniqueN(tmt_design$Mixture)
  num_conditions = data.table::uniqueN(tmt_design$Condition)
  num_subjects = uniqueN(tmt_design$NestedBioRep)

  mixture_effects = rnorm(num_mixtures, 0, sd_mix)
  names(mixture_effects) = unique(tmt_design$Mixture)
  subject_effects = rnorm(num_subjects, 0, sd_sub)
  names(subject_effects) = unique(tmt_design$NestedBioRep)
  condition_effects = rnorm(num_conditions * num_mixtures, 0, sd_cond_mix)
  names(condition_effects) = unique(tmt_design$MixCond)

  random_error = rnorm(seq_len(nrow(tmt_design)), 0, sd_error)
  mix_cond_error = condition_effects[tmt_design$MixCond]
  mix_error = mixture_effects[tmt_design$Mixture]
  sub_error = subject_effects[tmt_design$NestedBioRep]

  if (length(log2FC) == 1) {
    conditions = unique(tmt_design$Condition)
    log2_fcs = c(0, log2FC)
    names(log2_fcs) = conditions
    change = log2_fcs[tmt_design$Condition]
  } else {
    prot_conds = unique(tmt_design$ProtCond)
    log2_fcs = lapply(log2FC, function(x) c(x, 0))
    log2_fcs = unlist(log2_fcs, TRUE, FALSE)
    log2_fcs = c(log2_fcs, rep(0, length(prot_conds) - length(log2_fcs)))
    names(log2_fcs) = prot_conds
    change = log2_fcs[tmt_design$ProtCond]
  }

  abundances = baseline + change * tmt_design$IsSignificant +
    mix_error + mix_cond_error + sub_error + random_error
  tmt_design = tmt_design[, lapply(.SD, as.factor)]
  tmt_design$Abundance = abundances
  tmt_design
}
