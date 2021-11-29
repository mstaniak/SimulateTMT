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
  simulated_data[, ProtCond := paste(Protein, Condition, sep = "__")]
  simulated_data[, RunChannel := paste(Run, Channel, sep = "_")]
  significant_proteins = unique(simulated_data$Protein)[seq_len(num_significant)]
  simulated_data[, IsSignificant := Protein %in% significant_proteins]
  simulated_data[, NestedBioRep := paste(Mixture, Condition, BioReplicate, sep = "_")]
  simulated_data
}

add_variability = function(labels, sd) {
  num_objects = data.table::uniqueN(labels)
  effects = rnorm(num_objects, 0, sd)
  names(effects) = unique(labels)
  effects[labels]
}

add_condition_variability = function(input, sd_cond_mix) {
  tmt_by_log2fc = split(input, input[["change"]])
  tmt_by_log2fc = lapply(tmt_by_log2fc, function(x) {
    log2fc = unique(x$change)
    sd_cond = sd_cond_mix[[as.character(log2fc)]]
    x[, ConditionVariability := add_variability(MixCond, sd_cond),
      by = "Protein"]
    x
  })
  data.table::rbindlist(tmt_by_log2fc)
}

#' @export
simulate_log_abundances = function(tmt_design, baseline, log2FC, sd_mix,
                                   sd_cond_mix, sd_sub, error_type, sd_error = NULL) {

  tmt_design[, MixtureVariability := add_variability(Mixture, sd_mix),
             by = "Protein"]
  tmt_design[, SubjectVariability := add_variability(NestedBioRep, sd_sub),
             by = "Protein"]
  num_conditions = data.table::uniqueN(tmt_design$Condition)

  if (length(log2FC) == 1) {
    conditions = unique(tmt_design$Condition)
    log2_fcs = c(log2FC, 0)
    names(log2_fcs) = conditions
    change_values = log2_fcs[tmt_design$Condition]
    tmt_design[, change := change_values]
  } else {
    log2_fcs = lapply(log2FC, function(x) {
      if (length(x) < num_conditions) {
        c(x, rep(0, num_conditions - length(x)))
      } else {
        x
      }
    })
    log2_fcs = unlist(log2_fcs, TRUE, FALSE)
    log2FC_dt = data.table::data.table(Protein = rep(sort(unique(tmt_design$Protein))[seq_along(log2FC)],
                                                     each = num_conditions),
                                       Condition = rep(reorder(unique(tmt_design$Condition),
                                                               as.numeric(as.character(unique(tmt_design$Condition)))),
                                                       times = length(log2FC)),
                                       change = log2_fcs)
    tmt_design = merge(tmt_design, log2FC_dt, by = c("Protein", "Condition"), all.x = TRUE)
    tmt_design[, change := ifelse(is.na(change), 0, change)]
  }
  tmt_design = add_condition_variability(tmt_design, sd_cond_mix)

  if (error_type == "constant") {
    tmt_design[, RandomError := rnorm(.N, 0, sd_error)]
  } else {
    tmt_design[, RandomErrorSigma := ifelse(IsSignificant,  0.1 * sqrt(max(change)), rlnorm(1, -2.5, 0.7)),
               by = c("Protein")]
    tmt_design[, RandomError := rnorm(.N, 0, unique(RandomErrorSigma)),
               by = c("Protein")]
    tmt_design[, RandomErrorSigma := NULL]
  }

  tmt_design[, Abundance := baseline + change * IsSignificant + MixtureVariability +
               SubjectVariability + ConditionVariability + RandomError]
  tmt_design[, MixtureVariability := NULL]
  tmt_design[, SubjectVariability := NULL]
  tmt_design[, ConditionVariability := NULL]
  tmt_design[, RandomError := NULL]
  tmt_design[, change := NULL]
  tmt_design = tmt_design[, lapply(.SD, function(x) if (is.character(x)) as.factor(x) else x)]
  tmt_design
}
