calc_test_tail_dep_coef <- function(tail_dep_coef_type,
                                    q_lin_pred,
                                    flux_threshold,
                                    block_size,
                                    obs_preds) {
  if (is.null(obs_preds)) return(NA_real_)
  if ("lin_pred" %in% names(obs_preds)) {
    test_tail_dep_coef <- calc_tail_dep_coef(
      obs_preds$lin_pred, obs_preds$y,
      tail_dep_coef_type,
      q_lin_pred, flux_threshold,
      block_size
    )
  } else {
    test_tail_dep_coef <- NA_real_
  }
  test_tail_dep_coef
}

make_confusion_mat <- function(obs_preds) {
  if (is.null(obs_preds)) return(NULL)
  # Use factor() to prevent zero rows and columns from being dropped
  with(
    filter(obs_preds, !was_y_imputed),
    table( # Specify the levels so that they'll be sorted properly
      pred = factor(pred, levels = c("TRUE", "FALSE")),
      obs = factor(obs, levels = c("TRUE", "FALSE"))
    )
  )
}

calc_event_rate <- function(confusion_mat) {
  if (is.null(confusion_mat)) return(NA_real_)
  sum(confusion_mat[, "TRUE"]) / sum(confusion_mat)
}

calc_alarm_rate <- function(confusion_mat) {
  if (is.null(confusion_mat)) return(NA_real_)
  sum(confusion_mat["TRUE", ]) / sum(confusion_mat)
}

# The true positive rate (TPR) is the probability of an alarm given an event
calc_tpr <- function(confusion_mat) {
  if (is.null(confusion_mat)) return(NA_real_)
  # If there are no events, return 1, the best possible TPR
  coalesce(confusion_mat["TRUE", "TRUE"] / sum(confusion_mat[, "TRUE"]), 1)
}

# The false positive rate (FPR) is the probability of an alarm given a non-event
calc_fpr <- function(confusion_mat) {
  if (is.null(confusion_mat)) return(NA_real_)
  # If there are no non-events, return 0, the best possible FPR
  coalesce(confusion_mat["TRUE", "FALSE"] / sum(confusion_mat[, "FALSE"]), 0)
}

# The precision is the probability of an event given an alarm
calc_precision <- function(confusion_mat) {
  if (is.null(confusion_mat)) return(NA_real_)
  confusion_mat["TRUE", "TRUE"] / sum(confusion_mat["TRUE", ])
}

# The True Skill Statistic (TSS) is defined as the difference between the TPR
# and the FPR
calc_tss <- function(confusion_mat) {
  if (is.null(confusion_mat)) return(NA_real_)
  calc_tpr(confusion_mat) - calc_fpr(confusion_mat)
}