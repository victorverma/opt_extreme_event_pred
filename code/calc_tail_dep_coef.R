calc_tail_dep_coef <- function(x, y,
                               type = c("empirical", "modified CFG"),
                               x_threshold = NULL,
                               y_threshold = NULL,
                               block_size = 1) {
  match.arg(type)
  if (type == "empirical") {
    stopifnot(!is.null(x_threshold), !is.null(y_threshold))
    x_vals <- x >= x_threshold
    y_vals <- y >= y_threshold
    levels <- c("TRUE", "FALSE")
    confusion_mat <- table(factor(x_vals, levels), factor(y_vals, levels))
    return(confusion_mat["TRUE", "TRUE"] / sum(confusion_mat["TRUE", ]))
  } else if (type == "modified CFG") {
    n <- length(x)
    if (block_size > 1) {
      group_nums <- (seq_len(n) - 1) %/% block_size
      x <- tapply(x, group_nums, max)
      y <- tapply(y, group_nums, max)
      n <- length(x)
    }
    x_ranks <- rank(x) / (n + 1)
    y_ranks <- rank(y) / (n + 1)
    m <- pmax(x_ranks, y_ranks)
    return(2 - 1 / mean(log(1 / m)))
  }
}