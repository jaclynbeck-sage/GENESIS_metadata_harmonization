print_unique_column_vals <- function(df, numeric_columns) {
  for (cn in colnames(df)) {
    if (!(cn %in% c("individualID", "projid"))) {
      vals <- unique(df[, cn])
      if (cn %in% numeric_columns) {
        vals <- sort(vals)
      }
      cat(cn, "(", paste(vals, collapse = ",  "), ")\n\n")
    }
  }
}

print_columns_with_nas <- function(df) {
  vals <- c()
  for (cn in colnames(df)) {
    if (any(is.na(df[, cn]))) {
      vals <- c(vals, cn)
    }
  }

  cat(paste(vals, collapse = ",  "), "\n")
}
