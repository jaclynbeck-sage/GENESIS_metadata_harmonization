# This file contains functions related to QC and validation of data sets pre-
# and post-harmonization.

library(dplyr)

# Print summary of dataset quality ---------------------------------------------

# Prints summaries of all expected columns for visual inspection:
#   1. Reports which columns are missing from the data frame
#   2. Reports how many NA values are in each column
#   3. Prints tables of unique values in each column vs the count of each value,
#      excluding `ageDeath` and `PMI`.
#   4. Checks if there are any `ageDeath` values of 90 or over that are not
#      censored as "90+".
#
# Arguments:
#   df - a data.frame where rows are individuals and columns are variables
#   <X>_col - the real name of the column in `df` that corresponds to <X>.
#             Default values are provided, but many data sets have column names
#             that are different from the desired harmonized name (i.e. "CERAD"
#             instead of "amyCerad" or "ethnicity" instead of "isHispanic"). If
#             no column corresponding to <X> exists, `<X>_col` can be left as
#             the default value and will be reported as missing by this
#             function.
# Returns:
#   nothing
print_summary <- function(df,
                          ageDeath_col = "ageDeath",
                          isHispanic_col = "isHispanic",
                          pmi_col = "PMI",
                          race_col = "race",
                          sex_col = "sex",
                          apoe_col = "apoeGenotype",
                          apoe4status_col = "apoe4Status",
                          cerad_col = "amyCerad",
                          amyAny_col = "amyAny",
                          thal_col = "amyThal",
                          amyA_col = "amyA",
                          braak_nft_col = "Braak_NFT",
                          bscore_nft_col = "bScore_NFT",
                          braak_lb_col = "Braak_LB",
                          bscore_lb_col = "bScore_LB") {
  all_cols <- c(
    ageDeath_col, isHispanic_col, pmi_col, race_col, sex_col, apoe_col,
    apoe4status_col, cerad_col, amyAny_col, thal_col, amyA_col, braak_nft_col,
    bscore_nft_col, braak_lb_col, bscore_lb_col
  )

  missing <- setdiff(all_cols, colnames(df))
  cat("Missing columns:", paste(missing, collapse = ", "), "\n\n")

  cat("NA value check:\n\n")
  for (col_name in all_cols) {
    if (col_name %in% colnames(df)) {
      cat(col_name, "has", sum(is.na(df[, col_name])), "NA values\n")
    }
  }

  cat("\nUnique value check:\n\n")
  for (col_name in setdiff(all_cols, c(ageDeath_col, pmi_col))) {
    if (col_name %in% colnames(df)) {
      cat(col_name, "values:\n")
      tmp <- df |>
        group_by_at(col_name) |>
        count() |>
        ungroup() |>
        mutate_if(is.character, ~ paste0("\"", .x, "\"")) |>
        data.frame()

      print(tmp)
      cat("\n")
    }
  }

  if (ageDeath_col %in% colnames(df)) {
    cat("Age check:\n")
    ages_remove <- c("90+", "89+", "89+ ", "90_or_over", "Missing or unknown", "missing or unknown")
    tmp <- subset(df, !(df[, ageDeath_col] %in% ages_remove)) |>
      mutate(ageDeath = suppressWarnings(as.numeric(.data[[ageDeath_col]]))) |>
      subset(ageDeath >= 90)

    cat(ageDeath_col, "has", nrow(tmp), "uncensored ages.\n")
  }
}


# Validation of harmonized results ---------------------------------------------

# Given a data frame of harmonized metadata, this function validates the following:
#   1. There are no values in any harmonized field that aren't in the data dictionary
#   2. There are no `ageDeath` values above 89
#   3. The `ageDeath` and `PMI` columns only have numbers, NAs, or "90+"
#   4. Columns whose values are derived from other columns (`apoe4Status`,
#      `amyA`, `amyAny`, and `bScore`) have the correctly-derived values. This
#      check is needed to catch any accidental differences introduced by filling
#      in missing data from Diverse Cohorts / AMP-AD 1.0.
#
# If a column passes validation, the phrase "OK <column>" will print out if
# `verbose = TRUE`. If a column fails validation, the function will print out
# "X" plus a message describing the failure and failing values. Then an error
# will be raised after everything has been checked.
#
# Arguments:
#   metadata - a data frame of harmonized metadata, where rows are individuals
#          and columns are variables
#   spec - a `config` object describing the standardized values for each field,
#          as defined by this project's `GENESIS_harmonization.yml` file
#   verbose - if TRUE, the phrase "OK <column>" will print out if the column
#          passes validation. If FALSE, nothing will print out for columns that
#          pass. Columns that fail validation always have a print out, so
#          setting verbose = FALSE will result in only failures being printed.
#
# Returns:
#   nothing
validate_values <- function(metadata, spec, verbose = TRUE) {
  errors <- list()

  # ageDeath should have only NA, numbers, or "90+". No numbers should be above
  # 89.
  ageDeath <- na.omit(metadata$ageDeath) |>
    setdiff(spec$ageDeath$over90) |>
    as.numeric()

  if (any(ageDeath >= 90, na.rm = TRUE)) {
    errors[["ageDeath"]] <- ageDeath[ageDeath >= 90]
    cat("X  ageDeath\n")
  } else if (any(is.na(ageDeath))) {
    # NAs introduced by coercion to numeric, warn about initial value
    nas <- which(is.na(ageDeath))
    tmp <- na.omit(metadata$ageDeath) |>
      setdiff(spec$ageDeath$over90)

    errors[["ageDeath"]] <- tmp[nas]
    cat("X  ageDeath\n")
  } else if (verbose) {
    cat("OK ageDeath\n")
  }

  # PMI should have only NA or numbers
  pmi <- na.omit(metadata$PMI)
  if (length(pmi) > 0 && !is.numeric(pmi)) {
    errors[["PMI"]] <- "not numeric"
    cat("X  PMI\n")
  } else if (verbose) {
    cat("OK PMI\n")
  }

  # Other columns should have only string values that exist in the dictionary
  cols_check <- setdiff(spec$required_columns, c("individualID", "ageDeath", "PMI"))

  for (col_name in cols_check) {
    values <- metadata[, col_name]
    expected <- c(spec$missing, unlist(spec[[col_name]]))
    if (!all(values %in% expected)) {
      errors[[col_name]] <- setdiff(values, expected)
      cat("X ", col_name, "\n")
    } else if (verbose) {
      cat("OK", col_name, "\n")
    }
  }

  # Also check agreement between related columns
  if (!identical(metadata$apoe4Status, get_apoe4Status(metadata$apoeGenotype, spec))) {
    errors[["apoe4Status"]] <- "does not match apoeGenotype"
    cat("X apoe4Status does not match apoeGenotype\n")
  } else if (verbose) {
    cat("OK apoeGenotype vs apoe4Status\n")
  }

  if (!identical(metadata$amyA, get_amyA(metadata$amyThal, spec))) {
    errors[["amyA"]] <- "does not match amyThal"
    cat("X amyA does not match amyThal\n")
  } else if (verbose) {
    cat("OK amyThal vs amyA\n")
  }

  if (!identical(metadata$amyAny, get_amyAny(metadata$amyCerad, spec))) {
    errors[["amyAny"]] <- "does not match amyCerad"
    cat("X amyAny does not match amyCerad\n")
  } else if (verbose) {
    cat("OK amyCerad vs amyAny\n")
  }

  if (!identical(metadata$bScore_NFT, get_bScore(metadata$Braak_NFT, spec))) {
    errors[["bScore_NFT"]] <- "does not match Braak_NFT"
    cat("X bScore_NFT does not match Braak_NFT\n")
  } else if (verbose) {
    cat("OK Braak_NFT vs bScore_NFT\n")
  }

  if (!identical(metadata$bScore_LB, get_bScore(metadata$Braak_LB, spec))) {
    errors[["bScore_LB"]] <- "does not match Braak_LB"
    cat("X bScore_LB does not match Braak_LB\n")
  } else if (verbose) {
    cat("OK Braak_LB vs bScore_LB\n")
  }

  # Throw errors, if any
  if (length(errors) > 0) {
    cat("\n\n")
    for (col_name in names(errors)) {
      cat(col_name, "has invalid values:",
          paste(errors[[col_name]], collapse = ", "), "\n")
    }

    stop("Dataset failed validation.")
  }
}
