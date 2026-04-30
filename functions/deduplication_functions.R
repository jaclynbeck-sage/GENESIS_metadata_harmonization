# This file has functions that de-duplicate studies by comparing data from
# overlapping individuals and filling in missing information when one study has
# missing data and the other study has data.

library(dplyr)
library(stringr)

# De-duplicate metadata from different studies ---------------------------------

# This function takes a list of metadata data frames, concatenates them, and
# then attempts to resolve cases where different studies have different data for
# the same individual.
#
# For each individual ID that has duplicate rows, resolve duplicates:
#   1. For columns where some rows have NA and some have a unique non-NA value,
#      replace the NA value with that unique non-NA value.
#   2. For columns where some rows have "missing or unknown" and some have a
#      unique value other than that, replace "missing or unknown" with the
#      unique value.
#   3. For the ageDeath/PMI columns where rows have different numbers, report
#      the difference but leave the values as-is. If rows disagree only because
#      of precision, nothing is reported.
#   4. For cases where "MSBB" disagrees with another source, the value from
#      "MSBB" is always used except for the "apoeGenotype" and "apoe4Status"
#      columns, where the "NPS-AD" value is used instead if it is available.
#      MSBB data is preferentially used because it has been corrected by updated
#      data, while some NPS values have not.
#
# When duplicated data is un-resolvable, either because it is a special case
# that is intentionally flagged or because this function doesn't have anything
# implemented to handle it, the following information is printed out:
#   "<individualID> <column name> [<list of values for this individual/column>]"
#
# Note: Individual IDs from MSSM-related studies changed format in Diverse
# Cohorts. In order to be able to compare overlapping individuals, all IDs are
# temporarily converted to the Diverse Cohorts format, then reverted back to
# their original values for the returned data frame.
#
# Note: To shorten this function and make it more readable, some processing has
# been broken out into separate functions.
#
# Arguments:
#   df_list - a list of data frames, each of which must include an `individualID`
#       column as well as every column listed in `include_cols`
#   spec - a `config` object describing the standardized values for each field,
#       as defined by this project's `GENESIS_harmonization.yml` file
#   include_cols - a vector of column names to de-duplicate
#   exclude_cols - a vector of column names that should be excluded from
#       de-duplication even if they are listed in `include_cols`
#   verbose - if FALSE, only issues or un-resolvable data will be printed. If
#       TRUE, information on every column with duplicate values for each
#       individual will be reported even the duplication is resolved or ignored.
#
# Returns:
#   a single data frame containing all rows from all data frames in `df_list`,
#   with column values de-duplicated where possible. The data frame will contain
#   all columns present in any data frame in `df_list`, and values will be `NA`
#   for rows that come from data frames without that column.
deduplicate_studies <- function(df_list,
                                spec,
                                include_cols = spec$required_columns,
                                exclude_cols = c("study"),
                                verbose = TRUE) {
  # Make sure certain fields in each data frame are of the same type. Also
  # convert MSSM-style individual IDs to Diverse Cohorts-style IDs.
  df_list <- lapply(df_list, function(df_item) {
    df_item |>
      mutate(
        across(any_of(c("individualID", "apoeGenotype", "amyAny")), as.character),
        # Special case: MSBB/MSSM samples were re-named for Diverse Cohorts, this
        # makes them directly comparable
        original_individualID = individualID,
        individualID = str_replace(individualID, "AMPAD_MSSM_0+", ""),
        individualID = case_when(
          individualID == "29637" &
            dataContributionGroup == spec$dataContributionGroup$mssm ~ "29637_MSSM",
          individualID == "29582" &
            dataContributionGroup == spec$dataContributionGroup$mssm ~ "29582_MSSM",
          .default = as.character(individualID)
        ),
        # Special case: Some samples contributed by "Emory" are from "Mt Sinai
        # Brain Bank" and they need to be included with MSSM samples during
        # de-duplication
        original_dataContributionGroup = dataContributionGroup,
        dataContributionGroup = case_when(
          dataContributionGroup == spec$dataContributionGroup$emory &
            cohort == spec$cohort$msbb ~ spec$dataContributionGroup$mssm,
          .default = dataContributionGroup
        )
      )
  })

  meta_all <- purrr::list_rbind(df_list) |> as.data.frame()
  include_cols <- setdiff(include_cols, exclude_cols)

  # Find IDs that have multiple/duplicate rows. Adding dataContributionGroup
  # accounts for overlapping IDs between different studies that don't refer to
  # the same individual
  dupe_ids <- meta_all |>
    select(all_of(include_cols)) |>
    distinct() |>
    select(individualID, dataContributionGroup) |>
    mutate(
      group_id = paste(individualID, dataContributionGroup),
      duplicate = duplicated(group_id)
    ) |>
    subset(duplicate == TRUE) |>
    distinct()

  # Resolve duplicated data for each individual as best as possible
  for (row_id in 1:nrow(dupe_ids)) {
    ind_id <- dupe_ids$individualID[row_id]

    # This will be altered to resolve duplication, and will get added back to the
    # meta_all data frame
    meta_tmp <- subset(meta_all, individualID == ind_id &
                         dataContributionGroup == dupe_ids$dataContributionGroup[row_id])

    for (col_name in include_cols) {
      unique_vals <- unique(meta_tmp[, col_name])

      if (length(unique_vals) > 1) {
        # For reporting un-resolved mismatches
        report_string <- paste(
          ind_id, col_name, "[", paste(unique_vals, collapse = ", "), "]\n"
        )

        # The ageDeath/PMI validation function will print out the string if
        # there's a mismatch. All other columns should print out here.
        if (verbose) {
          cat(report_string)
        }

        # Remove NA, "", and "missing or unknown" values and see what's left
        leftover <- setdiff(unique_vals, c(spec$missing, "")) |>
          na.omit()

        # If one unique value left, replace all values with that value
        if (length(leftover) == 1) {
          meta_tmp[, col_name] <- leftover
        } else if (length(leftover) == 0) {
          # Nothing left, use "missing or unknown" if it's there, otherwise set
          # to NA.
          if (any(unique_vals == spec$missing)) {
            meta_tmp[, col_name] <- spec$missing
          } else {
            meta_tmp[, col_name] <- NA
          }
        } else {
          # More than 1 value left after omission of NA and "missing or unknown".
          # Do some column-specific handling of duplicates
          if (col_name %in% c("apoeGenotype", "apoe4Status") &&
              "NPS-AD" %in% meta_tmp$study &&
              meta_tmp[meta_tmp$study == "NPS-AD", col_name] != spec$missing) {
            # Use the NPS-AD value for apoe genotype / status where it disagrees
            # with other data sets
            meta_tmp[, col_name] <- meta_tmp[meta_tmp$study == "NPS-AD", col_name]
          } else if ("MSBB" %in% meta_tmp$study) {
            if (col_name %in% c("race", "isHispanic")) {
              # Do not correct NPS-AD race or isHispanic values
              meta_tmp[meta_tmp$study != "NPS-AD", col_name] <-
                meta_tmp[meta_tmp$study == "MSBB", col_name]
            } else {
              meta_tmp[, col_name] <- meta_tmp[meta_tmp$study == "MSBB", col_name]
            }
          } else if (col_name %in% c("ageDeath", "PMI")) {
            meta_tmp <- deduplicate_ageDeath_pmi(meta_tmp, leftover, col_name, report_string)

          } else {
            # Column is something else that we don't have specific handling for,
            # so we report it but don't try to resolve duplication.
            cat(report_string)
          }
        }
      }
    }

    # Replace original rows with de-duplicated data
    rows_replace <- meta_all$individualID == ind_id &
      meta_all$dataContributionGroup == unique(meta_tmp$dataContributionGroup)
    meta_all[rows_replace, ] <- meta_tmp
  }

  # Revert back to original individualIDs and dataContributionGroups
  meta_all <- meta_all |>
    mutate(individualID = original_individualID,
           dataContributionGroup = original_dataContributionGroup) |>
    select(-original_individualID, -original_dataContributionGroup)

  return(meta_all)
}


# Age- and PMI-specific handling for de-duplication ----------------------------

# This function is used inside `deduplicate_studies` to resolve duplication of
# age and PMI data for a single individual. If this function is called, then
# there are at least 2 distinct, non-NA values in the `ageDeath` or `PMI`
# column that are assigned to this individual. This function checks whether the
# numbers differ only by precision (e.g. 3.5 vs 3.547), or if the numbers are
# completely different (e.g. 4 vs 10). The former case is ignored, and the
# latter case will be reported to the console. Currently, no modification is
# done to the values themselves, as NPS-AD wishes to keep all of their values
# as-is even where they differ from Diverse Cohorts / AMP-AD 1.0.
#
# Arguments:
#   meta_tmp - a data frame with 2 or more rows, where all rows have the same
#     individual ID and there are 2 or more distinct values in the `ageDeath`
#     and/or `PMI` column
#   leftover - a vector of unique values from <col_name> for this individual,
#     which has had `NA` values removed
#   col_name - the name of the column being handled, either "ageDeath" or "PMI"
#   report_string - the string that gets printed out if there is a real
#     difference between values that is not due to precision. The string
#     contains the `individualID`, column name, and unique values in the column.
#
# Returns:
#   meta_tmp unaltered (giving us the option to alter it in a future update)
deduplicate_ageDeath_pmi <- function(meta_tmp, leftover, col_name, report_string) {
  # By this point there is more than one unique, non-NA age value in `leftover`.
  # Check if all values are roughly equal to make sure there isn't an actual
  # mis-match. There may be more than 2 values so this is the quickest way to
  # check.
  num_vals <- suppressWarnings(as.numeric(leftover)) |>
    na.omit()

  equivalent <- sapply(num_vals, all.equal, num_vals[1], tolerance = 1e-3)

  # If there's a real mismatch, try using the MSBB value if it exists. If not,
  # then use the NPS-AD if it exists. If not, report the mismatch. `num_vals`
  # should have one unique value if all numbers are roughly equal, AND
  # `num_vals` should be the same length as `leftover`. If it's not, that means
  # not all values in `leftover` are numeric (i.e. one may be 90+). We report it
  # but don't try and resolve the duplication. Note that `all.equal` returns a
  # string with the difference between two numbers if they are not equal, rather
  # than FALSE, so we have to check for != TRUE instead of == FALSE.
  if (any(equivalent != TRUE)) {
    if ("MSBB" %in% meta_tmp$study) {
      meta_tmp[, col_name] <- meta_tmp[meta_tmp$study == "MSBB", col_name]
    }
    else if ("NPS-AD" %in% meta_tmp$study) {
      meta_tmp[, col_name] <- meta_tmp[meta_tmp$study == "NPS-AD", col_name]
    } else {
      cat(report_string)
    }
  } else if (length(leftover) != ncol(meta_tmp)) {
    # No real mismatch but at least one value was NA and there are multiple
    # close-enough values. Use MSBB value first if it exists, then NPS-AD if it
    # exists, then Diverse Cohorts if it exists. Otherwise print.
    na_vals <- which(is.na(meta_tmp[, col_name]))

    if ("MSBB" %in% meta_tmp$study &&
        !is.na(meta_tmp[meta_tmp$study == "MSBB", col_name])) {
      meta_tmp[na_vals, col_name] <- meta_tmp[meta_tmp$study == "MSBB", col_name]
    }
    else if ("NPS-AD" %in% meta_tmp$study &&
             !is.na(meta_tmp[meta_tmp$study == "NPS-AD", col_name])) {
      meta_tmp[na_vals, col_name] <- meta_tmp[meta_tmp$study == "NPS-AD", col_name]
    } else if ("AMP-AD_DiverseCohorts" %in% meta_tmp$study &&
               !is.na(meta_tmp[meta_tmp$study == "AMP-AD_DiverseCohorts", col_name])) {
      meta_tmp[na_vals, col_name] <- meta_tmp[meta_tmp$study == "AMP-AD_DiverseCohorts", col_name]
    } else {
      cat("Unresolved NA fill: ", report_string)
    }
  }

  # Otherwise we leave the ages as-is even if there's difference in precision

  return(meta_tmp)
}
