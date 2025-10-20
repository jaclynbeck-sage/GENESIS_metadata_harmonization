# This file contains utility functions and variables that are used in multiple
# steps/scripts and are not dataset-specific. Functions include things like
# Synapse upload/download, quality control printouts, converting raw values to
# harmonized values, harmonized value validation, and de-duplication of data
# across multiple studies.

# Global variable -- all columns that are harmonized and should be in each
# metadata table
expectedColumns <- c(
  "individualID", "dataContributionGroup", "cohort", "sex", "race",
  "isHispanic", "ageDeath", "PMI", "apoeGenotype", "apoe4Status", "amyCerad",
  "amyAny", "amyThal", "amyA", "Braak_NFT", "bScore_NFT", "Braak_LB",
  "bScore_LB", "study"
)


# Quality control printouts
#
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
print_qc <- function(df,
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
    # "89+" applies to NPS-AD only, all other studies use "90+" or "90_or_over"
    ages_remove <- c("90+", "89+", "90_or_over", "Missing or unknown", "missing or unknown")
    tmp <- subset(df, !(df[, ageDeath_col] %in% ages_remove)) |>
      mutate(ageDeath = suppressWarnings(as.numeric(.data[[ageDeath_col]]))) |>
      subset(ageDeath >= 90)

    cat(ageDeath_col, "has", nrow(tmp), "uncensored ages.\n")
  }
}


# Validation of harmonized results
#
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
# "X" plus a message describing the failure and failing values.
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
  # ageDeath should have only NA, numbers, or "90+". No numbers should be above
  # 89.
  ageDeath <- na.omit(metadata$ageDeath) |>
    setdiff(spec$over90) |>
    as.numeric()

  if (any(ageDeath >= 90, na.rm = TRUE)) {
    cat("X  ageDeath has uncensored ages above 90.\n")
  } else if (any(is.na(ageDeath))) {
    # NAs introduced by coercion to numeric, warn about initial value
    nas <- which(is.na(ageDeath))
    tmp <- na.omit(metadata$ageDeath) |>
      setdiff(spec$over90)
    cat(
      "X  ageDeath has invalid age values:",
      paste(tmp[nas], collapse = ", "), "\n"
    )
  } else if (verbose) {
    cat("OK ageDeath\n")
  }

  # PMI should have only NA or numbers
  pmi <- na.omit(metadata$PMI)
  if (length(pmi) > 0 && !is.numeric(pmi)) {
    cat("X  PMI is not numeric\n")
  } else if (verbose) {
    cat("OK PMI\n")
  }

  # Other columns should have only string values that exist in the dictionary
  cols_check <- setdiff(expectedColumns, c("individualID", "ageDeath", "PMI"))

  for (col_name in cols_check) {
    values <- metadata[, col_name]
    expected <- c(spec$missing, unlist(spec[[col_name]]))
    if (!all(values %in% expected)) {
      cat(
        "X ", col_name, "has unexpected values:",
        paste(setdiff(values, expected), collapse = ", "), "\n"
      )
    } else if (verbose) {
      cat("OK", col_name, "\n")
    }
  }

  # Also check agreement between related columns
  if (!identical(metadata$apoe4Status, get_apoe4Status(metadata$apoeGenotype, spec))) {
    cat("X apoe4Status does not match apoeGenotype\n")
  } else if (verbose) {
    cat("OK apoeGenotype vs apoe4Status\n")
  }

  if (!identical(metadata$amyA, get_amyA(metadata$amyThal, spec))) {
    cat("X amyA does not match amyThal\n")
  } else if (verbose) {
    cat("OK amyThal vs amyA\n")
  }

  if (!identical(metadata$amyAny, get_amyAny(metadata$amyCerad, spec))) {
    cat("X amyAny does not match amyCerad\n")
  } else if (verbose) {
    cat("OK amyCerad vs amyAny\n")
  }

  if (!identical(metadata$bScore_NFT, get_bScore(metadata$Braak_NFT, spec))) {
    cat("X bScore_NFT does not match Braak_NFT\n")
  } else if (verbose) {
    cat("OK Braak_NFT vs bScore_NFT\n")
  }

  if (!identical(metadata$bScore_LB, get_bScore(metadata$Braak_LB, spec))) {
    cat("X bScore_LB does not match Braak_LB\n")
  } else if (verbose) {
    cat("OK Braak_LB vs bScore_LB\n")
  }
}


# Convert numbers to Braak stage values
#
# This function converts numerical Braak values (0-6) to "Stage " + the Roman
# numeral version of the number, as defined by the GENESIS data dictionary.
#
# Arguments:
#   num - a vector containing numerical Braak stages from 0 to 6 or `NA`
#   spec - a `config` object describing the standardized values for each field,
#          as defined by this project's `GENESIS_harmonization.yml` file
#
# Returns:
#   a vector with all Braak values converted to "Stage " + Roman numeral, and
#   all `NA` values converted to "missing or unknown". This function always
#   defaults to returning the original value as a character string if it doesn't
#   meet any of the criteria in the `case_match` statement, so the value will
#   fail validation and can be examined.
to_Braak_stage <- function(num, spec) {
  return(case_match(num,
    0 ~ spec$Braak_NFT$none,
    1 ~ spec$Braak_NFT$stage1,
    2 ~ spec$Braak_NFT$stage2,
    3 ~ spec$Braak_NFT$stage3,
    4 ~ spec$Braak_NFT$stage4,
    5 ~ spec$Braak_NFT$stage5,
    6 ~ spec$Braak_NFT$stage6,
    NA ~ spec$missing,
    .default = as.character(num)
  ))
}


# Get bScore values based on Braak
#
# This function turns (harmonized) Braak scores into the corresponding values
# for bScore in the data dictionary:
#   "None" => "None"
#   "Stage I", "Stage II" => "Braak Stage I-II"
#   "Stage III", "Stage IV" => "Braak Stage III-IV"
#   "Stage V", "Stage VI" => "Braak Stage V-VI"
#   "missing or unknown" => "missing or unknown"
#
# Arguments:
#   Braak - a vector containing harmonized, data dictionary-compliant Braak
#           scores, which are either "Stage " + a Roman numeral or "missing or
#           unknown".
#   spec - a `config` object describing the standardized values for each field,
#          as defined by this project's `GENESIS_harmonization.yml` file
#
# Returns:
#   a vector of strings with bScore values derived from Braak. Values should be
#   as described above.
#
#   This function always defaults to returning the original Braak value as a
#   character string if it doesn't meet any of the criteria in the `case_match`
#   statement, so the value will fail validation and can be examined.
get_bScore <- function(Braak, spec) {
  return(case_match(Braak,
    spec$Braak_NFT$none ~ spec$bScore_NFT$none,
    c(spec$Braak_NFT$stage1, spec$Braak_NFT$stage2) ~ spec$bScore_NFT$stage1_2,
    c(spec$Braak_NFT$stage3, spec$Braak_NFT$stage4) ~ spec$bScore_NFT$stage3_4,
    c(spec$Braak_NFT$stage5, spec$Braak_NFT$stage6) ~ spec$bScore_NFT$stage5_6,
    NA ~ spec$missing,
    .default = as.character(Braak)
  ))
}


# Get amyAny values based on amyCerad
#
# This function turns (harmonized) amyCerad scores into the corresponding values
# for amyAny in the data dictionary:
#   "None/No AD/C0" => "0"
#   "Sparse/Possible/C1", "Moderate/Probable/C2", "Frequent/Definite/C3" => "1"
#   "missing or unknown" => "missing or unknown"
#
# Arguments:
#   amyCerad - a vector containing harmonized, data dictionary-compliant
#          amyCerad scores, which should all be strings
#   spec - a `config` object describing the standardized values for each field,
#          as defined by this project's `GENESIS_harmonization.yml` file
#
# Returns:
#   a vector of strings with amyAny values derived from amyCerad. Values should
#   be as described above.
#
#   This function always defaults to returning the original amyCerad value as a
#   character string if it doesn't meet any of the criteria in the `case_match`
#   statement, so the value will fail validation and can be examined.
get_amyAny <- function(amyCerad, spec) {
  return(case_match(amyCerad,
    spec$amyCerad$none ~ spec$amyAny$amyNo,
    spec$amyCerad$sparse ~ spec$amyAny$amyYes,
    spec$amyCerad$moderate ~ spec$amyAny$amyYes,
    spec$amyCerad$frequent ~ spec$amyAny$amyYes,
    NA ~ spec$missing,
    .default = as.character(amyCerad)
  ))
}


# Get amyA values based on amyThal
#
# This function turns (harmonized) amyThal scores into the corresponding values
# for amyA in the data dictionary:
#   "None" => "None"
#   "Phase 1", "Phase 2" => "Thal Phase 1 or 2"
#   "Phase 3" => "Thal Phase 3"
#   "Phase 4", "Phase 5" => "Thal Phase 4 or 5"
#   "missing or unknown" => "missing or unknown"
#
# Arguments:
#   amyThal - a vector containing harmonized, data dictionary-compliant
#          amyThal scores, which should all be strings
#   spec - a `config` object describing the standardized values for each field,
#          as defined by this project's `GENESIS_harmonization.yml` file
#
# Returns:
#   a vector of strings with amyA values derived from amyThal, with values as
#   described above.
#
#   This function always defaults to returning the original amyThal value as a
#   character string if it doesn't meet any of the criteria in the `case_match`
#   statement, so the value will fail validation and can be examined.
get_amyA <- function(amyThal, spec) {
  return(case_match(amyThal,
    spec$amyThal$none ~ spec$amyA$none,
    c(spec$amyThal$phase1, spec$amyThal$phase2) ~ spec$amyA$phase1_2,
    spec$amyThal$phase3 ~ spec$amyA$phase3,
    c(spec$amyThal$phase4, spec$amyThal$phase5) ~ spec$amyA$phase4_5,
    NA ~ spec$missing,
    .default = as.character(amyThal)
  ))
}


# Get APOE4 status based on genotype
#
# This function turns (harmonized) apoeGenotype scores into the corresponding
# values for apoe4Status in the data dictionary:
#   "22", "23", "33" => "no"
#   "24", "34", "44" => "yes"
#   "missing or unknown" => "missing or unknown"
#
# Arguments:
#   apoeGenotype - a vector containing harmonized, data dictionary-compliant
#          apoeGenotype scores, which should all be strings
#   spec - a `config` object describing the standardized values for each field,
#          as defined by this project's `GENESIS_harmonization.yml` file
#
# Returns:
#   a vector of strings with apoe4Status values derived from apoeGenotype, with
#   values as described above.
#
#   This function always defaults to returning the original apoeGenotype value
#   as a character string if it doesn't meet any of the criteria in the
#   `case_match` statement, so the value will fail validation and can be
#   examined.
get_apoe4Status <- function(apoeGenotype, spec) {
  return(case_match(apoeGenotype,
    c(
      spec$apoeGenotype$e2e4,
      spec$apoeGenotype$e3e4,
      spec$apoeGenotype$e4e4
    ) ~ spec$apoe4Status$e4yes,
    c(
      spec$apoeGenotype$e2e2,
      spec$apoeGenotype$e2e3,
      spec$apoeGenotype$e3e3
    ) ~ spec$apoe4Status$e4no,
    NA ~ spec$missing,
    .default = as.character(apoeGenotype)
  ))
}


# Censor ages 90 or above
#
# This function censors any ages 90 or above by replacing them with "90+". It
# also converts some studies' versions of this from "90_or_over" or "89+" to
# "90+". Empty strings and "missing or unknown" values should be replaced with
# `NA`.
#
# Arguments:
#   ages - a vector of ages, which may be strings or numerical
#   spec - a `config` object describing the standardized values for each field,
#          as defined by this project's `GENESIS_harmonization.yml` file
#
# Returns:
#   a vector of strings with age values properly censored
censor_ages <- function(ages, spec) {
  return(case_when(
    ages %in% c("90+", "90_or_over", "89+") ~ spec$over90,
    ages == "" ~ NA,
    ages == spec$missing ~ NA,
    suppressWarnings(as.numeric(ages)) >= 90 ~ spec$over90,
    .default = as.character(ages)
  ))
}


# Write a metadata data frame to a file
#
# This is a wrapper around `write.csv` that has some extra handling for values
# that contain commas and end of line characters (\n). Values with commas
# need to be escaped with quotes for a CSV file, and some data sets have them
# escaped already and some don't. \n characters are removed entirely, as quote
# escaping doesn't affect them.
#
# Arguments:
#   metadata - a data frame of harmonized metadata where rows are individuals
#              and columns are variables
#   filename - the full path and name of the file to be written. This function
#              automatically inserts "_harmonized" just before ".csv" in the
#              file name.
#
# Returns:
#   the new file name that was written, which should have "_harmonized" added
#   before ".csv"
write_metadata <- function(metadata, filename) {
  # Put quotes around values with commas
  for (column in colnames(metadata)) {
    if (is.character(metadata[, column])) {
      # Columns that contain commas and aren't already escaped with quotes
      commas <- grepl(",", metadata[, column]) & !grepl("\"", metadata[, column])
      metadata[commas, column] <- paste0("\"", metadata[commas, column], "\"")

      # Remove any "\n" characters
      metadata[, column] <- str_replace_all(metadata[, column], "\n", "")
    }
  }

  new_filename <- file.path(
    "data", "output", str_replace(filename, "\\.(csv|txt)", "_harmonized.csv")
  )
  write.csv(metadata, new_filename,
    row.names = FALSE, quote = FALSE
  )

  return(new_filename)
}


# Upload a file to Synapse
#
# This is a wrapper around `synStore` to shorten code slightly. It uploads a
# file to Synapse but only if the contents are different than what is currently
# on Synapse. This check is done because if the file is the same, synStore will
# erase any version comments in Synapse even if the file itself doesn't get a
# new version number. This function also makes sure that synStore doesn't erase
# any annotations on the file in Synapse.
#
# Arguments:
#   filename - the full path and name of the file to upload
#   folder_id - the Synapse ID of the folder on Synapse where the file should be
#              uploaded.
#
# Returns:
#   a Synapse `File` object containing information about the uploaded file
synapse_upload <- function(filename, folder_id) {
  syn_info <- synapse_get_info(filename, folder_id)

  if (!is.null(syn_info)) {
    md5 <- tools::md5sum(filename)

    # Don't actually update the file.
    if (md5 == syn_info$get("_file_handle")$contentMd5) {
      message(str_glue("\"{basename(filename)}\" matches the file on Synapse ",
                       "and will not be re-uploaded."))
      return(syn_info)
    }
  }

  syn_file <- File(filename, parent = folder_id)
  syn_file <- synStore(syn_file, forceVersion = FALSE, set_annotations = FALSE)
  return(syn_file)
}


# Search a folder on Synapse to see if a file is already there. If it is,
# return info on the file without downloading it. Otherwise, return NULL.
synapse_get_info <- function(filename, folder_id) {
  id <- synFindEntityId(basename(filename), folder_id)
  if (is.null(id)) {
    return(NULL)
  }
  synGet(id, downloadFile = FALSE)
}


# Download a file from Synapse
#
# This is a wrapper around `synGet` to shorten code slightly. All downloads go
# into "data/downloads", and if a file with that name already exists, the old
# file is overwritten with the new one to avoid making multiple copies.
#
# Arguments:
#   syn_id - the Synapse ID of the file on Synapse to download
#
# Returns:
#   a Synapse `File` object containing information about the downloaded file
synapse_download <- function(syn_id) {
  synGet(syn_id,
    downloadLocation = file.path("data", "downloads"),
    ifcollision = "overwrite.local"
  )
}


# Check Synapse for new file versions
#
# All Synapse IDs used for this code have the file's version number included for
# reproducibility. This function checks to see if there are newer versions
# available on Synapse than what is specified in the code, and prints a message
# if that's the case.
#
# Arguments:
#   syn_id_list - a named list where each item is a Synapse ID of the format
#         "syn123" or "syn123.5", where the optional number after the decimal is
#         the file version on Synapse. If no version is specified, a warning
#         is printed stating that the latest version of the file on Synapse will
#         be used.
#
# Returns:
#   nothing
check_new_versions <- function(syn_id_list) {
  for (dataset_name in names(syn_id_list)) {
    syn_id <- syn_id_list[[dataset_name]]
    vals <- str_split_1(syn_id, pattern = "\\.")

    if (length(vals) != 2 || is.na(suppressWarnings(as.numeric(vals[2])))) {
      cat(
        str_glue(
          "WARNING: no valid version specified for '{syn_id}' ({dataset_name}). ",
          "The latest version will be used for harmonization."
        ),
        "\n"
      )
    } else {
      id <- vals[1]
      version <- vals[2]
      syn_file <- synGet(id, downloadFile = FALSE)

      if (syn_file$versionNumber != version) {
        cat(
          str_glue(
            "WARNING: there is a new version of {id} ({dataset_name}): ",
            "{version} => {syn_file$versionNumber}. Version {version} will be ",
            "used for harmonization."
          ),
          "\n"
        )
      }
    }
  }
}


# De-duplicate metadata from different studies
#
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
#   4. For cases where cohort values disagree, resolve as follows:
#       a) Replace "ROSMAP" with the cohort value from Diverse Cohorts /
#          AMP-AD 1.0,
#       b) Ignore cases where the disagreement is "Mayo Clinic" vs "Banner",
#       c) Otherwise, report the difference but don't change any values
#   5. For cases where "MSBB_corrections" disagrees with another source,
#      the value from "MSBB_corrections" is always used except for the
#      "apoeGenotype" and "apoe4Status" columns, where the "NPS-AD" value is
#      used instead if it is available.
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
                                include_cols = expectedColumns,
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

  meta_all <- purrr::list_rbind(df_list)
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
          } else if ("MSBB_corrections" %in% meta_tmp$study) {
            if (col_name %in% c("race", "isHispanic")) {
              # Do not correct NPS-AD race or isHispanic values
              meta_tmp[meta_tmp$study != "NPS-AD", col_name] <-
                meta_tmp[meta_tmp$study == "MSBB_corrections", col_name]
            } else {
              meta_tmp[, col_name] <- meta_tmp[meta_tmp$study == "MSBB_corrections", col_name]
            }
          } else if (col_name %in% c("ageDeath", "PMI")) {
            meta_tmp <- deduplicate_ageDeath_pmi(meta_tmp, leftover, col_name, report_string)
          } else if (col_name == "cohort") {
            meta_tmp <- deduplicate_cohort(meta_tmp, leftover, col_name, spec, report_string)
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


# Age- and PMI-specific handling for de-duplication
#
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

  # If there's a real mismatch, try using the NPS-AD value if it exists. If not,
  # report the mismatch. `num_vals` should have one unique value if all numbers
  # are roughly equal, AND `num_vals` should be the same length as `leftover`.
  # If it's not, that means not all values in `leftover` are numeric (i.e. one
  # may be 90+). We report it but don't try and resolve the duplication. Note
  # that `all.equal` returns a string with the difference between two numbers if
  # they are not equal, rather than FALSE, so we have to check for != TRUE
  # instead of == FALSE.
  if (any(equivalent != TRUE)) {
    if ("NPS-AD" %in% meta_tmp$study) {
      meta_tmp[, col_name] <- meta_tmp[meta_tmp$study == "NPS-AD", col_name]
    } else {
      cat(report_string)
    }
  } else if (length(leftover) != ncol(meta_tmp)) {
    # No real mismatch but at least one value was NA and there are multiple
    # close-enough values. Use NPS-AD value first if it exists, then Diverse
    # Cohorts if it exists. Otherwise print.
    na_vals <- which(is.na(meta_tmp[, col_name]))

    if ("NPS-AD" %in% meta_tmp$study) {
      meta_tmp[na_vals, col_name] <- meta_tmp[meta_tmp$study == "NPS-AD", col_name]
    } else if ("AMP-AD_DiverseCohorts" %in% meta_tmp$study) {
      meta_tmp[na_vals, col_name] <- meta_tmp[meta_tmp$study == "AMP-AD_DiverseCohorts", col_name]
    } else {
      cat("Unresolved NA fill: ", report_string)
    }
  }

  # Otherwise we leave the ages as-is even if there's difference in precision

  return(meta_tmp)
}


# Cohort-specific handling for de-duplication
#
# This function is used inside `deduplicate_studies` to resolve duplication of
# cohort values for a single individual. If this function is called, then
# there are at least 2 distinct, non-NA values in the `cohort` column that are
# assigned to this individual. This function handles two special cases and
# defaults to just reporting the duplication if neither special case applies:
#   1. NPS-AD reports the `cohort` of all ROSMAP samples as "ROSMAP" rather than
#      identifying them as "ROS" or "MAP" separately. All of these samples exist
#      in Diverse Cohorts or ROSMAP 1.0 metadata, which has the correct
#      separation into "ROS" and "MAP", so we replace NPS-AD's ROSMAP `cohort`
#      values with the correct values from DC/ROSMAP 1.0.
#   2. Mayo Clinic 1.0 metadata reports cohort as "Mayo Clinic" while the same
#      samples in Diverse Cohorts are reported as "Banner". This difference
#      is ignored as only the Diverse Cohorts metadata is used directly.
#
# Arguments:
#   meta_tmp - a data frame with 2 or more rows, where all rows have the same
#     individual ID and there are 2 or more distinct values in the `cohort`
#     column
#   leftover - a vector of unique values from <col_name> for this individual,
#     which has had `NA` values removed
#   col_name - the name of the column being handled ("cohort")
#   report_string - the string that gets printed out if there is a real
#     difference between values that is not due to precision. The string
#     contains the `individualID`, column name, and unique values in the column.
#
# Returns:
#   meta_tmp, which will have some `cohort` values replaced if the original
#   value was "ROSMAP". Other values and columns are left as-is.
deduplicate_cohort <- function(meta_tmp, leftover, col_name, spec, report_string) {
  # If there is more than one left over value and the values don't meet the two
  # special cases below, report it but don't try to resolve duplication.

  # Special case: NPS-AD reports cohort on some samples as "ROSMAP", which needs
  # to instead use the cohort value from Diverse Cohorts or ROSMAP 1.0 data.
  if ("ROSMAP" %in% leftover) {
    cohort_val <- setdiff(leftover, "ROSMAP")

    # If there is still more than one value for cohort, or no remaining values,
    # report it but don't resolve the de-duplication
    if (length(cohort_val) != 1) {
      cat(report_string)
    } else {
      # Otherwise resolve duplication
      meta_tmp[, col_name] <- cohort_val
    }
  } else {
    # Special case: Mayo data may be labeled as "Mayo Clinic" in the original
    # Mayo metadata or "Banner" in Diverse Cohorts, but we don't need to change
    # this in either metadata file or print it out.
    is_mayo_banner <- identical(sort(leftover), c(spec$cohort$banner, spec$cohort$mayo))

    if (!is_mayo_banner) {
      cat(report_string)
    }
  }

  return(meta_tmp)
}


# Fill missing AMP-AD 1.0 values in Diverse Cohorts
#
# The Diverse Cohorts metadata has an `individualID_AMPAD_1.0` column that is
# supposed to have the corresponding AMP-AD 1.0 ID if that sample exists in 1.0
# data. However this field has a lot of NAs for samples that do exist in 1.0
# data, whether intentional or not. This function fills in the matching ID from
# 1.0 data in these cases.
#
# Arguments:
#   meta_all - a data.frame of all Diverse Cohorts, AMP-AD 1.0, and NPS-AD data
#     as returned by `deduplicate_studies()`.
#
# Returns:
#   meta_all but with appropriate `individualID_AMPAD_1.0` values filled in
fill_missing_ampad1.0_ids <- function(meta_all, spec) {
  dc <- subset(meta_all, study == "AMP-AD_DiverseCohorts") |>
    select(individualID, cohort, individualID_AMPAD_1.0, study)

  ampad_1.0 <- subset(meta_all, study %in% c("MayoRNAseq", "MSBB", "ROSMAP")) |>
    mutate(
      original_individualID = individualID,
      individualID = str_replace(individualID, "AMPAD_MSSM_[0]+", ""),
      individualID = case_when(
        individualID == "29637" &
          dataContributionGroup == spec$dataContributionGroup$mssm ~ "29637_MSSM",
        individualID == "29582" &
          dataContributionGroup == spec$dataContributionGroup$mssm ~ "29582_MSSM",
        .default = as.character(individualID)
      )
    ) |>
    select(individualID, original_individualID, cohort)

  matches_1.0 <- merge(dc, ampad_1.0)
  stopifnot(length(unique(matches_1.0$individualID)) == nrow(matches_1.0))

  cat(str_glue("Filling {sum(is.na(matches_1.0$individualID_AMPAD_1.0))} ",
               "missing AMPAD-1.0 IDs in Diverse Cohorts\n"))

  # This does nothing to the values that are already filled in, but replaces
  # NAs with the 1.0 ID
  matches_1.0$individualID_AMPAD_1.0 <- matches_1.0$original_individualID

  matches_1.0 <- select(matches_1.0, -original_individualID)

  col_order <- colnames(meta_all)

  meta_all |>
    # Replace column
    select(-individualID_AMPAD_1.0) |>
    merge(matches_1.0, all = TRUE, sort = FALSE) |>
    # Restore original column order
    select(all_of(col_order))
}


# Determine the value of ADoutcome for Diverse Cohorts data
#
# ADoutcome is not changed for data not from Diverse Cohorts, and data from
# Diverse Cohorts but contributed by Mayo.
#
# Arguments:
#   .data - a single row of a data frame, as from inside a rowwise() |> mutate() statement
#   spec - a `config` object describing the standardized values for each field,
#          as defined by this project's `GENESIS_harmonization.yml` file
#
# Returns:
#   a single value for ADoutcome, one of "Control", "AD", "Other", or "missing or unknown"
determineADoutcome <- function(.data, spec) {
  # Don't make an ADoutcome value for non-Diverse Cohorts data, and don't change
  # ADoutcome for samples coming from Mayo Clinic
  if (.data[["study"]] != "AMP-AD_DiverseCohorts" |
      .data[["derivedOutcomeBasedOnMayoDx"]] == TRUE) {
    return(.data[["ADoutcome"]])
  }

  high_Braak <- .data[["Braak_NFT"]] %in% c(spec$Braak_NFT$stage4,
                                            spec$Braak_NFT$stage5,
                                            spec$Braak_NFT$stage6)
  low_Braak <- .data[["Braak_NFT"]] %in% c(spec$Braak_NFT$none,
                                           spec$Braak_NFT$stage1,
                                           spec$Braak_NFT$stage2,
                                           spec$Braak_NFT$stage3)

  high_amyCerad <- .data[["amyCerad"]] %in% c(spec$amyCerad$moderate,
                                              spec$amyCerad$frequent)
  low_amyCerad <- .data[["amyCerad"]] %in% c(spec$amyCerad$none,
                                             spec$amyCerad$sparse)

  ADoutcome <- case_when(
    # AD: Braak IV-VI & Cerad Moderate or Frequent
    high_Braak & high_amyCerad ~ "AD",

    # Control: Braak 0-III & Cerad None or Sparse
    low_Braak & low_amyCerad ~ "Control",

    # Other: Braak IV-VI & Cerad None or Sparse
    high_Braak & low_amyCerad ~ "Other",

    # Other: Braak 0-III & Cerad Moderate or Frequent
    low_Braak & high_amyCerad ~ "Other",

    # Missing one or both of Braak and Cerad
    .data[["Braak_NFT"]] == spec$missing |
      .data[["amyCerad"]] == spec$missing ~ spec$missing,

    .default = .data[["ADoutcome"]]
  )

  return(ADoutcome)
}


# Update ADoutcome values for Diverse Cohorts data in case Braak or amyCerad
# values changed during de-duplication.
#
# ADoutcome is not changed for data not from Diverse Cohorts, and data from
# Diverse Cohorts but contributed by Mayo.
#
# Arguments:
#   meta_all - dataframe of harmonized data, which may contain data from any study
#   spec - a `config` object describing the standardized values for each field,
#          as defined by this project's `GENESIS_harmonization.yml` file
#
# Returns:
#   meta_all with ADoutcome values updated for Diverse Cohorts data only
updateADoutcome <- function(meta_all, spec) {
  meta_all |>
    rowwise() |>
    mutate(ADoutcome = determineADoutcome(.data, spec)) |>
    ungroup() |>
    as.data.frame()
}
