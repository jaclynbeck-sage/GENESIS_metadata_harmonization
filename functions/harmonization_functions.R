# This file contains generic harmonization-related functions that are not
# specific to any data set.

# Executed when sourced --------------------------------------------------------

library(dplyr)

# Load all the dataset-specific harmonization functions when this file is
# sourced. All functions are in the same folder: functions/dataset_functions.
dataset_functions <- list.files(file.path("functions", "dataset_functions"),
                                full.names = TRUE)
for (file in dataset_functions) {
  source(file)
}


# Generic harmonization function -----------------------------------------------

# Runs the appropriate dataset-specific function to rename and harmonize
# variables, fills in any NA values with "missing or unknown", and adds any
# missing columns to the data frame. It also de-duplicates studies that need it
# with AMP-AD 1.0 / Diverse Cohorts data.
#
# Arguments:
#   study_name - the name of the study
#   metadata - a `data.frame` of metadata from the source metadata file. Columns
#     are variables and rows are individuals.
#   spec - a `config` object describing the standardized values for each field,
#     as defined by this project's `GENESIS_harmonization.yml` file
#   harmonized_baseline - a `data.frame` of de-duplicated and harmonized
#     metadata from all AMP-AD 1.0 studies and Diverse Cohorts. If the study
#     does not need de-duplication, this should be NULL.
#   extra_metadata - a `data.frame` of extra metadata that is needed by several
#     studies, which has different information depending on study. If the study
#     does not need extra metadata, this should be NULL.
#
# Returns:
#   a `data.frame` with all relevant fields harmonized to the GENESIS data
#   dictionary. Columns not defined in the data dictionary are left as-is.
#
harmonize <- function(study_name, metadata, spec, nps_data = NULL) {
  # Study-specific harmonization
  metadata <- switch(
    study_name,
    "AMP-PD" = harmonize_AMP_PD(metadata, spec),
    "ASAP" = harmonize_ASAP(metadata, spec),
    "BD2" = harmonize_BD2(metadata, spec),
    "AMP-AD_DiverseCohorts" = harmonize_Diverse_Cohorts(metadata, spec),
    "MC-BrAD" = harmonize_MC_BrAD(metadata, spec),
    "MC_snRNA" = harmonize_MC_snRNA(metadata, spec),
    "McCarroll_SCZ" = harmonize_McCarroll_SCZ(metadata, spec),
    "McCarroll_HD" = harmonize_McCarroll_HD(metadata, spec),
    "MCMPS" = harmonize_MCMPS(metadata, spec),
    "MSBB" = harmonize_ADKP_studies(metadata, spec),
    "NPS-AD" = harmonize_NPS_AD(metadata, spec),
    "ROSMAP" = harmonize_ROSMAP(metadata, spec),
    "SEA-AD" = harmonize_SEA_AD(metadata, spec),
    "SMIB-AD" = harmonize_SMIB_AD(metadata, spec),
    .default = metadata
  )

  # Add any missing non-diagnosis fields
  missing_fields <- setdiff(spec$demographic_columns, colnames(metadata))
  for (field in missing_fields) {
    metadata[, field] <- spec$missing
  }

  # Add any missing diagnosis fields, which should have missing values set to
  # NA instead of "missing or unknown"
  missing_diagnosis <- setdiff(spec$diagnosis_columns, colnames(metadata))
  for (field in missing_diagnosis) {
    metadata[, field] <- NA
  }

  # Don't fill NA values with "missing" in the ageDeath or PMI columns.
  # Diagnosis columns are not included in this NA fill operation.
  cols_fill <- setdiff(spec$demographic_columns, c("ageDeath", "PMI"))

  metadata <- metadata |>
    mutate(
      # Fix fields that might be read in as numeric but should be characters
      across(any_of(cols_fill), as.character),
      # Fill NAs in character columns as "missing or unknown". cols_fill does
      # not include diagnosis columns
      across(any_of(cols_fill), ~ ifelse(is.na(.x), spec$missing, .x)),
      # Add or update derived columns
      apoe4Status = get_apoe4Status(apoeGenotype, spec),
      amyAny = get_amyAny(amyCerad, spec),
      amyA = get_amyA(amyThal, spec),
      bScore_NFT = get_bScore(Braak_NFT, spec),
      bScore_LB = get_bScore(Braak_LB, spec),
      # Add study name
      study = study_name
    )

  # If applicable, pull missing information from NPS-AD metadata
  if (!is.null(nps_data)) {
    metadata <- deduplicate_studies(
      list(metadata, nps_data),
      spec,
      verbose = FALSE
    ) |>
      subset(study == study_name) |>
      select(all_of(colnames(metadata)))

    # Extra step for BD2, which had individualID modified to match NPS-AD IDs
    # in the BD2 harmonization function
    if (study_name == spec$study$bd2) {
      metadata$individualID <- metadata$bd2_id
      metadata <- select(metadata, -bd2_id)
    }
  }

  # Put harmonized fields first in the data frame
  all_columns <- c(spec$demographic_columns, spec$diagnosis_columns)
  metadata <- metadata |>
    select(all_of(all_columns), !all_of(all_columns))

  return(metadata)
}


# Convert numbers to Braak stage values ----------------------------------------

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


# Get bScore values based on Braak ---------------------------------------------

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
  return(
    case_match(
      Braak,
      spec$Braak_NFT$none ~ spec$bScore_NFT$none,
      c(spec$Braak_NFT$stage1, spec$Braak_NFT$stage2) ~ spec$bScore_NFT$stage1_2,
      c(spec$Braak_NFT$stage3, spec$Braak_NFT$stage4) ~ spec$bScore_NFT$stage3_4,
      c(spec$Braak_NFT$stage5, spec$Braak_NFT$stage6) ~ spec$bScore_NFT$stage5_6,
      NA ~ spec$missing,
      .default = as.character(Braak)
    )
  )
}


# Get amyAny values based on amyCerad ------------------------------------------

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
  return(
    case_match(
      amyCerad,
      spec$amyCerad$none ~ spec$amyAny$amyNo,
      spec$amyCerad$sparse ~ spec$amyAny$amyYes,
      spec$amyCerad$moderate ~ spec$amyAny$amyYes,
      spec$amyCerad$frequent ~ spec$amyAny$amyYes,
      NA ~ spec$missing,
      .default = as.character(amyCerad)
    )
  )
}


# Get amyA values based on amyThal ---------------------------------------------

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
  return(
    case_match(
      amyThal,
      spec$amyThal$none ~ spec$amyA$none,
      c(spec$amyThal$phase1, spec$amyThal$phase2) ~ spec$amyA$phase1_2,
      spec$amyThal$phase3 ~ spec$amyA$phase3,
      c(spec$amyThal$phase4, spec$amyThal$phase5) ~ spec$amyA$phase4_5,
      NA ~ spec$missing,
      .default = as.character(amyThal)
    )
  )
}


# Get APOE4 status based on genotype -------------------------------------------

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
  return(
    case_match(
      apoeGenotype,
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
    )
  )
}


# Censor ages 90 or above ------------------------------------------------------

# This function censors any ages 90 or above by replacing them with "90+". It
# also converts some studies' versions of this from "90_or_over", "89+", or
# "89+ " to "90+". Empty strings and "missing or unknown" values should be
# replaced with `NA`.
#
# Arguments:
#   ages - a vector of ages, which may be strings or numerical
#   spec - a `config` object describing the standardized values for each field,
#          as defined by this project's `GENESIS_harmonization.yml` file
#
# Returns:
#   a vector of strings with age values properly censored
censor_ages <- function(ages, spec) {
  return(
    case_when(
      ages %in% c("90+", "90_or_over", "89+", "89+ ") ~ spec$ageDeath$over90,
      ages == "" ~ NA,
      ages == spec$missing ~ NA,
      suppressWarnings(as.numeric(ages)) >= 90 ~ spec$ageDeath$over90,
      .default = as.character(ages)
    )
  )
}


# Determine the value of ADoutcome ---------------------------------------------

# This function uses Braak_NFT and amyCerad to determine ADoutcome (AD, Other,
# or Control) with the following criteria:
#   AD = Braak_NFT >= Stage IV, amyCerad = "Moderate" or "Frequent"
#   Control = Braak_NFT <= Stage III, amyCerad = "None" or "Sparse"
#   Other = other combinations (high Braak/low Cerad, low Braak/high Cerad)
#
# Exception: In the AMP-AD_DiverseCohorts study, samples from MayoClinic do not
# have ADoutcome re-computed, as Mayo samples do not have amyCerad values and
# the existing ADoutcome values have already been carefully standardized based
# on pathology for those samples.
#
# This function is intended to be used inside a mutate() function and requires
# multiple columns from the data frame.
#
# Arguments:
#   .data - a data frame, as from inside mutate() statement
#   spec - a `config` object describing the standardized values for each field,
#          as defined by this project's `GENESIS_harmonization.yml` file
#
# Returns:
#   a column of ADoutcome values, which are "Control", "AD", "Other", or
#   "missing or unknown"
determineADoutcome <- function(.data, spec) {
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

  # Don't change ADoutcome for samples coming from Mayo Clinic (Diverse Cohorts
  # only). Defining these variables allows the function to be used for other
  # studies which don't have a MayoDx or ADoutcome field already and which would
  # crash the case_when statement if we didn't adjust for that.
  is_Mayo_dx <- rep(FALSE, length(high_Braak))
  ADoutcome <- rep(spec$missing, length(high_Braak))

  if ("derivedOutcomeBasedOnMayoDx" %in% colnames(.data)) {
    is_Mayo_dx <- .data[["derivedOutcomeBasedOnMayoDx"]]
    ADoutcome <- .data[["ADoutcome"]]
  }

  ADoutcome <- case_when(
    # Leave Mayo diagnoses alone
    is_Mayo_dx ~ ADoutcome,

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

    .default = ADoutcome
  )

  return(ADoutcome)
}


# Turn a text column into a binary 1/0 diagnosis column ------------------------

# This function takes a column of text values and assigns binary values based on
# whether they match a target value. For every value in the column, the binary
# version becomes:
#   1 if the value matches the target value
#   0 if the value doesn't match and is non-missing
#   NA if the value is missing
#
# Arguments:
#   column_values - a column of values from a data frame
#   true_value - the text value in the column that should be considered "1"
#   spec - a `config` object describing the standardized values for each field,
#          as defined by this project's `GENESIS_harmonization.yml` file
#
# Returns:
#   a column of values that have been standardized to 1, 0, and NA based on the
#   contents of column_values
make_binary_column <- function(column_values, true_value, spec) {
  case_match(column_values,
             true_value ~ 1,
             c(spec$missing, NA) ~ NA,
             .default = 0)
}
