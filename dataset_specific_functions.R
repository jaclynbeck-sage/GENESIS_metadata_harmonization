# This file contains one function per GENESIS data set that performs all
# harmonization operations. Specific changes to each data set are documented
# in the comments of the corresponding functions. These functions are tied to
# specific versions of each source file and may need to be updated if the source
# metadata file is updated.

# Generic harmonization function
#
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
harmonize <- function(study_name, metadata, spec, harmonized_baseline = NULL) {
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

  # Add any missing fields
  missing_fields <- setdiff(spec$required_columns, colnames(metadata))
  for (field in missing_fields) {
    metadata[, field] <- spec$missing
  }

  # Don't fill NA values with "missing" in the ageDeath or PMI columns
  cols_fill <- setdiff(spec$required_columns, c("ageDeath", "PMI"))

  metadata <- metadata |>
    mutate(
      # Fix fields that might be read in as numeric but should be characters
      across(any_of(cols_fill), as.character),
      # Fill NAs in character columns as "missing or unknown"
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

  # If applicable, pull missing information from AMP-AD 1.0 and Diverse Cohorts metadata
  if (!is.null(harmonized_baseline)) {
    metadata <- deduplicate_studies(
      list(metadata, harmonized_baseline),
      spec,
      verbose = FALSE
    ) |>
      subset(study == study_name) |>
      select(all_of(colnames(metadata)))

    # Extra step for BD2, which had individualID modified to match AMP-AD 1.0 IDs
    # in the BD2 harmonization function
    if (study_name == spec$study$bd2) {
      metadata$individualID <- metadata$bd2_id
      metadata <- select(metadata, -bd2_id)
    }
  }

  # Put harmonized fields first in the data frame
  metadata <- metadata |>
    select(all_of(spec$required_columns), !all_of(spec$required_columns))

  return(metadata)
}
