# Generic harmonize function for PEC-related studies
#
# Source files: locally downloaded files for CommonMind Consortium (CMC) and
#   SZBDMulti-Seq provided by Jaro on May 12, 2026.
#
# Modifications needed for version May 12, 2026:
#   * Rename columns:
#     * `reportedGender` => `sex`
#     * `ethnicity` => `isHispanic`
#     * `individualIdSource` => `cohort`
#     * `Braak` => `Braak_NFT`
#   * Convert Braak_NFT from numeric to Roman numerals
#   * Create binary diagnosis columns for BD, SCZ, Control
#   * Update `cohort` values to conform to the data dictionary
#   * Create `dataContributionGroup` based on `cohort`
harmonize_PEC <- function(metadata, spec) {
  metadata |>
    dplyr::rename(sex = reportedGender,
                  isHispanic = ethnicity,
                  cohort = individualIdSource,
                  Braak_NFT = Braak) |>
    dplyr::mutate(
      # Per Jaro, assume all "89+" ages are "90+"
      ageDeath = ifelse(ageDeath == "89+", spec$age$over90, ageDeath),
      # Braak_NFT needs to be converted to numeric for the case where the column
      # is all NAs
      Braak_NFT = to_Braak_stage(as.numeric(Braak_NFT), spec),
      # The HBTRC cohort is already the right value
      cohort = ifelse(cohort == "MSSM", spec$cohort$msbb, cohort),
      dataContributionGroup = case_match(
        cohort,
        spec$cohort$msbb ~ spec$dataContributionGroup$mssm,
        spec$cohort$harvard ~ spec$dataContributionGroup$harvard,
        .default = cohort
      ),
      BD = make_binary_column(primaryDiagnosis,
                              c("Bipolar disorder", "Bipolar Disorder"),
                              spec),
      SCZ = make_binary_column(primaryDiagnosis, "Schizophrenia", spec),
      Control = (BD + SCZ == 0)
    )
}

# Harmonize CommonMind Consortium (CMC) data
#
# Harmonize using the generic PEC function but also fix CMC-specific issues
# in the data:
#   * Update all `isHispanic` values to conform to the data dictionary
harmonize_PEC_CMC <- function(metadata, spec) {
  harmonize_PEC(metadata, spec) |>
    dplyr::mutate(
      isHispanic = case_match(
        isHispanic,
        "Hispanic" ~ spec$isHispanic$true_val,
        c("6", "6.85") ~ spec$missing, # error in CMC file for two rows
        .default = isHispanic
      )
    )
}


# Harmonize SZBD data
#
# Harmonize using the generic PEC function but also fix some SZBD-specific
# issues in the data:
#   * Update all `race` and `isHispanic` values to conform to the data dictionary
harmonize_PEC_SZBD <- function(metadata, spec) {
  harmonize_PEC(metadata, spec) |>
    dplyr::mutate(
      race = case_match(
        race,
        "w" ~ spec$race$White,
        "u" ~ spec$missing,
        .default = race
      ),
      isHispanic = ifelse(isHispanic %in% c("nr", "u"),
                          spec$missing, isHispanic)
    )
}
