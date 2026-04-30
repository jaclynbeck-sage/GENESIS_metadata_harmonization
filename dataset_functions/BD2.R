# Harmonize BD2 / GEN-B8 metadata
#
# Modifies the BD2 donor metadata file to conform to the GENESIS data
# dictionary. There are 14 samples that overlap with Diverse Cohorts/AMP-AD 1.0,
# so missing Braak/amyCerad/amyThal values for those samples are filled in from
# that data in the main harmonization function. Source metadata file: local
# download provided by Jaroslav Bendl on July 16, 2025.
#
# Modifications needed for version July 16, 2025:
#   * Rename columns:
#     * `SubNum` => `individualID`
#     * `Brain_bank` => `dataContributionGroup`
#     * `Age` => `ageDeath`
#     * `Sex` => `sex`
#     * `Race` => `race`
#   * Update `race` values to conform to the data dictionary
#   * Change `sex` values to all lower case
#   * Samples with `race` = "Hispanic" changed to `race` = "Other",
#     `isHispanic` = "True"
#   * Add a `cohort` column with either "Mount Sinai Brain Bank" or "UPitt"
#   * Change `dataContributionGroup` "UPitt" to "University of Pittsburgh"
#   * Harmonize the `Dx` column into a binary `BD` column
#   * Harmonize the `ClinConPrimDxText` into a binary `SCZ` column
#   * Create a binary `Control` column based on `Dx` and `ClinConPrimDxText`
#
# Arguments:
#   metadata - a `data.frame` of metadata from the source metadata file. Columns
#     are variables and rows are individuals.
#   spec - a `config` object describing the standardized values for each field,
#     as defined by this project's `GENESIS_harmonization.yml` file
#
# Returns:
#   a `data.frame` with all relevant fields harmonized to the GENESIS data
#   dictionary. Columns not defined in the data dictionary are left as-is.
#
harmonize_BD2 <- function(metadata, spec) {
  meta_new <- metadata |>
    dplyr::rename(
      individualID = SubNum,
      dataContributionGroup = Brain_bank,
      ageDeath = Age,
      sex = Sex,
      race = Race
    ) |>
    dplyr::mutate(
      # Shouldn't need censoring but just in case
      ageDeath = censor_ages(ageDeath, spec),
      isHispanic = ifelse(race == "Hispanic",
                          spec$isHispanic$true_val,
                          spec$isHispanic$false_val),
      race = case_match(
        race,
        "Black" ~ spec$race$Black,
        "Hispanic" ~ spec$race$other,
        .default = race
      ),
      sex = tolower(sex),
      cohort = case_match(
        dataContributionGroup,
        "MSSM" ~ spec$cohort$msbb,
        "UPitt" ~ spec$cohort$upitt,
        .default = dataContributionGroup
      ),
      dataContributionGroup = case_match(
        dataContributionGroup,
        "MSSM" ~ spec$dataContributionGroup$mssm,
        "UPitt" ~ spec$dataContributionGroup$upitt,
        .default = dataContributionGroup
      ),
      # Not enough Braak or Cerad values to create an AD diagnosis column.
      # There are no missing values for Dx. Most of the values in
      # ClinConPrimDxText and ClinConSecDxText are missing/NA, so we only use
      # their non-missing values to pull out diagnoses for SCZ, MDD, and OTHER
      # and assume missing values mean no diagnosis for those diseases.
      BD = make_binary_column(Dx, "BD", spec),
      SCZ = ifelse(grepl("Schizophrenia", ClinConPrimDxText),
                   1, 0),
      MDD = ifelse(grepl("Major depressive disorder", ClinConPrimDxText) |
                     grepl("Major depressive disorder", ClinConSecDxText),
                   1, 0),
      Other = ifelse(grepl("Schizoaffective", ClinConPrimDxText),
                     1, 0),
      Control = ifelse(Dx == "Control" & is.na(ClinConPrimDxText) &
                         is.na(ClinConSecDxText),
                       1, 0)
      # TODO this study has a "BD_type" column
    )

  # Modification for de-duplication in the main harmonization function
  meta_new <- meta_new |>
    mutate(
      bd2_id = individualID,
      individualID = case_when(
        nchar(Synapse_GENESIS) > 0 ~ Synapse_GENESIS,
        nchar(Synapse_MSSM) > 0 ~ Synapse_MSSM,
        .default = as.character(individualID)
      )
    )

  return(meta_new)
}
