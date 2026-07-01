# Harmonize BD2 / GEN-B8 metadata
#
# Modifies the BD2 donor metadata file to conform to the GENESIS data
# dictionary. There are 14 samples that overlap with Diverse Cohorts/AMP-AD 1.0,
# so missing Braak/amyCerad/amyThal values for those samples are filled in from
# that data in the main harmonization function. Source metadata file: local
# download provided by Jaroslav Bendl on June 29, 2026.
#
# Modifications needed for version June 29, 2026:
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
#   * Create binary `MS` and `SCZ` columns based on `ClinConPrimDxText`
#   * Create binary `MDD` column based on `ClinConPrimDxText` and `ClinConSecDxText`
#   * Create binary `Other` column based on `ClinConPrimDxText`
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
  metadata |>
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
      isHispanic = case_match(
        race,
        "Hispanic" ~ spec$isHispanic$true_val,
        NA ~ spec$missing,
        .default = spec$isHispanic$false_val
      ),
      race = case_match(
        race,
        "Black" ~ spec$race$Black,
        "Hispanic" ~ spec$race$other,
        "Asiansian" ~ spec$race$Asian, # Fix typo in source data
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
      SCZ = grep_to_binary_column(ClinConPrimDxText, "Schizophrenia"),

      # Major depressive disorder might be in either column so we paste them
      # together and pattern match
      MDD = grep_to_binary_column(paste(ClinConPrimDxText, ClinConSecDxText),
                                  "Major depressive disorder"),
      MS = grep_to_binary_column(ClinConPrimDxText, "Multiple sclerosis"),
      Other = grep_to_binary_column(ClinConPrimDxText,
                                    "Schizoaffective|Encephalopathy")
      # TODO this study has a "BD_type" column
    )
}
