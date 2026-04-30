# Harmonize AMP-PD / GEN-A3 metadata
#
# Modifies the AMP-PD donor metadata file to conform to the GENESIS data
# dictionary. Source metadata files: local download provided by Jaroslav Bendl
# on Aug 18, 2025 and the AMP-PD portal clinical Demographics.csv file (version
# 2023_v4release_1027).
#
# Modifications needed for version Aug 18, 2025 / 2023_v4release_1027:
#   * Rename columns:
#     * `participant_id` => `individualID`
#     * `Info.Institution` => `dataContributionGroup`
#     * `Demographics.age_at_baseline` => `ageDeath`
#     * `Demographics.sex` => `sex`
#     * `ethnicity` => `isHispanic`
#     * `Demographics.ethnicity` => geneticAncestry`
#     * `PMI.PMI_hours` => `PMI`
#     * `LBD_Cohort_Path_Data.path_braak_nft` => `Braak_NFT`
#     * `LBD_Cohort_Path_Data.path_cerad` => `amyCerad`
#     * `Info.BraakGroup` => `bScore_LB`
#   * Update `race` values to conform to the data dictionary
#   * Change `sex` values to all lower case
#   * Add a `cohort` column based on dataContributionGroup values
#   * Change `dataContributionGroup` values to conform to data dictionary
#   * Harmonize the `Info.Diagnosis` field into binary `PD` and `Control` fields
#
# Note on diagnosis:
#   There are several possible diagnosis fields: Info.Diagnosis,
#   PD_Medical_History.diagnosis, and PD_Medical_History.most_recent_diagnosis.
#   The first two agree with each other, but PD_Medical_History.diagnosis has 4
#   samples with values of "Dementia With Lewy Bodies" or "Parkinsonism", and
#   Info.Diagnosis groups these in with "PD".
#   PD_Medical_History.most_recent_diagnosis appears to be missing a lot of data
#   that exists in the other two fields. Given that, I've chosen to use
#   Info.Diagnosis as the ground truth field for the harmonized PD column.
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
harmonize_AMP_PD <- function(metadata, spec) {
  metadata |>
    dplyr::rename(
      individualID = participant_id,
      dataContributionGroup = Info.Institution,
      ageDeath = Demographics.age_at_baseline,
      sex = Demographics.sex,
      isHispanic = ethnicity,
      geneticAncestry = Demographics.ethnicity,
      PMI = PMI.PMI_hours,
      Braak_NFT = LBD_Cohort_Path_Data.path_braak_nft,
      Braak_LB = LBD_Cohort_Path_Data.path_braak_lb,
      amyCerad = LBD_Cohort_Path_Data.path_cerad,
      bScore_LB = Info.BraakGroup
    ) |>
    dplyr::mutate(
      ageDeath = censor_ages(ageDeath, spec),
      isHispanic = case_match(
        isHispanic,
        "Unknown" ~ spec$missing,
        "Hispanic or Latino" ~ spec$isHispanic$true_val,
        "Not Hispanic or Latino" ~ spec$isHispanic$false_val,
        .default = isHispanic
      ),
      race = ifelse(race == "Unknown", spec$missing, race),
      sex = tolower(sex),
      Braak_NFT = to_Braak_stage(Braak_NFT, spec),
      Braak_LB = to_Braak_stage(Braak_LB, spec),
      # For lack of better information, cohort and dataContributionGroup are
      # effectively the same
      cohort = case_match(
        dataContributionGroup,
        "HA" ~ spec$cohort$harvard,
        "MS" ~ spec$cohort$msbb,
        "UD" ~ spec$cohort$udall,
        "UM" ~ spec$cohort$umiami,
        .default = dataContributionGroup
      ),
      dataContributionGroup = case_match(
        dataContributionGroup,
        "HA" ~ spec$dataContributionGroup$harvard,
        "MS" ~ spec$dataContributionGroup$mssm,
        "UD" ~ spec$dataContributionGroup$udall,
        "UM" ~ spec$dataContributionGroup$umiami,
        .default = dataContributionGroup
      ),
      # Binary diagnosis columns
      PD = make_binary_column(Info.Diagnosis, "PD", spec),
      Control = make_binary_column(Info.Diagnosis, "Control", spec)
    )
}
