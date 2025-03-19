# Harmonize Diverse Cohorts metadata
#
# Makes minor edits to the Diverse Cohorts individual metadata file, which was
# previously harmonized using a very similar data dictionary to what is needed
# for GENESIS. Source metadata file: syn51757646 (version 20) on Synapse.
#
# Modifications needed for version 20:
#   * Rename `PMI` column to `pmi`
#   * Change `ageDeath` and `pmi` value `missing or unknown` to `NA`
#   * Change `pmi` to a numeric column
#   * Rename `isHispanic` values from [`TRUE`, `FALSE`] to [`True`, `False`]
#   * Rename `race` value `Black` to `Black or African American`
#   * Add `apoeStatus` column
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
harmonize_Diverse_Cohorts <- function(metadata, spec) {
  metadata |>
    dplyr::rename(pmi = PMI) |>
    dplyr::mutate(
      ageDeath = case_when(
        ageDeath == spec$missing ~ NA,
        .default = ageDeath
      ),
      pmi = case_when(
        pmi == spec$missing ~ NA,
        .default = suppressWarnings(as.numeric(pmi))
      ),
      isHispanic = case_when(
        isHispanic == "TRUE" ~ spec$isHispanic$hisp_true,
        isHispanic == "FALSE" ~ spec$isHispanic$hisp_false,
        .default = isHispanic
      ),
      race = case_when(
        race == "Black" ~ spec$race$Black,
        .default = race
      ),
      apoe4Status = get_apoe4Status(apoeGenotype, spec)
    )
}


# Harmonize ROSMAP metadata
#
# Modifies the ROSMAP individual metadata file to conform to the GENESIS data
# dictionary. Source metadata file: syn3191087 (version 11) on Synapse. The
# clinical codebook describing the meanings of values in each column can be
# found at syn3191090 on Synapse.
#
# Modifications needed for version 11:
#   * Rename multiple columns:
#     * `spanish` => `isHispanic`
#     * `age_death` => `ageDeath`
#     * `msex` => `sex`
#     * `apoe_genotype` => `apoeGenotype`
#     * `ceradsc` => `amyCerad`
#     * `braaksc` => `Braak`
#     * `Study` => `cohort`
#   * Convert `ageDeath` empty string values to `NA`
#   * Convert `NA` values in the `sex`, `isHispanic`, `race`, `apoeGenotype`,
#       `Braak`, and `amyCerad` columns to `missing or unknown`
#   * Convert `sex` values [0, 1] to [`female`, `male`]
#   * Convert `race` numerical values to values in data dictionary
#     * NOTE: "Hawaiian / Pacific Islanders" are categorized as "other" for this
#       effort as there are only 2 of them and the GENESIS dictionary does not
#       have a separate category for them.
#   * Convert `isHispanic` values [1, 2] to [`True`, `False`]
#   * Convert `Braak` numerical values to "None" or "Stage " + Roman numeral
#   * Convert `amyCerad` numerical values to values in data dictionary
#   * Add `apoe4status`, `bScore`, and `amyAny` columns
#   * Add `amyThal` and `amyA` columns with all `missing or unknown` values
#   * Add `dataContributionGroup` with value = "Rush"
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
harmonize_ROSMAP <- function(metadata, spec) {
  metadata |>
    dplyr::rename(
      isHispanic = spanish,
      ageDeath = age_death,
      sex = msex,
      apoeGenotype = apoe_genotype,
      amyCerad = ceradsc,
      Braak = braaksc,
      cohort = Study
    ) |>
    dplyr::mutate(
      sex = case_when(
        sex == 1 ~ spec$sex$male,
        sex == 0 ~ spec$sex$female,
        .default = spec$missing
      ),
      ageDeath = censor_ages(ageDeath, spec),
      race = case_when(
        race == 1 ~ spec$race$White,
        race == 2 ~ spec$race$Black,
        race == 3 ~ spec$race$Amer_Ind,
        race == 4 ~ spec$race$other, # Hawaiian / Pacific Islanders
        race == 5 ~ spec$race$Asian,
        race == 6 ~ spec$race$other,
        race == 7 ~ spec$missing,
        .default = spec$missing
      ),
      isHispanic = case_when(
        isHispanic == 1 ~ spec$isHispanic$hisp_true,
        isHispanic == 2 ~ spec$isHispanic$hisp_false,
        .default = spec$missing
      ),
      apoeGenotype = case_when(
        is.na(apoeGenotype) ~ spec$missing,
        .default = as.character(apoeGenotype)
      ),
      apoe4Status = get_apoe4Status(apoeGenotype, spec),
      Braak = case_when(
        Braak >= 0 ~ to_Braak_stage(Braak, spec),
        .default = spec$missing
      ),
      bScore = get_bScore(Braak, spec),
      amyCerad = case_when(
        amyCerad == 1 ~ spec$amyCerad$frequent,
        amyCerad == 2 ~ spec$amyCerad$moderate,
        amyCerad == 3 ~ spec$amyCerad$sparse,
        amyCerad == 4 ~ spec$amyCerad$none,
        .default = spec$missing
      ),
      amyAny = get_amyAny(amyCerad, spec),
      amyThal = spec$missing,
      amyA = spec$missing,
      dataContributionGroup = spec$dataContributionGroup$rush
    )
}


# Harmonize SEA-AD metadata
#
# Modifies the SEA-AD individual metadata file to conform to the GENESIS data
# dictionary. The version of SEA-AD that is on Synapse is missing
# Hispanic/Latino information that is present in the version released by the
# Allen Institute on brain-map.org. We use the version on Synapse, because it
# has been curated and approved for release in the AD Knowledge Portal, but pull
# in the missing Hispanic/Latino information from the AI version.
#
# Source metadata files:
#   * syn31149116 (version 7) on Synapse
#   * the "Donor Metadata" file downloaded from https://portal.brain-map.org.
#     Full URL: https://brainmapportal-live-4cc80a57cd6e400d854-f7fdcae.divio-media.net/filer_public/b4/c7/b4c727e1-ede1-4c61-b2ee-bf1ae4a3ef68/sea-ad_cohort_donor_metadata_072524.xlsx
#
# Modifications needed for version 7:
#   * Rename several columns:
#     * `Hispanic.Latino` => `isHispanic`
#     * `CERAD` => `amyCerad`
#     * `Thal phase` => `amyThal`
#   * Convert `isHispanic` values [`Yes`, `No`] to [`True`, `False`]
#   * Individuals may have multiple races in the `race` column, separated by a
#       semi-colon. Any value with `American Indian` in it is assigned to
#       `American Indian or Alaska Native`.
#   * Convert `race` value `Other` to `other`
#   * Convert `amyCerad` numerical values to values in data dictionary
#   * Convert `Braak` numerical values to "None" or "Stage " + Roman numeral
#   * Convert `amyThal` values from `Thal #` to `None` or `Phase #`
#   * Convert `NA` values in the `isHispanic`, `race`, `apoeGenotype`, `Braak`,
#       `amyThal`, and `amyCerad` columns to `missing or unknown`
#   * Update `apoe4Status` column with values in data dictionary
#   * Add `bScore`, and `amyAny` columns
#   * Add `cohort` ("SEA-AD") and `dataContributionGroup` ("Allen Institute")
#       columns.
#
# Arguments:
#   metadata_synapse - a `data.frame` of metadata from the Synapse metadata
#     file. Columns are variables and rows are individuals.
#   metadata_allen - a `data.frame` of metadata from brain-map.org. Columns are
#     variables and rows are individuals.
#   spec - a `config` object describing the standardized values for each field,
#     as defined by this project's `GENESIS_harmonization.yml` file
#
# Returns:
#   a `data.frame` with all relevant fields harmonized to the GENESIS data
#   dictionary. Columns not defined in the data dictionary are left as-is.
#
harmonize_SEA_AD <- function(metadata_synapse, metadata_allen, spec) {
  # Keep only the Hispanic/Latino column and IDs from the Allen Institute data
  metadata_allen <- metadata_allen %>%
    select("Donor ID", "Hispanic/Latino")

  meta <- merge(metadata_synapse, metadata_allen,
    by.x = "individualID",
    by.y = "Donor ID",
    all = TRUE
  )

  meta %>%
    dplyr::rename(
      isHispanic = "Hispanic/Latino",
      amyCerad = CERAD,
      amyThal = "Thal phase"
    ) %>%
    dplyr::mutate(
      #ageDeath = censor_ages(ageDeath, spec),
      isHispanic = case_when(
        isHispanic == "No" ~ spec$isHispanic$hisp_false,
        isHispanic == "Yes" ~ spec$isHispanic$hisp_true,
        .default = spec$missing
      ),
      race = case_when(
        race == "Other" ~ spec$race$other,
        grepl("American Indian", race) ~ spec$race$Amer_Ind,
        is.na(race) ~ spec$missing,
        .default = race
      ),
      apoeGenotype = case_when(
        is.na(apoeGenotype) ~ spec$missing,
        .default = as.character(apoeGenotype)
      ),
      apoe4Status = get_apoe4Status(apoeGenotype, spec),
      amyCerad = case_when(
        amyCerad == 0 ~ spec$amyCerad$none,
        amyCerad == 1 ~ spec$amyCerad$sparse,
        amyCerad == 2 ~ spec$amyCerad$moderate,
        amyCerad == 3 ~ spec$amyCerad$frequent,
        .default = spec$missing
      ),
      amyAny = get_amyAny(amyCerad, spec),
      amyThal = case_when(
        is.na(amyThal) ~ spec$missing,
        amyThal == "Thal 0" ~ spec$amyThal$none,
        .default = str_replace(amyThal, "Thal", "Phase")
      ),
      amyA = get_amyA(amyThal, spec),
      Braak = case_when(
        is.na(Braak) ~ spec$missing,
        .default = to_Braak_stage(Braak, spec)
      ),
      bScore = get_bScore(Braak, spec),
      cohort = spec$cohort$sea_ad,
      dataContributionGroup = spec$dataContributionGroup$allen_institute
    )
}
