# Harmonize Diverse Cohorts metadata
#
# Makes minor edits to the Diverse Cohorts individual metadata file, which was
# previously harmonized using a very similar data dictionary to what is needed
# for GENESIS. Source metadata file: syn51757646 (version 20) on Synapse.
#
# Modifications needed for version 20:
#   * Rename `PMI` column to `pmi`
#   * Change `ageDeath` and `pmi` value "missing or unknown" to `NA`
#   * Change `pmi` to a numeric column
#   * Rename `isHispanic` values from ["TRUE", "FALSE"] to ["True", "False"]
#   * Rename `race` value "Black" to "Black or African American"
#   * Fixes some `amyA` values that don't match their `amyThal` values
#   * Add `apoeStatus` column
#
# NOTE: There are 8 individuals in v20 of the data with incorrect `amyCerad`
# values, which are corrected manually here based on comparison with ROSMAP
# metadata.
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
  cerad_fix_ids <- c(
    "R5309342", "R1753622", "R7767382", "R3519397", "R5307279",
    "R5553042", "R9594832", "R1951848"
  )

  metadata |>
    dplyr::rename(pmi = PMI) |>
    dplyr::mutate(
      ageDeath = censor_ages(ageDeath, spec),
      pmi = case_match(pmi,
        spec$missing ~ NA,
        .default = suppressWarnings(as.numeric(pmi))
      ),
      isHispanic = case_match(isHispanic,
        "TRUE" ~ spec$isHispanic$hisp_true,
        "FALSE" ~ spec$isHispanic$hisp_false,
        .default = isHispanic
      ),
      race = case_match(race,
        "Black" ~ spec$race$Black,
        .default = race
      ),
      apoe4Status = get_apoe4Status(apoeGenotype, spec),
      amyA = get_amyA(amyThal, spec),
      # Manual correction
      amyCerad = case_when(
        individualID %in% cerad_fix_ids ~ spec$amyCerad$frequent,
        .default = amyCerad
      ),
      # Update based on manual corrections
      amyAny = get_amyAny(amyCerad, spec)
    )
}


# Harmonize MayoRNASeq metadata
#
# Modifies the original MayoRNASeq individual metadata to conform to the GENESIS
# data dictionary. This file is not used directly by any GENESIS studies but is
# instead used to fill in missing information in GENESIS metadata that uses
# these samples. Source metadata file: syn23277389 (version 7) on Synapse.
#
# Modifications needed for version 7:
#   * Rename columns:
#     * `ethnicity` => `isHispanic`
#     * `CERAD` => `amyCerad`
#     * `Thal` => `amyThal`
#   * Change `ageDeath` value "90_or_over" to "90+"
#   * Rename `isHispanic` value "Caucasian" to "False"
#   * Replace `amyThal` value "0" with "None", and add "Phase" in front of other
#       `amyThal` numerical values
#   * Convert `Braak` value "0" with "None", and convert other numerical `Braak`
#       values to "Stage " + a Roman numeral
#   * Add `apoeStatus`, `amyA`, `amyAny`, and `bScore` columns
#   * Replace NA values in the `isHispanic`, `race`, `sex`, `apoeGenotype`,
#       `amyCerad`, `amyThal`, and `Braak` columns with "missing or unknown"
#   * Add columns `dataContributionGroup` = "Mayo" and `cohort` = "Mayo Clinic"
#
# NOTE: There is one individual with an incorrect `race` value, which is
# corrected manually here based on updated information from Diverse Cohorts.
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
harmonize_MayoRNASeq <- function(metadata, spec) {
  metadata |>
    rename(
      isHispanic = ethnicity,
      amyCerad = CERAD,
      amyThal = Thal
    ) |>
    mutate(
      ageDeath = censor_ages(ageDeath, spec),
      isHispanic = case_match(isHispanic,
        "Caucasian" ~ spec$isHispanic$hisp_false,
        NA ~ spec$missing,
        .default = isHispanic
      ),
      race = case_match(race,
        NA ~ spec$missing,
        .default = race
      ),
      # Manual correction
      race = case_when(
        individualID == "11387" ~ spec$race$other,
        .default = race
      ),
      sex = case_match(sex,
        NA ~ spec$missing,
        .default = sex
      ),
      apoeGenotype = case_match(apoeGenotype,
        NA ~ spec$missing,
        .default = as.character(apoeGenotype)
      ),
      apoe4Status = get_apoe4Status(apoeGenotype, spec),
      amyCerad = case_match(amyCerad,
        NA ~ spec$missing, # All values are NA in v7
        .default = amyCerad
      ),
      amyAny = get_amyAny(amyCerad, spec),
      amyThal = case_match(amyThal,
        NA ~ spec$missing,
        0 ~ spec$amyThal$none,
        .default = paste("Phase", amyThal)
      ),
      amyA = get_amyA(amyThal, spec),
      Braak = to_Braak_stage(floor(Braak), spec),
      bScore = get_bScore(Braak, spec),
      cohort = spec$cohort$mayo,
      dataContributionGroup = spec$dataContributionGroup$mayo
    )
}


# Harmonize MSBB metadata
#
# Modifies the original MSBB individual metadata to conform to the GENESIS
# data dictionary. This file is not used directly by any GENESIS studies but is
# instead used to fill in missing information in GENESIS metadata that uses
# these samples. Source metadata file: syn6101474 (version 9) on Synapse.
#
# Modifications needed for version 9:
#   * Rename columns:
#     * `ethnicity` => `isHispanic`
#     * `CERAD` => `amyCerad`
#     * `individualIdSource` => `dataContributionGroup`
#   * Convert `pmi` from minutes to hours
#   * Change `isHispanic` and `race` values to conform to the data dictionary
#   * Change `NA` values in the `isHispanic`, `race`, `apoeGenotype`,
#       `amyCerad`, and `Braak` columns to "missing or unknown"
#   * Add `apoe4Status`, `amyAny`, `amyA`, and `bScore` columns
#   * Add `cohort` = "Mt Sinai Brain Bank" column
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
harmonize_MSBB <- function(metadata, spec) {
  metadata |>
    rename(
      isHispanic = ethnicity,
      amyCerad = CERAD,
      dataContributionGroup = individualIdSource
    ) |>
    mutate(
      pmi = pmi / 60, # PMI is in minutes
      isHispanic = case_match(isHispanic,
        NA ~ spec$missing,
        c("A", "B", "W") ~ spec$isHispanic$hisp_false,
        "H" ~ spec$isHispanic$hisp_true,
        "U" ~ spec$missing,
        .default = isHispanic
      ),
      race = case_match(race,
        NA ~ spec$missing,
        "A" ~ spec$race$Asian,
        "B" ~ spec$race$Black,
        "H" ~ spec$race$other,
        "W" ~ spec$race$White,
        "U" ~ spec$missing,
        .default = race
      ),
      apoeGenotype = case_match(apoeGenotype,
        NA ~ spec$missing,
        .default = as.character(apoeGenotype)
      ),
      apoe4Status = get_apoe4Status(apoeGenotype, spec),
      amyCerad = case_match(amyCerad,
        NA ~ spec$missing,
        1 ~ spec$amyCerad$none,
        2 ~ spec$amyCerad$frequent,
        3 ~ spec$amyCerad$moderate,
        4 ~ spec$amyCerad$sparse,
        .default = as.character(amyCerad)
      ),
      amyAny = get_amyAny(amyCerad, spec),
      amyThal = spec$missing,
      amyA = get_amyA(amyThal, spec),
      Braak = to_Braak_stage(floor(Braak), spec),
      bScore = get_bScore(Braak, spec),
      cohort = spec$cohort$msbb
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
#       `Braak`, and `amyCerad` columns to "missing or unknown"
#   * Convert `sex` values [0, 1] to ["female", "male"]
#   * Convert `race` numerical values to values in data dictionary
#     * NOTE: "Hawaiian / Pacific Islanders" are categorized as "other" for this
#       effort as there are only 2 of them and the GENESIS dictionary does not
#       have a separate category for them.
#   * Convert `isHispanic` values [1, 2] to ["True", "False"]
#   * Convert `Braak` numerical values to "None" or "Stage " + Roman numeral
#   * Convert `amyCerad` numerical values to values in data dictionary
#   * Add `apoe4status`, `bScore`, `amyA`, and `amyAny` columns
#   * Add `amyThal` and `amyA` columns with all "missing or unknown" values
#   * Add `dataContributionGroup` with value = "Rush"
#
# NOTE: There is one individual with an incorrect value for `isHispanic`, which
# is manually corrected here based on updated data from Diverse Cohorts.
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
      sex = case_match(sex,
        1 ~ spec$sex$male,
        0 ~ spec$sex$female,
        NA ~ spec$missing,
        .default = as.character(sex)
      ),
      ageDeath = censor_ages(ageDeath, spec),
      race = case_match(race,
        1 ~ spec$race$White,
        2 ~ spec$race$Black,
        3 ~ spec$race$Amer_Ind,
        4 ~ spec$race$other, # Hawaiian / Pacific Islanders
        5 ~ spec$race$Asian,
        6 ~ spec$race$other,
        7 ~ spec$missing,
        NA ~ spec$missing,
        .default = as.character(race)
      ),
      isHispanic = case_match(isHispanic,
        1 ~ spec$isHispanic$hisp_true,
        2 ~ spec$isHispanic$hisp_false,
        NA ~ spec$missing,
        .default = as.character(isHispanic)
      ),
      # manual correction
      isHispanic = case_when(
        individualID == "R8412417" ~ spec$isHispanic$hisp_false,
        .default = isHispanic
      ),
      apoeGenotype = case_match(apoeGenotype,
        NA ~ spec$missing,
        .default = as.character(apoeGenotype)
      ),
      apoe4Status = get_apoe4Status(apoeGenotype, spec),
      Braak = to_Braak_stage(Braak, spec),
      bScore = get_bScore(Braak, spec),
      amyCerad = case_match(amyCerad,
        NA ~ spec$missing,
        1 ~ spec$amyCerad$frequent,
        2 ~ spec$amyCerad$moderate,
        3 ~ spec$amyCerad$sparse,
        4 ~ spec$amyCerad$none,
        .default = as.character(amyCerad)
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
#   * Convert `isHispanic` values ["Yes", "No"] to ["True", "False"]
#   * Individuals may have multiple races in the `race` column, separated by a
#       semi-colon. Any value with "American Indian" in it is assigned to
#       "American Indian or Alaska Native".
#   * Convert `race` value "Other" to "other"
#   * Convert `amyCerad` numerical values to values in data dictionary
#   * Convert `Braak` numerical values to "None" or "Stage " + Roman numeral
#   * Convert `amyThal` values from "Thal #" to "None" or "Phase #"
#   * Convert `NA` values in the `isHispanic`, `race`, `apoeGenotype`, `Braak`,
#       `amyThal`, and `amyCerad` columns to "missing or unknown"
#   * Update `apoe4Status` column with values in data dictionary
#   * Add `bScore`, `amyA`, and `amyAny` columns
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
  metadata_allen <- metadata_allen |>
    select("Donor ID", "Hispanic/Latino")

  meta <- merge(metadata_synapse, metadata_allen,
    by.x = "individualID",
    by.y = "Donor ID",
    all = TRUE
  )

  meta |>
    dplyr::rename(
      isHispanic = "Hispanic/Latino",
      amyCerad = CERAD,
      amyThal = "Thal phase"
    ) |>
    dplyr::mutate(
      isHispanic = case_match(isHispanic,
        "No" ~ spec$isHispanic$hisp_false,
        "Yes" ~ spec$isHispanic$hisp_true,
        "Unknown" ~ spec$missing,
        NA ~ spec$missing,
        .default = isHispanic
      ),
      race = case_when(
        race == "Other" ~ spec$race$other,
        grepl("American Indian", race) ~ spec$race$Amer_Ind,
        is.na(race) ~ spec$missing,
        .default = race
      ),
      apoeGenotype = case_match(apoeGenotype,
        NA ~ spec$missing,
        .default = as.character(apoeGenotype)
      ),
      apoe4Status = get_apoe4Status(apoeGenotype, spec),
      amyCerad = case_match(amyCerad,
        0 ~ spec$amyCerad$none,
        1 ~ spec$amyCerad$sparse,
        2 ~ spec$amyCerad$moderate,
        3 ~ spec$amyCerad$frequent,
        NA ~ spec$missing,
        .default = as.character(amyCerad)
      ),
      amyAny = get_amyAny(amyCerad, spec),
      amyThal = case_match(amyThal,
        NA ~ spec$missing,
        "Thal 0" ~ spec$amyThal$none,
        .default = str_replace(amyThal, "Thal", "Phase")
      ),
      amyA = get_amyA(amyThal, spec),
      Braak = to_Braak_stage(Braak, spec),
      bScore = get_bScore(Braak, spec),
      cohort = spec$cohort$sea_ad,
      dataContributionGroup = spec$dataContributionGroup$allen_institute
    )
}
