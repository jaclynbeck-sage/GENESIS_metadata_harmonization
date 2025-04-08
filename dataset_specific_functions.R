# This file contains one function per GENESIS data set that performs all
# harmonization operations. Specific changes to each data set are documented
# in the comments of the corresponding functions. These functions are tied to
# specific versions of each source file and may need to be updated if the source
# metadata file is updated.

# Harmonize Diverse Cohorts metadata
#
# Makes minor edits to the Diverse Cohorts individual metadata file, which was
# previously harmonized using a very similar data dictionary to what is needed
# for GENESIS. Source metadata file: syn51757646 (version 20) on Synapse.
#
# Modifications needed for version 20:
#   * Change `ageDeath` and `PMI` value "missing or unknown" to `NA`
#   * Change `PMI` to a numeric column
#   * Rename `isHispanic` values from ["TRUE", "FALSE"] to ["True", "False"]
#   * Rename `race` value "Black" to "Black or African American"
#   * Fixes some `amyA` values that don't match their `amyThal` values
#   * Add `apoeStatus` and `study` column
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
    dplyr::mutate(
      ageDeath = censor_ages(ageDeath, spec),
      PMI = case_match(PMI,
        spec$missing ~ NA,
        .default = suppressWarnings(as.numeric(PMI))
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
      ## Manual corrections
      amyCerad = case_when(
        individualID %in% cerad_fix_ids ~ spec$amyCerad$frequent,
        .default = amyCerad
      ),
      ##
      # Update based on manual corrections
      amyAny = get_amyAny(amyCerad, spec),
      study = spec$study$diverse_cohorts
    )
}


# Harmonize MayoRNAseq metadata
#
# Modifies the original MayoRNAseq individual metadata to conform to the GENESIS
# data dictionary. This file is not used directly by any GENESIS studies but is
# instead used to fill in missing information in GENESIS metadata that uses
# these samples. Source metadata file: syn23277389 (version 7) on Synapse.
#
# Modifications needed for version 7:
#   * Rename columns:
#     * `pmi` => `PMI`
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
#   * Add columns `dataContributionGroup` = "Mayo", `cohort` = "Mayo Clinic",
#     and `study` = "MayoRNAseq"
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
harmonize_MayoRNAseq <- function(metadata, spec) {
  metadata |>
    rename(
      PMI = pmi,
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
      ## Manual correction
      race = case_when(
        individualID == "11387" ~ spec$race$other,
        .default = race
      ),
      ##
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
      dataContributionGroup = spec$dataContributionGroup$mayo,
      study = spec$study$mayo
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
#     * `pmi` => `PMI`
#     * `ethnicity` => `isHispanic`
#     * `CERAD` => `amyCerad`
#     * `individualIdSource` => `dataContributionGroup`
#   * Convert `PMI` from minutes to hours
#   * Change `isHispanic` and `race` values to conform to the data dictionary
#   * Change `NA` values in the `isHispanic`, `race`, `apoeGenotype`,
#       `amyCerad`, and `Braak` columns to "missing or unknown"
#   * Add `apoe4Status`, `amyAny`, `amyA`, and `bScore` columns
#   * Add `cohort` = "Mt Sinai Brain Bank" and `study` = "MSBB"
#
# Arguments:
#   metadata - a `data.frame` of metadata from the source metadata file. Columns
#     are variables and rows are individuals.
#   spec - a `config` object describing the standardized values for each field,
#     as defined by this project's `GENESIS_harmonization.yml` file
#
# NOTE: There is one individual with incorrect `race` and `isHispanic` values
# and two individuals with incorrect `Braak` values. These values are corrected
# manually here based on Diverse Cohorts data.
#
# Returns:
#   a `data.frame` with all relevant fields harmonized to the GENESIS data
#   dictionary. Columns not defined in the data dictionary are left as-is.
#
harmonize_MSBB <- function(metadata, spec) {
  metadata |>
    rename(
      PMI = pmi,
      isHispanic = ethnicity,
      amyCerad = CERAD,
      dataContributionGroup = individualIdSource
    ) |>
    mutate(
      PMI = PMI / 60, # PMI is in minutes
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
      ## Manual corrections
      race = case_when(
        individualID == "AMPAD_MSSM_0000036634" ~ spec$race$other,
        .default = race
      ),
      isHispanic = case_when(
        individualID == "AMPAD_MSSM_0000036634" ~ spec$isHispanic$hisp_true,
        .default = isHispanic
      ),
      ##
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
      ## Manual corrections
      Braak = case_when(
        individualID == "AMPAD_MSSM_0000013078" ~ to_Braak_stage(4, spec),
        individualID == "AMPAD_MSSM_0000048247" ~ to_Braak_stage(6, spec),
        .default = Braak
      ),
      ##
      bScore = get_bScore(Braak, spec),
      cohort = spec$cohort$msbb,
      study = spec$study$msbb
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
#     * `pmi` => `PMI`
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
#   * Add `dataContributionGroup` = "Rush" and `study` = "ROSMAP"
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
      PMI = pmi,
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
      dataContributionGroup = spec$dataContributionGroup$rush,
      study = spec$study$rosmap
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
#     * `pmi` => `PMI`
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
#   * Add `cohort` ("SEA-AD"), `dataContributionGroup` ("Allen Institute")
#       and `study` ("SEA-AD")
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
      PMI = pmi,
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
      dataContributionGroup = spec$dataContributionGroup$allen_institute,
      study = spec$study$sea_ad
    )
}


# Harmonize NPS-AD / GEN-A1 metadata
#
# Modifies the NPS-AD individual metadata file to conform to the GENESIS data
# dictionary. A separate neuropathology file needs to be merged with the
# individual metadata before harmonization to obtain Braak scores. Source
# metadata files: syn55251012 (version 4, individual metadata) and syn55251003
# (version 1, neuropathology data) on Synapse.
#
# Note: Cerad mapping is defined in the NPS-AD data dictionary (syn57373364).
#
# Note: Ages are listed as "89+" in v4 of NPS-AD. All of these values should be
# replaced with "90+" or "89", which can mostly be filled in from Diverse
# Cohorts / AMP-AD 1.0 data. There are several individuals listed as "89+" with
# known discrepancies to the DC / 1.0 data, which are manually corrected here.
# The age column will be fixed in a future version of the NPS-AD metadata.
#
# Note: Even after de-duplication, most harmonized columns (ageDeath, PMI, sex,
# race, isHispanic, apoeGenotype, amyCerad, Braak) that have values that
# disagree with AMP-AD 1.0 and Diverse Cohorts data. NPS-AD determined these
# values algorithmically or by re-processing data, or assumes these are
# corrections to original metadata, and wishes these values to remain as-is, so
# we do not alter them.
#
# Modifications needed for version 4:
#   * Merge neuropathology data with individual metadata
#   * Rename columns:
#     * `ethnicity` => `isHispanic`
#     * `CERAD` => `amyCerad`
#     * `BRAAK_AD` => `Braak`
#   * Fix two values in `diverseCohortsIndividualIDFormat` to match the correct
#     format in Diverse Cohorts
#   * Manually correct `ageDeath` column as described above
#   * Convert `PMI` values from minutes to hours
#   * Convert `isHispanic` values to "True" or "False"
#   * Convert Braak numerical values to Roman numerals
#   * Convert `amyCerad` numerical values to values in data dictionary
#   * Convert `NA` values to "missing or unknown" in the `race`, `isHispanic`,
#     `apoeGenotype`, `Braak`, and `amyCerad` columns
#   * Add `apoe4Status`, `amyAny`, `amyThal`, `amyA`, and `bScore` columns
#   * Fix `cohort` values to match data dictionary
#   * Add the `dataContributionGroup` column with values appropriate to each
#     cohort.
#   * Use Diverse Cohorts and AMP-AD 1.0 data to get the correct `cohort` value
#     for samples marked as "ROSMAP"
#   * Add `study` = "NPS-AD"
#
# Arguments:
#   metadata - a `data.frame` of metadata from the source metadata file. Columns
#     are variables and rows are individuals.
#   neuropath - a `data.frame` of neuropathology data for each individual, which
#     can be matched to `metadata` by `individualID`. Columns are variables and
#     rows are individuals.
#   harmonized_baseline - a `data.frame` of de-duplicated and harmonized
#     metadata from all AMP-AD 1.0 studies and Diverse Cohorts
#   spec - a `config` object describing the standardized values for each field,
#     as defined by this project's `GENESIS_harmonization.yml` file
#
# Returns:
#   a `data.frame` with all relevant fields harmonized to the GENESIS data
#   dictionary. Columns not defined in the data dictionary are left as-is.
#
harmonize_NPS_AD <- function(metadata, neuropath, harmonized_baseline, spec) {
  metadata <- metadata |>
    # Braak is all NA in the individual metadata file but has values in the
    # neuropath file
    select(-Braak) |>
    merge(neuropath)

  meta_new <- metadata |>
    select(-Component) |>
    rename(
      isHispanic = ethnicity,
      amyCerad = CERAD,
      Braak = BRAAK_AD
    ) |>
    mutate(
      # Fix to allow comparison to Diverse Cohorts
      diverseCohortsIndividualIDFormat = case_match(diverseCohortsIndividualIDFormat,
        29637 ~ "29637_MSSM",
        29582 ~ "29582_MSSM",
        .default = as.character(diverseCohortsIndividualIDFormat)
      ),
      ## Manual corrections
      ageDeath = case_match(ageDeath,
        # Setting these to NA will force these values to be filled in by Diverse
        # Cohorts / AMP-AD 1.0 data during de-duplication
        "89+" ~ NA,
        .default = ageDeath),
      ageDeath = case_match(individualID,
        # Exceptions to the above: These three individuals have a different age
        # in NPS-AD than Diverse Cohorts / 1.0 data that should not be over-
        # written
        c("AMPAD_MSSM_0000011938", "AMPAD_MSSM_0000056009") ~ spec$over90,
        "AMPAD_MSSM_0000077061" ~ "89",
        .default = ageDeath
      ),
      ##
      PMI = PMI / 60,
      pmiUnits = "hours",
      race = case_match(race,
        NA ~ spec$missing,
        .default = race
      ),
      isHispanic = case_match(isHispanic,
        "Hispanic or Latino" ~ spec$isHispanic$hisp_true,
        "Not Hispanic or Latino" ~ spec$isHispanic$hisp_false,
        NA ~ spec$missing,
        .default = isHispanic
      ),
      apoeGenotype = case_match(apoeGenotype,
        NA ~ spec$missing,
        .default = as.character(apoeGenotype)
      ),
      apoe4Status = get_apoe4Status(apoeGenotype, spec),
      # Cerad mapping per NPS-AD data dictionary (syn57373364)
      amyCerad = case_match(amyCerad,
        1 ~ spec$amyCerad$none,
        2 ~ spec$amyCerad$sparse,
        3 ~ spec$amyCerad$moderate,
        4 ~ spec$amyCerad$frequent,
        NA ~ spec$missing,
        .default = as.character(amyCerad)
      ),
      amyAny = get_amyAny(amyCerad, spec),
      amyThal = spec$missing,
      amyA = spec$missing,
      Braak = to_Braak_stage(Braak, spec),
      bScore = get_bScore(Braak, spec),
      cohort = case_match(cohort,
        "MSBB" ~ spec$cohort$msbb,
        .default = cohort
      ),
      dataContributionGroup = case_match(cohort,
        spec$cohort$msbb ~ spec$dataContributionGroup$mssm,
        spec$cohort$hbcc ~ spec$dataContributionGroup$nimh,
        c(spec$cohort$ros, spec$cohort$map) ~ spec$dataContributionGroup$rush,
        "ROSMAP" ~ spec$dataContributionGroup$rush,
        .default = ""
      ),
      study = spec$study$nps_ad
    )

  # Pull missing information from AMP-AD 1.0 and Diverse Cohorts metadata,
  # including correct cohort information
  meta_new <- meta_new |>
    mutate(source = "GEN-A1")

  meta_new <- deduplicate_studies(
    list(meta_new, harmonized_baseline),
    spec,
    # NPS-AD did their own re-processing of these values, which should not be
    # over-written by any other data set.
    exclude_cols = c("amyCerad", "amyAny", "Braak", "bScore", "study"),
    verbose = FALSE
  ) |>
    subset(source == "GEN-A1") |>
    select(all_of(colnames(meta_new)), -source)

  return(meta_new)
}


# Harmonize SMIB-AD / GEN-A9 metadata
#
# Modifies the SMIB-AD individual metadata file to conform to the GENESIS data
# dictionary. Source metadata file: syn22432749 (version 1) on Synapse. Note:
# CERAD values are coded using ROSMAP's system.
#
# Modifications needed for version 1:
#   * Rename columns:
#     * `pmi` => `PMI`
#     * `ethnicity` => `isHispanic`
#     * `CERAD` => `amyCerad`
#   * Censor `ageDeath` values over 90
#   * Change `race` values from "European" to "White"
#   * Change `isHispanic` values from "European" to "False"
#   * Convert Braak numerical values to Roman numerals
#   * Convert `amyCerad` numerical values to values in data dictionary
#   * Add `apoe4Status`, `amyAny`, `amyThal`, `amyA`, and `bScore` columns
#   * Add a `cohort` column with either "SMRI" or "Banner"
#   * Add a `dataContributionGroup` column with values for Stanley and Banner
#   * Add `study` = "SMIB_AD"
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
harmonize_SMIB_AD <- function(metadata, spec) {
  metadata |>
    rename(
      PMI = pmi,
      isHispanic = ethnicity,
      amyCerad = CERAD
    ) |>
    mutate(
      ageDeath = censor_ages(ageDeath, spec),
      race = spec$race$White,
      isHispanic = spec$isHispanic$hisp_false,
      apoeGenotype = as.character(apoeGenotype),
      apoe4Status = get_apoe4Status(apoeGenotype, spec),
      Braak = to_Braak_stage(Braak, spec),
      bScore = get_bScore(Braak, spec),
      amyCerad = case_match(amyCerad,
        1 ~ spec$amyCerad$none,
        2 ~ spec$amyCerad$sparse,
        4 ~ spec$amyCerad$frequent,
        NA ~ spec$missing,
        .default = as.character(amyCerad)
      ),
      amyAny = get_amyAny(amyCerad, spec),
      amyThal = spec$missing,
      amyA = spec$missing,
      cohort = case_match(individualIdSource,
        NA ~ spec$cohort$smri,
        .default = spec$cohort$banner
      ),
      dataContributionGroup = case_match(cohort,
        spec$cohort$smri ~ spec$dataContributionGroup$stanley,
        .default = spec$dataContributionGroup$banner
      ),
      study = spec$study$smib_ad
    )
}


# Harmonize MCMPS / GEN-A10 metadata
#
# Modifies the MCMPS individual metadata file to conform to the GENESIS data
# dictionary. Source metadata file: syn25891193 (version 1) on Synapse. Note:
# these samples come from living tissue and do not have a value for `ageDeath`
# or `PMI`.
#
# Modifications needed for version 1:
#   * Rename columns:
#     * `pmi` => `PMI`
#     * `ethnicity` => `isHispanic`
#     * `CERAD` => `amyCerad`
#   * Trim extra spaces from `race` values
#   * Trim extra spaces from `isHispanic` values and convert to True/False
#   * Replace `NA` values with "missing or unknown" in the `apoeGenotype`,
#     `amyCerad`, and `Braak` columns
#   * Add `apoe4Status`, `amyAny`, `amyThal`, `amyA`, and `bScore` columns
#   * Add a `cohort` column containing "Mayo Clinic"
#   * Add a `dataContributionGroup` column containing "Mayo"
#   * Add `study` = "MCMPS"
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
harmonize_MCMPS <- function(metadata, spec) {
  metadata |>
    rename(
      PMI = pmi,
      isHispanic = ethnicity,
      amyCerad = CERAD
    ) |>
    mutate(
      race = str_trim(race),
      isHispanic = str_trim(isHispanic),
      isHispanic = case_match(isHispanic,
        "Hispanic or Latino" ~ spec$isHispanic$hisp_true,
        "Not Hispanic or Latino" ~ spec$isHispanic$hisp_false,
        "Middle Eastern" ~ spec$isHispanic$hisp_false,
        .default = isHispanic
      ),
      apoeGenotype = case_match(apoeGenotype,
        NA ~ spec$missing,
        .default = as.character(apoeGenotype)
      ),
      apoe4Status = get_apoe4Status(apoeGenotype, spec),
      amyCerad = spec$missing,
      amyAny = spec$missing,
      amyThal = spec$missing,
      amyA = spec$missing,
      Braak = spec$missing,
      bScore = spec$missing,
      cohort = spec$cohort$mayo,
      dataContributionGroup = spec$dataContributionGroup$mayo,
      study = spec$study$mcmps
    )
}


# Harmonize MC_snRNA / GEN-A11 metadata
#
# Modifies the MC_snRNA individual metadata file to conform to the GENESIS data
# dictionary. Source metadata file: syn31563038 (version 1) on Synapse. This
# data set has some sample overlap with AMP-AD 1.0 Mayo metadata and Diverse
# Cohorts metadata. There are some missing values in this data set that exist in
# these data sets, so we pull those values in.
#
# Modifications needed for version 1:
#   * Rename columns:
#     * `pmi` => `PMI`
#     * `ethnicity` => `isHispanic`
#     * `CERAD` => `amyCerad`
#   * Change `ageDeath` values of "90_or_over" to "90+"
#   * Change `isHispanic` values from "Caucasian" to "False"
#   * Replace `NA` values with "missing or unknown" in the `amyCerad` column
#   * Round numerical `Braak` values down and convert to Roman numerals
#   * Add `apoe4Status`, `amyAny`, `amyThal`, `amyA`, and `bScore` columns
#   * Add a `cohort` column containing "Mayo Clinic"
#   * Add a `dataContributionGroup` column containing "Mayo"
#   * Add `study` = "MC_snRNA
#
# Arguments:
#   metadata - a `data.frame` of metadata from the source metadata file. Columns
#     are variables and rows are individuals.
#   harmonized_baseline - a `data.frame` of de-duplicated and harmonized
#     metadata from all AMP-AD 1.0 studies and Diverse Cohorts
#   spec - a `config` object describing the standardized values for each field,
#     as defined by this project's `GENESIS_harmonization.yml` file
#
# Returns:
#   a `data.frame` with all relevant fields harmonized to the GENESIS data
#   dictionary. Columns not defined in the data dictionary are left as-is.
#
harmonize_MC_snRNA <- function(metadata, harmonized_baseline, spec) {
  meta_new <- metadata |>
    rename(
      PMI = pmi,
      isHispanic = ethnicity,
      amyCerad = CERAD
    ) |>
    mutate(
      ageDeath = censor_ages(ageDeath, spec),
      isHispanic = spec$isHispanic$hisp_false,
      apoeGenotype = as.character(apoeGenotype),
      apoe4Status = get_apoe4Status(apoeGenotype, spec),
      amyCerad = spec$missing,
      amyAny = spec$missing,
      amyThal = spec$missing,
      amyA = spec$missing,
      Braak = to_Braak_stage(floor(Braak), spec),
      bScore = get_bScore(Braak, spec),
      dataContributionGroup = spec$dataContributionGroup$mayo,
      cohort = spec$cohort$mayo,
      study = spec$study$mc_snrna
    )

  # Pull missing information from AMP-AD 1.0 and Diverse Cohorts metadata
  meta_new <- meta_new |>
    mutate(source = "GEN-A11")

  meta_new <- deduplicate_studies(
    list(meta_new, harmonized_baseline),
    spec,
    verbose = FALSE
  ) |>
    subset(source == "GEN-A11") |>
    select(all_of(colnames(meta_new)), -source)

  return(meta_new)
}


# Harmonize MC-BrAD / GEN-A12 metadata
#
# Modifies the MC-BrAD individual metadata file to conform to the GENESIS data
# dictionary. Source metadata file: syn51401700 (version 2) on Synapse. This
# data set has some sample overlap with AMP-AD 1.0 Mayo metadata and Diverse
# Cohorts metadata. There are some missing values in this data set that exist in
# these data sets, so we pull those values in.
#
# Modifications needed for version 2:
#   * Rename columns:
#     * `pmi` => `PMI`
#     * `ethnicity` => `isHispanic`
#     * `CERAD` => `amyCerad`
#     * `Thal` => `amyThal`
#   * Change `ageDeath` values of "90_or_over" to "90+"
#   * Replace `NA` values with "missing or unknown" in the `isHispanic`,
#   * `apoeGenotype`, `amyCerad`, and `amyThal` columns
#   * Round numerical `Braak` values down and convert to Roman numerals
#   * Convert `amyThal` values to conform to data dictionary
#   * Add `apoe4Status`, `amyAny`, `amyA`, and `bScore` columns
#   * Add a `cohort` column containing "Mayo Clinic"
#   * Add a `dataContributionGroup` column containing "Mayo"
#   * Add `study` = "MC-BrAD"
#
# Arguments:
#   metadata - a `data.frame` of metadata from the source metadata file. Columns
#     are variables and rows are individuals.
#   harmonized_baseline - a `data.frame` of de-duplicated and harmonized
#     metadata from all AMP-AD 1.0 studies and Diverse Cohorts
#   spec - a `config` object describing the standardized values for each field,
#     as defined by this project's `GENESIS_harmonization.yml` file
#
# Returns:
#   a `data.frame` with all relevant fields harmonized to the GENESIS data
#   dictionary. Columns not defined in the data dictionary are left as-is.
#
harmonize_MC_BrAD <- function(metadata, harmonized_baseline, spec) {
  meta_new <- metadata |>
    rename(
      PMI = pmi,
      isHispanic = ethnicity,
      amyCerad = CERAD,
      amyThal = Thal
    ) |>
    mutate(
      ageDeath = censor_ages(ageDeath, spec),
      isHispanic = spec$missing,
      apoeGenotype = case_match(apoeGenotype,
        NA ~ spec$missing,
        .default = as.character(apoeGenotype)
      ),
      apoe4Status = get_apoe4Status(apoeGenotype, spec),
      amyCerad = spec$missing,
      amyAny = spec$missing,
      amyThal = case_match(amyThal,
        0 ~ spec$amyThal$none,
        1 ~ spec$amyThal$phase1,
        NA ~ spec$missing,
        .default = as.character(amyThal)
      ),
      amyA = get_amyA(amyThal, spec),
      Braak = to_Braak_stage(floor(Braak), spec),
      bScore = get_bScore(Braak, spec),
      dataContributionGroup = spec$dataContributionGroup$mayo,
      cohort = spec$cohort$mayo,
      study = spec$study$mc_brad
    )

  # Pull missing information from AMP-AD 1.0 and Diverse Cohorts metadata
  meta_new <- deduplicate_studies(
    list(meta_new, harmonized_baseline),
    spec,
    verbose = FALSE
  ) |>
    subset(study == spec$study$mc_brad) |>
    select(all_of(colnames(meta_new)))

  return(meta_new)
}
