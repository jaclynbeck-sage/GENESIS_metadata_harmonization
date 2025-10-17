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
harmonize <- function(study_name, metadata, spec, harmonized_baseline = NULL,
                      extra_metadata = NULL) {
  # Study-specific harmonization
  metadata <- switch(
    study_name,
    "AMP-PD" = harmonize_AMP_PD(metadata, spec),
    "ASAP" = harmonize_ASAP(metadata, spec),
    "BD2" = harmonize_BD2(metadata, spec),
    "AMP-AD_DiverseCohorts" = harmonize_Diverse_Cohorts(metadata, spec),
    "MayoRNAseq" = harmonize_MayoRNAseq(metadata, spec),
    "MC-BrAD" = harmonize_MC_BrAD(metadata, spec),
    "MC_snRNA" = harmonize_MC_snRNA(metadata, spec),
    "McCarroll_SCZ" = harmonize_McCarroll_SCZ(metadata, spec),
    "MCMPS" = harmonize_MCMPS(metadata, spec),
    "MSBB" = harmonize_MSBB(metadata, spec),
    "NPS-AD" = harmonize_NPS_AD(metadata, extra_metadata, spec),
    "ROSMAP" = harmonize_ROSMAP(metadata, spec),
    "SEA-AD" = harmonize_SEA_AD(metadata, extra_metadata, spec),
    "SMIB-AD" = harmonize_SMIB_AD(metadata, spec),
    "MSBB_corrections" = harmonize_MSBB_corrections(metadata, spec), # TODO temporary
    .default = metadata
  )

  # Add any missing fields
  missing_fields <- setdiff(expectedColumns, colnames(metadata))
  for (field in missing_fields) {
    metadata[, field] <- spec$missing
  }

  # Don't fill NA values with "missing" in the ageDeath or PMI columns
  cols_fill <- setdiff(expectedColumns, c("ageDeath", "PMI"))

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
    select(all_of(expectedColumns), !all_of(expectedColumns))

  return(metadata)
}


# Harmonize Diverse Cohorts metadata
#
# Makes minor edits to the Diverse Cohorts individual metadata file, which was
# previously harmonized using a very similar data dictionary to what is needed
# for GENESIS. Source metadata file: syn51757646 (version 21) on Synapse.
#
# Modifications needed for version 21:
#   * Rename columns:
#     * `Braak` => `Braak_NFT`
#     * `bScore` => `bScore_NFT`
#   * Change `ageDeath` and `PMI` value "missing or unknown" to `NA`
#   * Change `PMI` to a numeric column
#   * Rename `isHispanic` values from ["TRUE", "FALSE"] to ["True", "False"]
#   * Update `dataContributionGroup` and `cohort` to conform to the data dictionary
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
    dplyr::rename(
      Braak_NFT = Braak,
      bScore_NFT = bScore
    ) |>
    dplyr::mutate(
      ageDeath = censor_ages(ageDeath, spec),
      PMI = ifelse(PMI == spec$missing, NA,
                   suppressWarnings(as.numeric(PMI))),
      isHispanic = case_match(
        as.character(isHispanic),
        "TRUE" ~ spec$isHispanic$hisp_true,
        "FALSE" ~ spec$isHispanic$hisp_false,
        .default = as.character(isHispanic)
      ),
      dataContributionGroup = case_match(
        dataContributionGroup,
        "Columbia" ~ spec$dataContributionGroup$columbia,
        "MSSM" ~ spec$dataContributionGroup$mssm,
        "Rush" ~ spec$dataContributionGroup$rush,
        "Emory" ~ spec$dataContributionGroup$emory,
        "Mayo" ~ spec$dataContributionGroup$mayo,
        .default = dataContributionGroup
      ),
      cohort = case_match(
        cohort,
        "Banner" ~ spec$cohort$banner,
        "Mt Sinai Brain Bank" ~ spec$cohort$msbb,
        .default = cohort
      )
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
#     * `Braak` => `Braak_NFT`
#   * Change `ageDeath` value "90_or_over" to "90+"
#   * Rename `isHispanic` value "Caucasian" to "False"
#   * Replace `amyThal` value "0" with "None", and add "Phase" in front of other
#       `amyThal` numerical values
#   * Convert `Braak_NFT` value "0" with "None", and convert other numerical
#       `Braak_NFT` values to "Stage " + a Roman numeral
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
harmonize_MayoRNAseq <- function(metadata, spec) {
  metadata |>
    rename(
      PMI = pmi,
      isHispanic = ethnicity,
      amyCerad = CERAD,
      amyThal = Thal,
      Braak_NFT = Braak
    ) |>
    mutate(
      ageDeath = censor_ages(ageDeath, spec),
      isHispanic = ifelse(isHispanic == "Caucasian",
                          spec$isHispanic$hisp_false,
                          isHispanic),
      ## Manual correction
      race = ifelse(individualID == "11387", spec$race$other, race),
      ##
      # Note: all amyCerad values are NA in this data set and will get filled
      # in with "missing or unknown" in the main harmonize() function
      amyThal = case_match(amyThal,
        NA ~ spec$missing,
        0 ~ spec$amyThal$none,
        .default = paste("Phase", amyThal)
      ),
      Braak_NFT = to_Braak_stage(floor(Braak_NFT), spec),
      cohort = spec$cohort$mayo,
      dataContributionGroup = spec$dataContributionGroup$mayo,
    )
}


# Harmonize MSBB metadata
#
# Modifies the original MSBB individual metadata to conform to the GENESIS
# data dictionary. This file is not used directly by any GENESIS studies but is
# instead used to fill in missing information in GENESIS metadata that uses
# these samples. Source metadata file: syn6101474 (version 10) on Synapse.
#
# Modifications needed for version 10:
#   * Rename columns:
#     * `pmi` => `PMI`
#     * `ethnicity` => `isHispanic`
#     * `CERAD` => `amyCerad`
#     * `Braak` => `Braak_NFT`
#   * Convert `PMI` from minutes to hours
#   * Change `isHispanic` and `race` values to conform to the data dictionary
#   * Add `dataContributionGroup` = "Mount Sinai School of Medicine",
#     `cohort` = "Mount Sinai Brain Bank", and `study` = "MSBB"
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
      PMI = pmi,
      isHispanic = ethnicity,
      amyCerad = CERAD,
      Braak_NFT = Braak
    ) |>
    mutate(
      PMI = PMI / 60, # PMI is in minutes
      isHispanic = case_match(isHispanic,
        c("A", "B", "O", "W") ~ spec$isHispanic$hisp_false,
        "H" ~ spec$isHispanic$hisp_true,
        "U" ~ spec$missing,
        .default = isHispanic
      ),
      race = case_match(race,
        "A" ~ spec$race$Asian,
        "B" ~ spec$race$Black,
        "H" ~ spec$race$other,
        "W" ~ spec$race$White,
        "O" ~ spec$race$other,
        "U" ~ spec$missing,
        .default = race
      ),
      amyCerad = case_match(amyCerad,
        1 ~ spec$amyCerad$none,
        2 ~ spec$amyCerad$frequent,
        3 ~ spec$amyCerad$moderate,
        4 ~ spec$amyCerad$sparse,
        .default = as.character(amyCerad)
      ),
      amyThal = spec$missing,
      Braak_NFT = to_Braak_stage(floor(Braak_NFT), spec),
      dataContributionGroup = spec$dataContributionGroup$mssm,
      cohort = spec$cohort$msbb
    )
}


# Note: individual "6160" is "White" in NPS-AD, MSBB 1.0, and Diverse cohorts,
# but "other" in MSBB corrections. Because of the consensus of the older data
# with NPS-AD, in this case we ignore the corrected MSBB data and leave the race
# as "White".
harmonize_MSBB_corrections <- function(metadata, spec) {
  metadata |>
    rename(
      ageDeath = Age,
      race = RaceLabel,
      sex = SexLabel,
      PMI = `PMI (min)`,
      amyCerad = CERAD_1,
      apoeGenotype = ApoE,
      Braak_NFT = `B&B Alz`
    ) |>
    mutate(
      ageDeath = censor_ages(ageDeath, spec),
      PMI = PMI / 60, # PMI is in minutes
      sex = tolower(sex),
      race = str_replace(race, " \\(nonHispanic\\)", ""),
      isHispanic = case_match(race,
        c("Asian", "Black", "White", "Other") ~ spec$isHispanic$hisp_false,
        "Hispanic" ~ spec$isHispanic$hisp_true,
        .default = race
      ),
      race = case_match(race,
        "Asian" ~ spec$race$Asian,
        "Black" ~ spec$race$Black,
        "Hispanic" ~ spec$race$other,
        "White" ~ spec$race$White,
        "Other" ~ spec$race$other,
        .default = race
      ),
      ## Manual change as noted above
      race = ifelse(individualID == "6160", "White", race),
      ##
      apoeGenotype = str_replace(apoeGenotype, "/", ""),
      amyCerad = case_match(amyCerad,
        0 ~ spec$amyCerad$none,
        1 ~ spec$amyCerad$sparse,
        2 ~ spec$amyCerad$moderate,
        3 ~ spec$amyCerad$frequent,
        .default = as.character(amyCerad)
      ),
      amyThal = spec$missing,
      Braak_NFT = to_Braak_stage(Braak_NFT, spec),
      cohort = spec$cohort$msbb,
      dataContributionGroup = spec$dataContributionGroup$mssm
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
#     * `braaksc` => `Braak_NFT`
#     * `Study` => `cohort`
#   * Convert `ageDeath` empty string values to `NA`
#   * Convert `sex` values [0, 1] to ["female", "male"]
#   * Convert `race` numerical values to values in data dictionary
#     * NOTE: "Hawaiian / Pacific Islanders" are categorized as "other" for this
#       effort as there are only 2 of them and the GENESIS dictionary does not
#       have a separate category for them.
#   * Convert `isHispanic` values [1, 2] to ["True", "False"]
#   * Convert `Braak_NFT` numerical values to "None" or "Stage " + Roman numeral
#   * Convert `amyCerad` numerical values to values in data dictionary
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
      Braak_NFT = braaksc,
      cohort = Study
    ) |>
    dplyr::mutate(
      sex = case_match(sex,
        1 ~ spec$sex$male,
        0 ~ spec$sex$female,
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
        .default = as.character(race)
      ),
      isHispanic = case_match(isHispanic,
        1 ~ spec$isHispanic$hisp_true,
        2 ~ spec$isHispanic$hisp_false,
        .default = as.character(isHispanic)
      ),
      # manual correction
      isHispanic = ifelse(individualID == "R8412417",
                          spec$isHispanic$hisp_false,
                          isHispanic),
      apoeGenotype = ifelse(is.na(apoeGenotype), spec$missing,
                            as.character(apoeGenotype)),
      Braak_NFT = to_Braak_stage(Braak_NFT, spec),
      amyCerad = case_match(amyCerad,
        1 ~ spec$amyCerad$frequent,
        2 ~ spec$amyCerad$moderate,
        3 ~ spec$amyCerad$sparse,
        4 ~ spec$amyCerad$none,
        .default = as.character(amyCerad)
      ),
      dataContributionGroup = spec$dataContributionGroup$rush,
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
#   * syn31149116 (version 8) on Synapse
#   * the "Donor Metadata" file downloaded from https://portal.brain-map.org.
#     Full URL: https://brainmapportal-live-4cc80a57cd6e400d854-f7fdcae.divio-media.net/filer_public/b4/c7/b4c727e1-ede1-4c61-b2ee-bf1ae4a3ef68/sea-ad_cohort_donor_metadata_072524.xlsx
#
# Modifications needed for version 8:
#   * Rename several columns:
#     * `pmi` => `PMI`
#     * `Hispanic.Latino` => `isHispanic`
#     * `CERAD` => `amyCerad`
#     * `Thal phase` => `amyThal`
#     * `Braak` => `Braak NFT`
#   * Convert `isHispanic` values ["Yes", "No"] to ["True", "False"]
#   * Individuals may have multiple races in the `race` column, separated by a
#       semi-colon. Any value with "American Indian" in it is assigned to
#       "American Indian or Alaska Native".
#   * Convert `race` value "Other" to "other"
#   * Convert `amyCerad` numerical values to values in data dictionary
#   * Convert `Braak_NFT` numerical values to "None" or "Stage " + Roman numeral
#   * Convert `amyThal` values from "Thal #" to "None" or "Phase #"
#   * Add `cohort` ("SEA-AD"), `dataContributionGroup` ("Allen Institute")
#       and `study` ("SEA-AD")
#   * Use Atherosclerosis values from the Allen metadata (even though these
#     aren't part of the harmonization), because v8 of the Synapse metadata
#     incorrectly combines "None" and NA values.
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
    select("Donor ID", "Hispanic/Latino", "Atherosclerosis")

  metadata_synapse <- metadata_synapse |> select(-Atherosclerosis) # v8 fix

  # v8 fix: match column names of v7
  colnames(metadata_synapse) <- colnames(metadata_synapse) |>
    str_replace_all("\\.", " ")

  meta <- merge(metadata_synapse, metadata_allen,
    by.x = "individualID",
    by.y = "Donor ID",
    all = TRUE
  ) |>
    dplyr::relocate(Atherosclerosis, .before = Arteriolosclerosis) # v8 fix

  # v8 fix: Set all empty strings to NA to match previous versions
  meta[meta == ""] <- NA

  meta |>
    dplyr::rename(
      PMI = pmi,
      isHispanic = "Hispanic/Latino",
      amyCerad = CERAD,
      amyThal = "Thal phase",
      Braak_NFT = Braak,
      `LATE-NC stage` = "LATE NC stage" # v8 fix to match v7
    ) |>
    dplyr::mutate(
      isHispanic = case_match(isHispanic,
        "No" ~ spec$isHispanic$hisp_false,
        "Yes" ~ spec$isHispanic$hisp_true,
        "Unknown" ~ spec$missing,
        .default = isHispanic
      ),
      race = case_when(
        race == "Other" ~ spec$race$other,
        grepl("American Indian", race) ~ spec$race$Amer_Ind,
        .default = race
      ),
      amyCerad = case_match(amyCerad,
        0 ~ spec$amyCerad$none,
        1 ~ spec$amyCerad$sparse,
        2 ~ spec$amyCerad$moderate,
        3 ~ spec$amyCerad$frequent,
        .default = as.character(amyCerad)
      ),
      amyThal = ifelse(amyThal == "Thal 0",
                       spec$amyThal$none,
                       str_replace(amyThal, "Thal", "Phase")),
      Braak_NFT = to_Braak_stage(Braak_NFT, spec),
      cohort = spec$cohort$sea_ad,
      dataContributionGroup = spec$dataContributionGroup$allen_institute,
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
# replaced with "90+" or "89", which can be filled in from Diverse Cohorts,
# AMP-AD 1.0 data, or the MSBB corrections data. The age column will be fixed in
# a future version of the NPS-AD metadata.
#
# Note: Even after de-duplication, the "race" and "isHispanic" columns will have
# values that disagree with AMP-AD 1.0 and Diverse Cohorts data. NPS-AD
# determined these values algorithmically or by re-processing data and wishes
# these values to remain as-is, so we do not alter them.
#
# Modifications needed for version 4:
#   * Merge neuropathology data with individual metadata
#   * Rename columns:
#     * `ethnicity` => `isHispanic`
#     * `CERAD` => `amyCerad`
#     * `BRAAK_AD` => `Braak_NFT`
#     * `race` => `geneticAncestry` # TODO
#   * Fix two values in `diverseCohortsIndividualIDFormat` to match the correct
#     format in Diverse Cohorts
#   * Manually correct `ageDeath` column as described above
#   * Convert `PMI` values from minutes to hours
#   * Convert `isHispanic` values to "True" or "False"
#   * Convert Braak_NFT numerical values to Roman numerals
#   * Convert `amyCerad` numerical values to values in data dictionary
#   * Fix `cohort` values to match data dictionary
#   * Add the `dataContributionGroup` column with values appropriate to each
#     cohort.
#
# Arguments:
#   metadata - a `data.frame` of metadata from the source metadata file. Columns
#     are variables and rows are individuals.
#   neuropath - a `data.frame` of neuropathology data for each individual, which
#     can be matched to `metadata` by `individualID`. Columns are variables and
#     rows are individuals.
#   spec - a `config` object describing the standardized values for each field,
#     as defined by this project's `GENESIS_harmonization.yml` file
#
# Returns:
#   a `data.frame` with all relevant fields harmonized to the GENESIS data
#   dictionary. Columns not defined in the data dictionary are left as-is.
#
harmonize_NPS_AD <- function(metadata, neuropath, spec) {
  metadata <- metadata |>
    # Braak is all NA in the individual metadata file but has values in the
    # neuropath file
    select(-Braak) |>
    merge(neuropath)

  metadata |>
    select(-Component) |>
    rename(
      #geneticAncestry = race, # TODO
      isHispanic = ethnicity,
      amyCerad = CERAD,
      Braak_NFT = BRAAK_AD
    ) |>
    mutate(
      # Fix to allow comparison to Diverse Cohorts
      diverseCohortsIndividualIDFormat = case_match(
        diverseCohortsIndividualIDFormat,
        29637 ~ "29637_MSSM",
        29582 ~ "29582_MSSM",
        .default = as.character(diverseCohortsIndividualIDFormat)
      ),
      ## Manual corrections
      # Setting these to NA will force these values to be filled in by Diverse
      # Cohorts / AMP-AD 1.0 data during de-duplication
      ageDeath = ifelse(ageDeath == "89+", NA, ageDeath),
      ##
      PMI = PMI / 60,
      pmiUnits = "hours",
      isHispanic = case_match(isHispanic,
        "Hispanic or Latino" ~ spec$isHispanic$hisp_true,
        "Not Hispanic or Latino" ~ spec$isHispanic$hisp_false,
        .default = isHispanic
      ),
      # Cerad mapping per NPS-AD data dictionary (syn57373364)
      amyCerad = case_match(amyCerad,
        1 ~ spec$amyCerad$none,
        2 ~ spec$amyCerad$sparse,
        3 ~ spec$amyCerad$moderate,
        4 ~ spec$amyCerad$frequent,
        .default = as.character(amyCerad)
      ),
      Braak_NFT = to_Braak_stage(Braak_NFT, spec),
      cohort = ifelse(cohort == "MSBB", spec$cohort$msbb, cohort),
      dataContributionGroup = case_match(cohort,
        spec$cohort$msbb ~ spec$dataContributionGroup$mssm,
        spec$cohort$hbcc ~ spec$dataContributionGroup$nimh,
        c(spec$cohort$ros, spec$cohort$map) ~ spec$dataContributionGroup$rush,
        "ROSMAP" ~ spec$dataContributionGroup$rush,
        .default = ""
      )
    )
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
#     * `Braak` => `Braak_NFT`
#   * Censor `ageDeath` values over 90
#   * Change `race` values from "European" to "White"
#   * Change `isHispanic` values from "European" to "False"
#   * Convert Braak_NFT numerical values to Roman numerals
#   * Convert `amyCerad` numerical values to values in data dictionary
#   * Add a `cohort` column with either "SMRI" or "Banner"
#   * Add a `dataContributionGroup` column with values for Stanley and Banner
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
      amyCerad = CERAD,
      Braak_NFT = Braak
    ) |>
    mutate(
      ageDeath = censor_ages(ageDeath, spec),
      race = spec$race$White,
      isHispanic = spec$isHispanic$hisp_false,
      Braak_NFT = to_Braak_stage(Braak_NFT, spec),
      amyCerad = case_match(amyCerad,
        1 ~ spec$amyCerad$none,
        2 ~ spec$amyCerad$sparse,
        4 ~ spec$amyCerad$frequent,
        .default = as.character(amyCerad)
      ),
      cohort = ifelse(is.na(individualIdSource),
                      spec$cohort$smri,
                      spec$cohort$banner),
      dataContributionGroup = ifelse(cohort == spec$cohort$smri,
                                     spec$dataContributionGroup$stanley,
                                     spec$dataContributionGroup$banner)
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
#     * `Braak` => `Braak_NFT`
#   * Trim extra spaces from `race` values
#   * Trim extra spaces from `isHispanic` values and convert to True/False
#   * Add a `cohort` column containing "Mayo Clinic"
#   * Add a `dataContributionGroup` column containing "Mayo"
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
      amyCerad = CERAD,
      Braak_NFT = Braak
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
      cohort = spec$cohort$mayo,
      dataContributionGroup = spec$dataContributionGroup$mayo
    )
}


# Harmonize MC_snRNA / GEN-A11 metadata
#
# Modifies the MC_snRNA individual metadata file to conform to the GENESIS data
# dictionary. Source metadata file: syn31563038 (version 1) on Synapse. This
# data set has some sample overlap with AMP-AD 1.0 Mayo metadata and Diverse
# Cohorts metadata. There are some missing values in this data set that exist in
# these data sets, so those get pulled in in the main harmonization function.
#
# Modifications needed for version 1:
#   * Rename columns:
#     * `pmi` => `PMI`
#     * `ethnicity` => `isHispanic`
#     * `CERAD` => `amyCerad`
#     * `Braak` => `Braak_NFT`
#   * Change `ageDeath` values of "90_or_over" to "90+"
#   * Change `isHispanic` values from "Caucasian" to "False"
#   * Round numerical `Braak_NFT` values down and convert to Roman numerals
#   * Add a `cohort` column containing "Mayo Clinic"
#   * Add a `dataContributionGroup` column containing "Mayo"
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
harmonize_MC_snRNA <- function(metadata, spec) {
  metadata |>
    rename(
      PMI = pmi,
      isHispanic = ethnicity,
      amyCerad = CERAD,
      Braak_NFT = Braak
    ) |>
    mutate(
      ageDeath = censor_ages(ageDeath, spec),
      isHispanic = spec$isHispanic$hisp_false,
      Braak_NFT = to_Braak_stage(floor(Braak_NFT), spec),
      dataContributionGroup = spec$dataContributionGroup$mayo,
      cohort = spec$cohort$mayo
    )
}


# Harmonize MC-BrAD / GEN-A12 metadata
#
# Modifies the MC-BrAD individual metadata file to conform to the GENESIS data
# dictionary. Source metadata file: syn51401700 (version 2) on Synapse. This
# data set has some sample overlap with AMP-AD 1.0 Mayo metadata and Diverse
# Cohorts metadata. There are some missing values in this data set that exist in
# these data sets, so those get pulled in in the main harmonization function.
#
# Modifications needed for version 2:
#   * Rename columns:
#     * `pmi` => `PMI`
#     * `ethnicity` => `isHispanic`
#     * `CERAD` => `amyCerad`
#     * `Thal` => `amyThal`
#     * `Braak` => `Braak_NFT`
#   * Change `ageDeath` values of "90_or_over" to "90+"
#   * Round numerical `Braak_NFT` values down and convert to Roman numerals
#   * Convert `amyThal` values to conform to data dictionary
#   * Add a `cohort` column containing "Mayo Clinic"
#   * Add a `dataContributionGroup` column containing "Mayo"
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
harmonize_MC_BrAD <- function(metadata, spec) {
  metadata |>
    rename(
      PMI = pmi,
      isHispanic = ethnicity,
      amyCerad = CERAD,
      amyThal = Thal,
      Braak_NFT = Braak
    ) |>
    mutate(
      ageDeath = censor_ages(ageDeath, spec),
      amyThal = case_match(amyThal,
        0 ~ spec$amyThal$none,
        1 ~ spec$amyThal$phase1,
        .default = as.character(amyThal)
      ),
      Braak_NFT = to_Braak_stage(floor(Braak_NFT), spec),
      dataContributionGroup = spec$dataContributionGroup$mayo,
      cohort = spec$cohort$mayo,
    )
}


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
                          spec$isHispanic$hisp_true,
                          spec$isHispanic$hisp_false),
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
      )
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
        "Hispanic or Latino" ~ spec$isHispanic$hisp_true,
        "Not Hispanic or Latino" ~ spec$isHispanic$hisp_false,
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
      )
    )
}


# Harmonize ASAP / GEN-A15 metadata
#
# Modifies the ASAP subject and clinical metadata files to conform to the
# GENESIS data dictionary. Data downloaded via command line on Oct 02, 2025,
# from https://cloud.parkinsonsroadmap.org
#
# Download commands used:
#   dnastack use cloud.parkinsonsroadmap.org
#   dnastack collections query -c prod-postmortem-derived-brain-sequencing-collection "SELECT * FROM \"collections\".\"prod_postmortem_derived_brain_sequencing_collection\".\"cohort_subject\"" -o csv > "data/downloads/ASAP_subject-export-2025-10-02.csv"
#   dnastack collections query -c prod-postmortem-derived-brain-sequencing-collection "SELECT * FROM \"collections\".\"prod_postmortem_derived_brain_sequencing_collection\".\"cohort_clinpath\"" -o csv > "data/downloads/ASAP_clinpath-export-2025-10-02.csv"
#
# Modifications needed for version Oct 02, 2025:
#   * Rename columns:
#     * `age_at_death` => `ageDeath`
#     * `duration_pmi` => `PMI`
#     * `ethnicity` => `isHispanic`
#     * `path_braak_nft` => `Braak_NFT`
#     * `path_braak_asyn` => `Braak_LB`
#     * `path_cerad` => `amyCerad`
#     * `path_thal` => `amyThal`
#     * `apoe_e4_status` => `apoeGenotype`
#     * `biobank_name` => `cohort`
#   * Add an `individualID` column that uses `asap_dataset_id` values. The
#     column `asap_dataset_id` is left in the data set despite being a duplicate
#     of individualID due to its ubiquitous use in other ASAP metadata files.
#     Leaving the column named as-is makes it easier to merge with other files.
#   * Censor ages over 90
#   * Set empty string ("") values in harmonized columns to NA so deduplication works
#   * Change `sex` values to all lower case
#   * Change `isHispanic` "Not Reported" values to "missing or unknown"
#   * Update `Braak_NFT`, `amyCerad`, `amyThal`, and `cohort` values to conform
#     to the data dictionary.
#   * Add `dataContributionGroup` values based on the `asap_team_id` values.
#   * Deduplicate rows with the same individual, which are identical except for
#     some columns where one row is missing a value and one has the value, or
#     where two different groups typed slightly different values for the same
#     column
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
harmonize_ASAP <- function(metadata, spec) {
  meta_new <- metadata |>
    dplyr::rename(
      ageDeath = age_at_death,
      PMI = duration_pmi,
      isHispanic = ethnicity,
      Braak_NFT = path_braak_nft,
      Braak_LB = path_braak_asyn,
      amyCerad = path_cerad,
      amyThal = path_thal,
      apoeGenotype = apoe_e4_status,
      cohort = biobank_name
    ) |>
    mutate(
      # This data set uses "" instead of NA, need to set to NA for some
      # harmonization functions to work
      across(any_of(expectedColumns), ~ ifelse(.x == "", NA, .x)),
      # Keep the field named "asap_subject_id" but add individualID column
      individualID = asap_subject_id,
      ageDeath = censor_ages(ageDeath, spec),
      sex = tolower(sex),
      # There aren't any non-missing or non "Not Reported" values
      isHispanic = ifelse(isHispanic == "Not Reported",
                          spec$missing,
                          isHispanic),
      # Braak_NFT is already in Roman numerals except for "0" and ""/NA
      Braak_NFT = case_match(
        Braak_NFT,
        NA ~ spec$missing,
        "0" ~ spec$Braak_NFT$none,
        .default = paste("Stage", Braak_NFT)
      ),
      Braak_LB = to_Braak_stage(Braak_LB, spec),
      # It's unclear whether missing values mean "None" or "missing" so they
      # have been left as NA so they are set to "missing".
      amyCerad = case_match(
        amyCerad,
        "Sparse" ~ spec$amyCerad$sparse,
        "Moderate" ~ spec$amyCerad$moderate,
        "Frequent" ~ spec$amyCerad$frequent,
        .default = amyCerad
      ),
      amyThal = case_match(
        amyThal,
        NA ~ spec$missing,
        "4/5" ~ spec$amyThal$phase4, # Round down to phase 4
        c("0", "0.0") ~ spec$amyThal$none,
        .default = paste("Phase", suppressWarnings(as.numeric(amyThal)))
      ),
      cohort = case_match(
        cohort,
        c("BSHRI", "Banner Sun Health Research Institute") ~ spec$cohort$banner,
        "QSBB_UK" ~ spec$cohort$qsbb,
        .default = cohort
      ),
      dataContributionGroup = spec$dataContributionGroup$asap,
    )

  # De-duplicate individuals within this data set, *NOT* with AMP-AD 1.0 / DivCo.
  meta_new <- deduplicate_studies(
    list(meta_new), spec,
    include_cols = colnames(meta_new),
    exclude_cols = c("asap_dataset_id", "subject_id", "asap_team_id"),
    verbose = FALSE
    ) |>
    distinct() |>
    # For any fields that had different non-missing values in them, concatenate
    # the values together, separated by "; ". This happens when two different
    # groups contributed metadata from the same individual and typed something
    # slightly different for the same column.
    group_by(individualID) |>
    summarize(
      across(
        everything(),
        ~ ifelse(is.character(.x) && !all(is.na(.x)),
                 paste(unique(.x), collapse = "; "),
                 unique(.x)) # Leave numeric or NA values alone but keep them in the data frame
      )
    ) |>
    distinct() |>
    as.data.frame()

  return(meta_new)
}


# Harmonize McCarroll SCZ / GEN-A16 metadata
#
# Modifies the McCarroll SCZ metadata file to conform to the GENESIS data
# dictionary. Data downloaded from the NeMo Archive on Oct 08, 2025 from:
# https://data.nemoarchive.org/other/grant/broad/mccarroll/transcriptome/sncell/10x_v3.1/human/processed/other/SZvillage_donorMetadata.txt
#
# Modifications needed for version Oct 08, 2025:
# * Rename columns:
#   * `Donor` => `individualID`
#   * `Sex` => `sex`
#   * `Age` => `ageDeath`
# * Censor ages above 90
# * Change `sex` values to all lower case
# * Add `dataContributionGroup` = ??
# * Add `cohort` = ??
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
harmonize_McCarroll_SCZ <- function(metadata, spec) {
  metadata |>
    dplyr::rename(
      individualID = Donor,
      sex = Sex,
      ageDeath = Age
    ) |>
    mutate(
      ageDeath = censor_ages(ageDeath, spec),
      sex = tolower(sex)#,
      #dataContributionGroup = ??,
      #cohort = ??
    )
}


# Harmonize McCarroll HD / GEN-A17 metadata
#
# Modifies the McCarroll HD metadata file to conform to the GENESIS data
# dictionary. Data downloaded from the NeMo Archive on Oct 08, 2025 from:
# https://data.nemoarchive.org/other/grant/broad_mccarroll_HD_2024/mccarroll/transcriptome/sncell/10x_v3.1/human/processed/other/donor_metadata.txt
#
# Modifications needed for version Oct 08, 2025:
# * Rename columns:
#   * `SID` => `individualID`
#   * `Sex` => `sex`
#   * `Age` => `ageDeath`
# * Change `ageDeath` values of ">89" to "90+"
# * Change `sex` values to all lower case
# * Add `dataContributionGroup` = ??
# * Add `cohort` = ??
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
harmonize_McCarroll_HD <- function(metadata, spec) {
  metadata |>
    dplyr::rename(
      individualID = SID,
      sex = Sex,
      ageDeath = Age
    ) |>
    mutate(
      ageDeath = ifelse(ageDeath == ">89", spec$over90, ageDeath),
      sex = tolower(sex)#,
      #dataContributionGroup = ??,
      #cohort = ??
    )
}
