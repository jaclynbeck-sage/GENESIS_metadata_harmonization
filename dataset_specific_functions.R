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
    "MC-BrAD" = harmonize_MC_BrAD(metadata, spec),
    "MC_snRNA" = harmonize_MC_snRNA(metadata, spec),
    "McCarroll_SCZ" = harmonize_McCarroll_SCZ(metadata, spec),
    "McCarroll_HD" = harmonize_McCarroll_HD(metadata, spec),
    "MCMPS" = harmonize_MCMPS(metadata, spec),
    "MSBB" = harmonize_ADKP_studies(metadata, spec),
    "NPS-AD" = harmonize_NPS_AD(metadata, spec),
    "ROSMAP" = harmonize_ROSMAP(metadata, spec),
    "SEA-AD" = harmonize_SEA_AD(metadata, extra_metadata, spec),
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


# Harmonize ADKP Harmonization Project studies (generic)
#
# Harmonize metadata from the ADKP Harmonization Project, which contains
# metadata files that have been harmonized with each other and de-duplicated.
# The data dictionary for these files is very similar to GENESIS, so only minor
# changes are needed.
#
# This function works on every file from the ADKP Harmonization project, but
# does not handle diagnosis because it is not a harmonized field in that study.
# Each individual study should have its own function that calls this one, and
# then harmonizes diagnosis.
#
# Sources: all files in the syn73713763 folder on Synapse
#   * Rename columns:
#     * `Braak` => `Braak_NFT`
#     * `bScore` => `bScore_NFT`
#   * Alter column values:
#     * `race` "Other" -> "other"
#     * `race` "Native Hawaiian or Other Pacific Islander" -> "other".
#         * There are only 2 individuals with this value and the GENESIS
#           dictionary does not have a separate category for them.
#     * `cohort` "Mayo Clinic Brain Bank" -> "Mayo Clinic"
#     * `cohort` "University of Kentucky" -> "Mayo Clinic"
#
# Note: individual "6160" is "White" in NPS-AD, un-harmonized MSBB 1.0, and
# un-harmonized Diverse cohorts, but "other" in corrected MSBB data. Because of
# the consensus of the older data with NPS-AD, for GENESIS we ignore the
# corrected MSBB data and leave the race as "White".
#
# Arguments:
#   metadata - a `data.frame` of metadata from the source metadata file. Columns
#     are variables and rows are individuals.
#   spec - a `config` object describing the standardized values for each field,
#      as defined by this project's `GENESIS_harmonization.yml` file
#
# Returns:
#   a `data.frame` with all relevant fields harmonized to the GENESIS data
#   dictionary. Columns not defined in the data dictionary are left as-is.
harmonize_ADKP_studies <- function(metadata, spec) {
  metadata |>
    dplyr::rename(
      Braak_NFT = Braak,
      bScore_NFT = bScore
    ) |>
    mutate(
      race = ifelse(race %in% c("Other", "Native Hawaiian or Other Pacific Islander"),
                    spec$race$other,
                    race),
      cohort = ifelse(cohort %in% c("Mayo Clinic Brain Bank", "University of Kentucky"),
                      spec$cohort$mayo,
                      cohort),
      ### Manual fix to individual 6160 ###
      race = ifelse(cohort == spec$cohort$msbb & grepl("6160", individualID),
                    "White", race)
      ###
    )
}


# Harmonize Diverse Cohorts metadata
#
# Calls harmonize_ADKP_studies and then creates an `AD` diagnosis column that is
# based on the ADoutcome values. "AD" is "True", "Control" and "Other" are
# "False", and missing values are "missing or unknown".
#
# Source metadata: syn73713769 (version 2) on Synapse
harmonize_Diverse_Cohorts <- function(metadata, spec) {
  harmonize_ADKP_studies(metadata, spec) |>
    mutate(AD = case_match(ADoutcome,
                           "AD" ~ spec$true_val,
                           c("Control", "Other") ~ spec$false_val,
                           .default = spec$missing),
           OTHER = case_match(ADoutcome,
                              "Other" ~ spec$true_val,
                              c("Control", "AD") ~ spec$false_val,
                              .default = spec$missing),
           Control = case_match(ADoutcome,
                                "Control" ~ spec$true_val,
                                c("AD", "Other") ~ spec$false_val,
                                .default = spec$missing))
}


# Harmonize MC_snRNA data
#
# Calls harmonize_ADKP_studies and then creates binary `AD` and `Control`
# diagnosis columns that are based on the `diagnosis` column. This study has no
# non-missing amyCerad values and only a few non-missing Thal values, so
# diagnosis cannot be determined based on pathology. The diagnosis supplied in
# the `diagnosis` field is used instead.
#
# Source file: syn73713776 (version 2) on Synapse
harmonize_MC_snRNA <- function(metadata, spec) {
  harmonize_ADKP_studies(metadata, spec) |>
    # There are no missing values in the diagnosis field
    mutate(AD = ifelse(diagnosis == "Alzheimer Disease",
                       spec$true_val, spec$false_val),
           Control = ifelse(diagnosis == "Control",
                            spec$true_val, spec$false_val))
}

# Harmonize MC-BrAD data
#
# Calls harmonize_ADKP_studies and then creates binary `PSP` and `Control`
# diagnosis columns that are based on the `diagnosis` field. The study
# documentation verified that control samples have no diagnosis of any of the
# common neuropathological disorders like AD, PD, HD, etc.
#
# Source file: syn73713775 (version 2) on Synapse
harmonize_MC_BrAD <- function(metadata, spec) {
  harmonize_ADKP_studies(metadata, spec) |>
    # There are no missing values in the diagnosis column
    mutate(PSP = ifelse(diagnosis == "progressive supranuclear palsy",
                       spec$true_val, spec$false_val),
           Control = ifelse(diagnosis == "control",
                            spec$true_val, spec$false_val))
}

# Harmonize MCMPS data
#
# Calls harmonize_ADKP_studies and then creates binary `Epilepsy` and `Tumor`
# columns based on the `surgeryReason` field. This study is unique in that the
# samples were taken from living tissue during surgery, so the surgury reason is
# used as the diagnosis.
#
# There are Braak scores but no Cerad and very few Thal, which is insufficient
# to add a pathology-based AD diagnosis. The `diagnosis` field is also empty,
# and the study description does not mention other disorders, so it is assumed
# that these patients did not have AD or other neuropathological diseases.
#
# Source file: syn73713777 (version 2) on Synapse
harmonize_MCMPS <- function(metadata, spec) {
  harmonize_ADKP_studies(metadata, spec) |>
    # There are no missing values in this field
    mutate(Epilepsy = ifelse(surgeryReason == "Epilepsy",
                             spec$true_val, spec$false_val),
           Tumor = ifelse(surgeryReason == "Tumor",
                          spec$true_val, spec$false_val))
}


# Harmonize NPS-AD metadata
#
# Harmonize NPS-AD metadata as with all the other ADKP studies, but fill `race`
# and `isHispanic` with `geneticAncestry` and `geneticAncestry_isHispanic` when
# those values are non-missing.
#
# Additionally, we convert the PD, FTD, MS, SCZ, BD, PTSD, ADHD, and OCD
# diagnosis columns to binary string columns, and create a binary `AD` diagnosis
# column based on ADoutcome.
#
# Source file: syn73713770 (version 1) on Synapse
harmonize_NPS_AD <- function(metadata, spec) {
  harmonize_ADKP_studies(metadata, spec) |>
    mutate(
      race = ifelse(is.na(geneticAncestry), race, geneticAncestry),
      isHispanic = ifelse(is.na(geneticAncestry_isHispanic),
                          isHispanic,
                          geneticAncestry_isHispanic),
      # Fix un-harmonized values from geneticAncestry_isHispanic
      isHispanic = case_match(
        isHispanic,
        "Hispanic or Latino" ~ spec$isHispanic$true_val,
        "Not Hispanic or Latino" ~ spec$isHispanic$false_val,
        .default = isHispanic
      ),
      ADoutcome = determineADoutcome(.data, spec),
      AD = case_match(ADoutcome,
                      "AD" ~ spec$true_val,
                      c("Control", "Other") ~ spec$false_val,
                      .default = spec$missing),
      OTHER = case_match(ADoutcome,
                         "Other" ~ spec$true_val,
                         c("Control", "AD") ~ spec$false_val,
                         .default = spec$missing),
      # There are a large number of NA values, which we convert to "missing or unknown"
      across(c(PD, FTD, MS, SCZ, BD, PTSD, ADHD, OCD),
             case_match(.x,
                        1 ~ spec$true_val,
                        0 ~ spec$false_val,
                        .default = spec$missing))
    )
}


# Harmonize ROSMAP metadata
#
# TODO
#
# Source file: syn73713768 (version 1) on Synapse
harmonize_ROSMAP <- function(metadata, spec) {
  harmonize_ADKP_studies(metadata, spec) |>
    mutate(
      ADoutcome = determineADoutcome(.data, spec), # TODO, look at cogdx, dcfdx_lv
      AD = case_match(ADoutcome,
                           "AD" ~ spec$true_val,
                           c("Control", "Other") ~ spec$false_val,
                           .default = spec$missing),
      OTHER = case_match(ADoutcome,
                         "Other" ~ spec$true_val,
                         c("Control", "AD") ~ spec$false_val,
                         .default = spec$missing),
      Control = case_match(ADoutcome,
                           "Control" ~ spec$true_val,
                           c("AD", "Other" ~ spec$false_val),
                           .default = spec$missing)
    )
}


# Harmonize SMIB-AD metadata
#
# Calls harmonize_ADKP_studies and then creates binary `AD` and `Control`
# diagnosis columns that are based on the `diagnosis` column. This study has no
# non-missing amyCerad values and only a few non-missing Thal values, so
# diagnosis cannot be determined based on pathology. The diagnosis supplied in
# the `diagnosis` field is used instead.
#
# Source file: syn73713779 (version 1) on Synapse
harmonize_SMIB_AD <- function(metadata, spec) {
  harmonize_ADKP_studies(metadata, spec) |>
    mutate(AD = ifelse(diagnosis == "Alzheimer Disease",
                       spec$true_val, spec$false_val),
           Control = ifelse(diagnosis == "no cognitive impairment",
                            spec$true_val, spec$false_val))
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
        "No" ~ spec$isHispanic$false_val,
        "Yes" ~ spec$isHispanic$true_val,
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
      ADoutcome = determineADoutcome(.data, spec), # TODO
      AD = case_match(
        ADoutcome,
        "AD" ~ spec$true_val,
        c("Control", "Other") ~ spec$false_val,
        spec$missing ~ spec$false_val, # All samples with missing Braak/Cerad are reference samples
        .default = spec$missing # TODO look at "diagnosis", "Cognitive.status", "Consensus.clinical.diagnosis", "Lewy.body.disease.pathology", "LATE.NC.stage"
      ),
      OTHER = case_match(
        ADoutcome,
        "Other" ~ spec$true_val,
        c("Control", "AD") ~ spec$false_val,
        spec$missing ~ spec$false_val, # All samples with missing Braak/Cerad are reference samples
        .default = spec$missing # TODO look at "diagnosis", "Cognitive.status", "Consensus.clinical.diagnosis", "Lewy.body.disease.pathology", "LATE.NC.stage"
      ),
      PD = ifelse(grepl("Parkinson", diagnosis), spec$true_val, spec$false_val),
      MCI = ifelse(grepl("mild cognitive impairment", diagnosis),
                   spec$true_val, spec$false_val),
      Dementia = ifelse(grepl("dementia", diagnosis),
                        spec$true_val, spec$false_val),
      LBD = ifelse(grepl("Lewy body", diagnosis), spec$true_val, spec$false_val),
      Tumor = ifelse(grepl("Cancer", diagnosis), spec$true_val, spec$false_val),
      Vascular = ifelse(grepl("vascular", diagnosis),
                        spec$true_val, spec$false_val)
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
      # Not enough Braak or Cerad values to use ADoutcome
      # There are no missing values for Dx. Most of the values in
      # ClinConPrimDxText and ClinConSecDxText are missing/NA, so we only use
      # their non-missing values to pull out diagnoses for SCZ, MDD, and OTHER
      # and assume missing values mean no diagnosis for those diseases.
      BD = ifelse(Dx == "BD", spec$true_val, spec$false_val),
      SCZ = ifelse(grepl("Schizophrenia", ClinConPrimDxText),
                   spec$true_val, spec$false_val),
      MDD = ifelse(grepl("Major depressive disorder", ClinConPrimDxText) |
                     grepl("Major depressive disorder", ClinConSecDxText),
                   spec$true_val, spec$false_val),
      OTHER = ifelse(grepl("Schizoaffective", ClinConPrimDxText),
                     spec$true_val, spec$false_val),
      Control = ifelse(Dx == "Control" & is.na(ClinConPrimDxText) &
                         is.na(ClinConSecDxText),
                       spec$true_val, spec$false_val)
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
      # There are no missing values in the Info.Diagnosis field
      PD = ifelse(Info.Diagnosis == "PD", spec$true_val, spec$false_val),
      Control = ifelse(Info.Diagnosis == "Control", spec$true_val, spec$false_val)
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
#   dnastack collections query -c prod-cohort-pmdbs-sc-rnaseq "SELECT * FROM \"collections\".\"prod_cohort_pmdbs_sc_rnaseq\".\"cohort_subject\"" -o csv > "data/downloads/ASAP_subject-export-2026-04-28.csv"
#   dnastack collections query -c prod-cohort-pmdbs-sc-rnaseq "SELECT * FROM \"collections\".\"prod_cohort_pmdbs_sc_rnaseq\".\"cohort_clinpath\"" -o csv > "data/downloads/ASAP_clinpath-export-2026-04-28.csv"
#
# Modifications needed for version April 28, 2026:
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
#   * Censor ages over 90. Some ages are already censored with values of "89+ "
#     (extra space included), and we assume those values should be "90+".
#   * Set empty string ("") and "Not Reported" values to NA so deduplication works
#   * Change `sex` values to all lower case
#   * Change `isHispanic` to True/False values
#   * Update `Braak_NFT`, `amyCerad`, `amyThal`, and `cohort` values to conform
#     to the data dictionary.
#   * Add `dataContributionGroup`
#   * Create a unified path_autopsy_dx field that combines path_autopsy_dx_main
#     and all of the path_autopsy_*_dx fields
#   * Manually assign some values for individuals with duplicate rows that have
#     discrepancies in diagnosis-related fields
#   * Deduplicate rows with the same individual, which should be identical
#     except for some columns where one row is missing a value and one has the
#     value, or where two different groups typed slightly different values for
#     the same column
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
  # First, combine main dx + secondary - eighth dx information for easy
  # searching of diagnosis names
  meta_new <- metadata |>
    rowwise() |>
    mutate(
      path_autopsy_dx = paste(
        c_across(starts_with("path_autopsy")) |>
          setdiff(c("", NA)) |> unique() |> sort(),
        collapse = "; ")
      ) |>
    ungroup() |>
    # Remove all path_autopsy_* variables except the one just created
    select(-(starts_with("path_autopsy") & !path_autopsy_dx))

  meta_new <- meta_new |>
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
      # This data set uses "" or "Not Reported" instead of NA, need to set to NA
      # for deduplication and some harmonization functions to work
      across(everything(), ~ ifelse(
        .x == "" | .x == "Not Reported",
        NA, .x)
      ),
      # Keep the field named "asap_subject_id" but add individualID column
      individualID = asap_subject_id,
      # Assume all 89+ are 90+ for this data set
      ageDeath = censor_ages(ageDeath, spec),
      sex = tolower(sex),
      isHispanic = case_match(
        isHispanic,
        "Hispanic or Latino" ~ spec$isHispanic$true_val,
        "Not Hispanic or Latino" ~ spec$isHispanic$false_val,
        .default = isHispanic
      ),
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
        "0" ~ spec$amyThal$none,
        .default = paste("Phase", suppressWarnings(as.numeric(amyThal)))
      ),
      cohort = case_match(
        cohort,
        c("Banner_Sun_Health_USA", "Banner Sun Health Research Institute") ~ spec$cohort$banner,
        "QSBB_UK" ~ spec$cohort$qsbb,
        .default = cohort
      ),
      dataContributionGroup = spec$dataContributionGroup$asap,

      ## Manual fixes to data with discrepancies

      # There are only 10 "Other" values, all from TEAM_SCHERZER data. We change
      # to "Control", i.e. not PD.
      gp2_phenotype = ifelse(gp2_phenotype == "Other", "Control", gp2_phenotype),

      # TEAM_HAFLER and TEAM_SCHERZER use "Healthy Control", while the other 4
      # teams use "No PD nor other neurological disorder", so we harmonize
      primary_diagnosis = ifelse(primary_diagnosis == "Healthy Control",
                                 "No PD nor other neurological disorder",
                                 primary_diagnosis),

      # ASAP_PMDBS_000003 and ASAP_PMDBS_000020 have two rows each with
      # different cognitive_status: MCI and Dementia. In both cases, they have a
      # diagnosis of PD with dementia so the cognitive_status is changed to
      # Dementia for both copies.
      # ASAP_PMDBS_000002 has two rows with Normal and MCI values, but the
      # primary_diagnosis_text for both has MCI in it and hx_dementia_mci = Yes,
      # so we set both rows to MCI
      # ASAP_PMDBS_000009 has two rows with Normal and MCI values, but is a
      # control sample with hx_dementia_mci = No, so we set both rows to Normal
      cognitive_status = case_match(
        asap_subject_id,
        c("ASAP_PMDBS_000003", "ASAP_PMDBS_000020") ~ "Dementia",
        "ASAP_PMDBS_000002" ~ "MCI",
        "ASAP_PMDBS_000009" ~ "Normal",
        .default = cognitive_status
      ),

      # This subject has two different PMIs listed in duplicate rows, so we use
      # the lowest one (chosen arbitrarily).
      PMI = ifelse(asap_subject_id == "ASAP_PMDBS_000225", 9.25, PMI),

      ## End manual fixes

      ADoutcome = determineADoutcome(.data, spec),
      AD = case_match(
        ADoutcome,
        "AD" ~ spec$true_val,
        c("Control", "Other") ~ spec$false_val,
        .default = spec$missing
      ),

      # gp2_phenotype is either PD or Control. There are no missing values.
      PD = ifelse(gp2_phenotype == "PD", spec$true_val, spec$false_val),

      # The last_diagnosis field also has 24 individuals with mention of Lewy
      # bodies, however 22 of these individuals have gp2_phenotype = PD, and
      # the other two are AD and Control, so we ignore this field.
      LBD = ifelse(grepl("Lewy body disease nos", path_autopsy_dx),
                   spec$true_val, spec$false_val),

      # Assume false if PSP not specified
      PSP = ifelse(grepl("Progressive supranuclear palsy", path_autopsy_dx),
                   spec$true_val, spec$false_val),

      # There are no missing primary_diagnosis values
      OTHER = case_when(
        primary_diagnosis == "Other neurological disorder" ~ spec$true_val,
        primary_diagnosis == "Alzheimer's disease" & ADoutcome == "Other" ~ spec$true_val,
        .default = spec$false_val
      ),

      # Values are either Normal, MCI, Dementia, or NA
      MCI = case_match(
        cognitive_status,
        "MCI" ~ spec$true_val,
        c("Normal", "Dementia") ~ spec$false_val,
        .default = spec$missing
      ),

      Dementia = case_match(
        cognitive_status,
        "Dementia" ~ spec$true_val,
        c("Normal", "MCI") ~ spec$false_val,
        .default = spec$missing
      ),

      # Assume false if Cerebrovascular disease is not listed
      Vascular = case_when(
        grepl("Cerebrovascular disease", path_autopsy_dx) ~ spec$true_val,
        .default = spec$false_val
      )
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
      path_autopsy_dx = str_split(path_autopsy_dx, "; ") |>
        unlist() |> unique() |> sort() |> paste(collapse = "; "),
      across(
        everything(),
        ~ ifelse(is.character(.x) && !all(is.na(.x)),
                 paste(sort(na.omit(unique(.x))), collapse = "; "),
                 unique(.x)) # Leave numeric or NA values alone but keep them in the data frame
      )
    ) |>
    distinct() |>
    as.data.frame()

  return(meta_new)
}

# TODO there is overlap between SCZ and HD, but 3 overlapping samples have slightly different PMIs

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
# * Add `dataContributionGroup` = "Broad Institute"
# * Add `cohort` = "HBTRC"
# * Harmonize `Schizophrenia` diagnosis into binary `SCZ` and `Control` columns
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
      sex = tolower(sex),
      dataContributionGroup = spec$dataContributionGroup$broad,
      cohort = spec$cohort$harvard,
      # There are no missing values in the Schizophrenia field
      SCZ = ifelse(Schizophrenia == "Affected",
                   spec$true_val, spec$false_val),
      Control = ifelse(Schizophrenia == "Unaffected",
                       spec$true_val, spec$false_val)
    ) |>
    select(-Schizophrenia)
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
# * Convert `PMI` to numeric values and replace "n/a" with NA
# * Change `sex` values to all lower case
# * Add `dataContributionGroup` = "Broad Institute"
# * Add `cohort` = "HBTRC"
# * Harmonize `Status` column into binary `HD` and `Control` columns
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
      ageDeath = ifelse(ageDeath == ">89", spec$ageDeath$over90, ageDeath),
      PMI = ifelse(PMI == "n/a", NA, suppressWarnings(as.numeric(PMI))),
      sex = tolower(sex),
      dataContributionGroup = spec$dataContributionGroup$broad,
      cohort = spec$cohort$harvard,
      # There are no missing values in the Status column
      HD = ifelse(Status == "Case", spec$true_val, spec$false_val),
      Control = ifelse(Status == "Control", spec$true_val, spec$false_val)
    )
}
