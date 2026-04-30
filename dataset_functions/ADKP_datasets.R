# This script contains a generic harmonization function for all studies that
# were harmonized for the ADKP Harmonization Project, plus individual functions
# for each study to handle diagnosis and any other special cases.
# Datasets with functions in this file:
#   AMP-AD_DiverseCohorts
#   MC-BrAD
#   MC_snRNA
#   MCMPS
#   NPS-AD
#   ROSMAP
#   SEA-AD
#   SMIB-AD

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
#
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
# Calls harmonize_ADKP_studies and then creates binary diagnosis columns based
# on the ADoutcome variable.
#
# Source metadata: syn73713769 (version 2) on Synapse
#
harmonize_Diverse_Cohorts <- function(metadata, spec) {
  harmonize_ADKP_studies(metadata, spec) |>
    mutate(AD = make_binary_column(ADoutcome, "AD", spec),
           Other = make_binary_column(ADoutcome, "Other", spec),
           Control = make_binary_column(ADoutcome, "Control", spec))
}


# Harmonize MC-BrAD data
#
# Calls harmonize_ADKP_studies and then creates binary `PSP` and `Control`
# diagnosis columns that are based on the `diagnosis` field. The study
# documentation verified that control samples have no diagnosis of any of the
# common neuropathological disorders like AD, PD, HD, etc.
#
# Source file: syn73713775 (version 2) on Synapse
#
harmonize_MC_BrAD <- function(metadata, spec) {
  harmonize_ADKP_studies(metadata, spec) |>
    mutate(PSP = make_binary_column(diagnosis, "progressive supranuclear palsy", spec),
           Control = make_binary_column(diagnosis, "control", spec))
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
#
harmonize_MC_snRNA <- function(metadata, spec) {
  harmonize_ADKP_studies(metadata, spec) |>
    mutate(AD = make_binary_column(diagnosis, "Alzheimer Disease", spec),
           Control = make_binary_column(diagnosis, "Control", spec))
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
#
harmonize_MCMPS <- function(metadata, spec) {
  harmonize_ADKP_studies(metadata, spec) |>
    mutate(Epilepsy = make_binary_column(surgeryReason, "Epilepsy", spec),
           Tumor = make_binary_column(surgeryReason, "Tumor", spec))
}


# Harmonize NPS-AD metadata
#
# Harmonize NPS-AD metadata as with all the other ADKP studies, but fill `race`
# and `isHispanic` with `geneticAncestry` and `geneticAncestry_isHispanic` when
# those values are non-missing.
#
# Additionally, we create binary `AD` and `Other` diagnosis columns based on
# ADoutcome. Binary `MCI` and `Dementia` columns are created based on the CDR
# column:
#   MCI: CDR = 0.5 (questionable dementia)
#   Dementia: CDR >= 1 (mild, moderate, severe, profound, or terminal dementia)
#
# Source file: syn73713770 (version 1) on Synapse
#
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
      AD = make_binary_column(ADoutcome, "AD", spec),
      Other = make_binary_column(ADoutcome, "Other", spec),
      MCI = case_when(
        CDR == 0.5 ~ 1,
        CDR == 0 | CDR >= 1 ~ 0,
        .default = NA
      ),
      Dementia = case_when(
        CDR >= 1 | FTD == 1 ~ 1,
        CDR <= 0.5 ~ 0,
        .default = NA
      )
      # TODO Control column?
    )
}


# Harmonize ROSMAP metadata
#
# Calls harmonize_ADKP_studies and then creates binary `AD` and `Other`
# diagnosis columns based on the ADoutcome variable.
# Binary `MCI` and `Dementia` columns are created based on the dcfdx_lv variable.
# The binary `Control` diagnosis column is based on all of the other diagnosis
# columns:
#   1 if all the other columns are 0,
#   0 if there is a 1 in any other column,
#   NA if any column has an NA value
#
# Note: There are two cognitive status columns (dcfdx_lv and cogdx) with the
# same coding:
#   1: no impairment
#   2: MCI
#   3: MCI + other cause of impairment
#   4: AD dementia
#   5: AD dementia + other dementia
#   6: other dementia
#
# The two columns mostly agree with each other, but have some discrepancies. As
# almost half of the cogdx values are missing, we use dcfdx_lv to determine
# MCI and Dementia status.
#
# Source file: syn73713768 (version 1) on Synapse
#
harmonize_ROSMAP <- function(metadata, spec) {
  harmonize_ADKP_studies(metadata, spec) |>
    mutate(
      ADoutcome = determineADoutcome(.data, spec),
      AD = make_binary_column(ADoutcome, "AD", spec),
      Other = make_binary_column(ADoutcome, "Other", spec),
      MCI = case_match(dcfdx_lv,
                       c(2, 3) ~ 1, # 2 and 3 = MCI
                       c(1, 4, 5, 6) ~ 0,
                       .default = NA),
      Dementia = case_match(dcfdx_lv,
                            c(4, 5, 6) ~ 1, # 4, 5, and 6 = Dementia
                            c(1, 2, 3) ~ 0,
                            .default = NA),
      # Only mark someone as Control if they don't have a 1 in any of the above
      # diagnosis columns. NA values in any of these columns propagate as NA in
      # the Control column.
      tmp_sum = AD + Other + MCI + Dementia,
      Control = case_when(
        tmp_sum == 0 ~ 1,
        tmp_sum > 0 ~ 0,
        .default = NA
      )
    ) |>
    select(-tmp_sum)
}


# Harmonize SEA-AD metadata
#
# Calls harmonize_ADKP_studies and then creates binary diagnosis variables for
#
# Source file: syn73713778 (version 1) on Synapse
#
harmonize_SEA_AD <- function(metadata, spec) {
  harmonize_ADKP_studies(metadata, spec) |>
    mutate(
      ADoutcome = determineADoutcome(.data, spec),
      # Reference samples have missing Braak/Cerad but should have a "Control" ADoutcome
      ADoutcome = ifelse(dataset == "Reference", "Control", ADoutcome),
      # Binary diagnosis columns
      AD = make_binary_column(ADoutcome, "AD", spec),
      PD = ifelse(grepl("Parkinsons Disease", Consensus.clinical.diagnosis), 1, 0),
      MCI = ifelse(grepl("MCI", Consensus.clinical.diagnosis), 1, 0),

      # Reference samples end up with NA Dementia values, change to 0
      Dementia = make_binary_column(Cognitive.status, "Dementia", spec),
      Dementia = ifelse(dataset == "Reference", 0, Dementia),

      LBD = ifelse(grepl("Lewy body disease", Consensus.clinical.diagnosis), 1, 0),

      # Other is 1 if the diagnosis has "Other" or "Multiple System Atrophy", or
      # if ADoutcome is "Other"
      Other = ifelse(grepl("Other|Multiple System Atrophy", Consensus.clinical.diagnosis), 1, 0),
      Other = ifelse(ADoutcome == "Other", 1, Other),

      Tumor = ifelse(grepl("brain mets", Consensus.clinical.diagnosis), 1, 0),
      Vascular = ifelse(grepl("Vascular Dementia", Consensus.clinical.diagnosis), 1, 0),

      # Control is 1 if Consensus.clinical.diagnosis = "Control" AND
      # ADoutcome = "Control".
      # Reference samples end up with NA Control values, change to 1.
      Control = make_binary_column(Consensus.clinical.diagnosis, "Control", spec),
      Control = ifelse(ADoutcome != "Control", 0, Control),
      Control = ifelse(dataset == "Reference", 1, Control)
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
#
harmonize_SMIB_AD <- function(metadata, spec) {
  harmonize_ADKP_studies(metadata, spec) |>
    mutate(AD = make_binary_column(diagnosis, "Alzheimer Disease", spec),
           Control = make_binary_column(diagnosis, "no cognitive impairment", spec))
}

