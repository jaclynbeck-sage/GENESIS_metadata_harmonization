# Harmonize ASAP / GEN-A15 metadata
#
# Modifies the ASAP subject and clinical metadata files to conform to the
# GENESIS data dictionary. Data downloaded via command line on April 28, 2026
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
      # Braak_NFT is already in Roman numerals except for "0" and NA
      Braak_NFT = case_match(
        Braak_NFT,
        NA ~ spec$missing,
        "0" ~ spec$Braak_NFT$none,
        .default = paste("Stage", Braak_NFT)
      ),
      # Braak_LB is numbered 0-6, convert to Roman numerals
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
        "Banner_Sun_Health_USA" ~ spec$cohort$banner,
        "Banner Sun Health Research Institute" ~ spec$cohort$banner,
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
      # shortcut function from ADKP_datasets.R
      AD = make_binary_column(ADoutcome, "AD", spec),

      # gp2_phenotype is either PD or Control. There are no missing values.
      PD = make_binary_column(gp2_phenotype, "PD", spec),

      # The last_diagnosis field also has 24 individuals with mention of Lewy
      # bodies, however 22 of these individuals have gp2_phenotype = PD, and
      # the other two are AD and Control, so we ignore this field.
      LBD = ifelse(grepl("Lewy body disease nos", path_autopsy_dx),
                   1, 0),

      # Assume false if PSP not specified
      PSP = ifelse(grepl("Progressive supranuclear palsy", path_autopsy_dx),
                   1, 0),

      # There are no missing primary_diagnosis values
      Other = case_when(
        primary_diagnosis == "Other neurological disorder" ~ 1,
        primary_diagnosis == "Alzheimer's disease" & ADoutcome == "Other" ~ 1,
        .default = 0
      ),

      # cognitive_status values are either Normal, MCI, Dementia, or NA
      MCI = make_binary_column(cognitive_status, "MCI", spec),
      Dementia = make_binary_column(cognitive_status, "Dementia", spec),

      # Assume false if Cerebrovascular disease is not listed
      Vascular = case_when(
        grepl("Cerebrovascular disease", path_autopsy_dx) ~ 1,
        .default = 0
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
      # Recombine path_autopsy_dx using only unique values
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
