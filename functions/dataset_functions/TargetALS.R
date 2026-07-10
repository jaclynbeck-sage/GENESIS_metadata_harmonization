# Harmonize TargetALS data
#
# Source file: locally downloaded file provided by Jaro, dated Jan 16, 2026
#
# Modifications needed for version Jan 16, 2026:
#   * Rename columns:
#     * `Externalsubjectid` => `individualID`
#     * `PMI (hrs)` => `PMI`
#     * `Sex` => `sex`
#     * `Race` => `race`
#     * `Age.At.Death` => `ageDeath`
#     * `Site.Specimen.Collected` => `cohort`
#   * Trim whitespace from individualIDs
#   * Create an `isHispanic` colum based on `Ethnicity`. We keep the Ethnicity
#     column because it has additional information in it beyond Hispanic/Latino
#     status.
#   * Update `race` to conform to the data dictionary
#   * Censor ages over 90
#   * Make PMI numeric
#   * Update sex to conform to the data dictionary
#   * Extract Braak_NFT and amyThal values from the `Comorbidities` column
#   * Add `dataContributionGroup` = "Mount Sinai School of Medicine"
#   * Update `cohort` to conform to the data dictionary
#   * Create binary `AD`, `ALS`, and `FTD` diagnosis columns based on the
#     `selection_group` column
#   * Create binary `PD`, `Tumor`, and `Vascular` diagnosis columns based on the
#     `Subject Group Subcategory` column
#   * Create binary `DLBD` and `Dementia` diagnosis columns based on the
#     `Mnd.With.Dementia` column
#   * Create binary `MDD`, `MS`, `Epilepsy`, and `Other` diagnosis columns based
#     on the `Comorbidities` column
#   * Add to the `DLBD` and `Vascular` columns with data from the
#     `Comorbidities` column
#   * Strip quote characters from the `Comorbidities` column so the file writes
#     correctly
#   * Manually fix the ageDeath of one individual with duplicate rows. This
#     individual has two different ageDeath values, so we use the largest value.
harmonize_TargetALS <- function(metadata, spec) {
  metadata |>
    dplyr::rename(
      individualID = Externalsubjectid,
      PMI = `PMI (hrs)`,
      sex = Sex,
      race = Race,
      ageDeath = Age.At.Death,
      cohort = Site.Specimen.Collected
    ) |>
    dplyr::mutate(
      individualID = str_trim(individualID),
      ageDeath = censor_ages(ageDeath, spec),
      PMI = as.numeric(PMI),
      isHispanic = case_match(
        Ethnicity,
        "Hispanic or Latino" ~ spec$isHispanic$true_val,
        "Not Hispanic or Latino" ~ spec$isHispanic$false_val,
        "Not Hispanic or Latino, Ashkenazi Jewish" ~ spec$isHispanic$false_val,
        c("Unknown", "Not Reported") ~ spec$missing,
        .default = Ethnicity
      ),
      race = case_match(
        race,
        "Indian/Asian" ~ spec$race$Asian,
        "White/Scottish" ~ spec$race$White,
        c(NA, "Unknown", "Not reported", "Not Reported") ~ spec$missing,
        .default = race
      ),
      sex = tolower(sex),
      cohort = ifelse(cohort == "Edinburgh", spec$cohort$edinburgh, cohort),
      dataContributionGroup = spec$dataContributionGroup$mssm,

      # Extract Braak_NFT -- values are written as "Braak tangle stage <value>",
      # where <value> is either 1 (number) or Roman numerals
      Braak_NFT = str_extract(Comorbidities, "Braak tangle stage [0-9A-Z]*") |>
        str_replace("Braak tangle stage", "Stage") |>
        str_replace("1", "I"),

      # Extract amyThal - values are written as "Thal amyloid phase <value>,
      # where <value> is a number from 1-4
      amyThal = str_extract(Comorbidities, "Thal amyloid phase [0-9]") |>
        str_replace("Thal amyloid phase", "Phase"),

      # Diagnosis columns
      AD = pmax(
        make_binary_column(selection_group, "ALS, Alzheimers", spec),
        grep_to_binary_column(`Subject Group Subcategory`, "Alzheimer's Disease")
      ),
      ALS = grep_to_binary_column(selection_group, "^ALS"),
      FTD = make_binary_column(selection_group, "ALS/FTD", spec),

      PD = pmax(
        grep_to_binary_column(`Subject Group Subcategory`, "Parkinson's Disease"),
        grep_to_binary_column(Comorbidities, "Parkinsons")
      ),
      Tumor = grep_to_binary_column(`Subject Group Subcategory`, "Metastatic Carcinoma"),
      Vascular = pmax(
        grep_to_binary_column(`Subject Group Subcategory`, "Cerebrovascular disease"),
        grep_to_binary_column(Comorbidities,
                              "Cerebrovascular|CAA|non(-| )amyloid SVD|arteriosclerosis")
      ),

      DLBD = pmax(
        make_binary_column(Mnd.With.Dementia, "Lewy Body Dementia", spec),
        grep_to_binary_column(Comorbidities, "Lewy Body Disease")
      ),
      Dementia = make_binary_column(Mnd.With.Dementia,
                                    c("Yes", "AD", "Lewy Body Dementia"),
                                    spec),

      MDD = grep_to_binary_column(Comorbidities, "Depression"),
      MS = grep_to_binary_column(Comorbidities, "Multiple sclerosis"),
      Epilepsy = grep_to_binary_column(Comorbidities, "Epilepsy"),
      Other = grep_to_binary_column(Comorbidities,
                                    "stroke|encephalopathy|meningitis"),

      # Strip quotes from Comorbidities column
      Comorbidities = str_replace_all(Comorbidities, "\"", ""),

      ### Manual fix to duplicate record
      ageDeath = ifelse(individualID == "NEUNP654LLE",
                        max(ageDeath[individualID == "NEUNP654LLE"]),
                        ageDeath)
      ### End manual fix
    ) |>
    as.data.frame()
}
