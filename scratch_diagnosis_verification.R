spec <- config::get(file = "GENESIS_column_spec.yml")
output_dir <- file.path("data", "output")

# Helper functions

print_filled_diagnoses <- function(df, spec) {
  diag_sum <- colSums(!is.na(df[, spec$diagnosis_columns]))
  cols_filled <- names(diag_sum[diag_sum > 0])

  diags <- tidyr::pivot_longer(df,
                               cols = all_of(cols_filled),
                               names_to = "diagnosis_field",
                               values_to = "value")

  cat("\n", "Filled diagnoses:")
  print(table(diags$diagnosis_field, diags$value))
}

validate_ad_v_adoutcome <- function(df){
  tb <- table(df$ADoutcome, df$AD, useNA = "ifany")
  cat("\n", "ADoutcome vs AD")
  print(tb)

  # No ADoutcome = AD are missing an AD diagnosis
  validate_all_ones(df$AD, df$ADoutcome, "AD")

  # No ADoutcome other than AD has an AD diagnosis.
  validate_all_zeros(df$AD, df$ADoutcome, "AD")

  # All missing or unknown values have an NA diagnosis
  if (spec$missing %in% rownames(tb)) {
    stopifnot(all(tb[spec$missing, c("0", "1")] == 0))
    stopifnot(tb[spec$missing, 3] == sum(df$ADoutcome == spec$missing))
  }
  cat("AD vs ADoutcome: PASS", "\n")
}

validate_other_v_adoutcome <- function(df){
  tb <- table(df$ADoutcome, df$Other, useNA = "ifany")
  cat("\n", "ADoutcome vs Other")
  print(tb)

  # No ADoutcome = Other are missing an Other diagnosis
  validate_all_ones(df$Other, df$ADoutcome, "Other")

  # No ADoutcome other than Other has an Other diagnosis.
  validate_all_zeros(df$Other, df$ADoutcome, "Other")

  # All missing or unknown values have an NA diagnosis
  if (spec$missing %in% rownames(tb)) {
    stopifnot(all(tb[spec$missing, c("0", "1")] == 0))
    stopifnot(tb[spec$missing, 3] == sum(df$ADoutcome == spec$missing))
  }
  cat("Other vs ADoutcome: PASS", "\n")
}

validate_no_mci_dementia_overlap <- function(df) {
  tb <- table(df$MCI, df$Dementia) # ignore NAs
  print(tb)

  stopifnot(tb["1", "1"] == 0)
  cat("MCI vs Dementia: PASS", "\n")
}

validate_dementia_ftd_overlap <- function(df) {
  tb <- table(df$Dementia, df$FTD, useNA = "ifany")
  print(tb)

  if (1 %in% df$FTD) {
    # All FTD = 1 should have Dementia = 1
    validate_all_ones(df$Dementia, df$FTD, "1")

    # No FTD = 1 should have Dementia = 0. Order of columns are swapped here because
    # FTD = 0 can have Dementia = 1.
    validate_all_zeros(df$FTD, df$Dementia, "1")
  }
}

validate_all_ones <- function(diagnosis_col, compare_col, compare_match_val) {
  tb <- table(compare_col, diagnosis_col, useNA = "ifany")

  # All compare_col column values that equal `compare_match_val` should be 1
  stopifnot(tb[compare_match_val, "1"] == sum(tb[compare_match_val, ]))
}

validate_all_zeros <- function(diagnosis_col, compare_col, compare_match_val) {
  # Change NA to strings for easier indexing and comparison
  diagnosis_col[is.na(diagnosis_col)] <- "NA"
  compare_col[is.na(compare_col)] <- "NA"

  tb <- table(compare_col, diagnosis_col, useNA = "ifany")

  # All column values that *don't* equal `compare_match_val` should not be 1
  other_rows <- setdiff(rownames(tb), compare_match_val)
  stopifnot(all(tb[other_rows, "1"] == 0))

  # All column values that *don't* equal `compare_match_val` should be 0 or NA
  for (value in other_rows) {
    other_total <- tb[value, "0"]

    # If NA column exists, add to the total
    if ("NA" %in% colnames(tb)) {
      other_total <- other_total + tb[value, "NA"]
    }

    stopifnot(other_total == sum(compare_col == value, na.rm = TRUE))
  }
}




# DivCo

df <- read.csv(file.path(output_dir, "AMP-AD_DiverseCohorts_individual_metadata_harmonized.csv"))

print_filled_diagnoses(df, spec)

validate_ad_v_adoutcome(df)
validate_other_v_adoutcome(df)
validate_no_mci_dementia_overlap(df)
validate_dementia_ftd_overlap(df)


# AMP-PD

df <- read.csv(file.path(output_dir, "AMP-PD_donor_metadata_harmonized.csv"))

print_filled_diagnoses(df, spec)

tb <- table(df$Info.Diagnosis, df$PD, useNA = "ifany")
print(tb)

# All PD samples are 1 and all Control samples are 0
validate_all_ones(df$PD, df$Info.Diagnosis, "PD")
validate_all_zeros(df$PD, df$Info.Diagnosis, "PD")

# No missing values
stopifnot(all(!is.na(df$PD)))

# These mostly agree, displaying for reference
print(table(df$PD_Medical_History.most_recent_diagnosis, df$PD))
print(table(df$PD_Medical_History.diagnosis, df$PD))

cat("PD vs Info.Diagnosis: PASS", "\n")


# ASAP

df <- read.csv(file.path(output_dir, "ASAP_PMDBS_metadata_harmonized.csv"))

print_filled_diagnoses(df, spec)

tb <- table(df$ADoutcome, df$AD, useNA = "ifany")
print(tb)

validate_ad_v_adoutcome(df)

validate_all_ones(df$PD, df$gp2_phenotype, "PD")
validate_all_zeros(df$PD, df$gp2_phenotype, "PD")

validate_all_ones(df$PD, df$primary_diagnosis, "Idiopathic PD")
validate_all_zeros(df$PD, df$primary_diagnosis,
                   c("Alzheimer's disease", "Idiopathic PD", "Other neurological disorder"))

validate_all_ones(df$DLBD, grepl("Lewy", df$path_autopsy_dx), "TRUE")
validate_all_zeros(df$DLBD, grepl("Lewy", df$path_autopsy_dx), "TRUE")

validate_all_ones(df$PSP, grepl("supranuclear", df$path_autopsy_dx), "TRUE")
validate_all_zeros(df$PSP, grepl("supranuclear", df$path_autopsy_dx), "TRUE")

validate_all_ones(df$Other, df$primary_diagnosis, "Other neurological disorder")
validate_all_ones(df$Other, df$ADoutcome, "Other")
validate_all_zeros(df$Other,
                   df$primary_diagnosis == "Other neurological disorder" |
                     df$ADoutcome == "Other",
                   "TRUE")

validate_all_ones(df$MCI, df$cognitive_status, "MCI")
validate_all_zeros(df$MCI, df$cognitive_status, "MCI")

validate_all_ones(df$Dementia, df$cognitive_status, "Dementia")
validate_all_zeros(df$Dementia, df$cognitive_status, "Dementia")

validate_no_mci_dementia_overlap(df)

validate_all_ones(df$Vascular, grepl("Cerebrovascular", df$path_autopsy_dx), "TRUE")
validate_all_zeros(df$Vascular, grepl("Cerebrovascular", df$path_autopsy_dx), "TRUE")


# BD2

df <- read.csv(file.path(output_dir, "BD2_metadata_harmonized.csv"))

print_filled_diagnoses(df, spec)

tb <- table(df$Dx, df$BD, useNA = "ifany")
cat("Dx vs BD", "\n")
print(tb)

# All Dx = BD values should be 1, other values should be 0
validate_all_ones(df$BD, df$Dx, "BD")
validate_all_zeros(df$BD, df$Dx, "BD")

# Internal consistency check: All rows with a valid BD_type value should
# have BD = 1
validate_all_ones(df$BD, grepl("Bipolar", df$BD_type), "TRUE")
validate_all_ones(df$BD, grepl("Bipolar", df$BD_type), "TRUE")

# There should be no missing BD values
stopifnot(all(!is.na(df$BD)))

# All ClinConPrimDxText or ClinConSecDxText with Major depressive disorder should
# have MDD = 1
combined_cols <- paste(df$ClinConPrimDxText, df$ClinConSecDxText)
validate_all_ones(df$MDD, grepl("Major depressive", combined_cols), "TRUE")
validate_all_zeros(df$MDD, grepl("Major depressive", combined_cols), "TRUE")

# There should be no missing MDD values
stopifnot(all(!is.na(df$MDD)))

# All ClinConPrimDxText or ClinConSecDxText with Multiple sclerosis should
# have MS = 1
validate_all_ones(df$MS, grepl("Multiple sclerosis", combined_cols), "TRUE")
validate_all_zeros(df$MS, grepl("Multiple sclerosis", combined_cols), "TRUE")

# There should be no missing MS values
stopifnot(all(!is.na(df$MS)))

# All Dx = Schizoaffective should have an Other diagnosis of 1
# Other Dx values can have a diagnosis of Other so we don't check all_zeros
validate_all_ones(df$Other, df$Dx, "Schizoaffective")

# All ClinConPrimDxText or ClinConSecDxText with Encephalopathy or Schizoaffective
# should have an Other diagnosis of 1. Other samples can have a diagnosis of
# Other too so we don't check all_zeros
validate_all_ones(df$Other, grepl("Enceph|Schizoaffective", combined_cols), "TRUE")

# There should be no missing Other values
stopifnot(all(!is.na(df$Other)))

# All ClinConPrimDxText with Schizophrenia should have SCZ = 1
validate_all_ones(df$SCZ, grepl("Schizophrenia", combined_cols), "TRUE")
validate_all_zeros(df$SCZ, grepl("Schizophrenia", combined_cols), "TRUE")

# There should be no missing SCZ values
stopifnot(all(!is.na(df$SCZ)))

# One cross-filled AD = 1 value, should match ClinConPrimDxText
validate_all_ones(df$AD, grepl("Alz", combined_cols), "TRUE")
validate_all_zeros(df$AD, grepl("Alz", combined_cols), "TRUE")

validate_no_mci_dementia_overlap(df)
validate_dementia_ftd_overlap(df)


# MC_snRNA

df <- read.csv(file.path(output_dir, "MC_snRNA_individual_metadata_harmonized.csv"))

print_filled_diagnoses(df, spec)

tb <- table(df$diagnosis, df$AD, useNA = "ifany")
cat("diagnosis vs AD", "\n")
print(tb)

validate_all_ones(df$AD, df$diagnosis, "Alzheimer Disease")
validate_all_zeros(df$AD, df$diagnosis, "Alzheimer Disease")

# There should be no missing AD values
stopifnot(all(!is.na(df$AD)))


# MC-BrAD

df <- read.csv(file.path(output_dir, "MC-BrAD_individual_metadata_harmonized.csv"))

print_filled_diagnoses(df, spec)

tb <- table(df$diagnosis, df$PSP, useNA = "ifany")
cat("diagnosis vs PSP", "\n")
print(tb)

validate_all_ones(df$PSP, df$diagnosis, "progressive supranuclear palsy")
validate_all_zeros(df$PSP, df$diagnosis, "progressive supranuclear palsy")

# There should be no missing PSP values
stopifnot(all(!is.na(df$PSP)))


# McCarroll_HD

df <- read.csv(file.path(output_dir, "McCarroll_HD_donor_metadata_harmonized.csv"))

print_filled_diagnoses(df, spec)

tb <- table(df$Status, df$HD, useNA = "ifany")
cat("Status vs HD", "\n")
print(tb)

validate_all_ones(df$HD, df$Status, "Case")
validate_all_zeros(df$HD, df$Status, "Case")

# There should be no missing HD values
stopifnot(all(!is.na(df$HD)))


# McCarroll_SCZ

df <- read.csv(file.path(output_dir, "McCarroll_SZvillage_donorMetadata_harmonized.csv"))

print_filled_diagnoses(df, spec)

tb <- table(df$Schizophrenia, df$SCZ, useNA = "ifany")
cat("Schizophrenia vs SCZ", "\n")
print(tb)

validate_all_ones(df$SCZ, df$Schizophrenia, "Affected")
validate_all_zeros(df$SCZ, df$Schizophrenia, "Affected")

# There should be no missing SCZ values
stopifnot(all(!is.na(df$SCZ)))


# MCMPS

df <- read.csv(file.path(output_dir, "MCMPS_individual_metadata_harmonized.csv"))

print_filled_diagnoses(df, spec)

tb <- table(df$surgeryReason, df$Epilepsy, useNA = "ifany")
cat("surgeryReason vs Epilepsy", "\n")
print(tb)

validate_all_ones(df$Epilepsy, df$surgeryReason, "Epilepsy")
validate_all_zeros(df$Epilepsy, df$surgeryReason, "Epilepsy")

validate_all_ones(df$Tumor, df$surgeryReason, "Tumor")
validate_all_zeros(df$Tumor, df$surgeryReason, "Tumor")

# There should be no missing Epilepsy or Tumor values
stopifnot(all(!is.na(df$Epilepsy)) & all(!is.na(df$Tumor)))


# NPS-AD

df <- read.csv(file.path(output_dir, "NPS-AD_individual_metadata_harmonized.csv"))

print_filled_diagnoses(df, spec)

tb <- table(df$ADoutcome, df$AD, useNA = "ifany")
cat("ADoutcome vs AD", "\n")
print(tb)

# AD vs ADoutcome
validate_ad_v_adoutcome(df)

# CDR vs Dementia
validate_all_ones(df$Dementia, df$CDR >= 1, "TRUE")

# Some NA CDRs have Dementia = 1, so we exclude NA CDRs from the all zeros check
validate_all_zeros(df$Dementia, df$CDR >= 1 | is.na(df$CDR), "TRUE")

# CDR vs MCI
validate_all_ones(df$MCI, df$CDR == 0.5, "TRUE")

# Some NA CDRs have MCI = 1, so we exclude NA CDRs from the all zeros check
validate_all_zeros(df$MCI, df$CDR == 0.5 | is.na(df$CDR), "TRUE")

# MCI vs Dementia
validate_no_mci_dementia_overlap(df)
validate_dementia_ftd_overlap(df)

# All ADoutcome = Other should have Other = 1. Two "Other = 1" values are
# cross-filled from BD2 so we don't check for all_zeros
validate_all_ones(df$Other, df$ADoutcome, "Other")


# PEC_CMC

df <- read.csv(file.path(output_dir, "snRNAseq_clinical_metadata_A5_PEC_Roussos_harmonized.csv"))

print_filled_diagnoses(df, spec)

tb <- table(df$primaryDiagnosis, df$BD, useNA = "ifany")
cat("primaryDiagnosis vs BD", "\n")
print(tb)

validate_all_ones(df$BD, df$primaryDiagnosis, "Bipolar disorder")
validate_all_zeros(df$BD, df$primaryDiagnosis, "Bipolar disorder")

validate_all_ones(df$SCZ, df$primaryDiagnosis, "Schizophrenia")
validate_all_zeros(df$SCZ, df$primaryDiagnosis, "Schizophrenia")

# There should be no missing BD or SCZ values
stopifnot(all(!is.na(df$BD)) & all(!is.na(df$SCZ)))


# PEC_SZBD

df <- read.csv(file.path(output_dir, "snRNAseq_clinical_metadata_A6_PEC_RuzickaKellis_harmonized.csv"))

print_filled_diagnoses(df, spec)

tb <- table(df$primaryDiagnosis, df$BD, useNA = "ifany")
cat("primaryDiagnosis vs BD", "\n")
print(tb)

validate_all_ones(df$BD, df$primaryDiagnosis, "Bipolar Disorder")
validate_all_zeros(df$BD, df$primaryDiagnosis, "Bipolar Disorder")

validate_all_ones(df$SCZ, df$primaryDiagnosis, "Schizophrenia")
validate_all_zeros(df$SCZ, df$primaryDiagnosis, "Schizophrenia")

# There should be no missing BD or SCZ values
stopifnot(all(!is.na(df$BD)) & all(!is.na(df$SCZ)))


# ROSMAP

df <- read.csv(file.path(output_dir, "ROSMAP_clinical_harmonized.csv"))

print_filled_diagnoses(df, spec)

tb <- table(df$ADoutcome, df$AD, useNA = "ifany")
cat("ADoutcome vs AD", "\n")
print(tb)

validate_ad_v_adoutcome(df)
validate_other_v_adoutcome(df)

validate_all_ones(df$Dementia, df$dcfdx_lv >= 4, "TRUE")
validate_all_zeros(df$Dementia, df$dcfdx_lv >= 4, "TRUE")

validate_all_ones(df$MCI, df$dcfdx_lv %in% c(2, 3), "TRUE")
validate_all_zeros(df$MCI, df$dcfdx_lv %in% c(2, 3), "TRUE")

validate_no_mci_dementia_overlap(df)
validate_dementia_ftd_overlap(df)


# SEA-AD

df <- read.csv(file.path(output_dir, "SEA-AD_individual_metadata_harmonized.csv"))

print_filled_diagnoses(df, spec)

tb <- table(df$ADoutcome, df$AD, useNA = "ifany")
cat("ADoutcome vs AD", "\n")
print(tb)

validate_ad_v_adoutcome(df)

validate_all_ones(df$Dementia, df$Cognitive.status, "Dementia")
validate_all_zeros(df$Dementia, df$Cognitive.status, "Dementia")

validate_all_ones(df$Dementia, grepl("dementia", df$diagnosis), "TRUE")
validate_all_zeros(df$Dementia, grepl("dementia", df$diagnosis), "TRUE")

validate_all_ones(df$DLBD, grepl("Lewy body", df$diagnosis), "TRUE")
validate_all_zeros(df$DLBD, grepl("Lewy body", df$diagnosis), "TRUE")

validate_all_ones(df$DLBD, grepl("Lewy body", df$Consensus.clinical.diagnosis, ignore.case = TRUE), "TRUE")
validate_all_zeros(df$DLBD, grepl("Lewy body", df$Consensus.clinical.diagnosis, ignore.case = TRUE), "TRUE")

validate_all_ones(df$MCI, df$diagnosis, "mild cognitive impairment")
validate_all_zeros(df$MCI, df$diagnosis, "mild cognitive impairment")

validate_all_ones(df$MCI, grepl("MCI", df$Consensus.clinical.diagnosis), "TRUE")
validate_all_zeros(df$MCI, grepl("MCI", df$Consensus.clinical.diagnosis), "TRUE")

validate_no_mci_dementia_overlap(df)

# Note: There are 2 "diagnosis" values that contain "Parkinsons's disease" but only
# 1 "Consensus.clinical.diagnosis" value that contains "Parkinsons Disease Dementia".
# There is 1 "Consensus.clinical.diagnosis" value with "has Parkinsonian features" that
# has a "diagnosis" of "Parkinson's disease". The diagnosis field appears to be derived
# for the ADKP, so we use Consensus.clinical.diagnosis as ground truth since it's in the
# original data.
validate_all_ones(df$PD, grepl("Parkinsons", df$Consensus.clinical.diagnosis), "TRUE")
validate_all_zeros(df$PD, grepl("Parkinsons", df$Consensus.clinical.diagnosis), "TRUE")

validate_all_ones(df$Tumor, grepl("Brain Cancer", df$diagnosis), "TRUE")
validate_all_zeros(df$Tumor, grepl("Brain Cancer", df$diagnosis), "TRUE")

validate_all_ones(df$Tumor, grepl("brain mets", df$Consensus.clinical.diagnosis), "TRUE")
validate_all_zeros(df$Tumor, grepl("brain mets", df$Consensus.clinical.diagnosis), "TRUE")

validate_all_ones(df$Vascular, grepl("Vascular", df$Consensus.clinical.diagnosis), "TRUE")
validate_all_zeros(df$Vascular, grepl("Vascular", df$Consensus.clinical.diagnosis), "TRUE")

# All ADoutcome = "Other" should have a diagnosis of Other and all Consensus.clinical.diagnosis
# starting with "Other" or "Multiple System" should have a diagnosis of Other.
validate_all_ones(df$Other,
                  df$ADoutcome == "Other" |
                    grepl("^Other|Multiple System", df$Consensus.clinical.diagnosis),
                  "TRUE")
validate_all_zeros(df$Other,
                   df$ADoutcome == "Other" |
                     grepl("^Other|Multiple System", df$Consensus.clinical.diagnosis),
                   "TRUE")


# SMIB-AD

df <- read.csv(file.path(output_dir, "SMIB-AD_individual_human_metadata_harmonized.csv"))

print_filled_diagnoses(df, spec)

tb <- table(df$diagnosis, df$AD, useNA = "ifany")
cat("diagnosis vs AD", "\n")
print(tb)

validate_all_ones(df$AD, df$diagnosis, "Alzheimer Disease")
validate_all_zeros(df$AD, df$diagnosis, "Alzheimer Disease")

# There should be no missing values
stopifnot(all(!is.na(df$AD)))


# TargetALS

df <- read.csv(file.path(output_dir, "TargetALS_final_metadata_harmonized.csv"))

print_filled_diagnoses(df, spec)

tb <- table(df$selection_group, df$AD, useNA = "ifany")
cat("selection_group vs AD", "\n")
print(tb)

validate_all_ones(df$AD, df$selection_group, "ALS, Alzheimers")
validate_all_ones(df$AD, grepl("Alz", df$Subject.Group.Subcategory), "TRUE")
validate_all_zeros(df$AD,
                   df$selection_group == "ALS, Alzheimers" |
                     grepl("Alz", df$Subject.Group.Subcategory),
                   "TRUE")

validate_all_ones(df$ALS, grepl("ALS", df$selection_group), "TRUE")
validate_all_zeros(df$ALS, grepl("ALS", df$selection_group), "TRUE")

validate_all_ones(df$ALS, grepl("ALS", df$Subject.Group), "TRUE")
validate_all_zeros(df$ALS, grepl("ALS", df$Subject.Group), "TRUE")

validate_all_ones(df$ALS, grepl("ALS", df$Subject.Group.Subcategory), "TRUE")
validate_all_zeros(df$ALS, grepl("ALS", df$Subject.Group.Subcategory), "TRUE")

validate_all_ones(df$Dementia, df$Mnd.With.Dementia %in% c("AD", "Lewy Body Dementia", "Yes"), "TRUE")

validate_all_zeros(df$Dementia,
                   df$Mnd.With.Dementia %in% c("AD", "Lewy Body Dementia", "Yes") |
                     df$FTD == 1,
                   "TRUE")

validate_all_ones(df$DLBD,
                  df$Mnd.With.Dementia == "Lewy Body Dementia" |
                    grepl("Lewy Body", df$Comorbidities),
                  "TRUE")
validate_all_zeros(df$DLBD,
                   df$Mnd.With.Dementia == "Lewy Body Dementia" |
                     grepl("Lewy Body", df$Comorbidities),
                   "TRUE")

validate_all_ones(df$Epilepsy, grepl("Epilepsy", df$Comorbidities, ignore.case = TRUE), "TRUE")
validate_all_zeros(df$Epilepsy, grepl("Epilepsy", df$Comorbidities, ignore.case = TRUE), "TRUE")

validate_all_ones(df$FTD, df$selection_group, "ALS/FTD")
validate_all_zeros(df$FTD, df$selection_group, "ALS/FTD")

validate_all_ones(df$FTD, grepl("FTD", df$Subject.Group.Subcategory), "TRUE")
# One sample doesn't have FTD in Subject.Group.Subcategory but has it in selection_group
# so we don't test for all zeros here

validate_all_ones(df$MDD, grepl("depression", df$Comorbidities, ignore.case = TRUE), "TRUE")
validate_all_zeros(df$MDD, grepl("depression", df$Comorbidities, ignore.case = TRUE), "TRUE")

validate_all_ones(df$MS, grepl("Multiple Sclerosis", df$Comorbidities, ignore.case = TRUE), "TRUE")
validate_all_zeros(df$MS, grepl("Multiple Sclerosis", df$Comorbidities, ignore.case = TRUE), "TRUE")

validate_all_ones(df$PD, grepl("Parkinson", df$Subject.Group.Subcategory), "TRUE")
validate_all_ones(df$PD, grepl("Parkinson", df$Comorbidities, ignore.case = TRUE), "TRUE")
validate_all_zeros(df$PD,
                   grepl("Parkinson",
                         df$Subject.Group.Subcategory) |
                     grepl("Parkinson", df$Comorbidities, ignore.case = TRUE),
                   "TRUE")

validate_all_ones(df$Tumor, grepl("Metastatic", df$Subject.Group.Subcategory), "TRUE")
validate_all_zeros(df$Tumor, grepl("Metastatic", df$Subject.Group.Subcategory), "TRUE")

validate_all_ones(df$Vascular, grepl("CAA|Cerebrovascular", df$Comorbidities), "TRUE")
validate_all_ones(df$Vascular, grepl("arteriosclerosis", df$Comorbidities), "TRUE")
validate_all_ones(df$Vascular, grepl("SVD", df$Comorbidities), "TRUE")
validate_all_ones(df$Vascular, grepl("Cerebrovascular", df$Subject.Group.Subcategory), "TRUE")

validate_all_zeros(df$Vascular,
                   grepl("CAA|Cerebrovascular", df$Comorbidities) |
                     grepl("arteriosclerosis", df$Comorbidities) |
                     grepl("SVD", df$Comorbidities) |
                     grepl("Cerebrovascular", df$Subject.Group.Subcategory),
                   "TRUE")

validate_all_ones(df$Other, grepl("stroke|enceph|meningitis", df$Comorbidities), "TRUE")
validate_all_zeros(df$Other, grepl("stroke|enceph|meningitis", df$Comorbidities), "TRUE")

validate_dementia_ftd_overlap(df)


