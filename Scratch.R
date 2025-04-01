library(synapser)
library(dplyr)
library(purrr)
library(stringr)
library(readxl)

spec <- config::get(file = "GENESIS_harmonization.yml")
source("util_functions.R")
source("dataset_specific_functions.R")


syn_ids <- list(
  "GEN-A1" = "syn55251012.4", # NPS-AD
  "ROSMAP" = "syn64759878.4", # Harmonized file, for GEN-A2, GEN-A8, GEN-A13, GEN-B6
  "SEA-AD" = "syn31149116.7", # SEA-AD, for GEN-A4 and GEN-B5
  "GEN-A9" = "syn22432749.1", # SMIB-AD
  "GEN-A10" = "syn25891193.1", # MCMPS
  "GEN-A11" = "syn31563038.1", # MC_snRNA
  "GEN-A12" = "syn51401700.2", # MC-BrAD
  "GEN-A14" = "syn24610550.2", # HBI_scRNAseq
  "Diverse_Cohorts" = "syn64759872.4" # Harmonized file, for GEN-B4
  # "GEN-B1" = "TBD",
  # "GEN-B2" = "TBD",
  # "GEN-B3" = "TBD"
)

check_new_versions(syn_ids)

UPLOAD_SYNID <- "syn64759869"

manifest <- c()

synLogin()

harmonized_baseline <- readRDS(file.path("data", "tmp", "harmonized_originals.rds"))


# GEN-A1 -----------------------------------------------------------------------
# Has sample overlap with MSBB, ROSMAP, and Diverse Cohorts

meta_file <- synapse_download(syn_ids[["GEN-A1"]])
meta <- read.csv(meta_file$path)

colnames(meta)

print_qc(meta, isHispanic_col = "ethnicity", pmi_col = "PMI", cerad_col = "CERAD")

# TODO should pull metadata from Diverse Cohorts, MSBB, and ROSMAP for the non-NPS-AD samples,
# since some info like Braak is missing from this metadata. There is a clinical data
# file in the NPS-AD project with Braak etc
meta_new <- meta |>
  select(-Component) |>
  rename(
    isHispanic = ethnicity,
    pmi = PMI,
    amyCerad = CERAD
  ) |>
  mutate(
    # Fix to allow comparison to Diverse Cohorts
    diverseCohortsIndividualIDFormat = case_match(diverseCohortsIndividualIDFormat,
      29637 ~ "29637_MSSM",
      29582 ~ "29582_MSSM",
      .default = as.character(diverseCohortsIndividualIDFormat)
    ),
    pmi = pmi / 60,
    pmiUnits = "hours",
    # ageDeath = censor_ages(ageDeath, spec),
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
    Braak = spec$missing,
    bScore = spec$missing,
    cohort = case_match(cohort,
      "MSBB" ~ spec$cohort$msbb,
      .default = cohort
    ), # TODO need to differentiate ROS and MAP
    # TODO this might not be quite right
    dataContributionGroup = case_match(cohort,
      spec$cohort$msbb ~ spec$dataContributionGroup$mssm,
      spec$cohort$hbcc ~ spec$dataContributionGroup$nimh,
      c(spec$cohort$ros, spec$cohort$map) ~ spec$dataContributionGroup$rush,
      .default = ""
    )
  )


print_qc(meta_new)
validate_values(meta_new, spec)

new_filename <- write_metadata(meta_new, meta_file$name)
new_syn_id <- synapse_upload(new_filename, UPLOAD_SYNID)

manifest <- rbind(manifest, data.frame(
  GENESIS_study = "GEN-A1",
  ADKP_study = "NPS-AD",
  metadata_synid = paste0(new_syn_id$id, ".", new_syn_id$versionNumber)
))


# GEN-A2, GEN-A8, GEN-A13, GEN-B6 ----------------------------------------------
# All use harmonized ROSMAP metadata

manifest <- rbind(
  manifest,
  data.frame(
    GENESIS_study = c("GEN-A2", "GEN-A8", "GEN-A13", "GEN-B6"),
    ADKP_study = c("ROSMAP", "snRNAseqAD_TREM2", "snRNAseqPFC_BA10", "MIT_ROSMAP_Multiomics"),
    metadata_synid = syn_ids[["ROSMAP"]]
  )
)


# GEN-A4, GEN-B5 ---------------------------------------------------------------

# Uses SEA-AD metadata. The version of SEA-AD that is on Synapse is missing
# Hispanic/Latino information that is present in the version released by the
# Allen Institute on brain-map.org. We use the version on Synapse but pull in
# the missing information from the AI version.

meta_file <- synapse_download(syn_ids[["SEA-AD"]])

sea_ad_file <- file.path("data", "downloads", "sea-ad_cohort_donor_metadata_072524.xlsx")
download.file("https://brainmapportal-live-4cc80a57cd6e400d854-f7fdcae.divio-media.net/filer_public/b4/c7/b4c727e1-ede1-4c61-b2ee-bf1ae4a3ef68/sea-ad_cohort_donor_metadata_072524.xlsx",
              destfile = sea_ad_file
)

meta <- read_xlsx(meta_file$path)
meta_sea_ad <- read_xlsx(sea_ad_file)

colnames(meta)
colnames(meta_sea_ad)

print_qc(meta,
         cerad_col = "CERAD",
         thal_col = "Thal phase"
)

print_qc(meta_sea_ad,
         isHispanic_col = "Hispanic/Latino"
)

meta_new <- harmonize_SEA_AD(meta, meta_sea_ad, spec)

print_qc(meta_new)
validate_values(meta_new, spec)

# Original file is an Excel file, change filename to CSV file
file_write <- str_replace(meta_file$name, "xlsx", "csv")
new_filename <- write_metadata(meta_new, file_write)
new_syn_id <- synapse_upload(new_filename, UPLOAD_SYNID)

manifest <- rbind(
  manifest,
  data.frame(
    GENESIS_study = c("GEN-A4", "GEN-B5"),
    ADKP_study = c("SEA-AD", "SEA-AD"),
    metadata_synid = paste0(new_syn_id$id, ".", new_syn_id$versionNumber)
  )
)


# GEN-A9 -----------------------------------------------------------------------
# No overlap with other data sets.

# TODO unknown how they encode CERAD
# Best guess is using ROSMAP-style coding:
# 1 = "None"
# 2 = "Sparse"
# 4 = "Frequent"
# There are no samples with "3"

meta_file <- synapse_download(syn_ids[["GEN-A9"]])
meta <- read.csv(meta_file$path)

colnames(meta)

print_qc(meta, isHispanic_col = "ethnicity", cerad_col = "CERAD")

meta_new <- harmonize_SMIB_AD(meta, spec)

print_qc(meta_new)
validate_values(meta_new, spec)

new_filename <- write_metadata(meta_new, meta_file$name)
new_syn_id <- synapse_upload(new_filename, UPLOAD_SYNID)

manifest <- rbind(
  manifest,
  data.frame(
    GENESIS_study = "GEN-A9",
    ADKP_study = "SMIB-AD",
    metadata_synid = paste0(new_syn_id$id, ".", new_syn_id$versionNumber))
)


# GEN-A10 ----------------------------------------------------------------------

# Note: Samples come from Mayo Clinic but there is no overlap with AMP-AD 1.0 or
# Diverse Cohorts data.

meta_file <- synapse_download(syn_ids[["GEN-A10"]])
meta <- read.csv(meta_file$path)

colnames(meta)

print_qc(meta, isHispanic_col = "ethnicity", cerad_col = "CERAD")

meta_new <- harmonize_MCMPS(meta, spec)

print_qc(meta_new)
validate_values(meta_new, spec)

new_filename <- write_metadata(meta_new, meta_file$name)
new_syn_id <- synapse_upload(new_filename, UPLOAD_SYNID)

manifest <- rbind(
  manifest,
  data.frame(
    GENESIS_study = "GEN-A10",
    ADKP_study = "MCMPS",
    metadata_synid = paste0(new_syn_id$id, ".", new_syn_id$versionNumber))
)


# GEN-A11 ----------------------------------------------------------------------

# Has some sample overlap with AMP-AD 1.0 Mayo metadata and Diverse Cohorts
# metadata. After harmonization we pull in any missing values from 1.0 and
# Diverse Cohorts metadata.

meta_file <- synapse_download(syn_ids[["GEN-A11"]])
meta <- read.csv(meta_file$path)

colnames(meta)

print_qc(meta, isHispanic_col = "ethnicity", cerad_col = "CERAD")

meta_new <- harmonize_MC_snRNA(meta, harmonized_baseline, spec)

print_qc(meta_new)
validate_values(meta_new, spec)

new_filename <- write_metadata(meta_new, meta_file$name)
new_syn_id <- synapse_upload(new_filename, UPLOAD_SYNID)

manifest <- rbind(
  manifest,
  data.frame(
    GENESIS_study = "GEN-A11",
    ADKP_study = "MC_snRNA",
    metadata_synid = paste0(new_syn_id$id, ".", new_syn_id$versionNumber))
)


# GEN-A12 ----------------------------------------------------------------------

# Has some sample overlap with AMP-AD 1.0 Mayo metadata and Diverse Cohorts
# metadata. There are some missing values in this data set that exist in 1.0 or
# DC metadata, so we pull those values in.

meta_file <- synapse_download(syn_ids[["GEN-A12"]])
meta <- read.csv(meta_file$path)

colnames(meta)

print_qc(meta, isHispanic_col = "ethnicity", cerad_col = "CERAD", thal_col = "Thal")

meta_new <- harmonize_MC_BrAD(meta, harmonized_baseline, spec)

print_qc(meta_new)
validate_values(meta_new, spec)

new_filename <- write_metadata(meta_new, meta_file$name)
new_syn_id <- synapse_upload(new_filename, UPLOAD_SYNID)

manifest <- rbind(
  manifest,
  data.frame(
    GENESIS_study = "GEN-A12",
    ADKP_study = "MC-BrAD",
    metadata_synid = new_syn_id)
)


# GEN-A14 ----------------------------------------------------------------------

# This study uses data from MSSM but there is no sample overlap with AMP-AD 1.0
# MSBB or Diverse Cohorts data.

# TODO unclear how they encode CERAD and I can't guess based on correlation with
# Braak.
# TODO we might not be using this data set

meta_file <- synapse_download(syn_ids[["GEN-A14"]])
meta <- read.csv(meta_file$path)

colnames(meta)

print_qc(meta, isHispanic_col = "ethnicity", cerad_col = "CERAD")

meta_new <- meta |>
  rename(
    isHispanic = ethnicity,
    amyCerad = CERAD
  ) |>
  mutate(
    pmi = pmi / 60,
    ageDeath = censor_ages(ageDeath, spec),
    # "Hispanic" status is encoded in the race column
    isHispanic = case_match(race,
      "Hispanic" ~ spec$isHispanic$hisp_true,
      NA ~ spec$missing,
      .default = isHispanic
    ),
    race = case_match(race,
      "Black" ~ spec$race$Black,
      "Hispanic" ~ spec$missing,
      .default = race
    ),
    apoeGenotype = spec$missing,
    apoe4Status = get_apoe4Status(apoeGenotype, spec),
    Braak = to_Braak_stage(Braak, spec),
    bScore = get_bScore(Braak, spec),
    amyCerad = spec$missing, # TODO
    amyAny = get_amyAny(amyCerad, spec),
    amyThal = spec$missing,
    amyA = spec$missing,
    dataContributionGroup = spec$dataContributionGroup$mssm,
    cohort = spec$cohort$msbb
  )

print_qc(meta_new)
validate_values(meta_new, spec)

new_filename <- write_metadata(meta_new, meta_file$name)
new_syn_id <- synapse_upload(new_filename, UPLOAD_SYNID)

manifest <- rbind(
  manifest,
  data.frame(
    GENESIS_study = "GEN-A14",
    ADKP_study = "HBI_scRNAseq",
    metadata_synid = new_syn_id)
)


# GEN-B4 -----------------------------------------------------------------------
# Uses Diverse Cohorts metadata

manifest <- rbind(
  manifest,
  data.frame(
    GENESIS_study = "GEN-B4",
    ADKP_study = "AMP-AD_DiverseCohorts",
    metadata_synid = syn_ids[["Diverse_Cohorts"]])
)


# Upload manifest file ---------------------------------------------------------

manifest <- manifest |> arrange(GENESIS_study)
manifest_file <- file.path("data", "output", "metadata_manifest.csv")
write.csv(manifest, manifest_file,
  row.names = FALSE, quote = FALSE
)
synapse_upload(manifest_file, UPLOAD_SYNID)


# Combine all harmonized data into one data frame

df_list <- apply(manifest, 1, function(m_row) {
  m_file <- synapse_download(m_row[["metadata_synid"]])
  meta <- read.csv(m_file$path) |>
    mutate(
      individualID = as.character(individualID),
      apoeGenotype = as.character(apoeGenotype),
      amyAny = as.character(amyAny),
      GENESIS_study = m_row[["GENESIS_study"]],
      ADKP_study = m_row[["ADKP_study"]]
    )
  return(meta)
}, simplify = FALSE)

df_all <- purrr::list_rbind(df_list)

validate_values(df_all, spec)

# df_all <- df_all |>
#  group_by(individualID) |>
#  mutate(genesis_study = paste(sort(genesis_study), collapse = "; "))

new_file <- write_metadata(df_all, "GENESIS_metadata_combined.csv")
synapse_upload(new_file, UPLOAD_SYNID)
