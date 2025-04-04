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
  "GEN-A1_neuropath" = "syn55251003.1", # NPS-AD neuropathology
  "ROSMAP" = "syn64759878.4", # Harmonized file, for GEN-A2, GEN-A8, GEN-A13, GEN-B6
  "SEA-AD" = "syn31149116.7", # SEA-AD, for GEN-A4 and GEN-B5
  "GEN-A9" = "syn22432749.1", # SMIB-AD
  "GEN-A10" = "syn25891193.1", # MCMPS
  "GEN-A11" = "syn31563038.1", # MC_snRNA
  "GEN-A12" = "syn51401700.2", # MC-BrAD
  "Diverse_Cohorts" = "syn64759872.5" # Harmonized file, for GEN-B4
  # "GEN-B1" = "TBD",
  # "GEN-B2" = "TBD",
  # "GEN-B3" = "TBD"
)

synLogin()
check_new_versions(syn_ids)

UPLOAD_SYNID <- "syn64759869"

manifest <- c()

# Assume this file has been generated by Step 1.
harmonized_baseline <- readRDS(file.path("data", "tmp", "AMP1.0_DiverseCohorts_harmonized.rds"))


# GEN-A1 -----------------------------------------------------------------------
# Has sample overlap with MSBB, ROSMAP, and Diverse Cohorts.
# Special case: There is an additional file with neuropathology data that should
# be pulled into the individual metadata.

meta_file <- synapse_download(syn_ids[["GEN-A1"]])
meta <- read.csv(meta_file$path)

neuro_file <- synapse_download(syn_ids[["GEN-A1_neuropath"]])
neuropath <- read.csv(neuro_file$path) |>
  dplyr::rename(individualID = IndividualID)

colnames(meta)

print_qc(meta, isHispanic_col = "ethnicity", pmi_col = "PMI", cerad_col = "CERAD")
print_qc(neuropath, cerad_col = "CERAD", braak_col = "BRAAK_AD")

meta_new <- harmonize_NPS_AD(meta, neuropath, harmonized_baseline, spec)

print_qc(meta_new)
validate_values(meta_new, spec)

new_filename <- write_metadata(meta_new, meta_file$name)
new_syn_id <- synapse_upload(new_filename, UPLOAD_SYNID)

manifest <- rbind(
  manifest,
  data.frame(
    GENESIS_study = "GEN-A1",
    ADKP_study = "NPS-AD",
    metadata_synid = paste0(new_syn_id$id, ".", new_syn_id$versionNumber)
  )
)


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

print_qc(meta, cerad_col = "CERAD", thal_col = "Thal phase")
print_qc(meta_sea_ad, isHispanic_col = "Hispanic/Latino")

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
    metadata_synid = paste0(new_syn_id$id, ".", new_syn_id$versionNumber)
  )
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
    metadata_synid = paste0(new_syn_id$id, ".", new_syn_id$versionNumber)
  )
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
    metadata_synid = paste0(new_syn_id$id, ".", new_syn_id$versionNumber)
  )
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
    metadata_synid = paste0(new_syn_id$id, ".", new_syn_id$versionNumber)
  )
)


# GEN-B4 -----------------------------------------------------------------------
# Uses Diverse Cohorts metadata

manifest <- rbind(
  manifest,
  data.frame(
    GENESIS_study = "GEN-B4",
    ADKP_study = "AMP-AD_DiverseCohorts",
    metadata_synid = syn_ids[["Diverse_Cohorts"]]
  )
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
