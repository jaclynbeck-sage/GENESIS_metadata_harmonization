library(synapser)
library(dplyr)
library(purrr)
library(stringr)
library(readxl)

spec <- config::get(file = "GENESIS_harmonization.yml")
source("util_functions.R")
source("dataset_specific_functions.R")

UPLOAD_SYNID <- "syn64759869"

syn_ids <- list(
  "Diverse_Cohorts" = "syn51757646.20",
  "MayoRNAseq" = "syn23277389.7",
  "MSBB" = "syn6101474.9",
  "ROSMAP" = "syn3191087.11",
  "SEA-AD" = "syn31149116.7"
)

check_new_versions(syn_ids)

harmonized_files <- c()

synLogin()

# MayoRNAseq -------------------------------------------------------------------

# GEN-A11 and GEN-A12 have > half their samples from the original MayoRNASeq
# metadata. They also share four individuals between them, all of which come
# from original Mayo metadata.
# GEN-A10 appears to have new samples from Mayo with no overlap to the original
# metadata.

meta_file <- synapse_download(syn_ids[["MayoRNAseq"]])
meta <- read.csv(meta_file$path)

colnames(meta)

print_qc(meta,
  isHispanic_col = "ethnicity",
  cerad_col = "CERAD",
  thal_col = "Thal"
)

meta_new <- harmonize_MayoRNASeq(meta, spec)

print_qc(meta_new)
validate_values(meta_new, spec)

new_filename <- write_metadata(meta_new, file.path("originals", meta_file$name))
harmonized_files <- c(harmonized_files, new_filename)


# MSBB -------------------------------------------------------------------------

# GEN-A1 has > 300 samples from the original MSBB metadata and 65 from ROSMAP.

meta_file <- synapse_download(syn_ids[["MSBB"]])
meta <- read.csv(meta_file$path)

colnames(meta)

print_qc(meta,
  isHispanic_col = "ethnicity",
  cerad_col = "CERAD"
)

meta_new <- harmonize_MSBB(meta, spec)

print_qc(meta_new)
validate_values(meta_new, spec)

new_filename <- write_metadata(meta_new, file.path("originals", meta_file$name))
harmonized_files <- c(harmonized_files, new_filename)


# ROSMAP -----------------------------------------------------------------------

# GEN-A2, GEN-A8, GEN-A13, and GEN-B6 all use original ROSMAP metadata

meta_file <- synapse_download(syn_ids[["ROSMAP"]])
meta <- read.csv(meta_file$path)

colnames(meta)

print_qc(meta,
  ageDeath_col = "age_death",
  sex_col = "msex",
  isHispanic_col = "spanish",
  apoe_col = "apoe_genotype",
  braak_col = "braaksc",
  cerad_col = "ceradsc"
)

meta_new <- harmonize_ROSMAP(meta, spec)

print_qc(meta_new)
validate_values(meta_new, spec)

new_filename <- write_metadata(meta_new, file.path("originals", meta_file$name))
harmonized_files <- c(harmonized_files, new_filename)


# Diverse Cohorts --------------------------------------------------------------

# Diverse Cohorts has data for samples from the original Mayo/MSBB/ROSMAP
# metadata but also has new samples. Metadata from Diverse Cohorts has already
# been harmonized using the same standards as GENESIS, minus some capitalization
# differences. Some of the old/original samples have additional information in
# the DC metadata that doesn't appear in the original files, but are missing
# information in DC that do appear in the original files, so we need to merge
# all four related metadata files at the end and do case-by-case harmonization
# where they disagree.

meta_file <- synapse_download(syn_ids[["Diverse_Cohorts"]])
meta <- read.csv(meta_file$path)

colnames(meta)

print_qc(meta, pmi_col = "PMI")

meta_new <- harmonize_Diverse_Cohorts(meta, spec)

print_qc(meta_new)
validate_values(meta_new, spec)

new_filename <- write_metadata(meta_new, file.path("originals", meta_file$name))
harmonized_files <- c(harmonized_files, new_filename)


# SEA-AD -----------------------------------------------------------------------

# GEN-A4 and GEN-B5 use SEA-AD data. The version of SEA-AD that is on Synapse is
# missing Hispanic/Latino information that is present in the version released by
# the Allen Institute on brain-map.org. We use the version on Synapse but pull
# in the missing information from the AI version.

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
new_filename <- write_metadata(meta_new, file.path("originals", file_write))
harmonized_files <- c(harmonized_files, new_filename)


# Merge all files into one data frame ------------------------------------------

meta_list <- lapply(harmonized_files, function(filename) {
  read.csv(filename) |>
    mutate(source_file = filename)
})

meta_all <- deduplicate_studies(meta_list, verbose = FALSE)

print_qc(meta_all)
validate_values(meta_all, spec)

# Re-write de-duplicated metadata files ----------------------------------------

# Use de-duplicated data but subset to only columns that exist in each
# individual metadata file
new_files <- sapply(meta_list, function(meta_old) {
  meta_new <- subset(meta_all, source_file == unique(meta_old$source_file)) |>
    select(all_of(colnames(meta_old))) |>
    # Remove the source file column we added
    select(-source_file)

  new_file <- str_replace(
    basename(unique(meta_old$source_file)),
    "_harmonized.csv",
    ".csv"
  )

  return(write_metadata(meta_new, new_file))
})

# Upload to GENESIS metadata space

# Temporary: Do not upload Mayo or MSBB
new_files <- new_files[!grepl("Mayo|MSBB", new_files)]

for (filename in new_files) {
  synapse_upload(filename, UPLOAD_SYNID)
}
