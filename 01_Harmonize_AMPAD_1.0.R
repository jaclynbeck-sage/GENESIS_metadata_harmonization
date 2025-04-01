# This script harmonizes and de-duplicates AMP-AD 1.0 metadata (MayoRNASeq,
# MSBB, and ROSMAP) and Diverse Cohorts metadata, filling in missing information
# where there is sample overlap between the four files. The de-duplicated data
# is then used as the standard to check quality of other overlapping GENESIS
# studies and fill in missing information in those studies as well.
#
# The harmonized ROSMAP 1.0 and Diverse Cohorts metadata files are also used
# directly by multiple GENESIS studies (GEN-A2, GEN-A8, GEN-A13, GEN-B6,
# GEN-B4), so they are uploaded to Synapse for GENESIS use.
#
# The MayoRNASeq 1.0 and MSBB 1.0 data sets are not used directly by any study
# so they are not uploaded to Synapse, to avoid confusion. However, several
# studies (GEN-A1, GEN-A11, GEN-A12) have samples from these data sets in their
# own metadata that are missing information. The missing information will be
# filled in by any available data from these 1.0 metadata files.
#
# TODO add provenance

library(synapser)
library(dplyr)
library(purrr)
library(stringr)

spec <- config::get(file = "GENESIS_harmonization.yml")
source("util_functions.R")
source("dataset_specific_functions.R")

UPLOAD_SYNID <- "syn64759869"

syn_ids <- list(
  "Diverse_Cohorts" = "syn51757646.20",
  "MayoRNAseq" = "syn23277389.7",
  "MSBB" = "syn6101474.9",
  "ROSMAP" = "syn3191087.11"
)

synLogin()
check_new_versions(syn_ids)

df_list <- list()

# MayoRNAseq -------------------------------------------------------------------

# GEN-A11 and GEN-A12 have > half their samples from the original MayoRNASeq
# metadata. They also share four individuals between them, all of which come
# from original Mayo metadata.

meta_file <- synapse_download(syn_ids[["MayoRNAseq"]])
meta <- read.csv(meta_file$path)

colnames(meta)

print_qc(meta,
  isHispanic_col = "ethnicity",
  cerad_col = "CERAD",
  thal_col = "Thal"
)

meta_new <- harmonize_MayoRNASeq(meta, spec) |>
  mutate(
    source = "MayoRNASeq",
    filename = meta_file$name
  )

print_qc(meta_new)
validate_values(meta_new, spec)

df_list[["MayoRNASeq"]] <- meta_new


# MSBB -------------------------------------------------------------------------

# GEN-A1 has > 300 samples from the original MSBB metadata

meta_file <- synapse_download(syn_ids[["MSBB"]])
meta <- read.csv(meta_file$path)

colnames(meta)

print_qc(meta,
  isHispanic_col = "ethnicity",
  cerad_col = "CERAD"
)

meta_new <- harmonize_MSBB(meta, spec) |>
  mutate(
    source = "MSBB",
    filename = meta_file$name
  )

print_qc(meta_new)
validate_values(meta_new, spec)

df_list[["MSBB"]] <- meta_new


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

meta_new <- harmonize_ROSMAP(meta, spec) |>
  mutate(
    source = "ROSMAP",
    filename = meta_file$name
  )

print_qc(meta_new)
validate_values(meta_new, spec)

df_list[["ROSMAP"]] <- meta_new


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

meta_new <- harmonize_Diverse_Cohorts(meta, spec) |>
  mutate(
    source = "Diverse Cohorts",
    filename = meta_file$name
  )

print_qc(meta_new)
validate_values(meta_new, spec)

df_list[["Diverse_Cohorts"]] <- meta_new


# Merge all files into one data frame ------------------------------------------

meta_all <- deduplicate_studies(df_list, verbose = FALSE)

print_qc(meta_all)
validate_values(meta_all, spec)

# Save de-duplicated data frame for step 2
saveRDS(meta_all, file.path("data", "tmp", "AMP1.0_DiverseCohorts_harmonized.rds"))


# Re-write de-duplicated metadata files ----------------------------------------

# Use de-duplicated data but subset to only columns that exist in each
# individual metadata file
new_files <- sapply(df_list, function(meta_old) {
  meta_new <- subset(meta_all, filename == unique(meta_old$filename)) |>
    select(all_of(colnames(meta_old))) |>
    # Remove the source file columns we added
    select(-source, -filename)

  new_file <- str_replace(
    basename(unique(meta_old$filename)),
    "_harmonized.csv",
    ".csv"
  )

  return(write_metadata(meta_new, new_file))
})

# Upload to GENESIS metadata space

# Note: Do not upload Mayo or MSBB
new_files <- new_files[!grepl("Mayo|MSBB", new_files)]

for (filename in new_files) {
  synapse_upload(filename, UPLOAD_SYNID)
}
