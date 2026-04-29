# This script harmonizes and de-duplicates AMP-AD 1.0 metadata (MSBB and ROSMAP)
# with Diverse Cohorts and NPS-AD metadata, filling in missing information where
# there is sample overlap between the four files. The de-duplicated data is then
# used as the standard to check quality of other overlapping GENESIS studies and
# fill in missing information in those studies as well.
#
# The harmonized ROSMAP 1.0, Diverse Cohorts, and NPS-AD metadata files are also
# used directly by multiple GENESIS studies (GEN-A2, GEN-A8, GEN-A13, GEN-B6,
# GEN-B4), so they are uploaded to Synapse for GENESIS use.
#
# The MSBB 1.0 data set is not used directly by any study, so it is not uploaded
# to Synapse, to avoid confusion. However, several studies (GEN-A1, GEN-A11,
# GEN-A12) have samples from this data set in their own metadata that are
# missing information. The missing information will be filled in by any
# available data from the MSBB 1.0 metadata file.
#
# TODO add provenance

library(synapser)
library(dplyr)
library(purrr)
library(stringr)

spec <- config::get(file = "GENESIS_harmonization.yml")
source("util_functions.R")
source("dataset_specific_functions.R")

# All files are from the ADKP Harmonization Project, so they have been
# harmonized and de-duplicated already and just need some minor edits to comply
# with the GENESIS data dictionary
syn_ids <- list(
  "Diverse_Cohorts" = "syn73713769.3",
  "MSBB" = "syn73713767.1",
  "ROSMAP" = "syn73713768.1",
  "NPS-AD" = "syn73713770.2"
)

synLogin()
check_new_versions(syn_ids)

df_list <- list()


# MSBB -------------------------------------------------------------------------

# NPS-AD has > 300 samples from the original MSBB metadata

meta_file <- synapse_download(syn_ids[["MSBB"]])
meta <- read.csv(meta_file$path)

colnames(meta)

print_qc(meta, braak_nft_col = "Braak", bscore_nft_col = "bScore")

meta_new <- harmonize(spec$study$msbb, meta, spec) |>
  mutate(filename = meta_file$name)

print_qc(meta_new)
validate_values(meta_new, spec)

df_list[["MSBB"]] <- meta_new


# ROSMAP -----------------------------------------------------------------------

# GEN-A2, GEN-A8, GEN-A13, and GEN-B6 all use original ROSMAP metadata

meta_file <- synapse_download(syn_ids[["ROSMAP"]])
meta <- read.csv(meta_file$path)

colnames(meta)

print_qc(meta, braak_nft_col = "Braak", bscore_nft_col = "bScore")

meta_new <- harmonize(spec$study$rosmap, meta, spec) |>
  mutate(filename = meta_file$name)

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

print_qc(meta, braak_nft_col = "Braak", bscore_nft_col = "bScore")

meta_new <- harmonize(spec$study$diverse_cohorts, meta, spec) |>
  mutate(filename = meta_file$name)

print_qc(meta_new)
validate_values(meta_new, spec)

df_list[["Diverse_Cohorts"]] <- meta_new


# NPS-AD -----------------------------------------------------------------------
# This dataset has overlap with MSBB, ROSMAP, and Diverse Cohorts. Some values
# from MSBB 1.0 have been corrected in this data, and some have not.

meta_file <- synapse_download(syn_ids[["NPS-AD"]])
meta <- read.csv(meta_file$path)

colnames(meta)

print_qc(meta, braak_nft_col = "Braak", bscore_nft_col = "bScore")

meta_new <- harmonize(spec$study$nps_ad, meta, spec) |>
  mutate(filename = meta_file$name)

print_qc(meta_new)
validate_values(meta_new, spec)

df_list[["NPS-AD"]] <- meta_new


# Merge all files into one data frame ------------------------------------------

meta_all <- deduplicate_studies(df_list, spec, verbose = FALSE)

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
    # Remove the source file column we added
    select(-filename)

  return(write_metadata(meta_new, basename(unique(meta_old$filename))))
})

# Upload to GENESIS metadata space

# Note: Do not upload MSBB
new_files <- new_files[!grepl("MSBB", new_files)]

for (filename in new_files) {
  synapse_upload(filename, spec$upload_synID)
}
