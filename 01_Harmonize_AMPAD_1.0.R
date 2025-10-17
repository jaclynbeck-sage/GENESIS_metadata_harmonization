# This script harmonizes and de-duplicates AMP-AD 1.0 metadata (MayoRNAseq,
# MSBB, and ROSMAP) with Diverse Cohorts and NPS-AD metadata, filling in missing
# information where there is sample overlap between the five files. The
# de-duplicated data is then used as the standard to check quality of other
# overlapping GENESIS studies and fill in missing information in those studies
# as well.
#
# Note: There are multiple metadata errors in the MSBB 1.0 data, which have
# been corrected either in NPS-AD data, Diverse Cohorts data, or an additional
# file of corrections. The correct data is propagated to all data using these
# samples (MSBB 1.0, Diverse Cohorts, and NPS-AD), which is why this script
# includes the NPS-AD data instead of putting it in the GENESIS harmonization
# script.
#
# The harmonized ROSMAP 1.0, Diverse Cohorts, and NPS-AD metadata files are also
# used directly by multiple GENESIS studies (GEN-A2, GEN-A8, GEN-A13, GEN-B6,
# GEN-B4), so they are uploaded to Synapse for GENESIS use.
#
# The MayoRNAseq 1.0 and MSBB 1.0 data sets are not used directly by any study
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

UPLOAD_SYNID <- "syn65931571"

syn_ids <- list(
  "Diverse_Cohorts" = "syn51757646.21",
  "MayoRNAseq" = "syn23277389.7",
  "MSBB" = "syn6101474.10",
  "ROSMAP" = "syn3191087.11",
  "NPS-AD" = "syn55251012.4", # intentionally one version behind due to v5 being a copy of the harmonized file generated here
  "NPS-AD_neuropath" = "syn55251003.1",
  # Access restricted to Sage internal
  "MSBB_corrections" = "syn66511661.2"
)

synLogin()
check_new_versions(syn_ids)

df_list <- list()

# MayoRNAseq -------------------------------------------------------------------

# GEN-A11 and GEN-A12 have > half their samples from the original MayoRNAseq
# metadata. They also share four individuals between them, all of which come
# from original Mayo metadata.

meta_file <- synapse_download(syn_ids[["MayoRNAseq"]])
meta <- read.csv(meta_file$path)

colnames(meta)

print_qc(meta,
  pmi_col = "pmi",
  isHispanic_col = "ethnicity",
  cerad_col = "CERAD",
  thal_col = "Thal",
  braak_nft_col = "Braak"
)

meta_new <- harmonize(spec$study$mayo, meta, spec) |>
  mutate(filename = meta_file$name)

print_qc(meta_new)
validate_values(meta_new, spec)

df_list[["MayoRNAseq"]] <- meta_new


# MSBB -------------------------------------------------------------------------

# NPS-AD has > 300 samples from the original MSBB metadata

meta_file <- synapse_download(syn_ids[["MSBB"]])
meta <- read.csv(meta_file$path)

colnames(meta)

print_qc(meta,
  pmi_col = "pmi",
  isHispanic_col = "ethnicity",
  cerad_col = "CERAD",
  braak_nft_col = "Braak"
)

meta_new <- harmonize(spec$study$msbb, meta, spec) |>
  mutate(filename = meta_file$name)

print_qc(meta_new)
validate_values(meta_new, spec)

df_list[["MSBB"]] <- meta_new


# Corrected MSBB samples -------------------------------------------------------

meta_file <- synapse_download(syn_ids[["MSBB_corrections"]])
meta <- readxl::read_xlsx(meta_file$path) |>
  as.data.frame()

colnames(meta)

print_qc(meta,
         ageDeath_col = "Age",
         race_col = "RaceLabel",
         sex_col = "SexLabel",
         pmi_col = "PMI (min)",
         apoe_col = "ApoE",
         cerad_col = "CERAD_1",
         braak_nft_col = "B&B Alz"
)

meta_new <- harmonize("MSBB_corrections", meta, spec)

print_qc(meta_new)

validate_values(meta_new, spec)

df_list[["MSBB_corrections"]] <- meta_new


# ROSMAP -----------------------------------------------------------------------

# GEN-A2, GEN-A8, GEN-A13, and GEN-B6 all use original ROSMAP metadata

meta_file <- synapse_download(syn_ids[["ROSMAP"]])
meta <- read.csv(meta_file$path)

colnames(meta)

print_qc(meta,
  ageDeath_col = "age_death",
  pmi_col = "pmi",
  sex_col = "msex",
  isHispanic_col = "spanish",
  apoe_col = "apoe_genotype",
  braak_nft_col = "braaksc",
  cerad_col = "ceradsc"
)

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
# Special case: There is an additional file with neuropathology data that should
# be pulled into the individual metadata.

meta_file <- synapse_download(syn_ids[["NPS-AD"]])
meta <- read.csv(meta_file$path)

neuro_file <- synapse_download(syn_ids[["NPS-AD_neuropath"]])
neuropath <- read.csv(neuro_file$path) |>
  dplyr::rename(individualID = IndividualID)

colnames(meta)

print_qc(meta, isHispanic_col = "ethnicity", cerad_col = "CERAD")
print_qc(neuropath, cerad_col = "CERAD", braak_nft_col = "BRAAK_AD")

meta_new <- harmonize(spec$study$nps_ad, meta, spec, extra_metadata = neuropath) |>
  mutate(filename = meta_file$name)

print_qc(meta_new)
validate_values(meta_new, spec) # Cohort will fail until after harmonization

df_list[["NPS-AD"]] <- meta_new


# Merge all files into one data frame ------------------------------------------

meta_all <- deduplicate_studies(df_list, spec, verbose = FALSE)

meta_all <- subset(meta_all, study != "MSBB_corrections") |>
  fill_missing_ampad1.0_ids(spec) |> # For Diverse Cohorts
  updateADoutcome(spec) # For Diverse Cohorts

print_qc(meta_all)
validate_values(meta_all, spec)

# Save de-duplicated data frame for step 2
saveRDS(meta_all, file.path("data", "tmp", "AMP1.0_DiverseCohorts_harmonized.rds"))


# Re-write de-duplicated metadata files ----------------------------------------

# Use de-duplicated data but subset to only columns that exist in each
# individual metadata file
df_list[["MSBB_corrections"]] <- NULL # Not needed at this point

new_files <- sapply(df_list, function(meta_old) {
  meta_new <- subset(meta_all, filename == unique(meta_old$filename)) |>
    select(all_of(colnames(meta_old))) |>
    # Remove the source file column we added
    select(-filename)

  new_file <- str_replace(
    basename(unique(meta_old$filename)),
    "_harmonized.csv",
    ".csv"
  )

  # Sorting for this data set gets changed during fill_missing_ampad1.0_ids() so
  # we sort it back. Also fill in missing "race" values with "geneticAncestry"
  # values.
  if (unique(meta_new$study) == "NPS-AD") {
    meta_new <- arrange(meta_new, individualID) #|>
      #mutate(race = ifelse(race == spec$missing, geneticAncestry, race)) # TODO
  }

  return(write_metadata(meta_new, new_file))
})

# Upload to GENESIS metadata space

# Note: Do not upload Mayo or MSBB
new_files <- new_files[!grepl("Mayo|MSBB", new_files)]

for (filename in new_files) {
  synapse_upload(filename, UPLOAD_SYNID)
}
