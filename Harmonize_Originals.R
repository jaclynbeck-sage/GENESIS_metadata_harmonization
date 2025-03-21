library(synapser)
library(dplyr)
library(purrr)
library(stringr)
library(readxl)

spec <- config::get(file = "GENESIS_harmonization.yml")
source("util_functions.R")
source("dataset_specific_functions.R")

syn_ids <- list(
  "Diverse_Cohorts" = "syn51757646.20",
  "MayoRNAseq" = "syn23277389.7",
  "MSBB" = "syn6101474.9",
  "ROSMAP" = "syn3191087.11",
  "SEA-AD" = "syn31149116.7"
)

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
    mutate(
      individualID = as.character(individualID),
      apoeGenotype = as.character(apoeGenotype),
      amyAny = as.character(amyAny),
      source_file = filename,
      # Special case: MSBB/MSSM samples were re-named for Diverse Cohorts, this
      # makes them directly comparable
      original_individualID = individualID,
      individualID = str_replace(individualID, "AMPAD_MSSM_0+", "")
    )
})
meta_all <- purrr::list_rbind(meta_list)

print_qc(meta_all)
validate_values(meta_all, spec)

# Resolve duplicates -----------------------------------------------------------
# Only for columns that are harmonized

# This accounts for overlapping IDs between different studies that don't refer
# to the same individual
dupe_ids <- meta_all |>
  select(individualID, dataContributionGroup) |>
  mutate(group_id = paste(individualID, dataContributionGroup),
         duplicate = duplicated(group_id)) |>
  subset(duplicate == TRUE) |>
  distinct()

# For each ID that has duplicate rows, resolve duplicates:
#   1. For columns where some rows have NA and some have a unique non-NA value,
#      replace the NA value with that unique non-NA value.
#   2. For columns where some rows have "missing or unknown" and some have a
#      unique value other than that, replace "missing or unknown" with the
#      unique value.
#   3. For the ageDeath/pmi columns where rows disagree because of precision,
#      use the most precise value.
for (row_id in 1:nrow(dupe_ids)) {
  ind_id <- dupe_ids$individualID[row_id]

  # This will be altered to resolve duplication, and will get added back to the
  # meta_all data frame
  meta_tmp <- subset(meta_all, individualID == ind_id &
                       dataContributionGroup == dupe_ids$dataContributionGroup[row_id])

  for (col_name in expectedColumns) {
    unique_vals <- unique(meta_tmp[, col_name])

    if (length(unique_vals) > 1) {
      #cat(ind_id, col_name, "[", paste(unique_vals, collapse = ", "), "]\n")

      # Use most precise ageDeath or pmi value
      if (col_name %in% c("ageDeath", "pmi")) {
        n_decimals <- str_replace(as.character(unique_vals), ".*\\.", "") |>
          nchar()
        meta_tmp[, col_name] <- unique_vals[which.max(n_decimals)]
      } else if (col_name == "cohort") {
        # If there is more than one left over value, report it but don't try to
        # resolve duplication. Special case: Mayo data may be labeled as "Mayo
        # Clinic" in the original Mayo metadata or "Banner" in Diverse Cohorts,
        # but we don't need to change this in either metadata file or print it
        # out.
        is_mayo_banner <- identical(sort(unique_vals),
                                    c("Banner", "Mayo Clinic"))
        if (!is_mayo_banner) {
          cat(ind_id, col_name, "[", paste(unique_vals, collapse = ", "), "]\n")
        }
      } else {
        # Remove NA and "missing or unknown" values to see what's left
        leftover <- unique_vals[!is.na(unique_vals) & (unique_vals != spec$missing)]

        if (length(leftover) == 1) {
          # One unique left over value
          meta_tmp[, col_name] <- leftover
        } else if (length(leftover) == 0) {
          # Nothing left, use "missing or unknown" if it's there, otherwise set
          # to NA.
          if (any(unique_vals == spec$missing)) {
            meta_tmp[, col_name] <- spec$missing
          } else {
            meta_tmp[, col_name] <- NA
          }
        } else {
          # If there is more than one left over value, report it but don't try
          # to resolve duplication.
          cat(ind_id, col_name, "[", paste(unique_vals, collapse = ", "), "]\n")
        }
      }
    }
  }

  # Replace original rows with de-duplicated data
  rows_replace <- meta_all$individualID == ind_id &
    meta_all$dataContributionGroup == unique(meta_tmp$dataContributionGroup)
  meta_all[rows_replace, ] <- meta_tmp
}

print_qc(meta_all)
validate_values(meta_all, spec)

# Re-write de-duplicated metadata files ----------------------------------------

# Use de-duplicated data but subset to only columns that exist in each
# individual metadata file
new_files <- lapply(meta_list, function(meta_old) {
  meta_new <- subset(meta_all, source_file == unique(meta_old$source_file)) |>
    select(all_of(colnames(meta_old))) |>
    # Undo any modifications we did to the individual ID (only matters for MSBB)
    mutate(individualID = original_individualID) |>
    select(-source_file, -original_individualID)

  new_file <- str_replace(unique(meta_old$source_file), ".csv",
                          "_deduplicated.csv")
  write.csv(meta_new, new_file)
  return(new_file)
})
