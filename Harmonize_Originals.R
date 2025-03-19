library(synapser)
library(dplyr)
library(purrr)
library(stringr)
library(readxl)

spec <- config::get(file = "GENESIS_harmonization.yml")
source("util_functions.R")
source("dataset_specific_functions.R")

syn_ids <- list("Diverse_Cohorts" = "syn51757646.20",
                "MayoRNAseq" = "syn23277389.7",
                "MSBB" = "syn6101474.9",
                "ROSMAP" = "syn3191087.11",
                "SEA-AD" = "syn31149116.7")

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
         thal_col = "Thal")

meta_new <- meta %>%
  rename(isHispanic = ethnicity,
         amyCerad = CERAD,
         amyThal = Thal) %>%
  mutate(
    ageDeath = censor_ages(ageDeath, spec),
    isHispanic = case_when(isHispanic == "Caucasian" ~ spec$isHispanic$hisp_false,
                           is.na(isHispanic) ~ spec$missing,
                           .default = isHispanic),
    race = case_when(is.na(race) ~ spec$missing,
                     .default = race),
    sex = case_when(is.na(sex) ~ spec$missing,
                    .default = sex),
    apoeGenotype = case_when(is.na(apoeGenotype) ~ spec$missing,
                             .default = as.character(apoeGenotype)),
    apoe4Status = get_apoe4Status(apoeGenotype, spec),
    amyCerad = case_when(is.na(amyCerad) ~ spec$missing,
                         .default = amyCerad),
    amyAny = get_amyAny(amyCerad, spec),
    amyThal = case_when(is.na(amyThal) ~ spec$missing,
                        amyThal == 0 ~ spec$amyThal$none,
                        .default = paste("Phase", amyThal)),
    amyA = get_amyA(amyThal, spec),
    Braak = case_when(is.na(Braak) ~ spec$missing,
                      Braak >= 0 ~ to_Braak_stage(floor(Braak), spec),
                      .default = as.character(Braak)),
    bScore = get_bScore(Braak, spec),
    cohort = "Mayo Clinic",
    dataContributionGroup = "Mayo"
  )

print_qc(meta_new)

new_filename <- write_metadata(meta_new, file.path("originals", meta_file$name))
harmonized_files <- c(harmonized_files, new_filename)


# MSBB -------------------------------------------------------------------------

# GEN-A1 has > 300 samples from the original MSBB metadata and 65 from ROSMAP.

meta_file <- synapse_download(syn_ids[["MSBB"]])
meta <- read.csv(meta_file$path)

colnames(meta)

print_qc(meta,
         isHispanic_col = "ethnicity",
         cerad_col = "CERAD")

meta_new <- meta %>%
  rename(isHispanic = ethnicity,
         amyCerad = CERAD,
         dataContributionGroup = individualIdSource) %>%
  mutate(
    ageDeath = censor_ages(ageDeath, spec),
    pmi = pmi / 60, # PMI is in minutes for MSBB
    isHispanic = case_when(is.na(isHispanic) ~ spec$missing,
                           isHispanic %in% c("A", "B", "W") ~ spec$isHispanic$hisp_false,
                           isHispanic == "H" ~ spec$isHispanic$hisp_true,
                           isHispanic == "U" ~ spec$missing,
                           .default = isHispanic),
    race = case_when(is.na(race) ~ spec$missing,
                     race == "A" ~ spec$race$Asian,
                     race == "B" ~ spec$race$Black,
                     race == "H" ~ spec$race$other,
                     race == "W" ~ spec$race$White,
                     race == "U" ~ spec$missing,
                     .default = race),
    apoeGenotype = case_when(is.na(apoeGenotype) ~ spec$missing,
                             .default = as.character(apoeGenotype)),
    apoe4Status = get_apoe4Status(apoeGenotype, spec),
    amyCerad = case_when(is.na(amyCerad) ~ spec$missing,
                         amyCerad == 1 ~ spec$amyCerad$none,
                         amyCerad == 2 ~ spec$amyCerad$frequent,
                         amyCerad == 3 ~ spec$amyCerad$moderate,
                         amyCerad == 4 ~ spec$amyCerad$sparse,
                         .default = as.character(amyCerad)),
    amyAny = get_amyAny(amyCerad, spec),
    amyThal = spec$missing,
    amyA = get_amyA(amyThal, spec),
    Braak = case_when(is.na(Braak) ~ spec$missing,
                      Braak >= 0 ~ to_Braak_stage(floor(Braak), spec),
                      .default = as.character(Braak)),
    bScore = get_bScore(Braak, spec),
    cohort = "Mt Sinai Brain Bank"
  )

print_qc(meta_new)

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
         cerad_col = "ceradsc")

meta_new <- harmonize_ROSMAP(meta, spec)

print_qc(meta_new)

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
              destfile = sea_ad_file)

meta <- read_xlsx(meta_file$path)
meta_sea_ad <- read_xlsx(sea_ad_file)

colnames(meta)
colnames(meta_sea_ad)

print_qc(meta,
         cerad_col = "CERAD",
         thal_col = "Thal phase")

print_qc(meta_sea_ad,
         isHispanic_col = "Hispanic/Latino")

meta_new <- harmonize_SEA_AD(meta, meta_sea_ad, spec)

print_qc(meta_new)

# Original file is an Excel file, change filename to CSV file
file_write <- str_replace(meta_file$name, "xlsx", "csv")
new_filename <- write_metadata(meta_new, file.path("originals", file_write))
harmonized_files <- c(harmonized_files, new_filename)


# Merge all files into one data frame ------------------------------------------

meta_all <- lapply(harmonized_files, function(filename) {
  read.csv(filename) %>%
    mutate(individualID = as.character(individualID),
           apoeGenotype = as.character(apoeGenotype),
           amyAny = as.character(amyAny))
})
meta_all <- purrr::list_rbind(meta_all) %>%
  distinct()

print_qc(meta_all)

dupe_ids <- meta_all$individualID[duplicated(meta_all$individualID)] %>%
  unique() %>%
  # Special case: These two IDs overlap between Mayo and MSBB and should not
  # be considered duplicates
  setdiff(c("1943", "12047"))

# For each ID that has duplicate rows, resolve duplicates:
#   1. For columns where some rows have NA and some have a unique non-NA value,
#      replace the NA value with that unique non-NA value.
#   2. For columns where some rows have "missing or unknown" and some have a
#      unique value other than that, replace "missing or unknown" with the
#      unique value.
#   3. For the ageDeath/pmi columns where rows disagree because of precision,
#      use the most precise value.
for (ind_id in dupe_ids) {
  meta_tmp <- subset(meta_all, individualID == ind_id)

  # This will be altered to resolve duplication, and will get added back to the
  # meta_all data frame
  meta_replace <- meta_tmp[1, ]

  for (col_name in colnames(meta_tmp)) {
    unique_vals <- unique(meta_tmp[, col_name])

    if (length(unique_vals) > 1) {
      #print(paste(ind_id, col_name, "[", paste(unique_vals, collapse = ", "), "]"))

      # Use most precise ageDeath or pmi value
      if (col_name %in% c("ageDeath", "pmi")) {
        n_decimals <- str_replace(as.character(unique_vals), ".*\\.", "") %>%
          nchar()
        meta_replace[, col_name] <- unique_vals[which.max(n_decimals)]

      } else if (col_name == "cohort") {
        # Special case: Mayo data may be labeled as "Mayo Clinic" or "Banner"
        # but be identical otherwise
        if (all(c("Mayo Clinic", "Banner") %in% unique_vals)) {
          meta_replace[, col_name] <- "Banner"
        } else {
          print(paste(ind_id, col_name, "[", paste(unique_vals, collapse = ", "), "]"))
        }
      } else {
        # Remove NA and "missing or unknown" values to see what's left
        leftover <- unique_vals[!is.na(unique_vals) & (unique_vals != spec$missing)]

        if (length(leftover) == 1) {
          # One unique left over value
          meta_replace[, col_name] <- leftover

        } else if (length(leftover) == 0) {
          # Nothing left, use "missing or unknown" if it's there, otherwise set
          # to NA.
          if (any(unique_vals == spec$missing)) {
            meta_replace[, col_name] <- spec$missing
          } else {
            meta_replace[, col_name] <- NA
          }

        } else {
          # There is more than one left over value. Need to inspect.
          print(paste(ind_id, col_name, "[", paste(unique_vals, collapse = ", "), "]"))
        }
      }
    }
  }

  meta_all[meta_all$individualID == ind_id, ] <- meta_replace
}
