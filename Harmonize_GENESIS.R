# This script harmonizes and de-duplicates multiple studies used for GENESIS,
# filling in missing information where there is sample overlap between data
# sets. Harmonized metadata is uploaded to Synapse for use.
#
# The following data sets are harmonized in this script:
#   GEN-A1 / NPS-AD
#   GEN-A2 / ROSMAP
#   GEN-A3 / AMP-PD
#   GEN-A4 / SEA-AD
#   GEN-A5 / CMC
#   GEN-A6 / SZBDMulti-Seq
#   GEN-A8 / snRNAseqAD_Trem2
#   GEN-A9 / SMIB-AD
#   GEN-A10 / MCMPS
#   GEN-A11 / MC_snRNA
#   GEN-A12 / MC-BrAD
#   GEN-A13 / snRNAseqPFC_BA10
#   GEN-A15 / ASAP
#   GEN-A16 / McCarroll_SCZ
#   GEN-A17 / McCarroll_HD
#   GEN-A18 / TargetALS
#   GEN-B4 / AMP-AD_DiverseCohorts
#   GEN-B5 / SEA-AD_multiome
#   GEN-B6 / MIT_ROSMAP_Multiomics
#   GEN-B8 / BD2
#   GEN-B11 / ROSMAP_CUIMC2

library(synapser)
library(dplyr)
library(purrr)
library(stringr)
library(readxl)

spec <- config::get(file = "GENESIS_column_spec.yml")
studies <- config::get(file = "GENESIS_study_spec.yml")

# Source all the helper function scripts
helper_functions <- list.files("functions", pattern = "\\.R", full.names = TRUE)
for (file in helper_functions) {
  source(file)
}

# TRUE will print column names of metadata + unique values for each dataset.
# FALSE will only print the status of the harmonized metadata for each dataset.
verbose <- FALSE

syn_ids <- sapply(studies, "[[", "syn_id") |> unlist() |> na.omit()
check_new_versions(syn_ids)

manifest <- c()
datasets <- list()


# GEN-A1 / NPS-AD --------------------------------------------------------------

nps_ad <- studies$nps_ad
meta_file <- synapse_download(nps_ad$syn_id)
meta <- read.csv(meta_file$path)

if (verbose) {
  colnames(meta)
  print_summary(meta, braak_nft_col = "Braak", bscore_nft_col = "bScore")
}

meta_new <- harmonize(nps_ad, meta, spec)

if (verbose) {
  print_summary(meta_new)
}

cat("\n", nps_ad$gen_name, "/", nps_ad$name, "\n")
validate_values(meta_new, spec)

datasets[[nps_ad$name]] <- list(
  data = meta_new,
  study = nps_ad,
  filename = meta_file$name
)


# GEN-A2, GEN-A8, GEN-A13, GEN-B6, GEN-B11 / ROSMAP ----------------------------
# These GENESIS studies all use original ROSMAP metadata

rosmap <- studies$rosmap
meta_file <- synapse_download(rosmap$syn_id)
meta <- read.csv(meta_file$path)

if (verbose) {
  colnames(meta)
  print_summary(meta, braak_nft_col = "Braak", bscore_nft_col = "bScore")
}

meta_new <- harmonize(rosmap, meta, spec)

if (verbose) {
  print_summary(meta_new)
}

cat("\n", rosmap$gen_name, "/", rosmap$name, "\n")
validate_values(meta_new, spec)

datasets[[rosmap$name]] <- list(
  data = meta_new,
  study = rosmap,
  assoc_studies = list(studies$snRNA_trem2, studies$snRNA_BA10,
                       studies$mit_rosmap, studies$cuimc2),
  filename = meta_file$name
)


# GEN-A3 / AMP-PD --------------------------------------------------------------

amp_pd <- studies$amp_pd

# AMP-PD uses locally-downloaded files that are not on Synapse
if (!file.exists(amp_pd$local_file$main)) {
  warning(str_glue("{amp_pd$name} file {amp_pd$local_file$main} doesn't ",
                   "exist! This dataset will be excluded from harmonization."))
} else if (!file.exists(amp_pd$local_file$demographics)) {
  warning(str_glue("{amp_pd$name} file {amp_pd$local_file$demographics} doesn't ",
                   "exist! This dataset will be excluded from harmonization."))
} else {
  meta <- read.csv(amp_pd$local_file$main)
  demographics <- read.csv(amp_pd$local_file$demographics) |>
    select(participant_id, race, ethnicity)

  meta <- merge(meta, demographics)

  if (verbose) {
    # colnames(meta) # Don't print, there are > 600 columns

    print_summary(
      meta,
      ageDeath_col = "Demographics.age_at_baseline",
      sex_col = "Demographics.sex",
      pmi_col = "PMI.PMI_hours",
      isHispanic_col = "ethnicity",
      braak_nft_col = "LBD_Cohort_Path_Data.path_braak_nft",
      braak_lb_col = "LBD_Cohort_Path_Data.path_braak_lb",
      cerad_col = "LBD_Cohort_Path_Data.path_cerad",
      bscore_lb_col = "Info.BraakGroup"
    )
  }

  meta_new <- harmonize(amp_pd, meta, spec)

  if (verbose) {
    print_summary(meta_new)
  }

  cat("\n", amp_pd$gen_name, "/", amp_pd$name, "\n")
  validate_values(meta_new, spec)

  datasets[[amp_pd$name]] <- list(
    data = meta_new,
    study = amp_pd,
    filename = basename(amp_pd$local_file$main)
  )
}


# GEN-A4, GEN-B5 / SEA-AD ------------------------------------------------------

sea_ad <- studies$sea_ad
meta_file <- synapse_download(sea_ad$syn_id)
meta <- read.csv(meta_file$path)

if (verbose) {
  colnames(meta)

  print_summary(meta, braak_nft_col = "Braak", bscore_nft_col = "bScore")
}

meta_new <- harmonize(sea_ad, meta, spec)

if (verbose) {
  print_summary(meta_new)
}

cat("\n", sea_ad$gen_name, "/", sea_ad$name, "\n")
validate_values(meta_new, spec)

datasets[[sea_ad$name]] <- list(
  data = meta_new,
  study = sea_ad,
  assoc_studies = list(studies$sea_ad_multi),
  filename = meta_file$name
)


# GEN-A5 / CMC (PsychENCODE) ---------------------------------------------------

# NOTE: There are overlaps between CMC and NPS-AD individuals, but the way they
# are named between the two studies is different. Per Jaro, we will wait until
# he compiles a list of overlapping samples to de-duplicate.

cmc <- studies$cmc
if (!file.exists(cmc$local_file)) {
  warning(str_glue("{cmc$name} file {cmc$local_file} doesn't exist! ",
                   "This dataset will be excluded from harmonization."))
} else {
  meta <- read.csv(cmc$local_file)

  if (verbose) {
    colnames(meta)

    print_summary(meta, sex_col = "reportedGender",
                  isHispanic_col = "ethnicity", braak_nft_col = "Braak")
  }

  meta_new <- harmonize(cmc, meta, spec)

  if (verbose) {
    print_summary(meta_new)
  }

  cat("\n", cmc$gen_name, "/", cmc$name, "\n")
  validate_values(meta_new, spec)

  datasets[[cmc$name]] <- list(
    data = meta_new,
    study = cmc,
    filename = basename(cmc$local_file)
  )
}


# GEN-A6 / SZBDMulti-Seq (PsychENCODE) -----------------------------------------

# NOTE: There are overlaps between CMC and NPS-AD individuals, but the way they
# are named between the two studies is different. Per Jaro, we will wait until
# he compiles a list of overlapping samples to de-duplicate.

szbd <- studies$szbd
if (!file.exists(szbd$local_file)) {
  warning(str_glue("{szbd$name} file {szbd$local_file} doesn't exist! ",
                   "This dataset will be excluded from harmonization."))
} else {
  meta <- read.csv(szbd$local_file)

  if (verbose) {
    colnames(meta)

    print_summary(meta, sex_col = "reportedGender",
                  isHispanic_col = "ethnicity", braak_nft_col = "Braak")
  }

  meta_new <- harmonize(szbd, meta, spec)

  if (verbose) {
    print_summary(meta_new)
  }

  cat("\n", szbd$gen_name, "/", szbd$name, "\n")
  validate_values(meta_new, spec)

  datasets[[szbd$name]] <- list(
    data = meta_new,
    study = szbd,
    filename = basename(szbd$local_file)
  )
}


# GEN-A9 / SMIB-AD -------------------------------------------------------------
# No overlap with other data sets.

smib_ad <- studies$smib_ad
meta_file <- synapse_download(smib_ad$syn_id)
meta <- read.csv(meta_file$path)

if (verbose) {
  colnames(meta)

  print_summary(meta, braak_nft_col = "Braak", bscore_nft_col = "bScore")
}

meta_new <- harmonize(smib_ad, meta, spec)

if (verbose) {
  print_summary(meta_new)
}

cat("\n", smib_ad$gen_name, "/", smib_ad$name, "\n")
validate_values(meta_new, spec)

datasets[[smib_ad$name]] <- list(
  data = meta_new,
  study = smib_ad,
  filename = meta_file$name
)


# GEN-A10 / MCMPS --------------------------------------------------------------

mcmps <- studies$mcmps
meta_file <- synapse_download(mcmps$syn_id)
meta <- read.csv(meta_file$path)

if (verbose) {
  colnames(meta)

  print_summary(meta, braak_nft_col = "Braak", bscore_nft_col = "bScore")
}

meta_new <- harmonize(mcmps, meta, spec)

if (verbose) {
  print_summary(meta_new)
}

cat("\n", mcmps$gen_name, "/", mcmps$name, "\n")
validate_values(meta_new, spec)

datasets[[mcmps$name]] <- list(
  data = meta_new,
  study = mcmps,
  filename = meta_file$name
)


# GEN-A11 / MC_snRNA -----------------------------------------------------------

mc_snrna <- studies$mc_snrna
meta_file <- synapse_download(mc_snrna$syn_id)
meta <- read.csv(meta_file$path)

if (verbose) {
  colnames(meta)

  print_summary(meta, braak_nft_col = "Braak", bscore_nft_col = "bScore")
}

meta_new <- harmonize(mc_snrna, meta, spec)

if (verbose) {
  print_summary(meta_new)
}

cat("\n", mc_snrna$gen_name, "/", mc_snrna$name, "\n")
validate_values(meta_new, spec)

datasets[[mc_snrna$name]] <- list(
  data = meta_new,
  study = mc_snrna,
  filename = meta_file$name
)


# GEN-A12 / MC-BrAD ------------------------------------------------------------

mc_brad <- studies$mc_brad
meta_file <- synapse_download(mc_brad$syn_id)
meta <- read.csv(meta_file$path)

if (verbose) {
  colnames(meta)

  print_summary(meta, braak_nft_col = "Braak", bscore_nft_col = "bScore")
}

meta_new <- harmonize(mc_brad, meta, spec)

if (verbose) {
  print_summary(meta_new)
}

cat("\n", mc_brad$gen_name, "/", mc_brad$name, "\n")
validate_values(meta_new, spec)

datasets[[mc_brad$name]] <- list(
  data = meta_new,
  study = mc_brad,
  filename = meta_file$name
)


# GEN-A15 / ASAP ---------------------------------------------------------------

asap <- studies$asap

# ASAP metadata is split across two files
if (!file.exists(asap$local_file$subject)) {
  warning(str_glue("{asap$name} file {asap$local_file$subject} doesn't exist! ",
                   "This dataset will be excluded from harmonization."))
} else if (!file.exists(asap$local_file$clinical)) {
  warning(str_glue("{asap$name} file {asap$local_file$clinical} doesn't exist! ",
                   "This dataset will be excluded from harmonization."))
} else {
  subj <- read.csv(asap$local_file$subject)
  clin <- read.csv(asap$local_file$clinical)

  meta <- merge(subj, clin)

  if (verbose) {
    colnames(meta)

    print_summary(
      meta, pmi_col = "duration_pmi", ageDeath_col = "age_at_death",
      braak_nft_col = "path_braak_nft", braak_lb_col = "path_braak_asyn",
      cerad_col = "path_cerad", thal_col = "path_thal",
      isHispanic_col = "ethnicity", apoe_col = "apoe_e4_status"
    )
  }

  meta_new <- harmonize(asap, meta, spec)

  if (verbose) {
    print_summary(meta_new)
  }

  cat("\n", asap$gen_name, "/", asap$name, "\n")
  validate_values(meta_new, spec)

  datasets[[asap$name]] <- list(
    data = meta_new,
    study = asap,
    filename = "ASAP_PMDBS_metadata.csv"
  )
}


# GEN-A16 / McCarroll_SCZ ------------------------------------------------------

mccarroll_scz <- studies$mccarroll_scz

if (!file.exists(mccarroll_scz$local_file)) {
  warning(str_glue("{mccarroll_scz$name} file {mccarroll_scz$local_file} doesn't exist! ",
                   "This dataset will be excluded from harmonization."))
} else {
  meta <- read.delim(mccarroll_scz$local_file)

  if (verbose) {
    colnames(meta)

    print_summary(meta, ageDeath_col = "Age", sex_col = "Sex")
  }

  meta_new <- harmonize(mccarroll_scz, meta, spec)

  if (verbose) {
    print_summary(meta_new)
  }

  cat("\n", mccarroll_scz$gen_name, "/", mccarroll_scz$name, "\n")
  validate_values(meta_new, spec)

  datasets[[mccarroll_scz$name]] <- list(
    data = meta_new,
    study = mccarroll_scz,
    filename = basename(mccarroll_scz$local_file)
  )
}


# GEN-A17 / McCarroll_HD -------------------------------------------------------

mccarroll_hd <- studies$mccarroll_hd

if (!file.exists(mccarroll_hd$local_file)) {
  warning(str_glue("{mccarroll_hd$name} file {mccarroll_hd$local_file} doesn't exist! ",
                   "This dataset will be excluded from harmonization."))
} else {
  meta <- read.delim(mccarroll_hd$local_file)

  if (verbose) {
    colnames(meta)

    print_summary(meta, ageDeath_col = "Age", sex_col = "Sex")
  }

  meta_new <- harmonize(mccarroll_hd, meta, spec)

  if (verbose) {
    print_summary(meta_new)
  }

  cat("\n", mccarroll_hd$gen_name, "/", mccarroll_hd$name, "\n")
  validate_values(meta_new, spec)

  datasets[[mccarroll_hd$name]] <- list(
    data = meta_new,
    study = mccarroll_hd,
    filename = basename(mccarroll_hd$local_file)
  )
}


# GEN-A18 / TargetALS ----------------------------------------------------------

target_als <- studies$target_als

if (!file.exists(target_als$local_file)) {
  warning(str_glue("{target_als$name} file {target_als$local_file} doesn't exist! ",
                   "This dataset will be excluded from harmonization."))
} else {
  meta <- read_xlsx(target_als$local_file)

  if (verbose) {
    colnames(meta)

    print_summary(meta, pmi_col = "PMI (hrs)", sex_col = "Sex",
                  isHispanic_col = "Ethnicity", race_col = "Race",
                  ageDeath_col = "Age.At.Death")
  }

  meta_new <- harmonize(target_als, meta, spec)

  if (verbose) {
    print_summary(meta_new)
  }

  cat("\n", target_als$gen_name, "/", target_als$name, "\n")
  validate_values(meta_new, spec)

  datasets[[target_als$name]] <- list(
    data = meta_new,
    study = target_als,
    filename = "TargetALS_final_metadata.csv"
  )
}


# GEN-B4 / AMP-AD_DiverseCohorts -----------------------------------------------

diverse_cohorts <- studies$diverse_cohorts
meta_file <- synapse_download(diverse_cohorts$syn_id)
meta <- read.csv(meta_file$path)

if (verbose) {
  colnames(meta)
  print_summary(meta, braak_nft_col = "Braak", bscore_nft_col = "bScore")
}

meta_new <- harmonize(diverse_cohorts, meta, spec)

if (verbose) {
  print_summary(meta_new)
}

cat("\n", diverse_cohorts$gen_name, "/", diverse_cohorts$name, "\n")
validate_values(meta_new, spec)

datasets[[diverse_cohorts$name]] <- list(
  data = meta_new,
  study = diverse_cohorts,
  filename = meta_file$name
)


# GEN-B8 / BD2 -----------------------------------------------------------------

bd2 <- studies$bd2

# BD2 data, uses locally-downloaded file that is not on Synapse
if (!file.exists(bd2$local_file)) {
  warning(str_glue("{bd2$name} file {bd2$local_file} doesn't exist! This dataset ",
                   "will be excluded from harmonization."))
} else {
  meta <- read.csv(bd2$local_file)

  if (verbose) {
    colnames(meta)
    print_summary(meta, ageDeath_col = "Age", race_col = "Race", sex_col = "Sex")
  }

  meta_new <- harmonize(bd2, meta, spec)

  if (verbose) {
    print_summary(meta_new)
  }

  cat("\n", bd2$gen_name, "/", bd2$name, "\n")
  validate_values(meta_new, spec)

  datasets[[bd2$name]] <- list(
    data = meta_new,
    study = bd2,
    filename = basename(bd2$local_file)
  )
}


# De-duplicate studies ---------------------------------------------------------

dedup_df <- deduplicate_studies(lapply(datasets, "[[", "data"),
                                spec, verbose = FALSE)

# Dementia should be 1 if FTD is 1
dedup_df$Dementia[dedup_df$FTD == 1] <- 1


# Write files and upload -------------------------------------------------------

for (data_item in datasets) {
  # subset deduplicated data to this study and to only the columns that existed
  # in the non-duplicated data
  new_data <- subset(dedup_df, study == data_item$study$name) |>
    select(all_of(colnames(data_item$data)))

  new_filename <- write_metadata(new_data, data_item$filename)
  new_syn_id <- synapse_upload(new_filename, spec$upload_synID)

  # Create manifest file
  manifest <- rbind(manifest, to_manifest_df(data_item$study, new_syn_id))

  if ("assoc_studies" %in% names(data_item)) {
    for (study in data_item$assoc_studies) {
      manifest <- rbind(manifest, to_manifest_df(study, new_syn_id))
    }
  }
}


# Upload manifest file ---------------------------------------------------------

manifest <- manifest |> arrange(GENESIS_study)
manifest_file <- file.path("data", "output", "metadata_manifest.csv")
write.csv(manifest, manifest_file,
  row.names = FALSE, quote = FALSE
)
synapse_upload(manifest_file, spec$upload_synID)


# Combine all harmonized data into one data frame

required_columns <- c(spec$demographic_columns, spec$diagnosis_columns)

combined_data <- dedup_df |>
  select(all_of(required_columns)) |>
  mutate(across(c(individualID, ageDeath, apoeGenotype), as.character)) |>
  merge(select(manifest, study, GENESIS_study))

cat("\nMerged metadata\n")
validate_values(combined_data, spec)

# df_all <- df_all |>
#  group_by(individualID) |>
#  mutate(genesis_study = paste(sort(genesis_study), collapse = "; "))

new_file <- write_metadata(combined_data, "GENESIS_metadata_combined.csv")
synapse_upload(new_file, spec$upload_synID)
