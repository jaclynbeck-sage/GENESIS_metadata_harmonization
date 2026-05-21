# This script harmonizes and de-duplicates multiple studies used for GENESIS,
# filling in missing information where there is sample overlap between a data
# set and the Diverse Cohorts / AMP-AD 1.0 / NPS-AD data. Harmonized metadata is
# uploaded to Synapse for use.
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
#   GEN-B4 / AMP-AD_DiverseCohorts
#   GEN-B5 / SEA-AD_multiome
#   GEN-B6 / MIT_ROSMAP_Multiomics
#   GEN-B8 / BD2
#   GEN-B11 / ROSMAP_CUIMC2

library(synapser)
library(dplyr)
library(purrr)
library(stringr)

spec <- config::get(file = "GENESIS_harmonization.yml")
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

# This dataset has overlap with BD2, ROSMAP, and Diverse Cohorts, and is used to
# fill in missing values in these studies.

meta_file <- synapse_download(studies$nps_ad$syn_id)
meta <- read.csv(meta_file$path)

if (verbose) {
  colnames(meta)
  print_summary(meta, braak_nft_col = "Braak", bscore_nft_col = "bScore")
}

meta_new <- harmonize(studies$nps_ad$name, meta, spec)

if (verbose) {
  print_summary(meta_new)
}

cat("\n", studies$nps_ad$gen_name, "/", studies$nps_ad$name, "\n")
validate_values(meta_new, spec)

datasets[[studies$nps_ad$name]] <- meta_new

new_filename <- write_metadata(meta_new, meta_file$name)
new_syn_id <- synapse_upload(new_filename, spec$upload_synID)

manifest <- rbind(
  manifest,
  data.frame(
    GENESIS_study = studies$nps_ad$gen_name,
    study = studies$nps_ad$name,
    metadata_synid = paste0(new_syn_id$id, ".", new_syn_id$versionNumber)
  )
)


# GEN-A2, GEN-A8, GEN-A13, GEN-B6, GEN-B11 / ROSMAP ----------------------------
# These GENESIS studies all use original ROSMAP metadata

meta_file <- synapse_download(studies$rosmap$syn_id)
meta <- read.csv(meta_file$path)

if (verbose) {
  colnames(meta)
  print_summary(meta, braak_nft_col = "Braak", bscore_nft_col = "bScore")
}

meta_new <- harmonize(studies$rosmap$name, meta, spec)

if (verbose) {
  print_summary(meta_new)
}

cat("\n", studies$rosmap$gen_name, "/", studies$rosmap$name, "\n")
validate_values(meta_new, spec)

datasets[[studies$rosmap$name]] <- meta_new

new_filename <- write_metadata(meta_new, meta_file$name)
new_syn_id <- synapse_upload(new_filename, spec$upload_synID)

ros_studies <- lapply(
  list(studies$rosmap, studies$snRNA_trem2, studies$snRNA_BA10,
       studies$mit_rosmap, studies$cuimc2),
  as.data.frame
) |>
  list_rbind() |>
  select(-syn_id) |>
  dplyr::rename(study = name, GENESIS_study = gen_name) |>
  mutate(metadata_synid = paste0(new_syn_id$id, ".", new_syn_id$versionNumber))

manifest <- rbind(manifest, ros_studies)


# GEN-A3 / AMP-PD --------------------------------------------------------------

# AMP-PD data, uses locally-downloaded files that are not on Synapse
if (!file.exists(studies$amp_pd$local_files$main)) {
  warning(str_glue("AMP-PD file {studies$amp_pd$local_files$main} doesn't ",
                   "exist! This dataset will be excluded from harmonization."))
} else if (!file.exists(studies$amp_pd$local_files$demographics)) {
  warning(str_glue("AMP-PD file {studies$amp_pd$local_files$demographics} doesn't ",
                   "exist! This dataset will be excluded from harmonization."))
} else {
  meta <- read.csv(studies$amp_pd$local_files$main)
  demographics <- read.csv(studies$amp_pd$local_files$demographics) |>
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

  meta_new <- harmonize(studies$amp_pd$name, meta, spec)

  if (verbose) {
    print_summary(meta_new)
  }

  cat("\n", studies$amp_pd$gen_name, "/", studies$amp_pd$name, "\n")
  validate_values(meta_new, spec)

  datasets[[studies$amp_pd$name]] <- meta_new

  new_filename <- write_metadata(meta_new, basename(studies$amp_pd$local_files$main))
  new_syn_id <- synapse_upload(new_filename, spec$upload_synID)

  manifest <- rbind(
    manifest,
    data.frame(
      GENESIS_study = studies$amp_pd$gen_name,
      study = studies$amp_pd$name,
      metadata_synid = paste0(new_syn_id$id, ".", new_syn_id$versionNumber)
    )
  )
}


# GEN-A4, GEN-B5 / SEA-AD ------------------------------------------------------

meta_file <- synapse_download(studies$sea_ad$syn_id)
meta <- read.csv(meta_file$path)

if (verbose) {
  colnames(meta)

  print_summary(meta, braak_nft_col = "Braak", bscore_nft_col = "bScore")
}

meta_new <- harmonize(studies$sea_ad$name, meta, spec)

if (verbose) {
  print_summary(meta_new)
}

cat("\n", studies$sea_ad$gen_name, "/", studies$sea_ad$name, "\n")
validate_values(meta_new, spec)

datasets[[studies$sea_ad$name]] <- meta_new

new_filename <- write_metadata(meta_new, meta_file$name)
new_syn_id <- synapse_upload(new_filename, spec$upload_synID)

manifest <- rbind(
  manifest,
  data.frame(
    GENESIS_study = c(studies$sea_ad$gen_name, studies$sea_ad_multi$gen_name),
    study = c(studies$sea_ad$name, studies$sea_ad_multi$name),
    metadata_synid = paste0(new_syn_id$id, ".", new_syn_id$versionNumber)
  )
)


# GEN-A5 / CMC (PsychENCODE) ---------------------------------------------------

# NOTE: There are overlaps between CMC and NPS-AD individuals, but the way they
# are named between the two studies is different. Per Jaro, we will wait until
# he compiles a list of overlapping samples to de-duplicate.

if (!file.exists(studies$cmc$local_files)) {
  warning(str_glue("PsychENCODE file {studies$cmc$local_files} doesn't exist! ",
                   "This dataset will be excluded from harmonization."))
} else {
  meta <- read.csv(studies$cmc$local_files)

  if (verbose) {
    colnames(meta)

    print_summary(meta, sex_col = "reportedGender",
                  isHispanic_col = "ethnicity", braak_nft_col = "Braak")
  }

  meta_new <- harmonize(studies$cmc$name, meta, spec)

  if (verbose) {
    print_summary(meta_new)
  }

  cat("\n", studies$cmc$gen_name, "/", studies$cmc$name, "\n")
  validate_values(meta_new, spec)

  datasets[[studies$cmc$name]] <- meta_new

  new_filename <- write_metadata(meta_new, basename(studies$cmc$local_files))
  new_syn_id <- synapse_upload(new_filename, spec$upload_synID)

  manifest <- rbind(
    manifest,
    data.frame(
      GENESIS_study = studies$cmc$gen_name,
      study = studies$cmc$name,
      metadata_synid = paste0(new_syn_id$id, ".", new_syn_id$versionNumber)
    )
  )
}


# GEN-A6 / SZBDMulti-Seq (PsychENCODE) -----------------------------------------

# NOTE: There are overlaps between CMC and NPS-AD individuals, but the way they
# are named between the two studies is different. Per Jaro, we will wait until
# he compiles a list of overlapping samples to de-duplicate.

if (!file.exists(studies$szbd$local_files)) {
  warning(str_glue("PsychENCODE file {studies$szbd$local_files} doesn't exist! ",
                   "This dataset will be excluded from harmonization."))
} else {
  meta <- read.csv(studies$szbd$local_files)

  if (verbose) {
    colnames(meta)

    print_summary(meta, sex_col = "reportedGender",
                  isHispanic_col = "ethnicity", braak_nft_col = "Braak")
  }

  meta_new <- harmonize(studies$szbd$name, meta, spec)

  if (verbose) {
    print_summary(meta_new)
  }

  cat("\n", studies$szbd$gen_name, "/", studies$szbd$name, "\n")
  validate_values(meta_new, spec)

  datasets[[studies$szbd$name]] <- meta_new

  new_filename <- write_metadata(meta_new, basename(studies$szbd$local_files))
  new_syn_id <- synapse_upload(new_filename, spec$upload_synID)

  manifest <- rbind(
    manifest,
    data.frame(
      GENESIS_study = studies$szbd$gen_name,
      study = studies$szbd$name,
      metadata_synid = paste0(new_syn_id$id, ".", new_syn_id$versionNumber)
    )
  )
}


# GEN-A9 / SMIB-AD -------------------------------------------------------------
# No overlap with other data sets.

meta_file <- synapse_download(studies$smib_ad$syn_id)
meta <- read.csv(meta_file$path)

if (verbose) {
  colnames(meta)

  print_summary(meta, braak_nft_col = "Braak", bscore_nft_col = "bScore")
}

meta_new <- harmonize(studies$smib_ad$name, meta, spec)

if (verbose) {
  print_summary(meta_new)
}

cat("\n", studies$smib_ad$gen_name, "/", studies$smib_ad$name, "\n")
validate_values(meta_new, spec)

datasets[[studies$smib_ad$name]] <- meta_new

new_filename <- write_metadata(meta_new, meta_file$name)
new_syn_id <- synapse_upload(new_filename, spec$upload_synID)

manifest <- rbind(
  manifest,
  data.frame(
    GENESIS_study = studies$smib_ad$gen_name,
    study = studies$smib_ad$name,
    metadata_synid = paste0(new_syn_id$id, ".", new_syn_id$versionNumber)
  )
)


# GEN-A10 / MCMPS --------------------------------------------------------------

meta_file <- synapse_download(studies$mcmps$syn_id)
meta <- read.csv(meta_file$path)

if (verbose) {
  colnames(meta)

  print_summary(meta, braak_nft_col = "Braak", bscore_nft_col = "bScore")
}

meta_new <- harmonize(studies$mcmps$name, meta, spec)

if (verbose) {
  print_summary(meta_new)
}

cat("\n", studies$mcmps$gen_name, "/", studies$mcmps$name, "\n")
validate_values(meta_new, spec)

datasets[[studies$mcmps$name]] <- meta_new

new_filename <- write_metadata(meta_new, meta_file$name)
new_syn_id <- synapse_upload(new_filename, spec$upload_synID)

manifest <- rbind(
  manifest,
  data.frame(
    GENESIS_study = studies$mcmps$gen_name,
    study = studies$mcmps$name,
    metadata_synid = paste0(new_syn_id$id, ".", new_syn_id$versionNumber)
  )
)


# GEN-A11 / MC_snRNA -----------------------------------------------------------

meta_file <- synapse_download(studies$mc_snrna$syn_id)
meta <- read.csv(meta_file$path)

if (verbose) {
  colnames(meta)

  print_summary(meta, braak_nft_col = "Braak", bscore_nft_col = "bScore")
}

meta_new <- harmonize(studies$mc_snrna$name, meta, spec)

if (verbose) {
  print_summary(meta_new)
}

cat("\n", studies$mc_snrna$gen_name, "/", studies$mc_snrna$name, "\n")
validate_values(meta_new, spec)

datasets[[studies$mc_snrna$name]] <- meta_new

new_filename <- write_metadata(meta_new, meta_file$name)
new_syn_id <- synapse_upload(new_filename, spec$upload_synID)

manifest <- rbind(
  manifest,
  data.frame(
    GENESIS_study = studies$mc_snrna$gen_name,
    study = studies$mc_snrna$name,
    metadata_synid = paste0(new_syn_id$id, ".", new_syn_id$versionNumber)
  )
)


# GEN-A12 / MC-BrAD ------------------------------------------------------------

meta_file <- synapse_download(studies$mc_brad$syn_id)
meta <- read.csv(meta_file$path)

if (verbose) {
  colnames(meta)

  print_summary(meta, braak_nft_col = "Braak", bscore_nft_col = "bScore")
}

meta_new <- harmonize(studies$mc_brad$name, meta, spec)

if (verbose) {
  print_summary(meta_new)
}

cat("\n", studies$mc_brad$gen_name, "/", studies$mc_brad$name, "\n")
validate_values(meta_new, spec)

datasets[[studies$mc_brad$name]] <- meta_new

new_filename <- write_metadata(meta_new, meta_file$name)
new_syn_id <- synapse_upload(new_filename, spec$upload_synID)

manifest <- rbind(
  manifest,
  data.frame(
    GENESIS_study = studies$mc_brad$gen_name,
    study = studies$mc_brad$name,
    metadata_synid = paste0(new_syn_id$id, ".", new_syn_id$versionNumber)
  )
)


# GEN-A15 / ASAP ---------------------------------------------------------------

# ASAP metadata is split across two files
if (!file.exists(studies$asap$local_files$subject)) {
  warning(str_glue("ASAP file {studies$asap$local_files$subject} doesn't exist! ",
                   "This dataset will be excluded from harmonization."))
} else if (!file.exists(studies$asap$local_files$clinical)) {
  warning(str_glue("ASAP file {studies$asap$local_files$clinical} doesn't exist! ",
                   "This dataset will be excluded from harmonization."))
} else {
  subj <- read.csv(studies$asap$local_files$subject)
  clin <- read.csv(studies$asap$local_files$clinical)

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

  meta_new <- harmonize(studies$asap$name, meta, spec)

  if (verbose) {
    print_summary(meta_new)
  }

  cat("\n", studies$asap$gen_name, "/", studies$asap$name, "\n")
  validate_values(meta_new, spec)

  datasets[[studies$asap$name]] <- meta_new

  new_filename <- write_metadata(meta_new, "ASAP_PMDBS_metadata.csv")
  new_syn_id <- synapse_upload(new_filename, spec$upload_synID)

  manifest <- rbind(
    manifest,
    data.frame(
      GENESIS_study = studies$asap$gen_name,
      study = studies$asap$name,
      metadata_synid = paste0(new_syn_id$id, ".", new_syn_id$versionNumber)
    )
  )
}


# GEN-A16 / McCarroll_SCZ ------------------------------------------------------

if (!file.exists(studies$mccarroll_scz$local_files)) {
  warning(str_glue("McCarroll SCZ file {studies$mccarroll_scz$local_files} doesn't exist! ",
                   "This dataset will be excluded from harmonization."))
} else {
  meta <- read.delim(studies$mccarroll_scz$local_files)

  if (verbose) {
    colnames(meta)

    print_summary(meta, ageDeath_col = "Age", sex_col = "Sex")
  }

  meta_new <- harmonize(studies$mccarroll_scz$name, meta, spec)

  if (verbose) {
    print_summary(meta_new)
  }

  cat("\n", studies$mccarroll_scz$gen_name, "/", studies$mccarroll_scz$name, "\n")
  validate_values(meta_new, spec)

  datasets[[studies$mccarroll_scz$name]] <- meta_new

  new_filename <- write_metadata(meta_new, basename(studies$mccarroll_scz$local_files))
  new_syn_id <- synapse_upload(new_filename, spec$upload_synID)

  manifest <- rbind(
    manifest,
    data.frame(
      GENESIS_study = studies$mccarroll_scz$gen_name,
      study = studies$mccarroll_scz$name,
      metadata_synid = paste0(new_syn_id$id, ".", new_syn_id$versionNumber)
    )
  )
}


# GEN-A17 / McCarroll_HD -------------------------------------------------------

if (!file.exists(studies$mccarroll_hd$local_files)) {
  warning(str_glue("McCarroll HD file {studies$mccarroll_hd$local_files} doesn't exist! ",
                   "This dataset will be excluded from harmonization."))
} else {
  meta <- read.delim(studies$mccarroll_hd$local_files)

  if (verbose) {
    colnames(meta)

    print_summary(meta, ageDeath_col = "Age", sex_col = "Sex")
  }

  meta_new <- harmonize(studies$mccarroll_hd$name, meta, spec)

  if (verbose) {
    print_summary(meta_new)
  }

  cat("\n", studies$mccarroll_hd$gen_name, "/", studies$mccarroll_hd$name, "\n")
  validate_values(meta_new, spec)

  datasets[[studies$mccarroll_hd$name]] <- meta_new

  new_filename <- write_metadata(meta_new, basename(studies$mccarroll_hd$local_files))
  new_syn_id <- synapse_upload(new_filename, spec$upload_synID)

  manifest <- rbind(
    manifest,
    data.frame(
      GENESIS_study = studies$mccarroll_hd$gen_name,
      study = studies$mccarroll_hd$name,
      metadata_synid = paste0(new_syn_id$id, ".", new_syn_id$versionNumber)
    )
  )
}


# GEN-B4 / AMP-AD_DiverseCohorts -----------------------------------------------
# Uses Diverse Cohorts metadata from step 1.

meta_file <- synapse_download(studies$diverse_cohorts$syn_id)
meta <- read.csv(meta_file$path)

if (verbose) {
  colnames(meta)
  print_summary(meta, braak_nft_col = "Braak", bscore_nft_col = "bScore")
}

meta_new <- harmonize(studies$diverse_cohorts$name, meta, spec)

if (verbose) {
  print_summary(meta_new)
}

cat("\n", studies$diverse_cohorts$gen_name, "/", studies$diverse_cohorts$name, "\n")
validate_values(meta_new, spec)

datasets[[studies$diverse_cohorts$name]] <- meta_new

new_filename <- write_metadata(meta_new, meta_file$name)
new_syn_id <- synapse_upload(new_filename, spec$upload_synID)

manifest <- rbind(
  manifest,
  data.frame(
    GENESIS_study = studies$diverse_cohorts$gen_name,
    study = studies$diverse_cohorts$name,
    metadata_synid = paste0(new_syn_id$id, ".", new_syn_id$versionNumber)
  )
)


# GEN-B8 / BD2 -----------------------------------------------------------------

# BD2 data, uses locally-downloaded file that is not on Synapse
if (!file.exists(studies$bd2$local_files)) {
  warning(str_glue("BD2 file {studies$bd2$local_files} doesn't exist! This dataset ",
                   "will be excluded from harmonization."))
} else {
  meta <- read.csv(studies$bd2$local_files)

  if (verbose) {
    colnames(meta)
    print_summary(meta, ageDeath_col = "Age", race_col = "Race", sex_col = "Sex")
  }

  meta_new <- harmonize(studies$bd2$name, meta, spec)

  if (verbose) {
    print_summary(meta_new)
  }

  cat("\n", studies$bd2$gen_name, "/", studies$bd2$name, "\n")
  validate_values(meta_new, spec)

  datasets[[studies$bd2$name]] <- meta_new

  new_filename <- write_metadata(meta_new, basename(studies$bd2$local_files))
  new_syn_id <- synapse_upload(new_filename, spec$upload_synID)

  manifest <- rbind(
    manifest,
    data.frame(
      GENESIS_study = studies$bd2$gen_name,
      study = studies$bd2$name,
      metadata_synid = paste0(new_syn_id$id, ".", new_syn_id$versionNumber)
    )
  )
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

df_list <- apply(manifest, 1, function(m_row) {
  m_file <- synapse_download(m_row[["metadata_synid"]])
  meta <- read.csv(m_file$path) |>
    mutate(
      across(c(individualID, ageDeath, apoeGenotype), as.character),
      GENESIS_study = m_row[["GENESIS_study"]],
      study = m_row[["study"]]
    ) |>
    select(all_of(required_columns), GENESIS_study)
  return(meta)
}, simplify = FALSE)

df_all <- purrr::list_rbind(df_list)

cat("\nMerged metadata\n")
validate_values(df_all, spec)

# df_all <- df_all |>
#  group_by(individualID) |>
#  mutate(genesis_study = paste(sort(genesis_study), collapse = "; "))

new_file <- write_metadata(df_all, "GENESIS_metadata_combined.csv")
synapse_upload(new_file, spec$upload_synID)


# WIP
dedup <- deduplicate_studies(datasets, spec, verbose = FALSE)

# The binary `Control` diagnosis column is based on all of the other diagnosis
# columns:
#   1 if all the other columns are 0 or NA,
#   0 if there is a 1 in any other column
dedup$Control <- (rowSums(dedup[, spec$diagnosis_columns], na.rm = TRUE) == 0) |>
  as.numeric()
