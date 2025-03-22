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
  "GEN-A2" = "syn3191087.11", # ROSMAP
  "GEN-A4" = "syn31149116.7", # SEA-AD
  "GEN-A8" = "syn3191087.11", # snRNAseqAD_TREM2, uses ROSMAP
  "GEN-A9" = "syn22432749.1", # SMIB-AD
  "GEN-A10" = "syn25891193.1", # MCMPS
  "GEN-A11" = "syn31563038.1", # MC_snRNA
  "GEN-A12" = "syn51401700.2", # MC-BrAD
  "GEN-A13" = "syn3191087.11", # snRNAseqPFC_BA10, uses ROSMAP
  "GEN-A14" = "syn24610550.2", # HBI_scRNAseq
  # "GEN-B1" = "TBD",
  # "GEN-B2" = "TBD",
  # "GEN-B3" = "TBD",
  "GEN-B4" = "syn51757646.20", # AMP-AD_DiverseCohorts
  "GEN-B5" = "syn31149116.7", # SEA-AD
  "GEN-B6" = "syn3191087.11" # MIT_ROSMAP_Multiomics metadata is syn52430346 but the study uses the ROSMAP file
)

UPLOAD_SYNID <- "syn64759869"

manifest <- c()

synLogin()


# GEN-A1 -----------------------------------------------------------------------
# Has sample overlap with MayoRNASeq, MSBB, ROSMAP, and Diverse Cohorts

meta_file <- synapse_download(syn_ids[["GEN-A1"]])
meta <- read.csv(meta_file$path)

colnames(meta)

print_qc(meta, isHispanic_col = "ethnicity", pmi_col = "PMI", cerad_col = "CERAD")

# TODO should pull metadata from Diverse Cohorts, MSBB, and ROSMAP for the non-NPS-AD samples,
# since some info like Braak is missing from this metadata
meta_new <- meta |>
  select(-Component) |>
  rename(
    isHispanic = ethnicity,
    pmi = PMI,
    amyCerad = CERAD
  ) |>
  mutate(
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
  study = "GEN-A1",
  metadata_synid = new_syn_id
))


# GEN-A2, GEN-A8, GEN-A13, GEN-B6 ----------------------------------------------
# All use ROSMAP metadata

# TODO get the harmonized file from Synapse
ros_file <- file.path(
  "data", "output",
  "ROSMAP_clinical_harmonized.csv"
)

new_syn_id <- synapse_upload(ros_file, UPLOAD_SYNID)

manifest <- rbind(
  manifest,
  data.frame(
    study = c("GEN-A2", "GEN-A8", "GEN-A13", "GEN-B6"),
    metadata_synid = new_syn_id
  )
)


# GEN-A4, GEN-B5 ---------------------------------------------------------------
# Uses SEA-AD metadata

# TODO get the harmonized file from Synapse
sea_ad_file <- file.path(
  "data", "output",
  "SEA-AD_individual_metadata_harmonized.csv"
)

new_syn_id <- synapse_upload(sea_ad_file, UPLOAD_SYNID)

manifest <- rbind(
  manifest,
  data.frame(
    study = c("GEN-A4", "GEN-B5"),
    metadata_synid = new_syn_id
  )
)


# GEN-A9 -----------------------------------------------------------------------
# TODO unknown how they encode CERAD

meta_file <- synapse_download(syn_ids[["GEN-A9"]])
meta <- read.csv(meta_file$path)

colnames(meta)

print_qc(meta, isHispanic_col = "ethnicity", cerad_col = "CERAD")

meta_new <- meta |>
  rename(
    isHispanic = ethnicity,
    amyCerad = CERAD,
    cohort = individualIdSource
  ) |>
  mutate(
    ageDeath = censor_ages(ageDeath, spec),
    race = spec$race$White,
    isHispanic = spec$isHispanic$hisp_false,
    apoeGenotype = as.character(apoeGenotype),
    apoe4Status = get_apoe4Status(apoeGenotype, spec),
    Braak = to_Braak_stage(Braak, spec),
    bScore = get_bScore(Braak, spec),
    amyCerad = spec$missing, # TODO unknown how they encode CERAD
    amyAny = get_amyAny(amyCerad, spec),
    amyThal = spec$missing,
    amyA = spec$missing,
    cohort = case_when(is.na(cohort) ~ spec$cohort$smri,
      .default = spec$cohort$banner
    ),
    dataContributionGroup = case_match(cohort,
      spec$cohort$smri ~ spec$dataContributionGroup$stanley,
      .default = spec$dataContributionGroup$banner
    )
  )


print_qc(meta_new)
validate_values(meta_new, spec)

new_filename <- write_metadata(meta_new, meta_file$name)
new_syn_id <- synapse_upload(new_filename, UPLOAD_SYNID)

manifest <- rbind(
  manifest,
  data.frame(study = "GEN-A9", metadata_synid = new_syn_id)
)


# GEN-A10 ----------------------------------------------------------------------

meta_file <- synapse_download(syn_ids[["GEN-A10"]])
meta <- read.csv(meta_file$path)

colnames(meta)

print_qc(meta, isHispanic_col = "ethnicity", cerad_col = "CERAD")

meta_new <- meta |>
  rename(
    isHispanic = ethnicity,
    amyCerad = CERAD,
    cohort = individualIdSource
  ) |>
  mutate(
    race = str_trim(race),
    isHispanic = str_trim(isHispanic),
    isHispanic = case_match(isHispanic,
      "Hispanic or Latino" ~ spec$isHispanic$hisp_true,
      "Not Hispanic or Latino" ~ spec$isHispanic$hisp_false,
      "Middle Eastern" ~ spec$isHispanic$hisp_false,
      .default = isHispanic
    ),
    apoeGenotype = case_match(apoeGenotype,
      NA ~ spec$missing,
      .default = as.character(apoeGenotype)
    ),
    apoe4Status = get_apoe4Status(apoeGenotype, spec),
    amyCerad = spec$missing,
    amyAny = spec$missing,
    amyThal = spec$missing,
    amyA = spec$missing,
    Braak = spec$missing,
    bScore = spec$missing,
    dataContributionGroup = spec$dataContributionGroup$mayo
  )

print_qc(meta_new)
validate_values(meta_new, spec)

new_filename <- write_metadata(meta_new, meta_file$name)
new_syn_id <- synapse_upload(new_filename, UPLOAD_SYNID)

manifest <- rbind(
  manifest,
  data.frame(study = "GEN-A10", metadata_synid = new_syn_id)
)


# GEN-A11 ----------------------------------------------------------------------
# Has some sample overlap with original Mayo metadata and Diverse Cohorts metadata

meta_file <- synapse_download(syn_ids[["GEN-A11"]])
meta <- read.csv(meta_file$path)

colnames(meta)

print_qc(meta, isHispanic_col = "ethnicity", cerad_col = "CERAD")

meta_new <- meta |>
  rename(
    isHispanic = ethnicity,
    amyCerad = CERAD
  ) |>
  mutate(
    ageDeath = censor_ages(ageDeath, spec),
    isHispanic = spec$isHispanic$hisp_false,
    apoeGenotype = as.character(apoeGenotype),
    apoe4Status = get_apoe4Status(apoeGenotype, spec),
    amyCerad = spec$missing,
    amyAny = spec$missing,
    amyThal = spec$missing,
    amyA = spec$missing,
    Braak = to_Braak_stage(floor(Braak), spec),
    bScore = get_bScore(Braak, spec),
    dataContributionGroup = spec$dataContributionGroup$mayo,
    cohort = spec$cohort$mayo
  )

print_qc(meta_new)
validate_values(meta_new, spec)

new_filename <- write_metadata(meta_new, meta_file$name)
new_syn_id <- synapse_upload(new_filename, UPLOAD_SYNID)

manifest <- rbind(
  manifest,
  data.frame(study = "GEN-A11", metadata_synid = new_syn_id)
)


# GEN-A12 ----------------------------------------------------------------------
# Has some sample overlap with original Mayo metadata and Diverse Cohorts metadata

meta_file <- synapse_download(syn_ids[["GEN-A12"]])
meta <- read.csv(meta_file$path)

colnames(meta)

print_qc(meta, isHispanic_col = "ethnicity", cerad_col = "CERAD", thal_col = "Thal")

meta_new <- meta |>
  rename(
    isHispanic = ethnicity,
    amyCerad = CERAD,
    amyThal = Thal
  ) |>
  mutate(
    ageDeath = censor_ages(ageDeath, spec),
    isHispanic = spec$missing,
    apoeGenotype = case_match(apoeGenotype,
      NA ~ spec$missing,
      .default = as.character(apoeGenotype)
    ),
    apoe4Status = get_apoe4Status(apoeGenotype, spec),
    amyCerad = spec$missing,
    amyAny = spec$missing,
    amyThal = case_match(amyThal,
      0 ~ spec$amyThal$none,
      1 ~ spec$amyThal$phase1,
      NA ~ spec$missing,
      .default = as.character(amyThal)
    ),
    amyA = get_amyA(amyThal, spec),
    Braak = to_Braak_stage(floor(Braak), spec),
    bScore = get_bScore(Braak, spec),
    dataContributionGroup = spec$dataContributionGroup$mayo,
    cohort = spec$cohort$mayo
  )

print_qc(meta_new)
validate_values(meta_new, spec)

new_filename <- write_metadata(meta_new, meta_file$name)
new_syn_id <- synapse_upload(new_filename, UPLOAD_SYNID)

manifest <- rbind(
  manifest,
  data.frame(study = "GEN-A12", metadata_synid = new_syn_id)
)


# GEN-A14 ----------------------------------------------------------------------

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
    isHispanic = case_when(race == "Hispanic" ~ spec$isHispanic$hisp_true,
      is.na(isHispanic) ~ spec$missing,
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
  data.frame(study = "GEN-A14", metadata_synid = new_syn_id)
)


# GEN-B4 -----------------------------------------------------------------------
# Uses Diverse Cohorts metadata

# TODO get the harmonized file from Synapse
dc_file <- file.path(
  "data", "output",
  "AMP-AD_DiverseCohorts_individual_metadata_harmonized.csv"
)

new_syn_id <- synapse_upload(dc_file, UPLOAD_SYNID)

manifest <- rbind(
  manifest,
  data.frame(study = "GEN-B4", metadata_synid = new_syn_id)
)


# Upload manifest file ---------------------------------------------------------

manifest <- manifest |> arrange(study)
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
      genesis_study = m_row[["study"]]
    )
  return(meta)
}, simplify = FALSE)

df_all <- purrr::list_rbind(df_list)

# df_all <- df_all |>
#  group_by(individualID) |>
#  mutate(genesis_study = paste(sort(genesis_study), collapse = "; "))

new_file <- write_metadata(df_all, "metadata_combined.csv")
synapse_upload(new_file, UPLOAD_SYNID)
