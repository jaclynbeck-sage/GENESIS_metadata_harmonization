library(synapser)
library(dplyr)
library(purrr)
library(stringr)
library(readxl)

spec <- config::get(file = "GENESIS_harmonization.yml")
source("util_functions.R")
source("dataset_specific_functions.R")

syn_ids <- list(
  "GEN-A1" = "syn55251012.4",  # NPS-AD
  "GEN-A2" = "syn3191087.11",  # ROSMAP
  "GEN-A4" = "syn31149116.7",  # SEA-AD
  "GEN-A8" = "syn3191087.11",  # snRNAseqAD_TREM2, same individual metadata as GEN-A2
  "GEN-A9" = "syn22432749.1",  # SMIB-AD
  "GEN-A10" = "syn25891193.1", # MCMPS
  "GEN-A11" = "syn31563038.1", # MC_snRNA
  "GEN-A12" = "syn51401700.2", # MC-BrAD
  "GEN-A13" = "syn3191087.11", # snRNAseqPFC_BA10, same individual metadata as GEN-A2
  "GEN-A14" = "syn24610550.2", # HBI_scRNAseq
  #"GEN-B1" = "TBD",
  #"GEN-B2" = "TBD",
  #"GEN-B3" = "TBD",
  "GEN-B4" = "syn51757646.20", # AMP-AD_DiverseCohorts
  "GEN-B5" = "syn31149116.7",  # SEA-AD, same individual metadata as GEN-A4
  "GEN-B6" = "syn3191087.11"   # MIT_ROSMAP_Multiomics metadata is syn52430346 but the study uses the ROSMAP file (like GEN-A2)
)

UPLOAD_SYNID <- "syn64759869"

manifest <- c()

synLogin()

# GEN-A1
meta_file <- synapse_download(syn_ids[["GEN-A1"]])
meta <- read.csv(meta_file$path)

colnames(meta)

print_qc(meta, isHispanic_col = "ethnicity", pmi_col = "PMI", cerad_col = "CERAD")

# TODO should pull metadata from Diverse Cohorts, MSBB, and ROSMAP for the non-NPS-AD samples,
# since some info like Braak is missing from this metadata
meta_new <- meta %>%
  select(-Component) %>%
  rename(isHispanic = ethnicity,
         pmi = PMI,
         amyCerad = CERAD) %>%
  mutate(pmi = pmi / 60,
         pmiUnits = "hours",
         ageDeath = censor_ages(ageDeath, spec),
         race = case_when(is.na(race) ~ spec$missing,
                          .default = race),
         isHispanic = case_when(isHispanic == "Hispanic or Latino" ~ spec$isHispanic$hisp_true,
                                isHispanic == "Not Hispanic or Latino" ~ spec$isHispanic$hisp_false,
                                is.na(isHispanic) ~ spec$missing),
         apoeGenotype = case_when(is.na(apoeGenotype) ~ spec$missing,
                                  .default = as.character(apoeGenotype)),
         apoe4Status = get_apoe4Status(apoeGenotype, spec),
         amyCerad = case_when(amyCerad == 1 ~ spec$amyCerad$none,
                              amyCerad == 2 ~ spec$amyCerad$sparse,
                              amyCerad == 3 ~ spec$amyCerad$moderate,
                              amyCerad == 4 ~ spec$amyCerad$frequent,
                              is.na(amyCerad) ~ spec$missing),
         amyAny = get_amyAny(amyCerad, spec),
         amyThal = spec$missing,
         amyA = spec$missing,
         Braak = spec$missing,
         bScore = spec$missing,
         # TODO this might not be quite right
         dataContributionGroup = case_when(cohort == "MSBB" ~ "MSSM",
                                           cohort == "HBCC" ~ "NIMH Human Brain Collection Core",
                                           cohort == "ROSMAP" ~ "Rush"),
         cohort = case_when(cohort == "MSBB" ~ "Mt Sinai Brain Bank",
                            .default = cohort)) # TODO need to differentiate ROS and MAP

print_qc(meta_new)
validate_values(meta_new, spec)

new_filename <- write_metadata(meta_new, meta_file$name)
new_syn_id <- synapse_upload(new_filename, UPLOAD_SYNID)

manifest <- rbind(manifest, data.frame(study = "GEN-A1",
                                       metadata_synid = new_syn_id))


# GEN-A2
meta_file <- synapse_download(syn_ids[["GEN-A2"]])
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
validate_values(meta_new, spec)

new_filename <- write_metadata(meta_new, meta_file$name)
new_syn_id <- synapse_upload(new_filename, UPLOAD_SYNID)

manifest <- rbind(manifest,
                  data.frame(study = c("GEN-A2", "GEN-A8", "GEN-A13", "GEN-B6"),
                             metadata_synid = new_syn_id))


# GEN-A4
meta_file <- synapse_download(syn_ids[["GEN-A4"]])

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
validate_values(meta_new, spec)

new_filename <- write_metadata(meta_new, meta_file$name)
new_syn_id <- synapse_upload(new_filename, UPLOAD_SYNID)

manifest <- rbind(manifest,
                  data.frame(study = c("GEN-A4", "GEN-B5"),
                             metadata_synid = new_syn_id))


# GEN-A9
meta_file <- synapse_download(syn_ids[["GEN-A9"]])
meta <- read.csv(meta_file$path)

colnames(meta)

print_qc(meta, isHispanic_col = "ethnicity", cerad_col = "CERAD")

meta_new <- meta %>%
  rename(isHispanic = ethnicity,
         amyCerad = CERAD,
         cohort = individualIdSource) %>%
  mutate(ageDeath = censor_ages(ageDeath, spec),
         race = "White",
         isHispanic = "False",
         apoeGenotype = as.character(apoeGenotype),
         apoe4Status = get_apoe4Status(apoeGenotype, spec),
         Braak = case_when(is.na(Braak) ~ spec$missing,
                           .default = to_Braak_stage(Braak, spec)),
         bScore = get_bScore(Braak, spec),
         amyCerad = spec$missing, # TODO
         amyAny = get_amyAny(amyCerad, spec),
         amyThal = spec$missing,
         amyA = spec$missing,
         cohort = case_when(is.na(cohort) ~ "SMRI",
                            .default = "Banner"),
         dataContributionGroup = case_when(cohort == "SMRI" ~ "Stanley Medical Research Institute",
                                           .default = "Banner Sun Health Research Institute"))


print_qc(meta_new)
validate_values(meta_new, spec)

new_filename <- write_metadata(meta_new, meta_file$name)
new_syn_id <- synapse_upload(new_filename, UPLOAD_SYNID)

manifest <- rbind(manifest,
                  data.frame(study = "GEN-A9", metadata_synid = new_syn_id))


# GEN-A10
meta_file <- synapse_download(syn_ids[["GEN-A10"]])
meta <- read.csv(meta_file$path)

colnames(meta)

print_qc(meta, isHispanic_col = "ethnicity", cerad_col = "CERAD")

meta_new <- meta %>%
  rename(isHispanic = ethnicity,
         amyCerad = CERAD,
         cohort = individualIdSource) %>%
  mutate(race = str_trim(race),
         isHispanic = case_when(str_trim(isHispanic) == "Hispanic or Latino" ~ spec$isHispanic$hisp_true,
                                isHispanic == "Not Hispanic or Latino" ~ spec$isHispanic$hisp_false,
                                isHispanic == "Middle Eastern" ~ spec$isHispanic$hisp_false),
         apoeGenotype = case_when(is.na(apoeGenotype) ~ spec$missing,
                                  .default = as.character(apoeGenotype)),
         apoe4Status = get_apoe4Status(apoeGenotype, spec),
         amyCerad = spec$missing,
         amyAny = spec$missing,
         amyThal = spec$missing,
         amyA = spec$missing,
         Braak = spec$missing,
         bScore = spec$missing,
         dataContributionGroup = "Mayo")

print_qc(meta_new)
validate_values(meta_new, spec)

new_filename <- write_metadata(meta_new, meta_file$name)
new_syn_id <- synapse_upload(new_filename, UPLOAD_SYNID)

manifest <- rbind(manifest,
                  data.frame(study = "GEN-A10", metadata_synid = new_syn_id))


# GEN-A11
meta_file <- synapse_download(syn_ids[["GEN-A11"]])
meta <- read.csv(meta_file$path)

colnames(meta)

print_qc(meta, isHispanic_col = "ethnicity", cerad_col = "CERAD")

meta_new <- meta %>%
  rename(isHispanic = ethnicity,
         amyCerad = CERAD) %>%
  mutate(ageDeath = censor_ages(ageDeath, spec),
         isHispanic = "False",
         apoeGenotype = as.character(apoeGenotype),
         apoe4Status = get_apoe4Status(apoeGenotype, spec),
         amyCerad = spec$missing,
         amyAny = spec$missing,
         amyThal = spec$missing,
         amyA = spec$missing,
         Braak = case_when(is.na(Braak) ~ spec$missing,
                           .default = to_Braak_stage(floor(Braak), spec)),
         bScore = get_bScore(Braak, spec),
         dataContributionGroup = "Mayo",
         cohort = "Mayo Clinic")

print_qc(meta_new)
validate_values(meta_new, spec)

new_filename <- write_metadata(meta_new, meta_file$name)
new_syn_id <- synapse_upload(new_filename, UPLOAD_SYNID)

manifest <- rbind(manifest,
                  data.frame(study = "GEN-A11", metadata_synid = new_syn_id))


# GEN-A12
meta_file <- synapse_download(syn_ids[["GEN-A12"]])
meta <- read.csv(meta_file$path)

colnames(meta)

print_qc(meta, isHispanic_col = "ethnicity", cerad_col = "CERAD", thal_col = "Thal")

meta_new <- meta %>%
  rename(isHispanic = ethnicity,
         amyCerad = CERAD,
         amyThal = Thal) %>%
  mutate(ageDeath = censor_ages(ageDeath, spec),
         isHispanic = spec$missing,
         apoeGenotype = case_when(is.na(apoeGenotype) ~ spec$missing,
                                  .default = as.character(apoeGenotype)),
         apoe4Status = get_apoe4Status(apoeGenotype, spec),
         amyCerad = spec$missing,
         amyAny = spec$missing,
         amyThal = case_when(amyThal == 0 ~ spec$amyThal$none,
                             amyThal == 1 ~ spec$amyThal$phase1,
                             is.na(amyThal) ~ spec$missing),
         amyA = get_amyA(amyThal, spec),
         Braak = case_when(is.na(Braak) ~ spec$missing,
                           Braak >= 0 ~ to_Braak_stage(floor(Braak), spec),
                           .default = as.character(Braak)),
         bScore = get_bScore(Braak, spec),
         dataContributionGroup = "Mayo",
         cohort = "Mayo Clinic")

print_qc(meta_new)
validate_values(meta_new, spec)

new_filename <- write_metadata(meta_new, meta_file$name)
new_syn_id <- synapse_upload(new_filename, UPLOAD_SYNID)

manifest <- rbind(manifest,
                  data.frame(study = "GEN-A12", metadata_synid = new_syn_id))


# GEN-A14
meta_file <- synapse_download(syn_ids[["GEN-A14"]])
meta <- read.csv(meta_file$path)

colnames(meta)

print_qc(meta, isHispanic_col = "ethnicity", cerad_col = "CERAD")

meta_new <- meta %>%
  rename(isHispanic = ethnicity,
         amyCerad = CERAD) %>%
  mutate(pmi = pmi / 60,
         ageDeath = censor_ages(ageDeath, spec),
         isHispanic = case_when(race == "Hispanic" ~ spec$isHispanic$hisp_true,
                                .default = spec$missing),
         race = case_when(race == "Black" ~ spec$race$Black,
                          race == "Hispanic" ~ spec$missing,
                          .default = race),
         apoeGenotype = spec$missing,
         apoe4Status = get_apoe4Status(apoeGenotype, spec),
         Braak = case_when(is.na(Braak) ~ spec$missing,
                           Braak >= 0 ~ to_Braak_stage(Braak, spec),
                           .default = as.character(Braak)),
         bScore = get_bScore(Braak, spec),
         amyCerad = spec$missing, # TODO
         amyAny = get_amyAny(amyCerad, spec),
         amyThal = spec$missing,
         amyA = spec$missing,
         dataContributionGroup = spec$dataContributionGroup$mssm,
         cohort = spec$cohort$msbb)

print_qc(meta_new)
validate_values(meta_new, spec)

new_filename <- write_metadata(meta_new, meta_file$name)
new_syn_id <- synapse_upload(new_filename, UPLOAD_SYNID)

manifest <- rbind(manifest,
                  data.frame(study = "GEN-A14", metadata_synid = new_syn_id))


# GEN-B4
meta_file <- synapse_download(syn_ids[["GEN-B4"]])
meta <- read.csv(meta_file$path)

colnames(meta)

print_qc(meta, pmi_col = "PMI")

meta_new <- harmonize_Diverse_Cohorts(meta, spec)

print_qc(meta_new)
validate_values(meta_new, spec)

new_filename <- write_metadata(meta_new, meta_file$name)
new_syn_id <- synapse_upload(new_filename, UPLOAD_SYNID)

manifest <- rbind(manifest,
                  data.frame(study = "GEN-B4", metadata_synid = new_syn_id))


# Upload manifest file

manifest <- manifest %>% arrange(study)
manifest_file <- file.path("data", "output", "metadata_manifest.csv")
write.csv(manifest, manifest_file,
          row.names = FALSE, quote = FALSE)
synapse_upload(manifest_file, UPLOAD_SYNID)


# Combine all harmonized data into one data frame

df_list <- apply(manifest, 1, function(m_row) {
  m_file <- synapse_download(m_row[["metadata_synid"]])
  meta <- read.csv(m_file$path) %>%
    mutate(individualID = as.character(individualID),
           apoeGenotype = as.character(apoeGenotype),
           amyAny = as.character(amyAny),
           genesis_study = m_row[["study"]])
  return(meta)
}, simplify = FALSE)

df_all <- purrr::list_rbind(df_list)

df_all <- df_all %>%
  group_by(individualID) %>%
  mutate(genesis_study = paste(sort(genesis_study), collapse = "; "))

new_file <- write_metadata(df_all, "metadata_combined.csv")
synapse_upload(new_file, UPLOAD_SYNID)
