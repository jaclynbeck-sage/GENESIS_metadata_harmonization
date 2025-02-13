library(synapser)
library(dplyr)
library(purrr)
library(stringr)
library(readxl)

spec <- config::get(file = "GENESIS_harmonization.yml")
source("util_functions.R")

# TODO after harmonization, confirm that all harmonized fields have all correct values

syn_ids <- list(
  "GEN-A1" = "syn55251012",  # NPS-AD
  "GEN-A2" = "syn3191087",   # ROSMAP
  "GEN-A4" = "syn31149116",  # SEA-AD
  "GEN-A8" = "syn3191087",   # snRNAseqAD_TREM2, same individual metadata as GEN-A2
  "GEN-A9" = "syn22432749",  # SMIB-AD
  "GEN-A10" = "syn25891193", # MCMPS
  "GEN-A11" = "syn31563038", # MC_snRNA
  "GEN-A12" = "syn51401700", # MC-BrAD
  "GEN-A13" = "syn3191087",  # snRNAseqPFC_BA10, same individual metadata as GEN-A2
  "GEN-A14" = "syn24610550", # HBI_scRNAseq
  #"GEN-B1" = "TBD",
  #"GEN-B2" = "TBD",
  #"GEN-B3" = "TBD",
  "GEN-B4" = "syn51757646",  # AMP-AD_DiverseCohorts
  "GEN-B5" = "syn31149116",  # SEA-AD, same individual metadata as GEN-A4
  "GEN-B6" = "syn3191087"    # MIT_ROSMAP_Multiomics metadata is syn52430346 but the study uses the ROSMAP file (like GEN-A2)
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

new_filename <- write_metadata(meta_new, meta_file$name)
new_syn_id <- synapse_upload(new_filename, UPLOAD_SYNID)

manifest <- rbind(manifest, data.frame(study = "GEN-A1",
                                       metadata_synid = new_syn_id))


# GEN-A2
meta_file <- synapse_download(syn_ids[["GEN-A2"]])
meta <- read.csv(meta_file$path)

colnames(meta)

print_qc(meta, ageDeath_col = "age_death", sex_col = "msex", isHispanic_col = "spanish",
         apoe_col = "apoe_genotype", braak_col = "braaksc", cerad_col = "ceradsc")

meta_new <- meta %>%
  rename(isHispanic = spanish,
         ageDeath = age_death,
         sex = msex,
         apoeGenotype = apoe_genotype,
         amyCerad = ceradsc,
         Braak = braaksc,
         cohort = Study) %>%
  mutate(
    sex = case_when(sex == 1 ~ "male",
                    sex == 0 ~ "female",
                    is.na(sex) ~ spec$missing),
    ageDeath = censor_ages(ageDeath),
    race = case_when(is.na(race) ~ spec$missing,
                     race == 1 ~ "White",
                     race == 2 ~ "Black or African American",
                     race == 3 ~ "American Indian or Alaska Native",
                     race == 4 ~ "other", # Hawaiian / Pacific Islanders
                     race == 5 ~ "Asian",
                     race == 6 ~ "other",
                     race == 7 ~ spec$missing),
    isHispanic = case_when(isHispanic == 1 ~ "True",
                           isHispanic == 2 ~ "False",
                           is.na(isHispanic) ~ spec$missing),
    apoeGenotype = case_when(is.na(apoeGenotype) ~ spec$missing,
                             .default = as.character(apoeGenotype)),
    apoe4Status = get_apoe4Status(apoeGenotype),
    Braak = case_when(Braak == 0 ~ "None",
                      is.na(Braak) ~ spec$missing,
                      .default = paste("Stage", to_Roman_numerals(Braak))),
    bScore = get_bScore(Braak),
    amyCerad = case_when(amyCerad == 1 ~ "Frequent/Definite/C3",
                         amyCerad == 2 ~ "Moderate/Probable/C2",
                         amyCerad == 3 ~ "Sparse/Possible/C1",
                         amyCerad == 4 ~ "None/No AD/C0",
                         is.na(amyCerad) ~ spec$missing),
    amyAny = get_amyAny(amyCerad),
    amyThal = spec$missing,
    amyA = spec$missing,
    dataContributionGroup = "Rush"
  )

print_qc(meta_new)

new_filename <- write_metadata(meta_new, meta_file$name)
new_syn_id <- synapse_upload(new_filename, UPLOAD_SYNID)

manifest <- rbind(manifest,
                  data.frame(study = c("GEN-A2", "GEN-A8", "GEN-A13", "GEN-B6"),
                             metadata_synid = new_syn_id))


# GEN-A4
# This is SEA-AD data, which is more complicated than other data sets because
# the version on Synapse doesn't match the version distributed by the Allen
# Institute. We use the version on Synapse to determine what extra columns to
# keep because it's been curated/approved for the ADKP, but use the data that is
# in the AI file as it is more up to date.
meta_file <- synapse_download(syn_ids[["GEN-A4"]])
sea_ad_file <- file.path("data", "downloads", "sea-ad_cohort_donor_metadata_072524.xlsx")
download.file("https://brainmapportal-live-4cc80a57cd6e400d854-f7fdcae.divio-media.net/filer_public/b4/c7/b4c727e1-ede1-4c61-b2ee-bf1ae4a3ef68/sea-ad_cohort_donor_metadata_072524.xlsx",
              destfile = sea_ad_file)

meta <- read.csv(meta_file$path)
meta_sea_ad <- read_xlsx(sea_ad_file)
colnames(meta_sea_ad) <- make.names(colnames(meta_sea_ad), unique = TRUE)

colnames(meta)
colnames(meta_sea_ad)

# Remove columns from meta that will be added back by meta_sea_ad
meta <- meta %>%
  select(-sex, -race, -ethnicity, -ageDeath, -apoeGenotype, -apoe4Status, -pmi,
         -yearsEducation, -pH, -brainWeight, -CERAD:-LATE.NC.stage)

# Keep only columns from meta_sea_ad that are going into the new metadata
meta_sea_ad <- meta_sea_ad %>%
  select(-Primary.Study.Name, -Secondary.Study.Name, -specify.other.race,
         -Highest.level.of.education, -Age.of.onset.cognitive.symptoms,
         -Known.head.injury, -Have.they.had.neuroimaging,
         -Rapid.Frozen.Tissue.Type, -Ex.Vivo.Imaging,
         -Total.Microinfarcts..not.observed.grossly., -RIN,
         -Severely.Affected.Donor,
         -contains("Dx", ignore.case = TRUE)) # TODO eventually we want diagnosis

# SEA-AD has multiple race columns
race_cols <- grep("Race", colnames(meta_sea_ad), value = TRUE)

print_qc(meta_sea_ad,
         ageDeath_col = "Age.at.Death",
         isHispanic_col = "Hispanic.Latino",
         pmi_col = "PMI",
         sex_col = "Sex",
         race_col = race_cols,
         apoe_col = "APOE.Genotype",
         cerad_col = "CERAD.score",
         thal_col = "Thal")

# To make main mutate statement below more readable, collapse the race columns
# into one field
race_col_inds <- which(colnames(meta_sea_ad) %in% race_cols)
colnames(meta_sea_ad)[race_col_inds] <- colnames(meta_sea_ad)[race_col_inds] %>%
  str_replace("Race..choice.", "") %>%
  str_replace("\\.\\.", " or ") %>%
  str_replace_all("\\.", " ") %>%
  str_trim()

new_race_cols <- colnames(meta_sea_ad)[race_col_inds]
meta_sea_ad[, new_race_cols] <- meta_sea_ad[, new_race_cols] == "Checked"

meta_sea_ad$race <- apply(meta_sea_ad[, new_race_cols], 1, function(row) {
  paste(new_race_cols[row], collapse = ",")
})

meta_sea_ad <- select(meta_sea_ad, !all_of(new_race_cols))

meta_new <- merge(meta, meta_sea_ad, by.x = "individualID", by.y = "Donor.ID") %>%
  rename(ageDeath = Age.at.Death,
         isHispanic = Hispanic.Latino,
         sex = Sex,
         pmi = PMI,
         apoeGenotype = APOE.Genotype,
         amyCerad = CERAD.score,
         amyThal = Thal,
         yearsEducation = Years.of.education,
         pH = Brain.pH,
         brainWeight = Fresh.Brain.Weight,
         ADNC = Overall.AD.neuropathological.Change,
         Lewy.body.disease.pathology = Highest.Lewy.Body.Disease,
         LATE.NC.stage = LATE,
         Microinfarcts.in.screening.sections = Total.microinfarcts.in.screening.sections) %>%
  mutate(
    ageDeath = censor_ages(ageDeath),
    Age.of.Dementia.diagnosis = censor_ages(Age.of.Dementia.diagnosis, spec),
    isHispanic = case_when(isHispanic == "No" ~ spec$isHispanic$hisp_false,
                           isHispanic == "Yes" ~ spec$isHispanic$hisp_true,
                           .default = spec$missing),
    sex = str_to_lower(sex),
    race = case_when(race == "White,Other" ~ spec$race$other,
                     grepl("American Indian or Alaska Native", race) ~ spec$race$Amer_Ind,
                     .default = race),
    apoeGenotype = str_replace(apoeGenotype, "/", ""),
    apoe4Status = get_apoe4Status(apoeGenotype, spec),
    amyCerad = case_when(amyCerad == "Absent" ~ spec$amyCerad$none,
                         amyCerad == "Sparse" ~ spec$amyCerad$sparse,
                         amyCerad == "Moderate" ~ spec$amyCerad$moderate,
                         amyCerad == "Frequent" ~ spec$amyCerad$frequent),
    amyAny = get_amyAny(amyCerad, spec),
    amyThal = case_when(amyThal == "Thal 0" ~ spec$amyThal$none,
                        .default = str_replace(amyThal, "Thal", "Phase")),
    amyA = get_amyA(amyThal, spec),
    Braak = case_when(Braak == "Braak 0" ~ spec$Braak$none,
                      .default = str_replace(Braak, "Braak", "Stage")),
    bScore = get_bScore(Braak, spec),
    cohort = "SEA-AD",
    dataContributionGroup = "Allen Institute")

colnames(meta_new) <- str_replace(colnames(meta_new), "Last.", "") %>%
  str_replace(".in.months", "")

print_qc(meta_new)

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
                           .default = paste("Stage", to_Braak_stage(Braak, spec))),
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
         dataContributionGroup = "MSSM",
         cohort = "MSBB")

print_qc(meta_new)

new_filename <- write_metadata(meta_new, meta_file$name)
new_syn_id <- synapse_upload(new_filename, UPLOAD_SYNID)

manifest <- rbind(manifest,
                  data.frame(study = "GEN-A14", metadata_synid = new_syn_id))


# GEN-B4 - TODO this comes straight from Diverse Cohorts, consider using original DC file
meta_file <- synapse_download(syn_ids[["GEN-B4"]])
meta <- read.csv(meta_file$path)

colnames(meta)

print_qc(meta, pmi_col = "PMI")

meta_new <- meta %>%
  rename(pmi = PMI) %>%
  mutate(
    ageDeath = censor_ages(ageDeath, spec),
    pmi = case_when(pmi == spec$missing ~ NA,
                    .default = pmi),
    isHispanic = case_when(isHispanic == "TRUE" ~ spec$isHispanic$hisp_true,
                           isHispanic == "FALSE" ~ spec$isHispanic$hisp_false,
                           .default = spec$missing),
    race = case_when(race == "Black" ~ spec$race$Black,
                     .default = race),
    apoe4Status = get_apoe4Status(apoeGenotype, spec),
    bScore = get_bScore(Braak, spec))

print_qc(meta_new)

new_filename <- write_metadata(meta_new, meta_file$name)
new_syn_id <- synapse_upload(new_filename, UPLOAD_SYNID)

manifest <- rbind(manifest,
                  data.frame(study = "GEN-B4", metadata_synid = new_syn_id))


manifest <- manifest %>% arrange(study)
manifest_file <- file.path("data", "output", "metadata_manifest.csv")
write.csv(manifest, manifest_file,
          row.names = FALSE, quote = FALSE)
synapse_upload(manifest_file, UPLOAD_SYNID)


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
