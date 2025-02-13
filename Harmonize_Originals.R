library(synapser)
library(dplyr)
library(purrr)
library(stringr)
library(readxl)

spec <- config::get(file = "GENESIS_harmonization.yml")
source("util_functions.R")

syn_ids <- list("Diverse_Cohorts" = "syn51757646",
                "MayoRNAseq" = "syn23277389",
                "MSBB" = "syn6101474",
                "ROSMAP" = "syn3191087",
                "SEA-AD" = "syn31149116")

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

print_qc(meta, isHispanic_col = "ethnicity", cerad_col = "CERAD", thal_col = "Thal")

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


# MSBB -------------------------------------------------------------------------

# GEN-A1 has > 300 samples from the original MSBB metadata and 65 from ROSMAP.

meta_file <- synapse_download(syn_ids[["MSBB"]])
meta <- read.csv(meta_file$path)

colnames(meta)

print_qc(meta, isHispanic_col = "ethnicity", cerad_col = "CERAD")

meta_new <- meta %>%
  rename(isHispanic = ethnicity,
         amyCerad = CERAD,
         dataContributionGroup = individualIdSource) %>%
  mutate(
    ageDeath = censor_ages(ageDeath, spec),
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


# ROSMAP -----------------------------------------------------------------------

meta_file <- synapse_download(syn_ids[["ROSMAP"]])
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
    sex = case_when(sex == 1 ~ spec$sex$male,
                    sex == 0 ~ spec$sex$female,
                    is.na(sex) ~ spec$missing),
    ageDeath = censor_ages(ageDeath, spec),
    race = case_when(is.na(race) ~ spec$missing,
                     race == 1 ~ spec$race$White,
                     race == 2 ~ spec$race$Black,
                     race == 3 ~ spec$race$Amer_Ind,
                     race == 4 ~ spec$race$other, # Hawaiian / Pacific Islanders
                     race == 5 ~ spec$race$Asian,
                     race == 6 ~ spec$race$other,
                     race == 7 ~ spec$missing),
    isHispanic = case_when(isHispanic == 1 ~ spec$isHispanic$hisp_true,
                           isHispanic == 2 ~ spec$isHispanic$hisp_false,
                           is.na(isHispanic) ~ spec$missing),
    apoeGenotype = case_when(is.na(apoeGenotype) ~ spec$missing,
                             .default = as.character(apoeGenotype)),
    apoe4Status = get_apoe4Status(apoeGenotype, spec),
    Braak = case_when(is.na(Braak) ~ spec$missing,
                      Braak >= 0 ~ to_Braak_stage(Braak, spec),
                      .default = as.character(Braak)),
    bScore = get_bScore(Braak, spec),
    amyCerad = case_when(amyCerad == 1 ~ spec$amyCerad$frequent,
                         amyCerad == 2 ~ spec$amyCerad$moderate,
                         amyCerad == 3 ~ spec$amyCerad$sparse,
                         amyCerad == 4 ~ spec$amyCerad$none,
                         is.na(amyCerad) ~ spec$missing),
    amyAny = get_amyAny(amyCerad, spec),
    amyThal = spec$missing,
    amyA = spec$missing,
    dataContributionGroup = "Rush"
  )

print_qc(meta_new)

new_filename <- write_metadata(meta_new, file.path("originals", meta_file$name))


# Diverse Cohorts --------------------------------------------------------------

# Diverse Cohorts has data for samples from the original Mayo/MSBB/ROSMAP
# metadata but also has new samples. Some of the old/original samples have
# additional information in the DC metadata that doesn't appear in the original
# files, but are missing information in DC that do appear in the original files,
# so this metadata requires some case-by-case merging of data between the four
# metadata files.

meta_file <- synapse_download(syn_ids[["Diverse_Cohorts"]])
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
                           .default = isHispanic),
    race = case_when(race == "Black" ~ spec$race$Black,
                     .default = race),
    apoe4Status = get_apoe4Status(apoeGenotype, spec)
  )

print_qc(meta_new)

# TODO compare with Mayo/MSBB/ROSMAP
new_filename <- write_metadata(meta_new, file.path("originals", meta_file$name))
