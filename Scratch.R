library(synapser)
library(dplyr)
library(stringr)

source("util_functions.R")

# TODO after harmonization, confirm that all harmonized fields have all correct values

syn_ids <- list(
  "GEN-A1" = "syn55251012",  # NPS-AD
  "GEN-A2" = "syn3191087",   # ROSMAP
  "GEN-A4" = "syn31149116",  # SEA-AD
  "GEN-A8" = "syn3191087",   # snRNAseqAD_TREM2, same individual metadata as GEN-A2
  #"GEN-A8" = "syn22415705"  # This has a duplicate label and is mouse data?
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

expectedColumns <- c("individualID", "dataContributionGroup", "cohort", "sex",
                     "race", "isHispanic", "ageDeath", "pmi")
optionalColumns <- c("apoeGenotype", "amyCerad", "amyAny", "Braak", "bScore")

sex_values <- c("female", "male")
race_values <- c("American Indian or Alaska Native",
                 "Asian",
                 "Black or African American",
                 "White",
                 "other",
                 "Missing or unknown")
hispanic_values <- c("True", "False", "Missing or unknown")

synLogin()

# GEN-A1
meta_file <- synGet(syn_ids[["GEN-A1"]], 
                    downloadLocation = file.path("data", "downloads"),
                    ifcollision = "overwrite.local")
meta <- read.csv(meta_file$path)

colnames(meta)

setdiff(expectedColumns, colnames(meta))
setdiff(optionalColumns, colnames(meta))

print_qc(meta, isHispanic_col = "ethnicity", pmi_col = "PMI", cerad_col = "CERAD")

# TODO should pull metadata from Diverse Cohorts, MSBB, and ROSMAP for the non-NPS-AD samples,
# since some info like Braak is missing from this metadata
meta_new <- meta %>%
  select(-Component) %>%
  rename(isHispanic = ethnicity,
         pmi = PMI,
         amyCerad = CERAD,
         # TODO Check if this is correct, it might be things like "Rush", "MSBB", etc
         dataContributionGroup = individualMetadataSource) %>% 
  mutate(pmi = pmi / 60,
         pmiUnits = "hours",
         ageDeath = case_when(ageDeath == "89+" ~ "90+",
                              .default = ageDeath),
         race = case_when(is.na(race) ~ "Missing or unknown",
                          !is.na(race) ~ race),
         isHispanic = case_when(isHispanic == "Hispanic or Latino" ~ "True",
                                isHispanic == "Not Hispanic or Latino" ~ "False",
                                is.na(isHispanic) ~ "Missing or unknown"),
         apoeGenotype = case_when(is.na(apoeGenotype) ~ "Missing or unknown",
                                  .default = as.character(apoeGenotype)),
         amyCerad = case_when(amyCerad == 1 ~ "None/No AD/C0",
                              amyCerad == 2 ~ "Sparse/Possible/C1",
                              amyCerad == 3 ~ "Moderate/Probable/C2",
                              amyCerad == 4 ~ "Frequent/Definite/C3",
                              is.na(amyCerad) ~ "Missing or unknown"),
         amyAny = case_when(amyCerad == "None/No AD/C0" ~ "0",
                            amyCerad == "Missing or unknown" ~ amyCerad,
                            .default = "1"),
         Braak = "Missing or unknown",
         bScore = "Missing or unknown")

print_qc(meta_new)
print_columns_with_nas(meta_new)

write.csv(meta_new, file.path("data", "output", 
                              str_replace(meta_file$name, ".csv", "_harmonized.csv")))



# GEN-A2
meta_file <- synGet(syn_ids[["GEN-A2"]], 
                    downloadLocation = file.path("data", "downloads"),
                    ifcollision = "overwrite.local")
meta <- read.csv(meta_file$path)

colnames(meta)

setdiff(expectedColumns, colnames(meta))

print_qc(meta, ageDeath_col = "age_death", sex_col = "msex", isHispanic_col = "spanish",
         apoe_col = "apoe_genotype", braak_col = "braaksc", cerad_col = "ceradsc")

meta_new <- meta %>%
  rename(isHispanic = spanish,
         ageDeath = age_death,
         sex = msex,
         apoeGenotype = apoe_genotype,
         amyCerad = ceradsc,
         Braak = braaksc,
         # TODO Check if this is correct
         cohort = Study) %>% 
  mutate(
    sex = case_when(sex == 1 ~ "male",
                    sex == 0 ~ "female",
                    is.na(sex) ~ "Missing or unknown"),
    ageDeath = case_when(ageDeath == "" ~ "Missing or unknown",
                         .default = ageDeath),
    pmi = case_when(is.na(pmi) ~ "Missing or unknown",
                    .default = as.character(pmi)),
    race = case_when(is.na(race) ~ "Missing or unknown",
                     race == 1 ~ "White",
                     race == 2 ~ "Black or African American",
                     race == 3 ~ "American Indian or Alaska Native",
                     race == 4 ~ "??",
                     race == 5 ~ "Asian",
                     race == 6 ~ "other",
                     race == 7 ~ "Missing or unknown"),
    isHispanic = case_when(isHispanic == 1 ~ "True",
                           isHispanic == 2 ~ "False",
                           is.na(isHispanic) ~ "Missing or unknown"),
    apoeGenotype = case_when(is.na(apoeGenotype) ~ "Missing or unknown",
                             .default = as.character(apoeGenotype)),
    Braak = case_when(Braak == 0 ~ "None",
                      is.na(Braak) ~ "Missing or unknown",
                      .default = paste("Stage", to_Roman_numerals(Braak))),
    bScore = get_bScore(Braak),
    amyCerad = case_when(amyCerad == 1 ~ "Frequent/Definite/C3",
                         amyCerad == 2 ~ "Moderate/Probable/C2",
                         amyCerad == 3 ~ "Sparse/Possible/C1",
                         amyCerad == 4 ~ "None/No AD/C0",
                         is.na(amyCerad) ~ "Missing or unknown"),
    amyAny = get_amyAny(amyCerad),
    # TODO
    dataContributionGroup = "TBD"
  )

print_qc(meta_new)

write.csv(meta_new, file.path("data", "output", 
                              str_replace(meta_file$name, ".csv", "_harmonized.csv")))



# GEN-A4
meta_file <- synGet(syn_ids[["GEN-A4"]], 
                    downloadLocation = file.path("data", "downloads"),
                    ifcollision = "overwrite.local")
meta <- read.csv(meta_file$path)

colnames(meta)

setdiff(expectedColumns, colnames(meta))

print_qc(meta, isHispanic_col = "ethnicity", cerad_col = "CERAD")

tmp <- subset(meta, ageDeath != "90+")
any(as.numeric(tmp$ageDeath) >= 90)
any(is.na(meta$ageDeath))

meta_new <- meta %>%
  rename(isHispanic = TBD,
         # TODO Check if this is correct
         cohort = dataset,
         # TODO Check if this is correct
         dataContributionGroup = individualIdSource) %>% 
  mutate(pmi = case_when(is.na(pmi) ~ "Missing or unknown",
                         .default = pmi))
         # TODO isHispanic and race

print_qc(meta_new)



# GEN-A9
meta_file <- synGet(syn_ids[["GEN-A9"]], 
                    downloadLocation = file.path("data", "downloads"),
                    ifcollision = "overwrite.local")
meta <- read.csv(meta_file$path)

colnames(meta)

setdiff(expectedColumns, colnames(meta))

print_qc(meta, isHispanic_col = "ethnicity", cerad_col = "CERAD")

meta_new <- meta %>%
  rename(isHispanic = ethnicity,
         amyCerad = CERAD,
         dataContributionGroup = individualIdSource) %>% 
  mutate(pmi = case_when(is.na(pmi) ~ "Missing or unknown",
                         .default = as.character(pmi)),
         ageDeath = case_when(ageDeath == 90 ~ "90+",
                              .default = ageDeath),
         race = "White",
         isHispanic = "False",
         apoeGenotype = as.character(apoeGenotype),
         Braak = case_when(is.na(Braak) ~ "Missing or unknown",
                           .default = paste("Stage", to_Roman_numerals(Braak))),
         bScore = get_bScore(Braak),
         amyCerad = "TBD", # TODO
         amyAny = get_amyAny(amyCerad),
         dataContributionGroup = case_when(is.na(dataContributionGroup) ~ "Stanley Medical Research Institute",
                                           dataContributionGroup == "BannerSun" ~ "Banner Sun Health Research Institute"),
         cohort = dataContributionGroup)

print_qc(meta_new)

write.csv(meta_new, file.path("data", "output", 
                              str_replace(meta_file$name, ".csv", "_harmonized.csv")))


# GEN-A10
meta_file <- synGet(syn_ids[["GEN-A10"]], 
                    downloadLocation = file.path("data", "downloads"),
                    ifcollision = "overwrite.local")
meta <- read.csv(meta_file$path)

colnames(meta)

setdiff(expectedColumns, colnames(meta))

print_qc(meta, isHispanic_col = "ethnicity", cerad_col = "CERAD")

meta_new <- meta %>%
  rename(isHispanic = ethnicity,
         amyCerad = CERAD,
         dataContributionGroup = individualIdSource) %>% 
  mutate(race = str_trim(race),
         isHispanic = case_when(str_trim(isHispanic) == "Hispanic or Latino" ~ "True",
                                isHispanic == "Not Hispanic or Latino" ~ "False",
                                isHispanic == "Middle Eastern" ~ "False"),
         apoeGenotype = case_when(is.na(apoeGenotype) ~ "Missing or unknown",
                                  .default = as.character(apoeGenotype)),
         amyCerad = "Missing or unknown",
         amyAny = "Missing or unknown",
         Braak = "Missing or unknown",
         bScore = "Missing or unknown",
         cohort = dataContributionGroup)

print_qc(meta_new)

write.csv(meta_new, file.path("data", "output", 
                              str_replace(meta_file$name, ".csv", "_harmonized.csv")))


# GEN-A11
meta_file <- synGet(syn_ids[["GEN-A11"]], 
                    downloadLocation = file.path("data", "downloads"),
                    ifcollision = "overwrite.local")
meta <- read.csv(meta_file$path)

colnames(meta)

setdiff(expectedColumns, colnames(meta))

print_qc(meta, isHispanic_col = "ethnicity", cerad_col = "CERAD")

meta_new <- meta %>%
  rename(isHispanic = ethnicity,
         amyCerad = CERAD,
         dataContributionGroup = individualIdSource) %>% 
  mutate(ageDeath = case_when(ageDeath == "90_or_over" ~ "90+",
                              .default = ageDeath),
         isHispanic = "False",
         apoeGenotype = as.character(apoeGenotype),
         amyCerad = "Missing or unknown",
         amyAny = "Missing or unknown",
         Braak = case_when(is.na(Braak) ~ "Missing or unknown",
                           .default = paste("Stage", to_Roman_numerals(floor(Braak)))),
         bScore = get_bScore(Braak),
         cohort = dataContributionGroup)

print_qc(meta_new)

write.csv(meta_new, file.path("data", "output", 
                              str_replace(meta_file$name, ".csv", "_harmonized.csv")))


# GEN-A12
meta_file <- synGet(syn_ids[["GEN-A12"]], 
                    downloadLocation = file.path("data", "downloads"),
                    ifcollision = "overwrite.local")
meta <- read.csv(meta_file$path)

colnames(meta)

setdiff(expectedColumns, colnames(meta))

print_qc(meta, isHispanic_col = "ethnicity", cerad_col = "CERAD")

meta_new <- meta %>%
  rename(isHispanic = ethnicity,
         amyCerad = CERAD,
         dataContributionGroup = individualIdSource) %>% 
  mutate(ageDeath = case_when(ageDeath == "90_or_over" ~ "90+",
                              .default = ageDeath),
         isHispanic = "Missing or unknown",
         apoeGenotype = case_when(is.na(apoeGenotype) ~ "Missing or unknown",
                                  .default = as.character(apoeGenotype)),
         amyCerad = "Missing or unknown",
         amyAny = "Missing or unknown",
         Braak = case_when(is.na(Braak) ~ "Missing or unknown",
                           Braak == 0 ~ "None",
                           .default = paste("Stage", to_Roman_numerals(floor(Braak)))),
         bScore = get_bScore(Braak),
         cohort = dataContributionGroup)

print_qc(meta_new)

write.csv(meta_new, file.path("data", "output", 
                              str_replace(meta_file$name, ".csv", "_harmonized.csv")))


# GEN-A14
meta_file <- synGet(syn_ids[["GEN-A14"]], 
                    downloadLocation = file.path("data", "downloads"),
                    ifcollision = "overwrite.local")
meta <- read.csv(meta_file$path)

colnames(meta)

setdiff(expectedColumns, colnames(meta))

print_qc(meta, isHispanic_col = "ethnicity", cerad_col = "CERAD")

meta_new <- meta %>%
  rename(isHispanic = ethnicity,
         amyCerad = CERAD,
         dataContributionGroup = individualIdSource) %>% 
  mutate(pmi = pmi / 60,
         ageDeath = case_when(ageDeath == 90 ~ "90+",
                              .default = ageDeath),
         isHispanic = case_when(race == "Hispanic" ~ "True",
                                .default = "Missing or unknown"),
         race = case_when(race == "Black" ~ "Black or African American",
                          race == "Hispanic" ~ "Missing or unknown",
                          .default = race),
         apoeGenotype = "Missing or unknown",
         Braak = case_when(is.na(Braak) ~ "Missing or unknown",
                           Braak == 0 ~ "None",
                           .default = paste("Stage", to_Roman_numerals(Braak))),
         bScore = get_bScore(Braak),
         amyCerad = "TBD",
         amyAny = "TBD",
         cohort = dataContributionGroup)

print_qc(meta_new)

write.csv(meta_new, file.path("data", "output", 
                              str_replace(meta_file$name, ".csv", "_harmonized.csv")))


# GEN-B4 - TODO this comes straight from Diverse Cohorts, consider using original DC file
meta_file <- synGet(syn_ids[["GEN-B4"]], 
                    downloadLocation = file.path("data", "downloads"),
                    ifcollision = "overwrite.local")
meta <- read.csv(meta_file$path)

colnames(meta)

setdiff(expectedColumns, colnames(meta))

print_qc(meta, pmi_col = "PMI")

meta_new <- meta %>%
  rename(pmi = PMI) %>% 
  mutate(pmi = case_when(pmi == "missing or unknown" ~ NA,
                         .default = as.numeric(pmi)),
         ageDeath = case_when(ageDeath == "missing or unknown" ~ NA,
                              .default = ageDeath),
         isHispanic = case_when(isHispanic == "TRUE" ~ "True",
                                isHispanic == "FALSE" ~ "False",
                                .default = str_to_sentence(isHispanic)),
         race = case_when(race == "Black" ~ "Black or African American",
                          race == "missing or unknown" ~ "Missing or unknown",
                          .default = race),
         apoeGenotype = str_to_sentence(apoeGenotype),
         amyCerad = case_when(amyCerad == "missing or unknown" ~ "Missing or unknown",
                              .default = amyCerad),
         amyAny = get_amyAny(amyCerad),
         Braak = case_when(Braak == "missing or unknown" ~ "Missing or unknown",
                           .default = Braak),
         bScore = get_bScore(Braak))

print_qc(meta_new)

write.csv(meta_new, file.path("data", "output", 
                              str_replace(meta_file$name, ".csv", "_harmonized.csv")))

