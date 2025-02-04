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

print_unique_column_vals(meta, c("ageDeath", "PMI", "pmi"))
print_columns_with_nas(meta)

tmp <- subset(meta, ageDeath != "89+")
any(as.numeric(tmp$ageDeath) >= 89)

meta_new <- meta %>%
  rename(isHispanic = ethnicity,
         pmi = PMI,
         # TODO Check if this is correct, it might be things like "Rush", "MSBB", etc
         dataContributionGroup = individualMetadataSource) %>% 
  mutate(pmi = pmi / 60,
         pmiUnits = "hours",
         ageDeath = case_when(ageDeath == "89+" ~ "90+",
                              TRUE ~ ageDeath),
         race = case_when(is.na(race) ~ "Missing or unknown",
                          !is.na(race) ~ race),
         isHispanic = case_when(isHispanic == "Hispanic or Latino" ~ "True",
                                isHispanic == "Not Hispanic or Latino" ~ "False",
                                is.na(isHispanic) ~ "Missing or unknown"))

print_unique_column_vals(meta_new, c("ageDeath", "PMI", "pmi"))
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

print_unique_column_vals(meta, c("ageDeath", "pmi", "age_at_visit_max", 
                                 "age_first_ad_dx", "age_death", "cts_mmse30_lv", 
                                 "cts_mmse30_first_ad_dx"))
print_columns_with_nas(meta)

tmp <- subset(meta, age_death != "" & age_death != "90+")
any(as.numeric(tmp$age_death) >= 90)

meta_new <- meta %>%
  rename(isHispanic = spanish,
         ageDeath = age_death,
         sex = msex,
         # TODO Check if this is correct
         cohort = Study) %>% 
  mutate(
    sex = case_when(sex == 1 ~ "male",
                    sex == 0 ~ "female",
                    is.na(sex) ~ "Missing or unknown"),
    ageDeath = case_when(ageDeath == "" ~ "Missing or unknown",
                         TRUE ~ ageDeath),
    pmi = case_when(is.na(pmi) ~ "Missing or unknown",
                    TRUE ~ as.character(pmi)),
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
    # TODO
    dataContributionGroup = "TBD"
  )

print_unique_column_vals(meta_new, c("ageDeath", "pmi", "age_at_visit_max", 
                                     "age_first_ad_dx", "age_death", "cts_mmse30_lv", 
                                     "cts_mmse30_first_ad_dx"))
print_columns_with_nas(meta_new)

write.csv(meta_new, file.path("data", "output", 
                              str_replace(meta_file$name, ".csv", "_harmonized.csv")))



# GEN-A4
meta_file <- synGet(syn_ids[["GEN-A4"]], 
                    downloadLocation = file.path("data", "downloads"),
                    ifcollision = "overwrite.local")
meta <- read.csv(meta_file$path)

colnames(meta)

setdiff(expectedColumns, colnames(meta))

print_unique_column_vals(meta, c("ageDeath", "pmi"))
print_columns_with_nas(meta)

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
                         TRUE ~ pmi))
         # TODO isHispanic and race

print_unique_column_vals(meta_new, c("ageDeath", "pmi"))
print_columns_with_nas(meta_new)



# GEN-A9
meta_file <- synGet(syn_ids[["GEN-A9"]], 
                    downloadLocation = file.path("data", "downloads"),
                    ifcollision = "overwrite.local")
meta <- read.csv(meta_file$path)

colnames(meta)

setdiff(expectedColumns, colnames(meta))

print_unique_column_vals(meta, c("ageDeath", "pmi"))
print_columns_with_nas(meta)

tmp <- subset(meta, ageDeath != "90+")
any(as.numeric(tmp$ageDeath) >= 90)
tmp$ageDeath[tmp$ageDeath >= 90]

meta_new <- meta %>%
  rename(isHispanic = ethnicity,
         dataContributionGroup = individualIdSource) %>% 
  mutate(pmi = case_when(is.na(pmi) ~ "Missing or unknown",
                         TRUE ~ as.character(pmi)),
         ageDeath = case_when(ageDeath == 90 ~ "90+",
                              TRUE ~ ageDeath),
         race = "White",
         isHispanic = "False",
         dataContributionGroup = case_when(is.na(dataContributionGroup) ~ "Stanley Medical Research Institute",
                                           dataContributionGroup == "BannerSun" ~ "Banner Sun Health Research Institute"),
         cohort = dataContributionGroup)

print_unique_column_vals(meta_new, c("ageDeath", "pmi"))
print_columns_with_nas(meta_new)

write.csv(meta_new, file.path("data", "output", 
                              str_replace(meta_file$name, ".csv", "_harmonized.csv")))


# GEN-A10
meta_file <- synGet(syn_ids[["GEN-A10"]], 
                    downloadLocation = file.path("data", "downloads"),
                    ifcollision = "overwrite.local")
meta <- read.csv(meta_file$path)

colnames(meta)

setdiff(expectedColumns, colnames(meta))

print_unique_column_vals(meta, c("ageDeath", "pmi"))
print_columns_with_nas(meta)

meta_new <- meta %>%
  rename(isHispanic = ethnicity,
         dataContributionGroup = individualIdSource) %>% 
  mutate(pmi = "TBD",
         ageDeath = "TBD",
         race = str_trim(race),
         isHispanic = case_when(str_trim(isHispanic) == "Hispanic or Latino" ~ "True",
                                isHispanic == "Not Hispanic or Latino" ~ "False",
                                isHispanic == "Middle Eastern" ~ "False"),
         cohort = dataContributionGroup)

print_unique_column_vals(meta_new, c("ageDeath", "pmi"))
print_columns_with_nas(meta_new)

write.csv(meta_new, file.path("data", "output", 
                              str_replace(meta_file$name, ".csv", "_harmonized.csv")))


# GEN-A11
meta_file <- synGet(syn_ids[["GEN-A11"]], 
                    downloadLocation = file.path("data", "downloads"),
                    ifcollision = "overwrite.local")
meta <- read.csv(meta_file$path)

colnames(meta)

setdiff(expectedColumns, colnames(meta))

print_unique_column_vals(meta, c("ageDeath", "pmi"))
print_columns_with_nas(meta)

tmp <- subset(meta, ageDeath != "90_or_over")
any(as.numeric(tmp$ageDeath) >= 90)

meta_new <- meta %>%
  rename(isHispanic = ethnicity,
         dataContributionGroup = individualIdSource) %>% 
  mutate(pmi = "Missing or unknown",
         ageDeath = case_when(ageDeath == "90_or_over" ~ "90+",
                              TRUE ~ ageDeath),
         isHispanic = "False",
         cohort = dataContributionGroup)

print_unique_column_vals(meta_new, c("ageDeath", "pmi"))
print_columns_with_nas(meta_new)

write.csv(meta_new, file.path("data", "output", 
                              str_replace(meta_file$name, ".csv", "_harmonized.csv")))


# GEN-A12
meta_file <- synGet(syn_ids[["GEN-A12"]], 
                    downloadLocation = file.path("data", "downloads"),
                    ifcollision = "overwrite.local")
meta <- read.csv(meta_file$path)

colnames(meta)

setdiff(expectedColumns, colnames(meta))

print_unique_column_vals(meta, c("ageDeath", "pmi"))
print_columns_with_nas(meta)

tmp <- subset(meta, ageDeath != "90_or_over")
any(as.numeric(tmp$ageDeath) >= 90)

meta_new <- meta %>%
  rename(isHispanic = ethnicity,
         dataContributionGroup = individualIdSource) %>% 
  mutate(pmi = "Missing or unknown",
         ageDeath = case_when(ageDeath == "90_or_over" ~ "90+",
                              TRUE ~ ageDeath),
         isHispanic = "Missing or unknown",
         cohort = dataContributionGroup)

print_unique_column_vals(meta_new, c("ageDeath", "pmi"))
print_columns_with_nas(meta_new)

write.csv(meta_new, file.path("data", "output", 
                              str_replace(meta_file$name, ".csv", "_harmonized.csv")))


# GEN-A14
meta_file <- synGet(syn_ids[["GEN-A14"]], 
                    downloadLocation = file.path("data", "downloads"),
                    ifcollision = "overwrite.local")
meta <- read.csv(meta_file$path)

colnames(meta)

setdiff(expectedColumns, colnames(meta))

print_unique_column_vals(meta, c("ageDeath", "pmi"))
print_columns_with_nas(meta)

tmp <- subset(meta, ageDeath != "90+")
any(as.numeric(tmp$ageDeath) >= 90)
tmp$ageDeath[tmp$ageDeath >= 90]

meta_new <- meta %>%
  rename(isHispanic = ethnicity,
         dataContributionGroup = individualIdSource) %>% 
  mutate(pmi = pmi / 60,
         ageDeath = case_when(ageDeath == 90 ~ "90+",
                              TRUE ~ ageDeath),
         isHispanic = case_when(race == "Hispanic" ~ "True",
                                TRUE ~ "Missing or unknown"),
         race = case_when(race == "Black" ~ "Black or African American",
                          race == "Hispanic" ~ "Missing or unknown",
                          TRUE ~ race),
         cohort = dataContributionGroup)

print_unique_column_vals(meta_new, c("ageDeath", "pmi"))
print_columns_with_nas(meta_new)

write.csv(meta_new, file.path("data", "output", 
                              str_replace(meta_file$name, ".csv", "_harmonized.csv")))


# GEN-B4
meta_file <- synGet(syn_ids[["GEN-B4"]], 
                    downloadLocation = file.path("data", "downloads"),
                    ifcollision = "overwrite.local")
meta <- read.csv(meta_file$path)

colnames(meta)

setdiff(expectedColumns, colnames(meta))

print_unique_column_vals(meta, c("ageDeath", "PMI"))
print_columns_with_nas(meta)

tmp <- subset(meta, ageDeath != "90+" & ageDeath != "missing or unknown")
any(as.numeric(tmp$ageDeath) >= 90)

meta_new <- meta %>%
  rename(pmi = PMI) %>% 
  mutate(pmi = str_to_sentence(pmi),
         ageDeath = str_to_sentence(ageDeath),
         isHispanic = case_when(isHispanic == "TRUE" ~ "True",
                                isHispanic == "FALSE" ~ "False",
                                TRUE ~ str_to_sentence(isHispanic)),
         race = case_when(race == "Black" ~ "Black or African American",
                          race == "missing or unknown" ~ "Missing or unknown",
                          TRUE ~ race))

print_unique_column_vals(meta_new, c("ageDeath", "pmi"))
print_columns_with_nas(meta_new)

write.csv(meta_new, file.path("data", "output", 
                              str_replace(meta_file$name, ".csv", "_harmonized.csv")))

