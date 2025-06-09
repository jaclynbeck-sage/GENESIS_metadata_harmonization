library(synapser)
library(stringr)
library(dplyr)
library(purrr)

syn_ids <- list(
  "Diverse_Cohorts" = "syn51757646.20",
  "MSBB" = "syn6101474.9",
  "NPS_AD" = "syn55251012.4",
  "NPS_AD_neuropath" = "syn55251003.1",
  # Access restricted to Sage internal
  "MSBB_corrections" = "syn66511661.2"
)

synLogin()

# MSBB 1.0
meta_file <- synGet(syn_ids[["MSBB"]],
                    downloadLocation = file.path("data", "downloads"),
                    ifcollision = "overwrite.local")
meta_msbb <- read.csv(meta_file$path)

# MSBB corrections
meta_file <- synGet(syn_ids[["MSBB_corrections"]],
                    downloadLocation = file.path("data", "downloads"),
                    ifcollision = "overwrite.local")
meta_corrections <- readxl::read_xlsx(meta_file$path) |>
  as.data.frame()

# Diverse Cohorts
meta_file <- synGet(syn_ids[["Diverse_Cohorts"]],
                    downloadLocation = file.path("data", "downloads"),
                    ifcollision = "overwrite.local")
meta_dc <- read.csv(meta_file$path)

# NPS-AD
meta_file <- synGet(syn_ids[["NPS_AD"]],
                    downloadLocation = file.path("data", "downloads"),
                    ifcollision = "overwrite.local")
meta_nps <- read.csv(meta_file$path)

neuropath_file <- synGet(syn_ids[["NPS_AD_neuropath"]],
                         downloadLocation = file.path("data", "downloads"),
                         ifcollision = "overwrite.local")
path_nps <- read.csv(neuropath_file$path)

meta_nps <- merge(dplyr::select(meta_nps, -CERAD, -Braak),
                  path_nps,
                  by.x = "individualID", by.y = "IndividualID")


# Missing individualID_AMPAD_1.0 in Diverse Cohorts ----------------------------

# Must be done before re-naming MSBB 1.0 IDs or removing columns
has_1.0 <- subset(meta_dc, individualID_AMPAD_1.0 %in% meta_msbb$individualID)
matches_msbb <- subset(meta_dc, str_replace(individualID, "_MSSM", "") %in%
                         str_replace(meta_msbb$individualID, "AMPAD_MSSM_[0]+", "") &
                         cohort == "Mt Sinai Brain Bank")

missing_1.0_id <- setdiff(matches_msbb$individualID, has_1.0$individualID)
missing_1.0_values <- sapply(missing_1.0_id, function(m_id) {
  grep(m_id, meta_msbb$individualID, value = TRUE)
})


# Fix some columns to match the other data sets -----------------------------

# Helper function
to_Braak_stage <- function(num) {
  return(case_match(num,
                    0 ~ "None",
                    1 ~ "Stage I",
                    2 ~ "Stage II",
                    3 ~ "Stage III",
                    4 ~ "Stage IV",
                    5 ~ "Stage V",
                    6 ~ "Stage VI",
                    NA ~ "missing or unknown",
                    .default = as.character(num)
  ))
}

cols_include <- c("individualID", "sex", "race", "isHispanic",
                  "yearsEducation", "apoeGenotype", "ageDeath", "PMI",
                  "CERAD", "CERAD_original", "Braak", "Braak_original", "CDR")

meta_msbb <- meta_msbb |>
  dplyr::rename(PMI = pmi,
                isHispanic = ethnicity) |>
  mutate(
    individualID = str_replace(individualID, "AMPAD_MSSM_[0]+", ""),
    isHispanic = case_match(isHispanic,
                            c("A", "B", "W") ~ "FALSE",
                            "H" ~ "TRUE",
                            "U" ~ "missing or unknown"
    ),
    race = case_match(race,
                      "A" ~ "Asian",
                      "B" ~ "Black or African American",
                      "H" ~ "other",
                      "W" ~ "White",
                      "U" ~ "missing or unknown"
    ),
    CERAD_original = CERAD,
    # Turn CERAD into strings to avoid confusion between coding systems
    CERAD = case_match(
      CERAD,
      1 ~ "None/No AD/C0",
      2 ~ "Frequent/Definite/C3",
      3 ~ "Moderate/Probable/C2",
      4 ~ "Sparse/Possible/C1"
    ),
    Braak_original = Braak,
    Braak = to_Braak_stage(Braak)
  ) |>
  select(any_of(cols_include)) |>
  mutate(across(any_of(cols_include), as.character))

meta_corrections <- meta_corrections |>
  dplyr::rename(sex = SexLabel,
                race = RaceLabel,
                ageDeath = Age,
                apoeGenotype = ApoE,
                PMI = `PMI (min)`,
                Braak = `B&B Alz`,
                CERAD = CERAD_1) |>
  mutate(
    ageDeath = case_when(ageDeath >= 90 ~ "90+",
                         .default = as.character(ageDeath)),
    PMI = as.character(PMI),
    apoeGenotype = str_replace(apoeGenotype, "/", ""),
    sex = tolower(sex),
    race = str_replace(race, " (nonHispanic)", ""),
    isHispanic = case_match(race,
                            "Hispanic" ~ "TRUE",
                            .default = "FALSE"),
    race = case_match(race,
                      "Black" ~ "Black or African American",
                      "Other" ~ "other",
                      "Hispanic" ~ "other",
                      .default = race),
    CERAD_original = CERAD,
    # Turn CERAD into strings to avoid confusion between coding systems
    CERAD = case_match(
      CERAD,
      0 ~ "None/No AD/C0",
      1 ~ "Sparse/Possible/C1",
      2 ~ "Moderate/Probable/C2",
      3 ~ "Frequent/Definite/C3"
    ),
    Braak_original = Braak,
    Braak = to_Braak_stage(Braak)
  ) |>
  select(any_of(cols_include)) |>
  mutate(across(any_of(cols_include), as.character))

meta_dc <- meta_dc |>
  subset(cohort == "Mt Sinai Brain Bank") |>
  dplyr::rename(CERAD = amyCerad) |>
  # PMI is in hours, while every other data set is in minutes
  mutate(individualID = str_replace(individualID, "_MSSM", ""),
         PMI = suppressWarnings(as.numeric(PMI)) * 60) |>
  select(any_of(cols_include)) |>
  mutate(across(any_of(cols_include), as.character))

meta_nps <- meta_nps |>
  dplyr::rename(Braak = BRAAK_AD,
                isHispanic = ethnicity) |>
  mutate(
    individualID = str_replace(individualID, "AMPAD_MSSM_[0]+", ""),
    isHispanic = ifelse(isHispanic == "Hispanic or Latino", "TRUE", "FALSE"),
    CERAD_original = CERAD,
    CERAD = case_match(
      CERAD,
      NA ~ "missing or unknown",
      1 ~ "None/No AD/C0",
      2 ~ "Sparse/Possible/C1",
      3 ~ "Moderate/Probable/C2",
      4 ~ "Frequent/Definite/C3"
    ),
    Braak_original = Braak,
    Braak = to_Braak_stage(Braak)
  ) |>
  select(any_of(cols_include)) |>
  mutate(across(any_of(cols_include), as.character))


# Subset all data frames to overlapping samples only ---------------------------

all_ids <- c(meta_msbb$individualID, meta_corrections$individualID,
             meta_dc$individualID, meta_nps$individualID) |>
  table()

# Any samples that appear in at least 2 data sets, removing ROSMAP samples
overlaps <- names(all_ids)[which(all_ids >= 2)]
overlaps <- overlaps[!grepl("^R[0-9]+", overlaps)]

meta_msbb <- subset(meta_msbb, individualID %in% overlaps)
meta_corrections <- subset(meta_corrections, individualID %in% overlaps) # This should be all samples
meta_dc <- subset(meta_dc, individualID %in% overlaps)
meta_nps <- subset(meta_nps, individualID %in% overlaps)


# Compare overlapping samples --------------------------------------------------

meta_msbb$source <- "MSBB 1.0"
meta_msbb[is.na(meta_msbb)] <- "missing or unknown"

meta_corrections$source <- "MSBB Corrections"
meta_corrections[is.na(meta_corrections)] <- "missing or unknown"

meta_dc$source <- "Diverse Cohorts"
meta_dc[is.na(meta_dc)] <- "missing or unknown"

meta_nps$source <- "NPS-AD"
meta_nps[is.na(meta_nps)] <- "missing or unknown"

all_data <- purrr::list_rbind(list(meta_msbb, meta_corrections,
                                   meta_dc, meta_nps))

# This is still all samples...
dupes <- all_data |>
  select(-source, -CERAD_original, -Braak_original) |>
  distinct() |>
  pull(individualID) |>
  unique()

# helper function
mismatch_to_df <- function(col_name, mm_type, meta_tmp) {
  data.frame(individualID = unique(meta_tmp$individualID),
             column = col_name,
             type = mm_type,
             source = meta_tmp$source,
             value = meta_tmp[, col_name])
}

val_is_missing <- function(vals) {
  is.na(vals) | vals %in% "missing or unknown"
}

# Go sample by sample
mismatches <- lapply(dupes, function(dupe_id) {
  meta_tmp <- subset(all_data, individualID == dupe_id)

  mismatch_cols <- data.frame()

  for (col_name in setdiff(cols_include, c("individualID", "source", "CERAD_original", "Braak_original"))) {
    if (col_name == "CDR") {
      # Only MSBB 1.0 and NPS-AD have a CDR column
      vals <- unique(meta_tmp[meta_tmp$source %in% c("MSBB 1.0", "NPS-AD"), col_name])
    } else if (col_name %in% c("race", "isHispanic")) {
      # Ignore potential mismatches with NPS-AD race/ethnicity, which is not self-report
      vals <- unique(meta_tmp[meta_tmp$source != "NPS-AD", col_name])
    } else {
      vals <- unique(meta_tmp[, col_name])
    }

    if (all(val_is_missing(vals))) {
      next
    }

    if (length(vals) > 1 & any(val_is_missing(vals))) {
      mismatch_cols <- rbind(mismatch_cols,
                             mismatch_to_df(col_name, "missing", meta_tmp))
      meta_tmp <- subset(meta_tmp, !val_is_missing(meta_tmp[, col_name]))

      if (col_name == "CDR") {
        vals <- unique(meta_tmp[meta_tmp$source %in% c("MSBB 1.0", "NPS-AD"), col_name])
      } else {
        vals <- unique(meta_tmp[, col_name])
      }
    }

    if (length(vals) > 1) {
      if (col_name %in% c("ageDeath", "PMI")) {
        num_vals <- suppressWarnings(as.numeric(vals)) |>
          na.omit()

        equivalent <- sapply(num_vals, all.equal, num_vals[1], tolerance = 1e-3)
        if (any(equivalent != TRUE)) {
          mismatch_cols <- rbind(mismatch_cols,
                                 mismatch_to_df(col_name, "conflict", meta_tmp))
        }
      } else {
        mismatch_cols <- rbind(mismatch_cols,
                               mismatch_to_df(col_name, "conflict", meta_tmp))
      }
    }
  }

  return(mismatch_cols)
})

mismatches <- do.call(rbind, mismatches)

mismatches <- mismatches |>
  tidyr::pivot_wider(names_from = "source", values_from = "value")

table(mismatches$column, mismatches$type)


# Make a spreadsheet of missing values and conflicts

resolve_conflicts <- function(df) {
  df$dataset <- ""
  df$previous_value <- ""
  df$corrected_value <- ""
  df$correction_source <- ""

  for (ri in 1:nrow(df)) {
    row <- df[ri, ]
    vals <- unlist(row[1, c("MSBB 1.0", "NPS-AD", "Diverse Cohorts", "MSBB Corrections")])

    # ignore NPS-AD race values
    if (row$column %in% c("race", "isHispanic")) {
      vals <- vals[-2]
    }
    c_source <- names(vals)[which(!val_is_missing(vals))]

    if (row$type == "missing") {
      who <- names(vals)[which(vals %in% "missing or unknown")]
    } else {
      who <- c_source
    }

    if (length(unique(vals[c_source])) > 1) {
      # race and isHispanic only -- ignore NPS-AD values which are not self-reported
      if (row$column %in% c("race", "isHispanic") & "NPS-AD" %in% c_source) {
        c_source <- setdiff(c_source, "NPS-AD")
      }

      # Start dropping sources in order of lowest to highest priority
      if (length(unique(vals[c_source])) > 1 & "MSBB 1.0" %in% c_source) {
        c_source <- setdiff(c_source, "MSBB 1.0")
      }

      if (length(unique(vals[c_source])) > 1 & "Diverse Cohorts" %in% c_source) {
        c_source <- setdiff(c_source, "Diverse Cohorts")
      }

      # We assume NPS-AD and MSBB Corrections should both be correct. If not, we
      # have an unresolvable conflict that requires manual intervention.
      # Note from Laura: Assume MSBB Corrections is ground truth in case of
      # conflict with NPS-AD. Not implemented yet.
      if (length(unique(vals[c_source])) > 1) {
        print(paste0("Unresolvable conflict for ", row$column, ": ", row$individualID))
        print(vals[unique(c(who, c_source))])
        cat("\n")

        tmp_vals <- vals[c_source]
        c_source <- paste(c_source, collapse = "; ")
        vals[c_source] <- paste0("Unresolved conflict: ", paste(tmp_vals, collapse = " vs "))
      }
    }

    who <- setdiff(who, c_source)

    df$dataset[ri] <- paste(who, collapse = "; ")
    df$previous_value[ri] <- paste(unique(vals[who]), collapse = "; ")
    df$corrected_value[ri] <- unique(vals[c_source])
    df$correction_source[ri] <- paste(c_source, collapse = "; ")
  }

  df |>
    select(individualID, column, type, dataset, previous_value,
           corrected_value, correction_source) |>
    as.data.frame()
}

results <- resolve_conflicts(mismatches) |>
  arrange(column, type, individualID)

dc_missing_1.0_id <- data.frame(
  individualID = missing_1.0_id,
  column = "individualID_AMPAD_1.0",
  type = "missing",
  dataset = "Diverse Cohorts",
  previous_value = "NA",
  corrected_value = missing_1.0_values,
  correction_source = "MSBB 1.0"
)

results <- rbind(results, dc_missing_1.0_id)

res_stats <- results |>
  group_by(column, type) |>
  count() |>
  tidyr::pivot_wider(names_from = type, values_from = n, values_fill = 0) |>
  as.data.frame()

write.csv(results, "msbb_mismatch_summary.csv", row.names = FALSE, quote = FALSE)
write.csv(res_stats, "msbb_mismatch_stats.csv", row.names = FALSE, quote = FALSE)
