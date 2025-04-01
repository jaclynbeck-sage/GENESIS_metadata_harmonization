expectedColumns <- c(
  "individualID", "dataContributionGroup", "cohort", "sex", "race",
  "isHispanic", "ageDeath", "pmi", "apoeGenotype", "apoe4Status", "amyCerad",
  "amyAny", "amyThal", "amyA", "Braak", "bScore"
)


# TODO add all expected columns
print_qc <- function(df,
                     ageDeath_col = "ageDeath",
                     isHispanic_col = "isHispanic",
                     pmi_col = "pmi",
                     race_col = "race",
                     sex_col = "sex",
                     apoe_col = "apoeGenotype",
                     apoe4status_col = "apoe4Status",
                     cerad_col = "amyCerad",
                     amyAny_col = "amyAny",
                     thal_col = "amyThal",
                     amyA_col = "amyA",
                     braak_col = "Braak",
                     bscore_col = "bScore") {
  all_cols <- c(
    ageDeath_col, isHispanic_col, pmi_col, race_col, sex_col, apoe_col,
    apoe4status_col, cerad_col, amyAny_col, thal_col, amyA_col, braak_col,
    bscore_col
  )

  missing <- setdiff(all_cols, colnames(df))
  cat("Missing columns:", paste(missing, collapse = ", "), "\n\n")

  cat("NA value check:\n\n")
  for (col_name in all_cols) {
    if (col_name %in% colnames(df)) {
      cat(col_name, "has", sum(is.na(df[, col_name])), "NA values\n")
    }
  }

  cat("\nUnique value check:\n\n")
  for (col_name in setdiff(all_cols, c(ageDeath_col, pmi_col))) {
    if (col_name %in% colnames(df)) {
      cat(col_name, "values:\n")
      tmp <- df |>
        group_by_at(col_name) |>
        count() |>
        ungroup() |>
        mutate_if(is.character, ~ paste0("\"", .x, "\"")) |>
        data.frame()

      print(tmp)
      cat("\n")
    }
  }

  if (ageDeath_col %in% colnames(df)) {
    cat("Age check:\n")
    # "89+" applies to NPS-AD only, all other studies use "90+" or "90_or_over"
    ages_remove <- c("90+", "89+", "90_or_over", "Missing or unknown", "missing or unknown")
    tmp <- subset(df, !(df[, ageDeath_col] %in% ages_remove)) |>
      mutate(ageDeath = suppressWarnings(as.numeric(.data[[ageDeath_col]]))) |>
      subset(ageDeath >= 90)

    cat(ageDeath_col, "has", nrow(tmp), "uncensored ages.\n")
  }
}


validate_values <- function(metadata, spec, verbose = TRUE) {
  # ageDeath should have only NA, numbers, or "90+". No numbers should be above
  # 89.
  ageDeath <- na.omit(metadata$ageDeath) |>
    setdiff(spec$over90) |>
    as.numeric()

  if (any(ageDeath >= 90, na.rm = TRUE)) {
    cat("X  ageDeath has uncensored ages above 90.\n")
  } else if (any(is.na(ageDeath))) {
    # NAs introduced by coercion to numeric, warn about initial value
    nas <- which(is.na(ageDeath))
    tmp <- na.omit(metadata$ageDeath) |>
      setdiff(spec$over90)
    cat(
      "X  ageDeath has invalid age values:",
      paste(tmp[nas], collapse = ", "), "\n"
    )
  } else if (verbose) {
    cat("OK ageDeath\n")
  }

  # PMI should have only NA or numbers
  pmi <- na.omit(metadata$pmi)
  if (length(pmi) > 0 && !is.numeric(pmi)) {
    cat("X  pmi is not numeric\n")
  } else if (verbose) {
    cat("OK pmi\n")
  }

  # Other columns should have only string values that exist in the dictionary
  cols_check <- setdiff(expectedColumns, c("individualID", "ageDeath", "pmi"))

  for (col_name in cols_check) {
    values <- metadata[, col_name]
    expected <- c(spec$missing, unlist(spec[[col_name]]))
    if (!all(values %in% expected)) {
      cat(
        "X ", col_name, "has unexpected values:",
        paste(setdiff(values, expected), collapse = ", "), "\n"
      )
    } else if (verbose) {
      cat("OK", col_name, "\n")
    }
  }

  # Also check agreement between related columns
  if (!identical(metadata$apoe4Status, get_apoe4Status(metadata$apoeGenotype, spec))) {
    cat("X apoe4Status does not match apoeGenotype\n")
  } else if (verbose) {
    cat("OK apoeGenotype vs apoe4Status\n")
  }

  if (!identical(metadata$amyA, get_amyA(metadata$amyThal, spec))) {
    cat("X amyA does not match amyThal\n")
  } else if (verbose) {
    cat("OK amyThal vs amyA\n")
  }

  if (!identical(metadata$amyAny, get_amyAny(metadata$amyCerad, spec))) {
    cat("X amyAny does not match amyCerad\n")
  } else if (verbose) {
    cat("OK amyCerad vs amyAny\n")
  }

  if (!identical(metadata$bScore, get_bScore(metadata$Braak, spec))) {
    cat("X bScore does not match Braak\n")
  } else if (verbose) {
    cat("OK Braak vs bScore\n")
  }
}


# Always defaults to returning the original value if it doesn't meet any of the
# below criteria, so it will fail validation
to_Braak_stage <- function(num, spec) {
  return(case_match(num,
    0 ~ spec$Braak$none,
    1 ~ spec$Braak$stage1,
    2 ~ spec$Braak$stage2,
    3 ~ spec$Braak$stage3,
    4 ~ spec$Braak$stage4,
    5 ~ spec$Braak$stage5,
    6 ~ spec$Braak$stage6,
    NA ~ spec$missing,
    .default = as.character(num)
  ))
}

get_bScore <- function(Braak, spec) {
  return(case_match(Braak,
    spec$Braak$none ~ spec$bScore$none,
    c(spec$Braak$stage1, spec$Braak$stage2) ~ spec$bScore$stage1_2,
    c(spec$Braak$stage3, spec$Braak$stage4) ~ spec$bScore$stage3_4,
    c(spec$Braak$stage5, spec$Braak$stage6) ~ spec$bScore$stage5_6,
    NA ~ spec$missing,
    .default = as.character(Braak)
  ))
}

get_amyAny <- function(amyCerad, spec) {
  return(case_match(amyCerad,
    spec$amyCerad$none ~ spec$amyAny$zero,
    spec$amyCerad$sparse ~ spec$amyAny$one,
    spec$amyCerad$moderate ~ spec$amyAny$one,
    spec$amyCerad$frequent ~ spec$amyAny$one,
    NA ~ spec$missing,
    .default = as.character(amyCerad)
  ))
}

get_amyA <- function(amyThal, spec) {
  return(case_match(amyThal,
    spec$amyThal$none ~ spec$amyA$none,
    c(spec$amyThal$phase1, spec$amyThal$phase2) ~ spec$amyA$phase1_2,
    spec$amyThal$phase3 ~ spec$amyA$phase3,
    c(spec$amyThal$phase4, spec$amyThal$phase5) ~ spec$amyA$phase4_5,
    NA ~ spec$missing,
    .default = as.character(amyThal)
  ))
}

get_apoe4Status <- function(apoeGenotype, spec) {
  return(case_match(apoeGenotype,
    c(
      spec$apoeGenotype$e2e4,
      spec$apoeGenotype$e3e4,
      spec$apoeGenotype$e4e4
    ) ~ spec$apoe4Status$e4yes,
    c(
      spec$apoeGenotype$e2e2,
      spec$apoeGenotype$e2e3,
      spec$apoeGenotype$e3e3
    ) ~ spec$apoe4Status$e4no,
    NA ~ spec$missing,
    .default = as.character(apoeGenotype)
  ))
}


censor_ages <- function(ages, spec) {
  return(case_when(
    ages %in% c("90+", "90_or_over", "89+") ~ spec$over90,
    ages == "" ~ NA,
    ages == spec$missing ~ NA,
    suppressWarnings(as.numeric(ages)) >= 90 ~ spec$over90,
    .default = as.character(ages)
  ))
}


write_metadata <- function(metadata, filename) {
  # Put quotes around values with commas
  for (column in colnames(metadata)) {
    if (is.character(metadata[, column])) {
      # Columns that contain commas and aren't already escaped with quotes
      commas <- grepl(",", metadata[, column]) & !grepl("\"", metadata[, column])
      metadata[commas, column] <- paste0("\"", metadata[commas, column], "\"")
    }
  }

  new_filename <- file.path(
    "data", "output", str_replace(filename, ".csv", "_harmonized.csv")
  )
  write.csv(metadata, new_filename,
    row.names = FALSE, quote = FALSE
  )

  return(new_filename)
}


synapse_upload <- function(filename, folder_id) {
  syn_file <- File(filename, parent = folder_id)
  syn_file <- synStore(syn_file, forceVersion = FALSE)
  return(syn_file)
}


synapse_download <- function(syn_id) {
  synGet(syn_id,
    downloadLocation = file.path("data", "downloads"),
    ifcollision = "overwrite.local"
  )
}


check_new_versions <- function(syn_id_list) {
  for (dataset_name in names(syn_id_list)) {
    syn_id <- syn_id_list[[dataset_name]]
    vals <- str_split_1(syn_id, pattern = "\\.")

    if (length(vals) != 2 || is.na(suppressWarnings(as.numeric(vals[2])))) {
      cat(
        str_glue(
          "WARNING: no valid version specified for '{syn_id}' ({dataset_name}). ",
          "The latest version will be used for harmonization."
        ),
        "\n"
      )
    } else {
      id <- vals[1]
      version <- vals[2]
      syn_file <- synGet(id, downloadFile = FALSE)

      if (syn_file$versionNumber != version) {
        cat(
          str_glue(
            "WARNING: there is a new version of {id} ({dataset_name}): ",
            "{version} => {syn_file$versionNumber}. Version {version} will be ",
            "used for harmonization."
          ),
          "\n"
        )
      }
    }
  }
}


deduplicate_studies <- function(df_list,
                                include_cols = expectedColumns,
                                exclude_cols = c(),
                                verbose = TRUE) {
  # Make sure certain fields in each data frame are of the same type
  df_list <- lapply(df_list, function(df_item) {
    df_item |>
      mutate(
        individualID = as.character(individualID),
        apoeGenotype = as.character(apoeGenotype),
        amyAny = as.character(amyAny),
        # Special case: MSBB/MSSM samples were re-named for Diverse Cohorts, this
        # makes them directly comparable
        original_individualID = individualID,
        individualID = str_replace(individualID, "AMPAD_MSSM_0+", ""),
        individualID = case_when(
          individualID == 29637 & dataContributionGroup == "MSSM" ~ "29637_MSSM",
          individualID == 29582 & dataContributionGroup == "MSSM" ~ "29582_MSSM",
          .default = as.character(individualID)
        ),
      )
  })

  meta_all <- purrr::list_rbind(df_list)
  include_cols <- setdiff(include_cols, exclude_cols)

  # This accounts for overlapping IDs between different studies that don't refer
  # to the same individual
  dupe_ids <- meta_all |>
    select(individualID, dataContributionGroup) |>
    mutate(
      group_id = paste(individualID, dataContributionGroup),
      duplicate = duplicated(group_id)
    ) |>
    subset(duplicate == TRUE) |>
    distinct()

  # For each ID that has duplicate rows, resolve duplicates:
  #   1. For columns where some rows have NA and some have a unique non-NA value,
  #      replace the NA value with that unique non-NA value.
  #   2. For columns where some rows have "missing or unknown" and some have a
  #      unique value other than that, replace "missing or unknown" with the
  #      unique value.
  #   3. For the ageDeath/pmi columns where rows disagree because of precision,
  #      use the most precise value.
  for (row_id in 1:nrow(dupe_ids)) {
    ind_id <- dupe_ids$individualID[row_id]

    # This will be altered to resolve duplication, and will get added back to the
    # meta_all data frame
    meta_tmp <- subset(meta_all, individualID == ind_id &
      dataContributionGroup == dupe_ids$dataContributionGroup[row_id])

    for (col_name in include_cols) {
      unique_vals <- unique(meta_tmp[, col_name])

      if (length(unique_vals) > 1) {
        if (verbose) {
          cat(ind_id, col_name, "[", paste(unique_vals, collapse = ", "), "]\n")
        }

        # Use most precise ageDeath or pmi value
        if (col_name %in% c("ageDeath", "pmi")) {
          n_decimals <- str_replace(as.character(unique_vals), ".*\\.", "") |>
            nchar()
          meta_tmp[, col_name] <- unique_vals[which.max(n_decimals)]
        } else if (col_name == "cohort") {
          # If there is more than one left over value, report it but don't try to
          # resolve duplication. Special case: Mayo data may be labeled as "Mayo
          # Clinic" in the original Mayo metadata or "Banner" in Diverse Cohorts,
          # but we don't need to change this in either metadata file or print it
          # out.
          is_mayo_banner <- identical(
            sort(unique_vals),
            c("Banner", "Mayo Clinic")
          )
          if (!is_mayo_banner) {
            cat(ind_id, col_name, "[", paste(unique_vals, collapse = ", "), "]\n")
          }
        } else {
          # Remove NA and "missing or unknown" values to see what's left
          leftover <- unique_vals[!is.na(unique_vals) & (unique_vals != spec$missing)]

          if (length(leftover) == 1) {
            # One unique left over value
            meta_tmp[, col_name] <- leftover
          } else if (length(leftover) == 0) {
            # Nothing left, use "missing or unknown" if it's there, otherwise set
            # to NA.
            if (any(unique_vals == spec$missing)) {
              meta_tmp[, col_name] <- spec$missing
            } else {
              meta_tmp[, col_name] <- NA
            }
          } else {
            # If there is more than one left over value, report it but don't try
            # to resolve duplication.
            cat(ind_id, col_name, "[", paste(unique_vals, collapse = ", "), "]\n")
          }
        }
      }
    }

    # Replace original rows with de-duplicated data
    rows_replace <- meta_all$individualID == ind_id &
      meta_all$dataContributionGroup == unique(meta_tmp$dataContributionGroup)
    meta_all[rows_replace, ] <- meta_tmp
  }

  # Revert back to original individualIDs
  meta_all <- meta_all |>
    mutate(individualID = original_individualID) |>
    select(-original_individualID)

  return(meta_all)
}
