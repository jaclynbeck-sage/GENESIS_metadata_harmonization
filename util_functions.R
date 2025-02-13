expectedColumns <- c("individualID", "dataContributionGroup", "cohort", "sex",
                     "race", "isHispanic", "ageDeath", "pmi", "apoeGenotype",
                     "apoe4Status", "amyCerad", "amyAny", "amyThal", "amyA",
                     "Braak", "bScore")

print_unique_column_vals <- function(df, numeric_columns) {
  for (cn in colnames(df)) {
    if (!(cn %in% c("individualID", "projid"))) {
      vals <- unique(df[, cn])
      if (cn %in% numeric_columns) {
        vals <- sort(vals)
      }
      cat(cn, "(", paste(vals, collapse = ",  "), ")\n\n")
    }
  }
}

print_columns_with_nas <- function(df) {
  vals <- c()
  for (cn in colnames(df)) {
    if (any(is.na(df[, cn]))) {
      vals <- c(vals, cn)
    }
  }

  cat(paste(vals, collapse = ",  "), "\n")
}

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
  all_cols <- c(ageDeath_col, isHispanic_col, pmi_col, race_col, sex_col,
                apoe_col, apoe4status_col, cerad_col, amyAny_col, thal_col,
                amyA_col, braak_col, bscore_col)

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
      tmp <- df %>%
        group_by_at(col_name) %>%
        count() %>%
        ungroup() %>%
        mutate_if(is.character, ~paste0("\"", .x, "\"")) %>%
        data.frame()

      print(tmp)
      cat("\n")
    }
  }

  cat("Age check:\n")
  # "89+" applies to NPS-AD only, all other studies use "90+" or "90_or_over"
  ages_remove <- c("90+", "89+", "90_or_over", "Missing or unknown", "missing or unknown")
  tmp <- subset(df, !(df[, ageDeath_col] %in% ages_remove)) %>%
    mutate(ageDeath = as.numeric(.data[[ageDeath_col]])) %>%
    subset(ageDeath >= 90)

  cat(ageDeath_col, "has", nrow(tmp), "uncensored ages")
  if (nrow(tmp) > 0) {
    cat(":", tmp$ageDeath, "\n")
  } else {
    cat("\n")
  }
}


to_Braak_stage <- function(num, spec) {
  return(case_match(num,
                    0 ~ spec$Braak$none,
                    1 ~ spec$Braak$stage1,
                    2 ~ spec$Braak$stage2,
                    3 ~ spec$Braak$stage3,
                    4 ~ spec$Braak$stage4,
                    5 ~ spec$Braak$stage5,
                    6 ~ spec$Braak$stage6,
                    .default = spec$missing))
}

get_bScore <- function(Braak, spec) {
  return(case_when(Braak == spec$Braak$none ~ spec$bScore$none,
                   Braak %in% c(spec$Braak$stage1, spec$Braak$stage2) ~ spec$bScore$stage1_2,
                   Braak %in% c(spec$Braak$stage3, spec$Braak$stage4) ~ spec$bScore$stage3_4,
                   Braak %in% c(spec$Braak$stage5, spec$Braak$stage6) ~ spec$bScore$stage5_6,
                   .default = spec$missing))
}

get_amyAny <- function(amyCerad, spec) {
  return(case_when(amyCerad == spec$amyCerad$none ~ spec$amyAny$zero,
                   amyCerad %in% c(spec$amyCerad$sparse,
                                   spec$amyCerad$moderate,
                                   spec$amyCerad$frequent) ~ spec$amyAny$one,
                   .default = spec$missing))
}

get_amyA <- function(amyThal, spec) {
  return(case_when(amyThal == spec$amyThal$none ~ spec$amyA$none,
                   amyThal %in% c(spec$amyThal$phase1,
                                  spec$amyThal$phase2) ~ spec$amyA$phase1_2,
                   amyThal == spec$amyThal$phase3 ~ spec$amyA$phase3,
                   amyThal %in% c(spec$amyThal$phase4,
                                  spec$amyThal$phase5) ~ spec$amyA$phase4_5,
                   .default = spec$missing))
}

get_apoe4Status <- function(apoeGenotype, spec) {
  return(case_when(apoeGenotype %in% c(spec$apoeGenotype$e2e4,
                                       spec$apoeGenotype$e3e4,
                                       spec$apoeGenotype$e4e4) ~ spec$apoe4Status$e4yes,
                   apoeGenotype %in% c(spec$apoeGenotype$e2e2,
                                       spec$apoeGenotype$e2e3,
                                       spec$apoeGenotype$e3e3) ~ spec$apoe4Status$e4no,
                   .default = spec$missing))
}

censor_ages <- function(ages, spec) {
  return(case_when(ages %in% c("90+", "90_or_over", "89+") ~ spec$over90,
                   ages == "" ~ NA,
                   ages >= 90 ~ spec$over90,
                   ages < 90 ~ as.character(ages),
                   .default = NA))
}


write_metadata <- function(metadata, filename) {
  # Put quotes around values with commas
  for (column in colnames(metadata)) {
    if (is.character(metadata[, column])) {
      commas <- grepl(",", metadata[, column])
      metadata[commas, column] <- paste0("\"", metadata[commas, column], "\"")
    }
  }

  new_filename <- file.path("data", "output",
                            str_replace(filename, ".csv", "_harmonized.csv"))
  write.csv(metadata, new_filename,
            row.names = FALSE, quote = FALSE)

  return(new_filename)
}

synapse_upload <- function(filename, folder_id) {
  syn_file <- File(filename, parent = folder_id)
  syn_file <- synStore(syn_file, forceVersion = FALSE)
  return(syn_file$id)
}

synapse_download <- function(syn_id) {
  synGet(syn_id,
         downloadLocation = file.path("data", "downloads"),
         ifcollision = "overwrite.local")
}
