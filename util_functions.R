MISSING_STRING <- "missing or unknown"

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


to_Roman_numerals <- function(num) {
  return(case_match(num,
                    1 ~ "I",
                    2 ~ "II",
                    3 ~ "III",
                    4 ~ "IV",
                    5 ~ "V",
                    6 ~ "VI",
                    .default = MISSING_STRING))
}

get_bScore <- function(Braak) {
  return(case_when(Braak %in% c("None", MISSING_STRING) ~ Braak,
                   Braak %in% c("Stage I", "Stage II") ~ "Stage I-II",
                   Braak %in% c("Stage III", "Stage IV") ~ "Stage III-IV",
                   Braak %in% c("Stage V", "Stage VI") ~ "Stage V-VI",
                   .default = MISSING_STRING))
}

get_amyAny <- function(amyCerad) {
  return(case_when(amyCerad == "None/No AD/C0" ~ "0",
                   amyCerad == MISSING_STRING ~ MISSING_STRING,
                   .default = "1"))
}

get_amyA <- function(amyThal) {
  return(case_when(amyThal %in% c("None", MISSING_STRING) ~ amyThal,
                   amyThal %in% c("Phase 1", "Phase 2") ~ "Thal Phase 1 or 2",
                   amyThal == "Phase 3" ~ "Thal Phase 3",
                   amyThal %in% c("Phase 4", "Phase 5") ~ "Thal Phase 4 or 5",
                   .default = MISSING_STRING))
}

get_apoe4Status <- function(apoeGenotype) {
  return(case_when(str_detect(apoeGenotype, "4") ~ "Yes",
                   apoeGenotype == MISSING_STRING ~ MISSING_STRING,
                   .default = "No"))
}

censor_ages <- function(ages) {
  return(case_when(ages %in% c("90+", "90_or_over", "89+") ~ "90+",
                   ages == "" ~ NA,
                   ages >= 90 ~ "90+",
                   ages < 90 ~ as.character(ages),
                   .default = NA))
}


write_metadata <- function(metadata, filename) {
  write.csv(metadata,
            file.path("data", "output",
                      str_replace(filename, ".csv", "_harmonized.csv")),
            row.names = FALSE, quote = FALSE)
}
