# This file contains utility functions like writing to files and various Synapse
# operations.

library(synapser)
library(stringr)
synLogin()

# Write a metadata data frame to a file ----------------------------------------

# This is a wrapper around `write.csv` that has some extra handling for values
# that contain commas and end of line characters (\n). Values with commas
# need to be escaped with quotes for a CSV file, and some data sets have them
# escaped already and some don't. \n characters are removed entirely, as quote
# escaping doesn't affect them.
#
# Arguments:
#   metadata - a data frame of harmonized metadata where rows are individuals
#              and columns are variables
#   filename - the full path and name of the file to be written. This function
#              automatically inserts "_harmonized" just before ".csv" in the
#              file name.
#
# Returns:
#   the new file name that was written, which should have "_harmonized" added
#   before ".csv"
write_metadata <- function(metadata, filename) {
  # Put quotes around values with commas
  for (column in colnames(metadata)) {
    if (is.character(metadata[, column])) {
      # Columns that contain commas and aren't already escaped with quotes
      commas <- grepl(",", metadata[, column]) & !grepl("\"", metadata[, column])
      metadata[commas, column] <- paste0("\"", metadata[commas, column], "\"")

      # Remove any "\n" characters
      metadata[, column] <- str_replace_all(metadata[, column], "\n", "")
    }
  }

  # Some files already have "_harmonized" in the filename, so we remove it if it
  # exists to avoid duplicating "_harmonized". Also, always ensure the output
  # file is a csv, regardless of the input format.
  new_filename <- str_replace(filename,
                              "(_harmonized)?\\.(csv|txt)",
                              "_harmonized.csv")

  new_filename <- file.path("data", "output", new_filename)

  write.csv(metadata, new_filename,
    row.names = FALSE, quote = FALSE
  )

  return(new_filename)
}


# Upload a file to Synapse -----------------------------------------------------

# This is a wrapper around `synStore` to shorten code slightly. It uploads a
# file to Synapse but only if the contents are different than what is currently
# on Synapse. This check is done because if the file is the same, synStore will
# erase any version comments in Synapse even if the file itself doesn't get a
# new version number. This function also makes sure that synStore doesn't erase
# any annotations on the file in Synapse.
#
# Arguments:
#   filename - the full path and name of the file to upload
#   folder_id - the Synapse ID of the folder on Synapse where the file should be
#              uploaded.
#
# Returns:
#   a Synapse `File` object containing information about the uploaded file
synapse_upload <- function(filename, folder_id) {
  syn_info <- synapse_get_info(filename, folder_id)

  if (!is.null(syn_info)) {
    md5 <- tools::md5sum(filename)

    # Don't actually update the file.
    if (md5 == syn_info$get("_file_handle")$contentMd5) {
      message(str_glue("\"{basename(filename)}\" matches the file on Synapse ",
                       "and will not be re-uploaded."))
      return(syn_info)
    }
  }

  syn_file <- File(filename, parent = folder_id)
  syn_file <- synStore(syn_file, forceVersion = FALSE, set_annotations = FALSE)
  return(syn_file)
}


# Get annotations for a file on Synapse ----------------------------------------

# This function searches a folder on Synapse for a specific filename to see if
# the file exists. If the file exists, return info on the file without
# downloading it. Otherwise, return NULL.
synapse_get_info <- function(filename, folder_id) {
  id <- synFindEntityId(basename(filename), folder_id)
  if (is.null(id)) {
    return(NULL)
  }
  synGet(id, downloadFile = FALSE)
}


# Download a file from Synapse -------------------------------------------------

# This is a wrapper around `synGet` to shorten code slightly. All downloads go
# into "data/downloads", and if a file with that name already exists, the old
# file is overwritten with the new one to avoid making multiple copies.
#
# Arguments:
#   syn_id - the Synapse ID of the file on Synapse to download
#
# Returns:
#   a Synapse `File` object containing information about the downloaded file
synapse_download <- function(syn_id) {
  synGet(syn_id,
    downloadLocation = file.path("data", "downloads"),
    ifcollision = "overwrite.local"
  )
}


# Check Synapse for new file versions ------------------------------------------

# All Synapse IDs used for this code have the file's version number included for
# reproducibility. This function checks to see if there are newer versions
# available on Synapse than what is specified in the code, and prints a message
# if that's the case.
#
# Arguments:
#   syn_id_list - a named list where each item is a Synapse ID of the format
#         "syn123" or "syn123.5", where the optional number after the decimal is
#         the file version on Synapse. If no version is specified, a warning
#         is printed stating that the latest version of the file on Synapse will
#         be used.
#
# Returns:
#   nothing
check_new_versions <- function(syn_id_list) {
  for (dataset_name in names(syn_id_list)) {
    syn_id <- syn_id_list[[dataset_name]]
    vals <- str_split_1(syn_id, pattern = "\\.")

    if (length(vals) != 2 || is.na(suppressWarnings(as.numeric(vals[2])))) {
      warning(
        str_glue(
          "No valid version specified for '{syn_id}' ({dataset_name}). ",
          "The latest version will be used for harmonization."
        )
      )
    } else {
      id <- vals[1]
      version <- vals[2]
      syn_file <- synGet(id, downloadFile = FALSE)

      if (syn_file$versionNumber != version) {
        warning(
          str_glue(
            "There is a new version of {id} ({dataset_name}): ",
            "{version} => {syn_file$versionNumber}. Version {version} will be ",
            "used for harmonization."
          )
        )
      }
    }
  }
}
