# Harmonize McCarroll HD / GEN-A17 metadata
#
# Modifies the McCarroll HD metadata file to conform to the GENESIS data
# dictionary. Data downloaded from the NeMo Archive on Oct 08, 2025 from:
# https://data.nemoarchive.org/other/grant/broad_mccarroll_HD_2024/mccarroll/transcriptome/sncell/10x_v3.1/human/processed/other/donor_metadata.txt
#
# Modifications needed for version Oct 08, 2025:
# * Rename columns:
#   * `SID` => `individualID`
#   * `Sex` => `sex`
#   * `Age` => `ageDeath`
# * Change `ageDeath` values of ">89" to "90+"
# * Convert `PMI` to numeric values and replace "n/a" with NA
# * Change `sex` values to all lower case
# * Add `dataContributionGroup` = "Broad Institute"
# * Add `cohort` = "HBTRC"
# * Harmonize `Status` column into binary `HD` and `Control` columns
#
# Arguments:
#   metadata - a `data.frame` of metadata from the source metadata file. Columns
#     are variables and rows are individuals.
#   spec - a `config` object describing the standardized values for each field,
#     as defined by this project's `GENESIS_harmonization.yml` file
#
# Returns:
#   a `data.frame` with all relevant fields harmonized to the GENESIS data
#   dictionary. Columns not defined in the data dictionary are left as-is.
#
harmonize_McCarroll_HD <- function(metadata, spec) {
  metadata |>
    dplyr::rename(
      individualID = SID,
      sex = Sex,
      ageDeath = Age
    ) |>
    mutate(
      ageDeath = ifelse(ageDeath == ">89", spec$ageDeath$over90, ageDeath),
      PMI = ifelse(PMI == "n/a", NA, suppressWarnings(as.numeric(PMI))),
      sex = tolower(sex),
      dataContributionGroup = spec$dataContributionGroup$broad,
      cohort = spec$cohort$harvard,
      # Binary diagnosis columns
      HD = make_binary_column(Status, "Case", spec),
      Control = make_binary_column(Status, "Control", spec)
    )
}
