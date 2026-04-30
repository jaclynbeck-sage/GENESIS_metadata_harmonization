# Harmonize McCarroll SCZ / GEN-A16 metadata
#
# Modifies the McCarroll SCZ metadata file to conform to the GENESIS data
# dictionary. Data downloaded from the NeMo Archive on Oct 08, 2025 from:
# https://data.nemoarchive.org/other/grant/broad/mccarroll/transcriptome/sncell/10x_v3.1/human/processed/other/SZvillage_donorMetadata.txt
#
# Modifications needed for version Oct 08, 2025:
# * Rename columns:
#   * `Donor` => `individualID`
#   * `Sex` => `sex`
#   * `Age` => `ageDeath`
# * Censor ages above 90
# * Change `sex` values to all lower case
# * Add `dataContributionGroup` = "Broad Institute"
# * Add `cohort` = "HBTRC"
# * Harmonize `Schizophrenia` diagnosis into binary `SCZ` and `Control` columns
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
harmonize_McCarroll_SCZ <- function(metadata, spec) {
  metadata |>
    dplyr::rename(
      individualID = Donor,
      sex = Sex,
      ageDeath = Age
    ) |>
    mutate(
      ageDeath = censor_ages(ageDeath, spec),
      sex = tolower(sex),
      dataContributionGroup = spec$dataContributionGroup$broad,
      cohort = spec$cohort$harvard,
      # There are no missing values in the Schizophrenia field
      SCZ = make_binary_column(Schizophrenia, "Affected", spec),
      Control = make_binary_column(Schizophrenia, "Unaffected", spec)
    ) |>
    select(-Schizophrenia)
}
