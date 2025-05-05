library(DiAna)
library(dplyr)
library(kableExtra)
FAERS_version <- "24Q4"
# DiAna::setup_DiAna(quarter = FAERS_version)

# Load the data
import("DRUG")
import("INDI")
import("REAC")
import("THER")
import("DEMO")
import("OUTC")

# Selection of de-duplicated reports only -------------------------------------
Demo <- Demo[RB_duplicates_only_susp == FALSE]

Drug <- Drug[primaryid %in% Demo$primaryid]
Reac <- Reac[primaryid %in% Demo$primaryid]
Indi <- Indi[primaryid %in% Demo$primaryid]
Outc <- Outc[primaryid %in% Demo$primaryid]
Ther <- Ther[primaryid %in% Demo$primaryid]

# Selection of reports with at least one interacting drug
pids_inter <- Drug[role_cod %in% "I"]$primaryid

# Count substance repetitions for reports with interacting drugs
Drug_inter <- Drug[role_cod %in% "I"]
Drug_inter_count <- Drug_inter[, .N, by = .(substance)]
Drug_inter_count <- Drug_inter_count[!substance %in% NA]

kable(Drug_inter_count[order(-N)][1:10,], format = "html", caption = "Count of drugs reported as interacting") %>%
  kable_styling("striped", full_width = F) %>%
  column_spec(1, bold = T) %>%
  column_spec(2, color = "red") %>%
  add_header_above(c(" ", "Count")) %>%
  row_spec(0, bold = T, background = "#D9EAD3") %>%
  row_spec(1:10, background = "#F6F6F6")


# Descriptive for the ost reported interacting drugs













