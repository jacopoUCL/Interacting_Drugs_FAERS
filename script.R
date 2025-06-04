library(DiAna)
library(dplyr)
library(kableExtra)
library(knitr)
library(writexl)

FAERS_version <- "24Q4"
# DiAna::setup_DiAna(quarter = FAERS_version)

# Load the data
import("DRUG")
import("INDI")
import("REAC")
import("THER")
import("DEMO")
import("OUTC")
import_ATC()

# Selection of de-duplicated reports only -------------------------------------
Demo <- Demo[RB_duplicates_only_susp == FALSE]

Drug <- Drug[primaryid %in% Demo$primaryid]
Reac <- Reac[primaryid %in% Demo$primaryid]
Indi <- Indi[primaryid %in% Demo$primaryid]
Outc <- Outc[primaryid %in% Demo$primaryid]
Ther <- Ther[primaryid %in% Demo$primaryid]


# Selection of reports with at least one interacting drug ----------------------
pids_inter <- unique(Drug[role_cod == "I"]$primaryid) # 95.947

# Count substance repetitions for reports with interacting drugs ---------------
Drug_inter <- Drug[role_cod == "I"]
Drug_inter_count <- Drug_inter[!is.na(substance), .N, by= .(substance)][order(-N)]
transposed_df <- as.data.frame(t(Drug_inter_count[order(-N)][1:10,]))
colnames(transposed_df) <- paste0("Row", 1:ncol(transposed_df))
transposed_df <- cbind(Variable = rownames(transposed_df), transposed_df)
transposed_df <- transposed_df[,-1]
colnames(transposed_df) <- transposed_df[1,]
transposed_df <- transposed_df[-1,]

# Create styled kable
kable(transposed_df, format = "html", 
      caption = "<b><span style='font-size:18px;'>Count of drugs reported as interacting</span></b>") %>%
  kable_styling("striped", full_width = F) %>%
  row_spec(0, bold = T, background = "#D9EAD3") %>%
  column_spec(1, bold = T, background = "#F6F6F6")

# Descriptive for the most reported interacting drugs --------------------------
descriptive(pids_cases = pids_inter, drug = "warfarin", file_name = "Descriptive_warfarin.xlsx")
descriptive(pids_cases = pids_inter, drug = "tacrolimus", file_name = "Descriptive_tacrolimus.xlsx")
descriptive(pids_cases = pids_inter, drug = "acetylsalicylic acid", file_name = "Descriptive_acetylsalicylic_acid.xlsx")
descriptive(pids_cases = pids_inter, drug = "quetiapine", file_name = "Descriptive_quetiapine.xlsx")
descriptive(pids_cases = pids_inter, drug = "valproic acid", file_name = "Descriptive_valproic_acid.xlsx")

# Bivariate descriptive analysis of the selected reports (???) ----------------

# Convert the excel in a list of dataframes
des_war <- read_excel("results/Descriptive_warfarin.xlsx", col_types = c("text", "numeric", "numeric"))
des_tac <- read_excel("results/Descriptive_tacrolimus.xlsx", col_types = c("text", "numeric", "numeric"))
des_asa <- read_excel("results/Descriptive_acetylsalicylic_acid.xlsx", col_types = c("text", "numeric", "numeric"))
des_que <- read_excel("results/Descriptive_quetiapine.xlsx", col_types = c("text", "numeric", "numeric"))
des_val <- read_excel("results/Descriptive_valproic_acid.xlsx", col_types = c("text", "numeric", "numeric"))


# Can insert c("sex", "age_range", "Outcome", "continent", "role_cod", "Submission",  "Reporter",  "country",  "year",  
#                          "time_to_onset",  "age_in_years",  "wt_in_kgs",  "Reactions",  "Indications",  "Substances")
compute_bivariate_descriptives <- function(df = des_war, vars = c("sex", "age_range", "Outcome", "continent", "role_cod")) {
  # Transpose and reformat dataframe as in your code
  df <- as.data.frame(t(df))
  colnames(df) <- df[1, ]
  df <- df[-1, ]
  
  # Convert all columns to numeric
  for (i in seq_len(ncol(df))) {
    df[, i] <- as.numeric(df[, i])
  }
  
  # Compute N as sum of "Female", "Male", and "Unknown" in first row if those columns exist
  if (all(c("Female", "Male", "Unknown") %in% colnames(df))) {
    df$N[1] <- sum(df$Female[1], df$Male[1], df$Unknown[1], na.rm = TRUE)
  }
  
  # Identify categories by detecting NA in first row
  categories <- c()
  for (i in seq_len(ncol(df))) {
    if (is.na(df[1, i])) {
      categories <- c(categories, colnames(df)[i])
    }
  }
  
  # Keep only specific categories
  categories_yes <- categories[categories %in% vars]
  categories_no <- categories[!categories %in% vars]
  
  # Find category positions
  cat_positions_yes <- match(categories_yes, colnames(df))
  cat_positions_no <- match(categories_no, colnames(df))
  
  # remove non useful columns
  remove_indices <- integer()
  for (no_pos in cat_positions_no) {
    # Find the next yes position after the current no position
    next_yes <- cat_positions_yes[cat_positions_yes > no_pos]
    if (length(next_yes) > 0) {
      end_pos <- min(next_yes) - 1  # Up to the column before the next "yes"
    } else {
      end_pos <- ncol(df)  # If no "yes" after, go till end of df
    }
    remove_indices <- c(remove_indices, no_pos:end_pos)
  }
  remove_indices <- sort(unique(remove_indices))
  df <- df[, -remove_indices]
  
  cat_positions_yes <- match(categories_yes, colnames(df))
  
  # Split dataframe into separate dataframes by category
  war_dfs <- list()
  for (i in seq_along(categories_yes)) {
    start <- cat_positions_yes[i] + 1
    end <- if (i < length(categories_yes)) cat_positions_yes[i + 1] - 1 else ncol(df)
    if (start <= end) {
      war_dfs[[categories_yes[i]]] <- df[, start:end, drop = FALSE]
    }
  }
  
  # Compute bivariate descriptive statistics (cross tabulations)
  bivariate_des <- list()
  for (i in seq_along(war_dfs)) {
    for (j in seq_along(war_dfs)) {
      if (i != j) {
        df1 <- war_dfs[[i]]
        df2 <- war_dfs[[j]]
        tab <- t(as.matrix(df1)) %*% as.matrix(df2)
        tab_with_row_totals <- cbind(tab, Total = rowSums(tab, na.rm = TRUE))
        total_row <- colSums(tab_with_row_totals, na.rm = TRUE)
        tab_complete <- rbind(tab_with_row_totals, Total = total_row)
        name <- paste0(names(war_dfs)[i], "_", names(war_dfs)[j])
        bivariate_des[[name]] <- as.data.frame(tab_complete)
      }
    }
  }
  
  # Add category column and remove rownames for easier export
  bivariate_des_named <- lapply(bivariate_des, function(df) {
    df_out <- cbind(cat = rownames(df), df)
    rownames(df_out) <- NULL
    return(df_out)
  })
  
  return(bivariate_des_named)
}

des_war_biv <- compute_bivariate_descriptives(des_war)
des_tac_biv <- compute_bivariate_descriptives(des_tac)
des_asa_biv <- compute_bivariate_descriptives(des_asa)
des_que_biv <- compute_bivariate_descriptives(des_que)
des_val_biv <- compute_bivariate_descriptives(des_val)

# Export bivariate descriptive statistics to Excel files
write_xlsx(des_war_biv, "results/Bivariate_Descriptive_warfarin.xlsx")
write_xlsx(des_tac_biv, "results/Bivariate_Descriptive_tacrolimus.xlsx")
write_xlsx(des_asa_biv, "results/Bivariate_Descriptive_acetylsalicylic_acid.xlsx")
write_xlsx(des_que_biv, "results/Bivariate_Descriptive_quetiapine.xlsx")
write_xlsx(des_val_biv, "results/Bivariate_Descriptive_valproic_acid.xlsx")



# Ratio for most reported interacting drugs over non interacting reports --------------------------
# warfarin
pids_war_inter <- unique(Drug[role_cod == "I" & substance == "warfarin"]$primaryid)
pids_war_non_inter <- unique(Drug[role_cod != "I" & substance == "warfarin"]$primaryid)
ratio_war_inter <- length(pids_war_inter) / length(pids_war_non_inter)
# tacrolimus
pids_tac_inter <- unique(Drug[role_cod == "I" & substance == "tacrolimus"]$primaryid)
pids_tac_non_inter <- unique(Drug[role_cod != "I" & substance == "tacrolimus"]$primaryid)
ratio_tac_inter <- length(pids_tac_inter) / length(pids_tac_non_inter)
# acetylsalicylic acid
pids_asa_inter <- unique(Drug[role_cod == "I" & substance == "acetylsalicylic acid"]$primaryid)
pids_asa_non_inter <- unique(Drug[role_cod != "I" & substance == "acetylsalicylic acid"]$primaryid)
ratio_asa_inter <- length(pids_asa_inter) / length(pids_asa_non_inter)
# quetiapine
pids_que_inter <- unique(Drug[role_cod == "I" & substance == "quetiapine"]$primaryid)
pids_que_non_inter <- unique(Drug[role_cod != "I" & substance == "quetiapine"]$primaryid)
ratio_que_inter <- length(pids_que_inter) / length(pids_que_non_inter)
# valproic acid
pids_val_inter <- unique(Drug[role_cod == "I" & substance == "valproic acid"]$primaryid)
pids_val_non_inter <- unique(Drug[role_cod != "I" & substance == "valproic acid"]$primaryid)
ratio_val_inter <- length(pids_val_inter) / length(pids_val_non_inter)
# rivaroxaban
pids_riv_inter <- unique(Drug[role_cod == "I" & substance == "rivaroxaban"]$primaryid)
pids_riv_non_inter <- unique(Drug[role_cod != "I" & substance == "rivaroxaban"]$primaryid)
ratio_riv_inter <- length(pids_riv_inter) / length(pids_riv_non_inter)
# paracetamol
pids_par_inter <- unique(Drug[role_cod == "I" & substance == "paracetamol"]$primaryid)
pids_par_non_inter <- unique(Drug[role_cod != "I" & substance == "paracetamol"]$primaryid)
ratio_par_inter <- length(pids_par_inter) / length(pids_par_non_inter)
# ritonavir
pids_rit_inter <- unique(Drug[role_cod == "I" & substance == "ritonavir"]$primaryid)
pids_rit_non_inter <- unique(Drug[role_cod != "I" & substance == "ritonavir"]$primaryid)
ratio_rit_inter <- length(pids_rit_inter) / length(pids_rit_non_inter)
# clozapine
pids_clo_inter <- unique(Drug[role_cod == "I" & substance == "clozapine"]$primaryid)
pids_clo_non_inter <- unique(Drug[role_cod != "I" & substance == "clozapine"]$primaryid)
ratio_clo_inter <- length(pids_clo_inter) / length(pids_clo_non_inter)
# furosemide
pids_fur_inter <- unique(Drug[role_cod == "I" & substance == "furosemide"]$primaryid)
pids_fur_non_inter <- unique(Drug[role_cod != "I" & substance == "furosemide"]$primaryid)
ratio_fur_inter <- length(pids_fur_inter) / length(pids_fur_non_inter)

# summary table for the ratios
ratios_summary <- data.frame(
  Drug = c("Warfarin", "Tacrolimus", "Acetylsalicylic Acid", "Quetiapine", "Valproic Acid", 
           "Rivaroxaban", "Paracetamol", "Ritonavir", "Clozapine", "Furosemide"),
  Ratio = round((c(ratio_war_inter, ratio_tac_inter, ratio_asa_inter, ratio_que_inter, ratio_val_inter)*100),1)
)
ratios_summary <- ratios_summary[order(-ratios_summary$Ratio),]
ratios_summary$Ratio <- paste0(ratios_summary$Ratio, "%")
transposed_ratios <- as.data.frame(t(ratios_summary))
colnames(transposed_ratios) <- transposed_ratios[1,]
transposed_ratios <- transposed_ratios[-1,]

kable(transposed_ratios, format = "html", 
      caption = "<b><span style='font-size:18px;'>Ratio of Interacting Drugs over Non-Interacting Reports</span></b>") %>%
  kable_styling("striped", full_width = F) %>%
  row_spec(0, bold = T, background = "#D9EAD3") %>%
  column_spec(1, bold = T, background = "#F6F6F6")


# ADR -----------------------------------
# Count ADR repetitions for interacting drugs ---------------------
Reac_inter <- Reac[primaryid %in% pids_inter]
Reac_inter_count <- Reac_inter[!is.na(pt), .N, by= .(pt)][order(-N)]
transposed_reac <- as.data.frame(t(Reac_inter_count[order(-N)][1:10,]))
colnames(transposed_reac) <- transposed_reac[1,]
transposed_reac <- transposed_reac[-1,]

kable(transposed_reac, format = "html", 
      caption = "<b><span style='font-size:18px;'>Count of ADRs reported for interacting drugs</span></b>") %>%
  kable_styling("striped", full_width = F) %>%
  row_spec(0, bold = T, background = "#D9EAD3") %>%
  column_spec(1, bold = T, background = "#F6F6F6")

# Most reported ADRs for most reported interacting drugs --------------------------
# warfarin
war_react <- Reac[primaryid %in% pids_war_inter]
war_react_count <- war_react[!is.na(pt), .N, by = .(pt)][order(-N)]
transposed_war_react <- as.data.frame(t(war_react_count[order(-N)][1:10,]))
colnames(transposed_war_react) <- transposed_war_react[1,]
transposed_war_react <- transposed_war_react[-1,]

kable(transposed_war_react, format = "html", 
      caption = "<b><span style='font-size:18px;'>Count of ADRs reported for warfarin</span></b>") %>%
  kable_styling("striped", full_width = F) %>%
  row_spec(0, bold = T, background = "#D9EAD3") %>%
  column_spec(1, bold = T, background = "#F6F6F6")

# tacrolimus
tac_react <- Reac[primaryid %in% pids_tac_inter]
tac_react_count <- tac_react[!is.na(pt), .N, by = .(pt)][order(-N)]
transposed_tac_react <- as.data.frame(t(tac_react_count[order(-N)][1:10,]))
colnames(transposed_tac_react) <- transposed_tac_react[1,]
transposed_tac_react <- transposed_tac_react[-1,]

kable(transposed_tac_react, format = "html", 
      caption = "<b><span style='font-size:18px;'>Count of ADRs reported for tacrolimus</span></b>") %>%
  kable_styling("striped", full_width = F) %>%
  row_spec(0, bold = T, background = "#D9EAD3") %>%
  column_spec(1, bold = T, background = "#F6F6F6")

# acetylsalicylic acid
asa_react <- Reac[primaryid %in% pids_asa_inter]
asa_react_count <- asa_react[!is.na(pt), .N, by = .(pt)][order(-N)]
transposed_asa_react <- as.data.frame(t(asa_react_count[order(-N)][1:10,]))
colnames(transposed_asa_react) <- transposed_asa_react[1,]
transposed_asa_react <- transposed_asa_react[-1,]

kable(transposed_asa_react, format = "html", 
      caption = "<b><span style='font-size:18px;'>Count of ADRs reported for acetylsalicylic acid</span></b>") %>%
  kable_styling("striped", full_width = F) %>%
  row_spec(0, bold = T, background = "#D9EAD3") %>%
  column_spec(1, bold = T, background = "#F6F6F6")

# quetiapine
que_react <- Reac[primaryid %in% pids_que_inter]
que_react_count <- que_react[!is.na(pt), .N, by = .(pt)][order(-N)]
transposed_que_react <- as.data.frame(t(que_react_count[order(-N)][1:10,]))
colnames(transposed_que_react) <- transposed_que_react[1,]
transposed_que_react <- transposed_que_react[-1,]

kable(transposed_que_react, format = "html", 
      caption = "<b><span style='font-size:18px;'>Count of ADRs reported for quetiapine</span></b>") %>%
  kable_styling("striped", full_width = F) %>%
  row_spec(0, bold = T, background = "#D9EAD3") %>%
  column_spec(1, bold = T, background = "#F6F6F6")

# valproic acid
val_react <- Reac[primaryid %in% pids_val_inter]
val_react_count <- val_react[!is.na(pt), .N, by = .(pt)][order(-N)]
transposed_val_react <- as.data.frame(t(val_react_count[order(-N)][1:10,]))
colnames(transposed_val_react) <- transposed_val_react[1,]
transposed_val_react <- transposed_val_react[-1,]

kable(transposed_val_react, format = "html", 
      caption = "<b><span style='font-size:18px;'>Count of ADRs reported for valproic acid</span></b>") %>%
  kable_styling("striped", full_width = F) %>%
  row_spec(0, bold = T, background = "#D9EAD3") %>%
  column_spec(1, bold = T, background = "#F6F6F6")


# Non interacting drugs for most reported interacting drugs --------------
# warfarin
non_war_all_roles <- Drug[substance != "warfarin" & primaryid %in% pids_war_inter]
non_war_all_roles_count <- non_war_all_roles[!is.na(substance) & !is.na(role_cod), .N, by= .(substance, role_cod)][order(-N)]
transposed_non_war_all_roles <- as.data.frame(t(non_war_all_roles_count[order(-N)][1:10,]))
colnames(transposed_non_war_all_roles) <- transposed_non_war_all_roles[1,]
transposed_non_war_all_roles <- transposed_non_war_all_roles[-1,]

kable(transposed_non_war_all_roles, format = "html", 
      caption = "<b><span style='font-size:18px;'>Count of roles of Drugs for Warfarin</span></b>") %>%
  kable_styling("striped", full_width = F) %>%
  row_spec(0, bold = T, background = "#D9EAD3") %>%
  column_spec(1, bold = T, background = "#F6F6F6")

# tacrolimus
non_tac_all_roles <- Drug[substance != "tacrolimus" & primaryid %in% pids_tac_inter]
non_tac_all_roles_count <- non_tac_all_roles[!is.na(substance) & !is.na(role_cod), .N, by= .(substance, role_cod)][order(-N)]
transposed_non_tac_all_roles <- as.data.frame(t(non_tac_all_roles_count[order(-N)][1:10,]))
colnames(transposed_non_tac_all_roles) <- transposed_non_tac_all_roles[1,]
transposed_non_tac_all_roles <- transposed_non_tac_all_roles[-1,]

kable(transposed_non_tac_all_roles, format = "html", 
      caption = "<b><span style='font-size:18px;'>Count of roles of Drugs for Tacrolimus</span></b>") %>%
  kable_styling("striped", full_width = F) %>%
  row_spec(0, bold = T, background = "#D9EAD3") %>%
  column_spec(1, bold = T, background = "#F6F6F6")

# acetylsalicylic acid
non_asa_all_roles <- Drug[substance != "acetylsalicylic acid" & primaryid %in% pids_asa_inter]
non_asa_all_roles_count <- non_asa_all_roles[!is.na(substance) & !is.na(role_cod), .N, by= .(substance, role_cod)][order(-N)]
transposed_non_asa_all_roles <- as.data.frame(t(non_asa_all_roles_count[order(-N)][1:10,]))
colnames(transposed_non_asa_all_roles) <- transposed_non_asa_all_roles[1,]
transposed_non_asa_all_roles <- transposed_non_asa_all_roles[-1,]

kable(transposed_non_asa_all_roles, format = "html", 
      caption = "<b><span style='font-size:18px;'>Count of roles of Drugs for Acetylsalicylic Acid</span></b>") %>%
  kable_styling("striped", full_width = F) %>%
  row_spec(0, bold = T, background = "#D9EAD3") %>%
  column_spec(1, bold = T, background = "#F6F6F6")

# quetiapine
non_que_all_roles <- Drug[substance != "quetiapine" & primaryid %in% pids_que_inter]
non_que_all_roles_count <- non_que_all_roles[!is.na(substance) & !is.na(role_cod), .N, by= .(substance, role_cod)][order(-N)]
transposed_non_que_all_roles <- as.data.frame(t(non_que_all_roles_count[order(-N)][1:10,]))
colnames(transposed_non_que_all_roles) <- transposed_non_que_all_roles[1,]
transposed_non_que_all_roles <- transposed_non_que_all_roles[-1,]

kable(transposed_non_que_all_roles, format = "html", 
      caption = "<b><span style='font-size:18px;'>Count of roles of Drugs for Quetiapine</span></b>") %>%
  kable_styling("striped", full_width = F) %>%
  row_spec(0, bold = T, background = "#D9EAD3") %>%
  column_spec(1, bold = T, background = "#F6F6F6")

# valproic acid
non_val_all_roles <- Drug[substance != "valproic acid" & primaryid %in% pids_val_inter]
non_val_all_roles_count <- non_val_all_roles[!is.na(substance) & !is.na(role_cod), .N, by= .(substance, role_cod)][order(-N)]
transposed_non_val_all_roles <- as.data.frame(t(non_val_all_roles_count[order(-N)][1:10,]))
colnames(transposed_non_val_all_roles) <- transposed_non_val_all_roles[1,]
transposed_non_val_all_roles <- transposed_non_val_all_roles[-1,]

kable(transposed_non_val_all_roles, format = "html", 
      caption = "<b><span style='font-size:18px;'>Count of roles of Drugs for Valproic Acid</span></b>") %>%
  kable_styling("striped", full_width = F) %>%
  row_spec(0, bold = T, background = "#D9EAD3") %>%
  column_spec(1, bold = T, background = "#F6F6F6")


# Pltoot for interacting drugs by classes
atc_inter <- ATC[substance %in% Drug_inter_count$substance]
atc_inter_drug <- left_join(Drug_inter_count, atc_inter,  by = "substance")


library(ggplot2)
# Class 1
atc_inter_drug_c1 <- atc_inter_drug %>%
  group_by(Class1) %>%
  summarise(N = sum(N, na.rm = TRUE)) %>%
  ungroup() %>%
  as.data.frame() %>%
  filter(!is.na(Class1)) %>%
  arrange(desc(N))

ggplot(atc_inter_drug_c1, aes(x = reorder(Class1, -N), y = N)) +
  geom_bar(stat = "identity", fill = grDevices::rainbow(nrow(atc_inter_drug_c1))) +
  labs(title = "Count of Interacting Drugs by ATC Class 1",
       x = "ATC Class",
       y = "Count of Interacting Drugs") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# Class 2
# Nervous
cl2_ner <- unique(ATC[ATC$Class1 %in% atc_inter_drug_c1$Class1[1]]$Class2)

atc_inter_drug_c2_ner <- atc_inter_drug %>%
  group_by(Class2) %>%
  summarise(N = sum(N, na.rm = TRUE)) %>%
  ungroup() %>%
  as.data.frame() %>%
  filter(!is.na(Class2) & Class2 %in% cl2_ner) %>%
  arrange(desc(N))

ggplot(atc_inter_drug_c2_ner, aes(x = reorder(Class2, -N), y = N)) +
  geom_bar(stat = "identity", fill = grDevices::rainbow(9)[1]) +
  labs(title = "Count of Interacting Drugs by ATC Class 2 (Nervous System)",
       x = "ATC Class",
       y = "Count of Interacting Drugs") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# Cardiovascular
cl2_car <- unique(ATC[ATC$Class1 %in% atc_inter_drug_c1$Class1[2]]$Class2)

atc_inter_drug_c2_car <- atc_inter_drug %>%
  group_by(Class2) %>%
  summarise(N = sum(N, na.rm = TRUE)) %>%
  ungroup() %>%
  as.data.frame() %>%
  filter(!is.na(Class2) & Class2 %in% cl2_car) %>%
  arrange(desc(N))

ggplot(atc_inter_drug_c2_car, aes(x = reorder(Class2, -N), y = N)) +
  geom_bar(stat = "identity", fill = grDevices::rainbow(9)[2]) +
  labs(title = "Count of Interacting Drugs by ATC Class 2 (Cardiovascular System)",
       x = "ATC Class",
       y = "Count of Interacting Drugs") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Antiinfectives
cl2_ant <- unique(ATC[ATC$Class1 %in% atc_inter_drug_c1$Class1[3]]$Class2)

atc_inter_drug_c2_ant <- atc_inter_drug %>%
  group_by(Class2) %>%
  summarise(N = sum(N, na.rm = TRUE)) %>%
  ungroup() %>%
  as.data.frame() %>%
  filter(!is.na(Class2) & Class2 %in% cl2_ant) %>%
  arrange(desc(N))

ggplot(atc_inter_drug_c2_ant, aes(x = reorder(Class2, -N), y = N)) +
  geom_bar(stat = "identity", fill = grDevices::rainbow(15)[3]) +
  labs(title = "Count of Interacting Drugs by ATC Class 2 (Antiinfectives)",
       x = "ATC Class",
       y = "Count of Interacting Drugs") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# and so on ...


