# Authors: Jacopo Palombarini, Angela Boccia
# Last update: 9/12/2025

# Setup -----
# Packages  -----
library(DiAna)
library(dplyr)
library(kableExtra)
library(knitr)
library(writexl)
library(data.table)
library(ggplot2)
library(scales)
library(patchwork)
library(openxlsx)
library(readxl)
library(dplyr)
library(tidyr)
library(tibble)
library(RColorBrewer)

# DiAna setup -----
FAERS_version <- "24Q4"
# DiAna::setup_DiAna(quarter = FAERS_version)

# Load the data -----
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


#####

# Descriptive analyses -----
# Selection of reports with at least one interacting drug ----------------------
# C = Concomitant - PS = Principal suspect - SS = Secondary suspect - I = Interacting
pids_inter <- unique(Drug[role_cod == "I"]$primaryid) # 95.947
pids_non_inter <- unique(Drug[role_cod != "I"]$primaryid) # 14.739.289
subs_inter <- unique(Drug[role_cod != "I"]$substance) # 6317

# Unique combinations by primaryid and role_cod for all data
write.xlsx(
  Drug[, .SD[!duplicated(Drug, by = c("primaryid", "role_cod"))]]
  [, .N, by = role_cod]
  [, perc := round(N / sum(N) * 100, 2)]
  [order(-N)],
  overwrite = TRUE,
  file = "results/roles_count.xlsx"
)

# Unique combinations by primaryid and role_cod for pids in pids_inter
write.xlsx(
  Drug[primaryid %in% pids_inter][
    , .SD[!duplicated(.SD, by = c("primaryid", "role_cod"))]]
  [, .N, by = role_cod]
  [, perc := round(N / sum(N) * 100, 2)]
  [order(-N)],
  overwrite = TRUE,
  file = "results/roles_count2.xlsx"
)

# Count number of drugs per report within SS or I, with at least one I -----
drug_counts_all <- Drug[primaryid %in% pids_inter & role_cod %in% c("SS", "I"), .N, by = primaryid]
setnames(drug_counts_all, "N", "total_drugs")

# Table of reports by total_drugs for 95% cutoff
tab <- drug_counts_all[, .N, by = total_drugs][order(total_drugs)]
tab[, cumprop := cumsum(N) / sum(N)]
cutoffN <- tab[cumprop >= 0.95, min(total_drugs)]

# Join total_drugs to each drug row and count by role
bars <- Drug[primaryid %in% pids_inter & role_cod %in% c("SS", "I")][
  drug_counts_all, on = "primaryid"
][, .(count = .N), by = .(total_drugs, role_cod)]

# Tidy labels & percentages within each bar
bars[, role_cod := factor(role_cod, levels = c("I", "SS"), labels = c("I", "SS"))]
bars <- bars[, total := sum(count), by = total_drugs][, perc := round(100 * count / total, 1)][]

# Compute y-axis top
labels_df <- bars[, .(count = sum(count)), by = total_drugs]
max_y <- max(labels_df$count)
step_guess <- max(10, signif(max_y/5, 1))
y_top <- ceiling(max_y / step_guess) * step_guess

# Totals per bar (for top labels)
totals <- bars %>%
  dplyr::group_by(total_drugs) %>%
  dplyr::summarise(count_total = sum(count), .groups = "drop") %>%
  dplyr::mutate(perc_total = round(100 * count_total / sum(count_total), 1))
if (is.factor(bars$total_drugs)) {
  totals$total_drugs <- factor(totals$total_drugs, levels = levels(bars$total_drugs))
}

p1 <- ggplot(bars, aes(x = total_drugs, y = count, fill = role_cod)) +
  geom_col(position = "stack") +
  geom_segment(aes(x = cutoffN, xend = cutoffN, y = 0, yend = y_top - 1),
               linetype = "dashed", color = "black", size = 0.7, alpha = 0.6) +
  geom_text(aes(x = cutoffN + 0.2, y = y_top - 50, 
                label = "95% of reports"), 
            hjust = 0, size = 3, color = "black") +
  labs(x = "Number of drugs (SS or I)",
       y = "Number of reports (pids)", fill = "Drug role") +
  scale_x_continuous(breaks = 0:15, limits = c(0, 15),
                     expand = expansion(mult = c(0, 0.02))) +
  scale_y_continuous(breaks = seq(0, y_top, by = step_guess),
                     limits = c(0, y_top + step_guess * 0.25)) +
  scale_fill_manual(values = c("I" = "steelblue", "SS" = "grey70")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  # segment percentages per role
  geom_text(
    aes(label = paste0(perc, "%")),
    position = position_stack(vjust = 0.5),
    size = 3,
    color = "white"
  ) +
  # total percentage per bar (top label)
  geom_text(
    data = totals,
    mapping = aes(x = total_drugs, y = count_total, label = paste0(perc_total, "%")),
    vjust = -0.3,
    size = 3.5,
    fontface = "bold",
    color = "black",
    inherit.aes = FALSE
  ) +
  # Annotation with totals
  annotate("text", x = 10.5, y = y_top*0.75, 
           label = paste0("Within total reports with at least one SS or I: ", 
                          formatC(sum(totals$count_total), format = "d", big.mark = ","),
                          "\n- At least one SS: ", 
                          formatC(sum(bars$count[bars$role_cod == "SS"]), format = "d", big.mark = ","),
                          "\n- At least one I: ", 
                          formatC(sum(bars$count[bars$role_cod == "I"]), format = "d", big.mark = ",")),
           hjust = 0, size = 3.5, color = "black")
           
# Count number of drugs per report within PS, SS, C or I, with at least one I -----
Drug_susp <- copy(Drug)
Drug_susp[role_cod %in% c("PS", "SS"), role_cod := "S"]
Drug_susp[, role_cod := factor(role_cod, levels = c(C = "C", I = "I", S = "S"))]
drug_counts_all_2 <- Drug_susp[primaryid %in% pids_inter & !is.na(role_cod),.N, by = primaryid]
setnames(drug_counts_all_2, "N", "total_drugs")

# table of reports by total_drugs (for 95% cutoff on reports)
tab_2 <- drug_counts_all_2[, .N, by = total_drugs][order(total_drugs)]
tab_2[, cumprop := cumsum(N) / sum(N)]
cutoffN_2 <- tab_2[cumprop >= 0.95, min(total_drugs)]

# join total_drugs to each drug row, then count drugs by total_drugs and role_cod
bars_2 <- Drug_susp[primaryid %in% pids_inter & !is.na(role_cod)][
  drug_counts_all_2, on = "primaryid"
][, .(count = .N), by = .(total_drugs, role_cod)]

# tidy labels & percentages within each bar
bars_2[, role_cod := factor(role_cod, levels = c("S", "I", "C"), labels = c("S", "I", "C"))]
bars_2 <- bars_2[, total := sum(count), by = total_drugs][, perc := round(100 * count / total, 1)][]

# y-axis ceiling based on total number of drugs per bar
labels_df_2 <- bars_2[, .(count = sum(count)), by = total_drugs]
max_y_2 <- max(labels_df_2$count)
# choose a "nice" top (nearest 100 or 10, depending on scale)
step_guess <- max(10, signif(max_y_2/5, 1))
y_top_2 <- ceiling(max_y_2 / step_guess) * step_guess

# totals_2 as before
totals_2 <- bars_2 %>%
  dplyr::group_by(total_drugs) %>%
  dplyr::summarise(count_total = sum(count), .groups = "drop") %>%
  dplyr::mutate(perc_total = round(100 * count_total / sum(count_total), 1))
# if total_drugs is a factor in bars_2, keep the same levels
if (is.factor(bars_2$total_drugs)) {
  totals_2$total_drugs <- factor(totals_2$total_drugs, levels = levels(bars_2$total_drugs))
}

# Frequencies per role (exactly as in your example)
freq_roles <- Drug_susp[primaryid %in% pids_inter][
  , .SD[!duplicated(.SD, by = c("primaryid", "role_cod"))]][
    , .N, by = role_cod][
      , perc := round(N / sum(N) * 100, 2)][
        order(-N)]

# Convenience: pull N and perc into a named list/vector
getN   <- function(role) freq_roles[role_cod == role, N]
getPct <- function(role) freq_roles[role_cod == role, perc]

p2 <- ggplot(bars_2, aes(x = total_drugs, y = count, fill = role_cod)) +
  geom_col(position = "stack") +
  geom_segment(aes(x = cutoffN_2, xend = cutoffN_2, y = 0, yend = y_top_2 - 1),
               linetype = "dashed", color = "black", size = 0.7, alpha = 0.6) +
  geom_text(aes(x = cutoffN_2 + 0.2, y = y_top_2 - 50, 
                label = "95% of reports"), 
            hjust = 0, size = 3, color = "black") +
  labs(x = "Number of drugs (S (PS or SS) or I or C)",
       y = "Number of reports (pids)", fill = "Drug role") +
  scale_x_continuous(breaks = 0:25, limits = c(0, 25),
                     expand = expansion(mult = c(0, 0.02))) +
  scale_y_continuous(breaks = seq(0, y_top_2, by = step_guess),
                     limits = c(0, y_top_2 + step_guess * 0.25)) +
  scale_fill_manual(values = c("I" = "#BBDAF8", "S" = "#90C893", "C" = "#FFC899")) +
  geom_text(aes(label = paste0(perc, "%")),
            position = position_stack(vjust = 0.5),
            size = 3,
            color = "white") +
  geom_text(data = totals_2,
            mapping = aes(x = total_drugs, y = count_total, label = paste0(perc_total, "%")),
            vjust = -0.3,
            size = 3.5,
            fontface = "bold",
            color = "black",
            inherit.aes = FALSE) +
  theme_void() +
  theme(axis.title.y  = element_text(size = 10, face = "bold", angle = 90),
    axis.title.x  = element_text(size = 10, face = "bold"),
    axis.text   = element_text(size = 10, face = "bold"),
    axis.text.x = element_text(size = 10, face = "bold", angle = 45, hjust = 1, vjust = 1),
    axis.text.y = element_text(size = 10, face = "bold", angle = 45, hjust = 1.5),
    plot.title  = element_text(size = 10, face = "bold")) # +
  # annotate("text",
  #   x = 19.5, y = y_top_2 * 0.75,
  #   label = paste0(
  #     "Within total reports with at least one I: ",
  #     "\n- At least one I: ", 
  #     formatC(getN("I"), format = "d", big.mark = ","),
  #     "\n- At least one S: ",
  #     formatC(getN("S"), format = "d", big.mark = ","),
  #     "\n- At least one C: ",
  #     formatC(getN("C"), format = "d", big.mark = ",")), 
  #   hjust = 0, size = 3.5, color = "black")

p1 / p2 + plot_layout(heights = c(1, 1.2)) + 
  plot_annotation(title = "Reports by number of reported drugs", 
                  subtitle = "(1986-2024)")

## Plot only roles trend per numbers of reported drugs
p3 <- ggplot(bars_2, aes(x = total_drugs, y = perc, 
                   group = role_cod, color = role_cod)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 3) +
  scale_x_continuous(breaks = 1:35, limits = c(1, 35),
                     expand = expansion(mult = c(0, 0.02))) +
  scale_y_continuous(breaks = seq(0, 60, by = 5),
                     limits = c(0, 60)) +
  labs(x = "Number of drugs - S (PS or SS) or I or C -",
       y = "% by drug role", color = "Drug role") +
  scale_color_manual(values = c("I" = "#BBDAF8",
                                "S" = "#90C893", 
                                "C" = "#FFC899")) +
  theme_void() +
  theme(axis.title.y  = element_text(size = 10, face = "bold", angle = 90),
        axis.title.x  = element_text(size = 10, face = "bold"),
        axis.text   = element_text(size = 10, face = "bold"),
        axis.text.x = element_text(size = 10, face = "bold", angle = 45, hjust = 1, vjust = 1.5),
        axis.text.y = element_text(size = 10, face = "bold", angle = 45, hjust = 1.5),
        plot.title  = element_text(size = 10, face = "bold"))

p2 / p3 + plot_layout(heights = c(1, 1.2)) + 
  plot_annotation(title = "Reports by number of reported drugs and roles", 
                  subtitle = "(1986-2024, Reports with at least one interacting drug (I))")

ggsave("results/plots/plot_count.tiff", plot = p2 / p3, width = 14, height = 10, units = "in", dpi = 1000, compression = "lzw")

# Descriptive - all data -----
descriptive(pids_non_inter, file_name = "results/Descriptive_all_data_non_inter.xlsx")
descriptive(pids_inter, file_name = "results/Descriptive_all_data_inter.xlsx")
  
# Most reported interacting drugs -----
Drug_inter_count <- Drug[role_cod == "I" & !is.na(substance), .N, by= .(substance)][order(-N)]
Drug_inter_count$perc <- round(Drug_inter_count$N / sum(Drug_inter_count$N), 3) * 100
Drug_inter_count <- Drug_inter_count[1:10, ]
write.xlsx(Drug_inter_count, file = "results/Most_reported_interacting_drugs.xlsx")

# Bivariate descriptive analysis of the selected reports (???) ----------------

# Convert the excel in a list of dataframes
des_war <- read_excel("results/Descriptive_warfarin.xlsx", col_types = c("text", "numeric", "numeric"))
des_tac <- read_excel("results/Descriptive_tacrolimus.xlsx", col_types = c("text", "numeric", "numeric"))
des_asa <- read_excel("results/Descriptive_acetylsalicylic_acid.xlsx", col_types = c("text", "numeric", "numeric"))
des_que <- read_excel("results/Descriptive_quetiapine.xlsx", col_types = c("text", "numeric", "numeric"))
des_val <- read_excel("results/Descriptive_valproic_acid.xlsx", col_types = c("text", "numeric", "numeric"))

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
pids_1_inter <- unique(Drug[role_cod == "I" & substance == as.character(Drug_inter_count$substance[1])]$primaryid)
pids_1_non_inter <- unique(Drug[role_cod != "I" & substance == as.character(Drug_inter_count$substance[1])]$primaryid)
ratio_1_inter <- length(pids_1_inter) / length(pids_1_non_inter)
# tacrolimus
pids_2_inter <- unique(Drug[role_cod == "I" & substance == as.character(Drug_inter_count$substance[2])]$primaryid)
pids_2_non_inter <- unique(Drug[role_cod != "I" & substance == as.character(Drug_inter_count$substance[2])]$primaryid)
ratio_2_inter <- length(pids_2_inter) / length(pids_2_non_inter)
# acetylsalicylic acid
pids_3_inter <- unique(Drug[role_cod == "I" & substance == as.character(Drug_inter_count$substance[3])]$primaryid)
pids_3_non_inter <- unique(Drug[role_cod != "I" & substance == as.character(Drug_inter_count$substance[3])]$primaryid)
ratio_3_inter <- length(pids_3_inter) / length(pids_3_non_inter)
# quetiapine
pids_4_inter <- unique(Drug[role_cod == "I" & substance == as.character(Drug_inter_count$substance[4])]$primaryid)
pids_4_non_inter <- unique(Drug[role_cod != "I" & substance == as.character(Drug_inter_count$substance[4])]$primaryid)
ratio_4_inter <- length(pids_4_inter) / length(pids_4_non_inter)
# valproic acid
pids_5_inter <- unique(Drug[role_cod == "I" & substance == as.character(Drug_inter_count$substance[5])]$primaryid)
pids_5_non_inter <- unique(Drug[role_cod != "I" & substance == as.character(Drug_inter_count$substance[5])]$primaryid)
ratio_5_inter <- length(pids_5_inter) / length(pids_5_non_inter)
# rivaroxaban
pids_6_inter <- unique(Drug[role_cod == "I" & substance == as.character(Drug_inter_count$substance[6])]$primaryid)
pids_6_non_inter <- unique(Drug[role_cod != "I" & substance == as.character(Drug_inter_count$substance[6])]$primaryid)
ratio_6_inter <- length(pids_6_inter) / length(pids_6_non_inter)
# paracetamol
pids_7_inter <- unique(Drug[role_cod == "I" & substance == as.character(Drug_inter_count$substance[7])]$primaryid)
pids_7_non_inter <- unique(Drug[role_cod != "I" & substance == as.character(Drug_inter_count$substance[7])]$primaryid)
ratio_7_inter <- length(pids_7_inter) / length(pids_7_non_inter)
# ritonavir
pids_8_inter <- unique(Drug[role_cod == "I" & substance == as.character(Drug_inter_count$substance[8])]$primaryid)
pids_8_non_inter <- unique(Drug[role_cod != "I" & substance == as.character(Drug_inter_count$substance[8])]$primaryid)
ratio_8_inter <- length(pids_8_inter) / length(pids_8_non_inter)
# clozapine
pids_9_inter <- unique(Drug[role_cod == "I" & substance == as.character(Drug_inter_count$substance[9])]$primaryid)
pids_9_non_inter <- unique(Drug[role_cod != "I" & substance == as.character(Drug_inter_count$substance[9])]$primaryid)
ratio_9_inter <- length(pids_9_inter) / length(pids_9_non_inter)
# furosemide
pids_10_inter <- unique(Drug[role_cod == "I" & substance == as.character(Drug_inter_count$substance[10])]$primaryid)
pids_10_non_inter <- unique(Drug[role_cod != "I" & substance == as.character(Drug_inter_count$substance[10])]$primaryid)
ratio_10_inter <- length(pids_10_inter) / length(pids_10_non_inter)

Drug_ = vector()
for (i in 1:10) {
  Drug_[i] <- as.character(Drug_inter_count$substance[i])
}

# summary table for the ratios
ratios_summary <- data.frame(
  substance = Drug_,
  ratio = round((c(ratio_1_inter, ratio_2_inter, ratio_3_inter, ratio_4_inter, ratio_5_inter,
                   ratio_6_inter, ratio_7_inter, ratio_8_inter, ratio_9_inter, ratio_10_inter)*100),1)
)

most_reported <- merge(Drug_inter_count, ratios_summary, by = "substance")



# Descriptive for the most reported interacting drugs --------------------------
descriptive(pids_cases = pids_1_inter, drug = "warfarin", file_name = "results/Descriptive_warfarin.xlsx", 
            vars = c("sex", "Reporter", "age_range", "Outcome", "continent", "age_in_years", "Reactions", 
                     "Indications", "Substances", "role_cod", "time_to_onset", "year"))
descriptive(pids_cases = pids_2_inter, drug = "tacrolimus", file_name = "results/Descriptive_tacrolimus.xlsx", 
            vars = c("sex", "Reporter", "age_range", "Outcome", "continent", "age_in_years", "Reactions", 
                     "Indications", "Substances", "role_cod", "time_to_onset", "year"))
descriptive(pids_cases = pids_3_inter, drug = "acetylsalicylic acid", file_name = "results/Descriptive_acetylsalicylic_acid.xlsx", 
            vars = c("sex", "Reporter", "age_range", "Outcome", "continent", "age_in_years", "Reactions", 
                     "Indications", "Substances", "role_cod", "time_to_onset", "year"))
descriptive(pids_cases = pids_4_inter, drug = "quetiapine", file_name = "results/Descriptive_quetiapine.xlsx", 
            vars = c("sex", "Reporter", "age_range", "Outcome", "continent", "age_in_years", "Reactions", 
                     "Indications", "Substances", "role_cod", "time_to_onset", "year"))
descriptive(pids_cases = pids_5_inter, drug = "valproic acid", file_name = "results/Descriptive_valproic_acid.xlsx", 
            vars = c("sex", "Reporter", "age_range", "Outcome", "continent", "age_in_years", "Reactions", 
                     "Indications", "Substances", "role_cod", "time_to_onset", "year"))
descriptive(pids_cases = pids_6_inter, drug = "rivaroxaban", file_name = "results/Descriptive_rivaroxaban.xlsx", 
            vars = c("sex", "Reporter", "age_range", "Outcome", "continent", "age_in_years", "Reactions", 
                     "Indications", "Substances", "role_cod", "time_to_onset", "year"))
descriptive(pids_cases = pids_7_inter, drug = "paracetamol", file_name = "results/Descriptive_paracetamol.xlsx", 
            vars = c("sex", "Reporter", "age_range", "Outcome", "continent", "age_in_years", "Reactions", 
                     "Indications", "Substances", "role_cod", "time_to_onset", "year"))
descriptive(pids_cases = pids_8_inter, drug = "ritonavir", file_name = "results/Descriptive_ritonavir.xlsx", 
            vars = c("sex", "Reporter", "age_range", "Outcome", "continent", "age_in_years", "Reactions", 
                     "Indications", "Substances", "role_cod", "time_to_onset", "year"))
descriptive(pids_cases = pids_9_inter, drug = "clozapine", file_name = "results/Descriptive_clozapine.xlsx", 
            vars = c("sex", "Reporter", "age_range", "Outcome", "continent", "age_in_years", "Reactions", 
                     "Indications", "Substances", "role_cod", "time_to_onset", "year"))
descriptive(pids_cases = pids_10_inter, drug = "furosemide", file_name = "results/Descriptive_furosemide_acid.xlsx", 
            vars = c("sex", "Reporter", "age_range", "Outcome", "continent", "age_in_years", "Reactions", 
                     "Indications", "Substances", "role_cod", "time_to_onset", "year"))

# Chisq. test ----------------------
clean_desc <- function(path){
  read_excel(path) %>%
    select(`**Characteristic**`, N_cases) %>%
    mutate(N_cases = as.numeric(gsub(",", "", N_cases))) %>%
    pivot_wider(
      names_from = `**Characteristic**`,
      values_from = N_cases
    ) %>%
    select(-c(Unknown, N))
}

files <- c(
  int = "results/tables/univariate/Descriptive_all_data_inter.xlsx",
  non_int = "results/tables/univariate/Descriptive_all_data_non_inter.xlsx"
)

dfs <- lapply(files, clean_desc)
dfs$int <- dfs$int[,c(1:35, 167:173)]
dfs$non_int <- dfs$non_int[,c(1:35, 283:289)]

# sex 
tab_sex <- rbind(Int   = c(Male = as.numeric(dfs$int$Male[1]),   Female = as.numeric(dfs$int$Female[1])),
           Non_int = c(Male = as.numeric(dfs$non_int$Male[1]), Female = as.numeric(dfs$non_int$Female[1])))
tab_sex <- as.data.frame.matrix(t)
test_sex <- chisq.test(tab_sex)
test_sex$residuals
test_sex$p.value

# age range
tab_age <- rbind(Int = c(`Neonate (<28d)` = as.numeric(dfs$int$`Neonate (<28d)`[1]),
                         `Infant (28d-1y)` = as.numeric(dfs$int$`Infant (28d-1y)`[1]),
                         `Child (2y-11y)` = as.numeric(dfs$int$`Child (2y-11y)`[1]),
                         `Teenager (12y-17y)` = as.numeric(dfs$int$`Teenager (12y-17y)`[1]),
                         `Adult (18y-29y)` = as.numeric(dfs$int$`Adult (18y-29y)`[1]),
                         `Adult (30y-49y)` = as.numeric(dfs$int$`Adult (30y-49y)`[1]),
                         `Adult (50y-64y)` = as.numeric(dfs$int$`Adult (50y-64y)`[1]),
                         `Elderly (65y-74y)` = as.numeric(dfs$int$`Elderly (65y-74y)`[1]),
                         `Elderly (75y-84y)` = as.numeric(dfs$int$`Elderly (75y-84y)`[1]),
                         `Elderly (85y-99y)` = as.numeric(dfs$int$`Elderly (85y-99y)`[1]),
                         `Elderly (>99y)` = as.numeric(dfs$int$`Elderly (>99y)`[1])),
                 Non_int = c(`Neonate (<28d)` = as.numeric(dfs$non_int$`Neonate (<28d)`[1]),
                             `Infant (28d-1y)` = as.numeric(dfs$non_int$`Infant (28d-1y)`[1]),
                             `Child (2y-11y)` = as.numeric(dfs$non_int$`Child (2y-11y)`[1]),
                             `Teenager (12y-17y)` = as.numeric(dfs$non_int$`Teenager (12y-17y)`[1]),
                             `Adult (18y-29y)` = as.numeric(dfs$non_int$`Adult (18y-29y)`[1]),
                             `Adult (30y-49y)` = as.numeric(dfs$non_int$`Adult (30y-49y)`[1]),
                             `Adult (50y-64y)` = as.numeric(dfs$non_int$`Adult (50y-64y)`[1]),
                             `Elderly (65y-74y)` = as.numeric(dfs$non_int$`Elderly (65y-74y)`[1]),
                             `Elderly (75y-84y)` = as.numeric(dfs$non_int$`Elderly (75y-84y)`[1]),
                             `Elderly (85y-99y)` = as.numeric(dfs$non_int$`Elderly (85y-99y)`[1]),
                             `Elderly (>99y)` = as.numeric(dfs$non_int$`Elderly (>99y)`[1])))

tab_age <- as.data.frame.matrix(tab_age)

test_age <- chisq.test(tab_age)
test_age$residuals
test_age$p.value

# Outcome
tab_out <- rbind(Int = c(Death = as.numeric(dfs$int$Death[1]),
                         `Life threatening` = as.numeric(dfs$int$`Life threatening`[1]),
                         Disability = as.numeric(dfs$int$Disability[1]),
                         `Required intervention` = as.numeric(dfs$int$`Required intervention`[1]),
                         Hospitalization = as.numeric(dfs$int$Hospitalization[1]),
                         `Congenital anomaly` = as.numeric(dfs$int$`Congenital anomaly`[1]),
                         `Other serious` = as.numeric(dfs$int$`Other serious`[1]),
                         `Non Serious` = as.numeric(dfs$int$`Non Serious`[1])),
  
                 Non_int = c(Death = as.numeric(dfs$non_int$Death[1]),
                             `Life threatening` = as.numeric(dfs$non_int$`Life threatening`[1]),
                             Disability = as.numeric(dfs$non_int$Disability[1]),
                             `Required intervention` = as.numeric(dfs$non_int$`Required intervention`[1]),
                             Hospitalization = as.numeric(dfs$non_int$Hospitalization[1]),
                             `Congenital anomaly` = as.numeric(dfs$non_int$`Congenital anomaly`[1]),
                             `Other serious` = as.numeric(dfs$non_int$`Other serious`[1]),
                             `Non Serious` = as.numeric(dfs$non_int$`Non Serious`[1])))

tab_out <- as.data.frame.matrix(tab_out)

test_out <- chisq.test(tab_out)
test_out$residuals
test_out$p.value

# Continent
tab_con <- rbind(Int = c(`North America` = as.numeric(dfs$int$`North America`[1]),
                         Europe = as.numeric(dfs$int$Europe[1]),
                         Asia = as.numeric(dfs$int$Asia[1]),
                         `South America` = as.numeric(dfs$int$`South America`[1]),
                         Oceania = as.numeric(dfs$int$Oceania[1]),
                         Africa = as.numeric(dfs$int$Africa[1])),
  
                 Non_int = c(`North America` = as.numeric(dfs$non_int$`North America`[1]),
                             Europe = as.numeric(dfs$non_int$Europe[1]),
                             Asia = as.numeric(dfs$non_int$Asia[1]),
                             `South America` = as.numeric(dfs$non_int$`South America`[1]),
                             Oceania = as.numeric(dfs$non_int$Oceania[1]),
                             Africa = as.numeric(dfs$non_int$Africa[1])))

tab_con <- as.data.frame.matrix(tab_con)

test_con <- chisq.test(tab_con)
test_con$residuals
test_con$p.value

# summary
tab_pvals <- data.frame(Variable = c("Sex", "Age", "Outcome", "Continent"),
                        p_value  = c(test_sex$p.value, test_age$p.value,
                                     test_out$p.value, test_con$p.value))
tab_pvals$signif <- cut(tab_pvals$p_value,
                        breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
                        labels = c("***", "**", "*", "ns"))
# Plots
plot_residuals <- function(resid_mat, title) {
  df_resid <- as.data.frame(resid_mat)
  df_resid$Group <- rownames(df_resid)
  
  df_long <- df_resid |>
    pivot_longer(
      cols      = -Group,
      names_to  = "Category",
      values_to = "Residual"
    )
  
  ggplot(df_long, aes(x = Category, y = Group, fill = Residual)) +
    geom_tile() +
    geom_text(aes(label = round(Residual, 1))) +
    scale_fill_gradient2(
      midpoint = 0,
      low = "#BBDAF8",
      mid = "#FFFFFF",
      high = "#F4BED4",
      limits = c(min(df_long$Residual), max(df_long$Residual))) +
    labs(
      title = title,
      x = NULL,
      y = NULL,
      fill = "Std. residual") +
    theme_void() +
    theme(
      axis.title.y  = element_text(size = 10, face = "bold", angle = 90),
      axis.title.x  = element_text(size = 10, face = "bold"),
      axis.text   = element_text(size = 10, face = "bold"),
      axis.text.x = element_text(size = 10, face = "bold", angle = 45, hjust = 1, vjust = 1),
      axis.text.y = element_text(size = 10, face = "bold", angle = 45, hjust = 1),
      plot.title  = element_text(size = 10, face = "bold"))
  
}
fig_all <- plot_residuals(test_sex$residuals, "Sex") + 
  plot_residuals(test_age$residuals, "Age") +
  plot_residuals(test_out$residuals, "Outcome") +
  plot_residuals(test_con$residuals, "Continent") +
  plot_layout(ncol = 2)
p_chisq <- fig_all
ggsave("results/plots/plot_chisq.tiff", plot = p_chisq, width = 17, height = 13, units = "in", dpi = 1000, compression = "lzw")


# Cramer's V of imbalance
cramers_v <- function(tab) {
  chi       <- suppressWarnings(chisq.test(tab, correct = FALSE))
  N         <- sum(tab)
  k         <- min(nrow(tab), ncol(tab))
  V         <- sqrt(chi$statistic / (N * (k - 1)))
  as.numeric(V)
}
cramers_v_sex      <- cramers_v(tab_sex)
cramers_v_age      <- cramers_v(tab_age)
cramers_v_outcome  <- cramers_v(tab_out)
cramers_v_continent<- cramers_v(tab_con)

tab_effect <- data.frame(Variable  = c("Sex", "Age", "Outcome", "Continent"),
                         Cramers_V = c(cramers_v_sex, cramers_v_age, 
                                       cramers_v_outcome, cramers_v_continent))

final_tab <- cbind(tab_pvals[,-3], CramerV = tab_effect$Cramers_V)
write.xlsx(final_tab, file = "results/chisq_test.xlsx")

# save tables
wb <- createWorkbook()
add_chi_sheet <- function(wb, sheet_name, test_object, table_object) {
  addWorksheet(wb, sheet_name)
  writeData(wb, sheet_name, "Observed counts:", startRow = 1, startCol = 1)
  writeData(wb, sheet_name, as.data.frame.matrix(table_object), startRow = 2, startCol = 1)
  
  writeData(wb, sheet_name, "Expected counts:", startRow = 5 + nrow(table_object), startCol = 1)
  writeData(wb, sheet_name, as.data.frame.matrix(round(test_object$expected, 2)), 
            startRow = 6 + nrow(table_object), startCol = 1)
  
  writeData(wb, sheet_name, "Standardized residuals:", startRow = 9 + 2*nrow(table_object), startCol = 1)
  writeData(wb, sheet_name, as.data.frame.matrix(round(test_object$residuals, 2)), 
            startRow = 10 + 2*nrow(table_object), startCol = 1)
  
  writeData(wb, sheet_name, paste("p-value:", test_object$p.value), 
            startRow = 12 + 3*nrow(table_object), startCol = 1)
}
add_chi_sheet(wb, "Sex", test_sex, t(tab_sex))
add_chi_sheet(wb, "Age", test_age, t(tab_age))
add_chi_sheet(wb, "Outcome", test_out, t(tab_out))
add_chi_sheet(wb, "Continent", test_con, t(tab_con))
saveWorkbook(wb, "results/chi_tables.xlsx", overwrite = TRUE)

# Count ADR repetitions for interacting drugs ---------------------
Reac_inter_count <- Reac[primaryid %in% pids_inter & !is.na(pt), .N, by= .(pt)][order(-N)]
Reac_inter_count <- Reac_inter_count[1:10,]
write.xlsx(Reac_inter_count, file = "results/Most_reported_ADR_for_inter_drugs.xlsx")

# Most reported ADRs for most reported interacting drugs --------------------------
react_count_1 <- as.data.frame(Reac[primaryid %in% pids_1_inter & !is.na(pt), .N, by = .(pt)][order(-N)][1:10,])
react_count_2 <- as.data.frame(Reac[primaryid %in% pids_2_inter & !is.na(pt), .N, by = .(pt)][order(-N)][1:10,])
react_count_3 <- as.data.frame(Reac[primaryid %in% pids_3_inter & !is.na(pt), .N, by = .(pt)][order(-N)][1:10,])
react_count_4 <- as.data.frame(Reac[primaryid %in% pids_4_inter & !is.na(pt), .N, by = .(pt)][order(-N)][1:10,])
react_count_5 <- as.data.frame(Reac[primaryid %in% pids_5_inter & !is.na(pt), .N, by = .(pt)][order(-N)][1:10,])
react_count_6 <- as.data.frame(Reac[primaryid %in% pids_6_inter & !is.na(pt), .N, by = .(pt)][order(-N)][1:10,])
react_count_7 <- as.data.frame(Reac[primaryid %in% pids_7_inter & !is.na(pt), .N, by = .(pt)][order(-N)][1:10,])
react_count_8 <- as.data.frame(Reac[primaryid %in% pids_8_inter & !is.na(pt), .N, by = .(pt)][order(-N)][1:10,])
react_count_9 <- as.data.frame(Reac[primaryid %in% pids_9_inter & !is.na(pt), .N, by = .(pt)][order(-N)][1:10,])
react_count_10 <- as.data.frame(Reac[primaryid %in% pids_10_inter & !is.na(pt), .N, by = .(pt)][order(-N)][1:10,])
react_count <- list(react_count_1, react_count_2, react_count_3, react_count_4, react_count_5, 
                   react_count_6, react_count_7, react_count_8, react_count_9, react_count_10)

for (i in 1:10) { 
  colnames(react_count[[i]]) <- c("substance", "N")
  }

react_count_unique <- bind_rows(react_count)
write.xlsx(react_count_unique, file = "results/Most_reported_ADR_for_most_reported_inter_drugs.xlsx")


# Non interacting drugs for most reported interacting drugs --------------
# warfarin
non_int_1 <- Drug[substance != as.character(Drug_inter_count$substance[1]) & primaryid %in% pids_war_inter
                  ][!is.na(substance) & !is.na(role_cod), .N, by= .(substance, role_cod)
                    ][order(-N)]
non_int_1_C <- non_int_1[role_cod == "C"][1:10, -2]
non_int_1_SS <- non_int_1[role_cod == "SS"][1:10, -2]
non_int_1_PS <- non_int_1[role_cod == "PS"][1:10, -2]
non_int_1_list <- list(non_int_1_C, non_int_1_SS, non_int_1_PS)

non_int_2 <- Drug[substance != as.character(Drug_inter_count$substance[2]) & primaryid %in% pids_war_inter
                  ][!is.na(substance) & !is.na(role_cod), .N, by= .(substance, role_cod)
                    ][order(-N)]
non_int_2_C <- non_int_2[role_cod == "C"][1:10, -2]
non_int_2_SS <- non_int_2[role_cod == "SS"][1:10, -2]
non_int_2_PS <- non_int_2[role_cod == "PS"][1:10, -2]
non_int_2_list <- list(non_int_2_C, non_int_2_SS, non_int_2_PS)

non_int_3 <- Drug[substance != as.character(Drug_inter_count$substance[3]) & primaryid %in% pids_war_inter
                  ][!is.na(substance) & !is.na(role_cod), .N, by= .(substance, role_cod)
                    ][order(-N)]
non_int_3_C <- non_int_3[role_cod == "C"][1:10, -2]
non_int_3_SS <- non_int_3[role_cod == "SS"][1:10, -2]
non_int_3_PS <- non_int_3[role_cod == "PS"][1:10, -2]
non_int_3_list <- list(non_int_3_C, non_int_3_SS, non_int_3_PS)

non_int_4 <- Drug[substance != as.character(Drug_inter_count$substance[4]) & primaryid %in% pids_war_inter
                  ][!is.na(substance) & !is.na(role_cod), .N, by= .(substance, role_cod)
                    ][order(-N)]
non_int_4_C <- non_int_4[role_cod == "C"][1:10, -2]
non_int_4_SS <- non_int_4[role_cod == "SS"][1:10, -2]
non_int_4_PS <- non_int_4[role_cod == "PS"][1:10, -2]
non_int_4_list <- list(non_int_4_C, non_int_4_SS, non_int_4_PS)

non_int_5 <- Drug[substance != as.character(Drug_inter_count$substance[5]) & primaryid %in% pids_war_inter
                  ][!is.na(substance) & !is.na(role_cod), .N, by= .(substance, role_cod)
                    ][order(-N)]
non_int_5_C <- non_int_5[role_cod == "C"][1:10, -2]
non_int_5_SS <- non_int_5[role_cod == "SS"][1:10, -2]
non_int_5_PS <- non_int_5[role_cod == "PS"][1:10, -2]
non_int_5_list <- list(non_int_5_C, non_int_5_SS, non_int_5_PS)

non_int_6 <- Drug[substance != as.character(Drug_inter_count$substance[6]) & primaryid %in% pids_war_inter
                  ][!is.na(substance) & !is.na(role_cod), .N, by= .(substance, role_cod)
                    ][order(-N)]
non_int_6_C <- non_int_6[role_cod == "C"][1:10, -2]
non_int_6_SS <- non_int_6[role_cod == "SS"][1:10, -2]
non_int_6_PS <- non_int_6[role_cod == "PS"][1:10, -2]
non_int_6_list <- list(non_int_6_C, non_int_6_SS, non_int_6_PS)

non_int_7 <- Drug[substance != as.character(Drug_inter_count$substance[7]) & primaryid %in% pids_war_inter
                  ][!is.na(substance) & !is.na(role_cod), .N, by= .(substance, role_cod)
                    ][order(-N)]
non_int_7_C <- non_int_7[role_cod == "C"][1:10, -2]
non_int_7_SS <- non_int_7[role_cod == "SS"][1:10, -2]
non_int_7_PS <- non_int_7[role_cod == "PS"][1:10, -2]
non_int_7_list <- list(non_int_7_C, non_int_7_SS, non_int_7_PS)

non_int_8 <- Drug[substance != as.character(Drug_inter_count$substance[8]) & primaryid %in% pids_war_inter
                  ][!is.na(substance) & !is.na(role_cod), .N, by= .(substance, role_cod)
                    ][order(-N)]
non_int_8_C <- non_int_8[role_cod == "C"][1:10, -2]
non_int_8_SS <- non_int_8[role_cod == "SS"][1:10, -2]
non_int_8_PS <- non_int_8[role_cod == "PS"][1:10, -2]
non_int_8_list <- list(non_int_8_C, non_int_8_SS, non_int_8_PS)

non_int_9 <- Drug[substance != as.character(Drug_inter_count$substance[9]) & primaryid %in% pids_war_inter
                  ][!is.na(substance) & !is.na(role_cod), .N, by= .(substance, role_cod)
                    ][order(-N)]
non_int_9_C <- non_int_9[role_cod == "C"][1:10, -2]
non_int_9_SS <- non_int_9[role_cod == "SS"][1:10, -2]
non_int_9_PS <- non_int_9[role_cod == "PS"][1:10, -2]
non_int_9_list <- list(non_int_9_C, non_int_9_SS, non_int_9_PS)

non_int_10 <- Drug[substance != as.character(Drug_inter_count$substance[10]) & primaryid %in% pids_war_inter
                   ][!is.na(substance) & !is.na(role_cod), .N, by= .(substance, role_cod)
                     ][order(-N)]
non_int_10_C <- non_int_10[role_cod == "C"][1:10, -2]
non_int_10_SS <- non_int_10[role_cod == "SS"][1:10, -2]
non_int_10_PS <- non_int_10[role_cod == "PS"][1:10, -2]
non_int_10_list <- list(non_int_10_C, non_int_10_SS, non_int_10_PS)

non_int_list_all <- list(non_int_1_list, non_int_2_list, non_int_3_list, non_int_4_list, non_int_5_list, 
                         non_int_6_list, non_int_7_list, non_int_8_list, non_int_9_list, non_int_10_list)
# save
wb <- createWorkbook()
drug_names <- as.character(Drug_inter_count$substance[1:length(non_int_list_all)])
roles <- c("C", "SS", "PS")
addWorksheet(wb, "Non_interacting")
row_counter <- 1
for (i in seq_along(non_int_list_all)) {
  # Drug name
  writeData(wb, "Non_interacting", paste("Drug:", drug_names[i]), startRow = row_counter, startCol = 1)
  row_counter <- row_counter + 1
  
  for (j in seq_along(non_int_list_all[[i]])) {
    dt <- non_int_list_all[[i]][[j]]
    if (nrow(dt) == 0) next
    
    # Write role
    writeData(wb, "Non_interacting", paste("Role:", roles[j]), startRow = row_counter, startCol = 1)
    row_counter <- row_counter + 1
    
    # Write data.table
    writeData(wb, "Non_interacting", dt, startRow = row_counter, startCol = 1)
    row_counter <- row_counter + nrow(dt) + 2  # +2 to leave some space
  }
  
  row_counter <- row_counter + 1
}
saveWorkbook(wb, "results/non_interacting_drugs.xlsx", overwrite = TRUE)

# Most reported indications for most reported interacting drugs -----
Drug_u <- unique(Drug[primaryid %in% pids_inter], by = "primaryid")
Indi_u <- unique(Indi[primaryid %in% pids_inter], by = "primaryid")
temp <- merge(Drug_u, Indi_u, by = "primaryid", nomatch = 0L)

indi_war <- temp[role_cod == "I" & substance == as.character(int_data$substance[1])][,c(1,6)][, .N, by= .(indi_pt)][order(-N)][1:3,]
indi_tac <- temp[role_cod == "I" & substance == as.character(int_data$substance[2])][,c(1,6)][, .N, by= .(indi_pt)][order(-N)][1:3,]
indi_ace <- temp[role_cod == "I" & substance == as.character(int_data$substance[3])][,c(1,6)][, .N, by= .(indi_pt)][order(-N)][1:3,]
indi_que <- temp[role_cod == "I" & substance == as.character(int_data$substance[4])][,c(1,6)][, .N, by= .(indi_pt)][order(-N)][1:3,]
indi_val <- temp[role_cod == "I" & substance == as.character(int_data$substance[5])][,c(1,6)][, .N, by= .(indi_pt)][order(-N)][1:3,]
indi_riv <- temp[role_cod == "I" & substance == as.character(int_data$substance[6])][,c(1,6)][, .N, by= .(indi_pt)][order(-N)][1:3,]
indi_par <- temp[role_cod == "I" & substance == as.character(int_data$substance[7])][,c(1,6)][, .N, by= .(indi_pt)][order(-N)][1:3,]
indi_rit <- temp[role_cod == "I" & substance == as.character(int_data$substance[8])][,c(1,6)][, .N, by= .(indi_pt)][order(-N)][1:3,]
indi_clo <- temp[role_cod == "I" & substance == as.character(int_data$substance[9])][,c(1,6)][, .N, by= .(indi_pt)][order(-N)][1:3,]
indi_fur <- temp[role_cod == "I" & substance == as.character(int_data$substance[10])][,c(1,6)][, .N, by= .(indi_pt)][order(-N)][1:3,]

int_indi <- list(indi_war, indi_tac, indi_ace, indi_que, indi_val, 
                 indi_riv, indi_par, indi_rit, indi_clo, indi_fur)
names(int_indi) <- int_data$substance

# save
wb <- createWorkbook()
addWorksheet(wb, "Indications")
row_counter <- 1
drug_names <- names(int_indi)
for (i in seq_along(int_indi)) {
  writeData(wb, "Indications", paste("Drug:", drug_names[i]), startRow = row_counter, startCol = 1)
  row_counter <- row_counter + 1
  
  dt <- int_indi[[i]]
  if (is.null(dt) || nrow(dt) == 0) {
    row_counter <- row_counter + 1
    next
  }
  
  out <- as.data.frame(dt)
  out$indi_pt <- as.character(out$indi_pt)
  out <- out[order(out$N, decreasing = TRUE), , drop = FALSE]
  
  writeData(wb, "Indications", out, startRow = row_counter, startCol = 1, withFilter = TRUE)
  
  setColWidths(wb, "Indications", cols = 1:ncol(out), widths = "auto")
  
  row_counter <- row_counter + nrow(out) + 2
}
saveWorkbook(wb, "results/tables/rankings/Most_reported_indications_for_most_reported_interacting_drugs.xlsx", overwrite = TRUE)

# Plot for interacting drugs by classes -----
# Interacting drugs count per substance
Drug_inter_count <- Drug[role_cod == "I" & !is.na(substance), .N, by= .(substance)][order(-N)]
atc_inter <- ATC[substance %in% Drug_inter_count$substance]
atc_inter_drug <- left_join(Drug_inter_count, atc_inter,  by = "substance")

Drug_total_count <- Drug[!is.na(substance), .N, by= .(substance)]
atc_total <- ATC[substance %in% Drug_total_count$substance]
atc_total_drug <- left_join(Drug_total_count, atc_total,  by = "substance")

atc_inter_drug_c1 <- atc_inter_drug %>%
  group_by(Class1) %>%
  summarise(N = sum(N, na.rm = TRUE), .groups = "drop") %>%
  filter(!is.na(Class1)) %>%
  arrange(desc(N))

top_classes <- atc_inter_drug_c1$Class1[1:8]

drug_counts <- atc_inter_drug %>%
  filter(Class1 %in% top_classes) %>%
  group_by(Class1, substance) %>%
  summarise(N_sub = sum(N, na.rm = TRUE), .groups = "drop")

class_k <- drug_counts %>%
  group_by(Class1) %>%
  summarise(N_class = sum(N_sub), .groups = "drop") %>%
  mutate(
    k = pmax(1L, floor(N_class / 20000) * 3)
  )

drug_with_rank <- drug_counts %>%
  left_join(class_k, by = "Class1") %>%
  group_by(Class1) %>%
  arrange(desc(N_sub), .by_group = TRUE) %>%
  mutate(
    rank   = row_number(),
    is_top = rank <= k
  ) %>%
  ungroup()

top_df <- drug_with_rank %>%
  filter(is_top) %>%
  select(Class1, substance, N_sub)

others_df <- drug_with_rank %>%
  filter(!is_top) %>%
  group_by(Class1) %>%
  summarise(N_sub = sum(N_sub), .groups = "drop") %>%
  mutate(substance = "Others")

others_df <- others_df %>% filter(N_sub > 0)

strata_df <- bind_rows(top_df, others_df) %>%
  mutate(
    Class1   = factor(Class1, levels = top_classes),
    substance = as.factor(substance)
  )

class1_levels <- strata_df %>%
  group_by(substance) %>%
  summarise(N_tot = sum(N_sub), .groups = "drop") %>%
  arrange(N_tot) %>%
  pull(substance)
class1_levels <- c("Others", setdiff(class1_levels, "Others"))

strata_df <- strata_df %>%
  mutate(substance = factor(substance, levels = class1_levels))

pastel_palette <- c(
  "#E3F2FD", "#D0E6FB", "#BBDAF8", "#A6CEF4",
  "#90C2F0", "#7BB7EC", "#66ABE8", "#529FE4",
  "#D7F1FA", "#C3E9F7", "#AFE1F4", "#9BD9F0",
  "#87D1ED", "#73C9E9",
  "#E0F7FA", "#CCF1F5", "#B8EBF0", "#A4E5EB",
  "#90DFE6", "#7CD9E1",
  "#E8F5E9", "#D6ECD7", "#C5E3C6", "#B3DAB5",
  "#A2D1A4", "#90C893", "#7FBF82", "#6DB671",
  "#FFF8E1", "#FFF1CC", "#FFEAB7", "#FFE3A3",
  "#FFDC8E", "#FFD57A",
  "#FFEBD6", "#FFE0C2", "#FFD5AE", "#FFC899",
  "#FFBC85", "#FFAF70",
  "#FCE4EC", "#F8D1E0", "#F4BED4", "#F0ABC8", "#EC98BC",
  "#F3E5F5", "#E5CEF0", "#D7B8EB",
  "#ECEFF1", "#D6DCE0"
)

p0 <- ggplot(strata_df, aes(x = Class1, y = N_sub, fill = substance)) +
  geom_col(color = "black", linewidth = 0.15) +
  geom_text(
    aes(label = substance),
    position = position_stack(vjust = 0.5),
    size = 3,
    fontface = "bold",
    color = "black") +
  guides(fill = "none") +
  labs(
    title = "Interacting Drugs by ATC Class 1",
    x     = "ATC Class",
    y     = "Count of Interacting Drugs") +
  theme_void() +
  theme(
    axis.title.y  = element_text(size = 10, face = "bold", angle = 90),
    axis.title.x  = element_text(size = 10, face = "bold"),
    axis.text   = element_text(size = 10, face = "bold"),
    axis.text.x = element_text(size = 10, face = "bold", angle = 45, hjust = 1, vjust = 1),
    axis.text.y = element_text(size = 10, face = "bold", angle = 45, hjust = 1.5),
    plot.title  = element_text(size = 10, face = "bold")) +
  scale_fill_manual(values = pastel_palette)

ggsave("results/plots/plot0.tiff", plot = last_plot(), width = 12, height = 8, units = "in", dpi = 1000, compression = "lzw")

# Class 2
# Class 2 - Nervous
cl2_ner <- unique(ATC[ATC$Class1 %in% atc_inter_drug_c1$Class1[1]]$Class2)

drug_counts_ner <- atc_inter_drug %>%
  filter(Class2 %in% cl2_ner) %>%
  group_by(Class2, substance) %>%
  summarise(N_sub = sum(N, na.rm = TRUE), .groups = "drop")

class2_k_ner <- drug_counts_ner %>%
  group_by(Class2) %>%
  summarise(N_class = sum(N_sub), .groups = "drop") %>%
  mutate(
    k = pmax(1L, c(7, 1, 1, 6, 1, 11, 12))
  )

drug_with_rank_ner <- drug_counts_ner %>%
  left_join(class2_k_ner, by = "Class2") %>%
  group_by(Class2) %>%
  arrange(desc(N_sub), .by_group = TRUE) %>%
  mutate(
    rank   = row_number(),
    is_top = rank <= k
  ) %>%
  ungroup()

top_df_ner <- drug_with_rank_ner %>%
  filter(is_top) %>%
  select(Class2, substance, N_sub)

others_df_ner <- drug_with_rank_ner %>%
  filter(!is_top) %>%
  group_by(Class2) %>%
  summarise(N_sub = sum(N_sub), .groups = "drop") %>%
  mutate(substance = "Others") %>%
  filter(N_sub > 0)

strata_df_ner <- bind_rows(top_df_ner, others_df_ner) %>%
  left_join(class2_k_ner, by = "Class2") %>%
  mutate(
    Class2 = factor(Class2, levels = unique(Class2[order(-N_class)])),
    substance = as.factor(substance)
  ) %>%
  select(-N_class, -k)

ner_levels <- strata_df_ner %>%
  group_by(substance) %>%
  summarise(N_tot = sum(N_sub), .groups = "drop") %>%
  arrange(-1*N_tot) %>%
  pull(substance)
ner_levels <- c("Others", setdiff(ner_levels, "Others"))

strata_df_ner <- strata_df_ner %>%
  mutate(substance = factor(substance, levels = ner_levels))

p1 <- ggplot(strata_df_ner, aes(x = Class2, y = N_sub, fill = substance)) +
  geom_col(color = "black", linewidth = 0.15) +
  geom_text(
    aes(label = substance),
    position = position_stack(vjust = 0.5),
    size = 3,
    fontface = "bold",
    color = "black") +
  guides(fill = "none") +
  labs(
    title = "Interacting Drugs by ATC Class 2 (Nervous System)",
    x     = "ATC Class 2",
    y     = "Count of Interacting Drugs") +
  theme_void() +
  theme(
    axis.title.y  = element_text(size = 10, face = "bold", angle = 90),
    axis.title.x  = element_text(size = 10, face = "bold"),
    axis.text   = element_text(size = 10, face = "bold"),
    axis.text.x = element_text(size = 10, face = "bold", angle = 45, hjust = 1, vjust = 1),
    axis.text.y = element_text(size = 10, face = "bold", angle = 45, hjust = 1.5),
    plot.title  = element_text(size = 10, face = "bold")) +
  scale_fill_manual(values = pastel_palette)

# Class 2 - Cardiovascular
cl2_car <- unique(ATC[ATC$Class1 %in% atc_inter_drug_c1$Class1[2]]$Class2)

drug_counts_car <- atc_inter_drug %>%
  filter(Class2 %in% cl2_car) %>%
  group_by(Class2, substance) %>%
  summarise(N_sub = sum(N, na.rm = TRUE), .groups = "drop")

class2_k_car <- drug_counts_car %>%
  group_by(Class2) %>%
  summarise(N_class = sum(N_sub), .groups = "drop") %>%
  mutate(
    k = pmax(1L, c(9, 1, 7, 3, 6, 5, 5, NA, NA, NA))
  )

drug_with_rank_car <- drug_counts_car %>%
  left_join(class2_k_car, by = "Class2") %>%
  group_by(Class2) %>%
  arrange(desc(N_sub), .by_group = TRUE) %>%
  mutate(
    rank   = row_number(),
    is_top = rank <= k
  ) %>%
  ungroup()

top_df_car <- drug_with_rank_car %>%
  filter(is_top) %>%
  select(Class2, substance, N_sub)

others_df_car <- drug_with_rank_car %>%
  filter(!is_top) %>%
  group_by(Class2) %>%
  summarise(N_sub = sum(N_sub), .groups = "drop") %>%
  mutate(substance = "Others") %>%
  filter(N_sub > 0)

strata_df_car <- bind_rows(top_df_car, others_df_car) %>%
  left_join(class2_k_car, by = "Class2") %>%
  mutate(
    Class2 = factor(Class2, levels = unique(Class2[order(-N_class)])),
    substance = as.factor(substance)
  ) %>%
  select(-N_class, -k)

car_levels <- strata_df_car %>%
  group_by(substance) %>%
  summarise(N_tot = sum(N_sub), .groups = "drop") %>%
  arrange(-1*N_tot) %>%
  pull(substance)
car_levels <- c("Others", setdiff(car_levels, "Others"))

strata_df_car <- strata_df_car %>%
  mutate(substance = factor(substance, levels = car_levels))

p2 <- ggplot(strata_df_car, aes(x = Class2, y = N_sub, fill = substance)) +
  geom_col(color = "black", linewidth = 0.15) +
  geom_text(
    aes(label = substance),
    position = position_stack(vjust = 0.5),
    size = 3,
    fontface = "bold",
    color = "black") +
  guides(fill = "none") +
  labs(
    title = "Interacting Drugs by ATC Class 2 (Cardiovascular System)",
    x     = "ATC Class 2",
    y     = "Count of Interacting Drugs") +
  theme_void() +
  theme(
    axis.title.y  = element_text(size = 10, face = "bold", angle = 90),
    axis.title.x  = element_text(size = 10, face = "bold"),
    axis.text   = element_text(size = 10, face = "bold"),
    axis.text.x = element_text(size = 10, face = "bold", angle = 45, hjust = 1, vjust = 1),
    axis.text.y = element_text(size = 10, face = "bold", angle = 45, hjust = 1.5),
    plot.title  = element_text(size = 10, face = "bold")) +
  scale_fill_manual(values = pastel_palette)

# Class 2 - Antiinfectives
cl2_ant <- unique(ATC[ATC$Class1 %in% atc_inter_drug_c1$Class1[3]]$Class2)

drug_counts_ant <- atc_inter_drug %>%
  filter(Class2 %in% cl2_ant) %>%
  group_by(Class2, substance) %>%
  summarise(N_sub = sum(N, na.rm = TRUE), .groups = "drop")

class2_k_ant <- drug_counts_ant %>%
  group_by(Class2) %>%
  summarise(N_class = sum(N_sub), .groups = "drop") %>%
  mutate(
    k = pmax(1L, c(12, 2, 3, 12, NA, NA, NA))
  )

drug_with_rank_ant <- drug_counts_ant %>%
  left_join(class2_k_ant, by = "Class2") %>%
  group_by(Class2) %>%
  arrange(desc(N_sub), .by_group = TRUE) %>%
  mutate(
    rank   = row_number(),
    is_top = rank <= k
  ) %>%
  ungroup()

top_df_ant <- drug_with_rank_ant %>%
  filter(is_top) %>%
  select(Class2, substance, N_sub)

others_df_ant <- drug_with_rank_ant %>%
  filter(!is_top) %>%
  group_by(Class2) %>%
  summarise(N_sub = sum(N_sub), .groups = "drop") %>%
  mutate(substance = "Others") %>%
  filter(N_sub > 0)

strata_df_ant <- bind_rows(top_df_ant, others_df_ant) %>%
  left_join(class2_k_ant, by = "Class2") %>%
  mutate(
    Class2 = factor(Class2, levels = unique(Class2[order(-N_class)])),
    substance = as.factor(substance)
  ) %>%
  select(-N_class, -k)

ant_levels <- strata_df_ant %>%
  group_by(substance) %>%
  summarise(N_tot = sum(N_sub), .groups = "drop") %>%
  arrange(-1*N_tot) %>%
  pull(substance)
ant_levels <- c("Others", setdiff(ant_levels, "Others"))

strata_df_ant <- strata_df_ant %>%
  mutate(substance = factor(substance, levels = ant_levels))

p3 <- ggplot(strata_df_ant, aes(x = Class2, y = N_sub, fill = substance)) +
  geom_col(color = "black", linewidth = 0.15) +
  geom_text(
    aes(label = substance),
    position = position_stack(vjust = 0.5),
    size = 3,
    fontface = "bold",
    color = "black") +
  guides(fill = "none") +
  labs(
    title = "Interacting Drugs by ATC Class 2 (Antiinfectives)",
    x     = "ATC Class 2",
    y     = "Count of Interacting Drugs") +
  theme_void() +
  theme(
    axis.title.y  = element_text(size = 10, face = "bold", angle = 90),
    axis.title.x  = element_text(size = 10, face = "bold"),
    axis.text   = element_text(size = 10, face = "bold"),
    axis.text.x = element_text(size = 10, face = "bold", angle = 45, hjust = 1, vjust = 1),
    axis.text.y = element_text(size = 10, face = "bold", angle = 45, hjust = 1.5),
    plot.title  = element_text(size = 10, face = "bold")) +
  scale_fill_manual(values = pastel_palette)

# Class 2 - Immunomodulating
cl2_imm <- unique(ATC[ATC$Class1 %in% atc_inter_drug_c1$Class1[4]]$Class2)

drug_counts_imm <- atc_inter_drug %>%
  filter(Class2 %in% cl2_imm) %>%
  group_by(Class2, substance) %>%
  summarise(N_sub = sum(N, na.rm = TRUE), .groups = "drop")

class2_k_imm <- drug_counts_imm %>%
  group_by(Class2) %>%
  summarise(N_class = sum(N_sub), .groups = "drop") %>%
  mutate(
    k = pmax(1L, c(11, 1, NA, 8))
  )

drug_with_rank_imm <- drug_counts_imm %>%
  left_join(class2_k_imm, by = "Class2") %>%
  group_by(Class2) %>%
  arrange(desc(N_sub), .by_group = TRUE) %>%
  mutate(
    rank   = row_number(),
    is_top = rank <= k
  ) %>%
  ungroup()

top_df_imm <- drug_with_rank_imm %>%
  filter(is_top) %>%
  select(Class2, substance, N_sub)

others_df_imm <- drug_with_rank_imm %>%
  filter(!is_top) %>%
  group_by(Class2) %>%
  summarise(N_sub = sum(N_sub), .groups = "drop") %>%
  mutate(substance = "Others") %>%
  filter(N_sub > 0)

strata_df_imm <- bind_rows(top_df_imm, others_df_imm) %>%
  left_join(class2_k_imm, by = "Class2") %>%
  mutate(
    Class2 = factor(Class2, levels = unique(Class2[order(-N_class)])),
    substance = as.factor(substance)
  ) %>%
  select(-N_class, -k)

imm_levels <- strata_df_imm %>%
  group_by(substance) %>%
  summarise(N_tot = sum(N_sub), .groups = "drop") %>%
  arrange(-1*N_tot) %>%
  pull(substance)
imm_levels <- c("Others", setdiff(imm_levels, "Others"))

strata_df_imm <- strata_df_imm %>%
  mutate(substance = factor(substance, levels = imm_levels))

p4 <- ggplot(strata_df_imm, aes(x = Class2, y = N_sub, fill = substance)) +
  geom_col(color = "black", linewidth = 0.15) +
  geom_text(
    aes(label = substance),
    position = position_stack(vjust = 0.5),
    size = 3,
    fontface = "bold",
    color = "black") +
  guides(fill = "none") +
  labs(
    title = "Interacting Drugs by ATC Class 2 (Immunomodulating)",
    x     = "ATC Class 2",
    y     = "Count of Interacting Drugs") +
  theme_void() +
  theme(
    axis.title.y  = element_text(size = 10, face = "bold", angle = 90),
    axis.title.x  = element_text(size = 10, face = "bold"),
    axis.text   = element_text(size = 10, face = "bold"),
    axis.text.x = element_text(size = 10, face = "bold", angle = 45, hjust = 1, vjust = 1),
    axis.text.y = element_text(size = 10, face = "bold", angle = 45, hjust = 1.5),
    plot.title  = element_text(size = 10, face = "bold")) +
  scale_fill_manual(values = pastel_palette)

ggsave("results/plots/plot5.tiff", plot = (p1 | p2) / (p3 | p4), width = 22, height = 18, units = "in", dpi = 1000, compression = "lzw")

# Plot trend for most reported interacting drugs -----
int_data <- read_excel("results/tables/rankings/Most_reported_interacting_drugs.xlsx")
years_ <- function(path = NA) {
  d <- read_excel(path)
  d <- t(d)
  colnames(d) <- d[1,]
  d <- as.data.frame(d)[-c(1,3), (ncol(d)-15):ncol(d)]
  d <- as.data.frame(t(d))
  colnames(d) <- "N"
  return(d)
}
data_war <- years_("results/tables/univariate/Descriptive_warfarin.xlsx")
data_tac <- years_("results/tables/univariate/Descriptive_tacrolimus.xlsx")
data_ace <- years_("results/tables/univariate/Descriptive_acetylsalicylic_acid.xlsx")
data_que <- years_("results/tables/univariate/Descriptive_quetiapine.xlsx")
data_val <- years_("results/tables/univariate/Descriptive_valproic_acid.xlsx")
data_all_int <- years_("results/tables/univariate/Descriptive_all_data_inter.xlsx")
# data_riv <- years_("results/tables/univariate/Descriptive_rivaroxaban.xlsx")
# data_par <- years_("results/tables/univariate/Descriptive_paracetamol.xlsx")
# data_rit <- years_("results/tables/univariate/Descriptive_ritonavir.xlsx")
# data_clo <- years_("results/tables/univariate/Descriptive_clozapine.xlsx")
# data_fur <- years_("results/tables/univariate/Descriptive_furosemide_acid.xlsx")

df <- cbind(data_war, data_tac, data_ace, data_que, data_val) #, data_riv, data_par, data_rit, data_clo, data_fur
colnames(df) <- c(int_data$substance[1:5])
df <- df %>%
  rownames_to_column(var = "year") %>%
  mutate(
    across(
      -year,
      ~ as.numeric(gsub(",", "", as.character(.x)))
    )
  ) %>% 
  mutate(Avg. = rowSums(across(where(is.numeric))))

df$year <- as.numeric(df$year)
df$Avg. <- df$Avg./5

df_long <- df %>%
  pivot_longer(
    cols = -year,
    names_to  = "drug",
    values_to = "value"
  )

pastel_palette2 <- c("#529FE4", "red", "#9BD9F0", "#6DB671", "#FFAF70", "#EC98BC")

alphas <- c(rep(0.4, 5), 1)
names(alphas) <- unique(df_long$drug)

py <- ggplot(df_long, aes(x = year, y = value, group = drug)) +
  geom_line(aes(color = drug, alpha = drug), size = 1.3) +
  geom_point(
    aes(alpha = drug),
    size = 2.5,
    stroke = 1,
    fill = "grey") +
  geom_point(
    aes(alpha = drug, color = drug),
    size = 1.7,
    stroke = 1,
    fill = "white") +
  theme_minimal(base_size = 14) +
  labs(
    title = "Drug values over years",
    x = "Year",
    y = "Value",
    color = "Drug") +
  scale_x_continuous(breaks = unique(df_long$year)) +
  scale_color_manual(values = pastel_palette2) +
  scale_alpha_manual(values = alphas, guide = "none") +
  theme_void() +
  theme(
    axis.title.y  = element_text(size = 10, face = "bold", angle = 90),
    axis.title.x  = element_text(size = 10, face = "bold"),
    axis.text     = element_text(size = 10, face = "bold"),
    axis.text.x   = element_text(size = 10, face = "bold", angle = 45, hjust = 1, vjust = 1),
    axis.text.y   = element_text(size = 10, face = "bold", angle = 45, hjust = 1.5),
    plot.title    = element_text(size = 10, face = "bold"))

ggsave("results/plots/plot_years.tiff", plot = py, width = 12, height = 8, units = "in", dpi = 1000, compression = "lzw")

#####

# Omega -----
# Most reported couples of drugs -----
# For the most reported interacting drugs -----
pids_most_int <- list(pids_1_inter, pids_2_inter, pids_3_inter, pids_4_inter, pids_5_inter, 
                      pids_6_inter, pids_7_inter, pids_8_inter, pids_9_inter, pids_10_inter)
ggsave("results/plots/plot0.tiff", plot = last_plot(), width = 12, height = 8, units = "in", dpi = 1000, compression = "lzw")roles <- c("PS", "SS", "C")
tpl <- data.frame(PS = I(list(NULL)),
                  SS = I(list(NULL)),
                  C  = I(list(NULL)))
df_list <- replicate(10, tpl, simplify = FALSE)

for (i in seq_along(pids_most_int)) {
  for (j in seq_along(roles)) {
    res <- Drug[primaryid %in% pids_most_int[[i]] & substance != int_data$substance[i] & role_cod == roles[j], .N, by = .(substance)][order(-N)]
    if (nrow(res)) res <- res[seq_len(min(10, nrow(res)))]
    df_list[[i]][[roles[j]]] <- list(res)
  }
}

names(df_list) <- int_data$substance
names(react_count) <- int_data$substance
df_list$warfarin$PS
react_count$warfarin

# Omega application -----
library(komega)
# Select the mrid and then only the reports with the most frequent "number of reported drugs"
most_rep_comb <- function(top_drug = 1, top_reac = 10, top_indi = 10, 
                          top_c = 1, roles = NULL, verbose = FALSE, 
                          merge_reac = FALSE, merge_indi = FALSE) {
  
  stopifnot(exists("Drug"), exists("Indi"), exists("Reac"))
  stopifnot(data.table::is.data.table(Drug), data.table::is.data.table(Indi), data.table::is.data.table(Reac), 
            !is.na(top_reac), !is.na(top_indi), !is.na(top_drug), !is.na(top_c), 
            top_reac >= 0, top_indi >= 0, top_drug >= 0, top_c >= 0, 
            !is.null(top_reac), !is.null(top_indi), !is.null(top_drug), !is.null(top_c))
  
  # Setup dfs
  Drug_u <- Drug[role_cod == "I" & !is.na(substance)]
  
  Reac_u <- Reac[, .N, by = .(pt)][order(-N)]
  reac <- as.character(Reac_u$pt[1:min(10, nrow(Reac_u))])
  
  Indi_u <- Indi[, .N, by = .(indi_pt)][order(-N)]
  indi <- as.character(Indi_u$indi_pt[1:min(10, nrow(Indi_u))])
  
  # Filter to suspect drugs & optional indication / reaction ----------------
  if (merge_reac) {Drug_u <- merge(Drug_u, Reac,  by = "primaryid")[pt == reac[top_reac]]}
  if (merge_indi) {Drug_u <- merge(Drug_u, Indi, by = "primaryid")[indi_pt == indi[top_indi]]}
  
  # Count suspect substances and pick the requested top_drug ----------------------
  Drug_u <- Drug_u[, .N, by = .(substance)][order(-N)]
  if (nrow(Drug_u) == 0L) {
    if (verbose) message("No rows after filters; returning empty result.")
    return(data.table::data.table(combo = character(), N = integer(), perc = numeric()))
  }
  
  # Safe index for top_drug
  top_drug <- max(1L, min(top_drug, nrow(Drug_u)))
  top1_substance <- as.character(Drug_u$substance[top_drug])
  
  if (verbose) {
    print(Drug_u[1:min(5L, .N)])
    message("Top-1 substance selected: ", top1_substance)
  }
  
  # Get all primaryids having that top-1 substance as suspect (within pids_inter)
  pids_top1 <- Drug[role_cod == "I" & substance == top1_substance, unique(primaryid)]
  if (length(pids_top1) == 0L) {
    if (verbose) message("No primaryids for top-1 substance; returning empty result.")
    return(data.table::data.table(combo = character(), N = integer(), perc = numeric()))
  }
  
  # Build unique (primaryid, substance, role_cod) for those cases -----------
  DT_unique <- unique(
    Drug[primaryid %in% pids_top1 & !is.na(role_cod),
         .(primaryid, substance = as.character(substance), role_cod)]
  )
  
  # Optional role filter (now actually applied)
  if (!is.null(roles)) {
    DT_unique <- DT_unique[role_cod %in% roles]
  }
  
  # If everything got filtered out, exit early
  if (nrow(DT_unique) == 0L) {
    if (verbose) message("No drugs left after role filter; returning empty result.")
    return(data.table::data.table(combo = character(), N = integer(), perc = numeric()))
  }
  
  # Count unique drugs per primaryid, get the mode size ranking -------------
  drug_counts <- DT_unique[, .(total_drugs = data.table::uniqueN(substance)), by = primaryid]
  
  tab_top1 <- drug_counts[, .N, by = total_drugs][order(-N, -total_drugs)]
  tab_top1[, perc := round(100 * N / sum(N), 1)]
  if (verbose) print(tab_top1[1:min(5L, .N)])
  
  # Choose the top_c-th entry in that ranking safely
  top_c <- max(1L, min(top_c, nrow(tab_top1)))
  mode_size <- tab_top1[top_c, total_drugs]
  if (mode_size < 2) {
    top_c = top_c + 1
    mode_size <- tab_top1[top_c, total_drugs]
  }
  
  
  # Focus on cases with that size ------------------------------------------
  pids_mode <- drug_counts[total_drugs == mode_size, primaryid]
  if (length(pids_mode) == 0L) {
    if (verbose) message("No primaryids with the selected mode size; returning empty result.")
    return(data.table::data.table(combo = character(), N = integer(), perc = numeric()))
  }
  
  # Build the "combination" label (sorted, unique) that CONTAINS top1 -------
  # Since 'mode_size' equals the number of unique drugs in those cases,
  # the "combination of size k" is just the set itself; no need for combn().
  combo_labels <- DT_unique[primaryid %in% pids_mode, .(combo = {
    subs <- sort(unique(substance))
    if (length(subs) == mode_size && top1_substance %in% subs) {
      paste(subs, collapse = " + ")
    } else {NA_character_}}), by = primaryid][!is.na(combo)]
  
  if (nrow(combo_labels) == 0L) {
    if (verbose) message("No valid combinations found; returning empty result.")
    return(data.table::data.table(combo = character(), N = integer(), perc = numeric()))
  }
  
  # Count combinations and add percentages ---------------------------------
  combo_counts <- combo_labels[, .N, by = combo][order(-N, combo)]
  combo_counts[, perc := round(100 * N / sum(N), 1)]
  data.table::setcolorder(combo_counts, c("combo", "N", "perc"))
  
  # Optional summary messages ----------------------------------------------
  if (verbose) {
    msg <- paste0(
      if (!is.null(indi)) paste0("\nIndication: ", indi) else "",
      if (!is.null(reac)) paste0("\nADR: ", reac) else "",
      "\nTop-1 substance: ", top1_substance,
      "\nChosen regimen size (top_2): ", mode_size, "\n"
    )
    message(msg)
    print(combo_counts[1:min(10L, .N)])
  }
  
  # Return result with a couple of attributes for convenience
  attr(combo_counts, "top1_substance") <- top1_substance
  attr(combo_counts, "mode_size") <- mode_size
  return(list("top_comb" = combo_counts, "top_drugs" = Drug_u, 
              "top_quant" = tab_top1, "reac" = data.frame("reac" = reac), 
              "indi" = data.frame("indi" = indi)))
}
make_df <- function(comb, n_rows = NULL, top_drug = 1) {
  # anchor drug (first cell of top_drugs)
  drug1_ <- as.character(comb$top_drugs$substance[top_drug])
  d1_norm <- tolower(trimws(drug1_))
  
  # choose rows (all by default)
  combo_vector <- comb$top_comb$combo
  if (!is.null(n_rows)) combo_vector <- head(combo_vector, n_rows)
  
  # split combos on '+', trim, and drop empties
  split_list <- strsplit(combo_vector, "\\s*\\+\\s*")
  split_list <- lapply(split_list, function(v) {
    v <- trimws(v)
    v[nzchar(v)]
  })
  
  # remove **all** occurrences of drug1_, case-insensitive
  others <- lapply(split_list, function(v) v[tolower(trimws(v)) != d1_norm])
  
  # determine max number of remaining drugs to define columns drug2..drugk
  max_len <- if (length(others)) max(lengths(others)) else 0L
  
  # pad each row to max_len with NA, preserving order
  padded <- lapply(others, function(v) { length(v) <- max_len; v })
  if (max_len == 0L) {
    # no others anywhere: return a 1-col df with just drug1
    return(data.frame(drug1 = rep(drug1_, length(others)), stringsAsFactors = FALSE))
  }
  
  mat <- do.call(rbind, padded)
  df <- data.frame(drug1 = rep(drug1_, nrow(mat)), mat, stringsAsFactors = FALSE)
  colnames(df) <- c("drug1", paste0("drug", 2:(max_len + 1)))
  df
}

# Warfarin
## Atrial fibrillation                 
### Anaemia
#### PS
comb <- most_rep_comb(top_c = 2, top_drug = 2)
comb$top_comb
comb$top_drugs
comb$top_quant
reac <- comb$reac$reac[2]
indi <- comb$indi$indi[2]
drugs_df <- make_df(comb, n_rows = 10, top_drug = 2)

o_war_anemia_ps <- komega(drugs = drugs_df, reactions = reac, title_reac = "Nausea (PS)", indication = indi)

#### SS
comb <- most_rep_comb(merge_reac = T, merge_indi = T, top_reac = 2, top_indi = 2, roles = c("SS"))
reac <- comb$reac$reac[2]
indi <- comb$indi$indi[2]
drugs_df <- make_df(comb, n_rows = 10)

o_war_anemia_ss <- komega(drugs = drugs_df, reactions = reac, title_reac = "Nausea (SS)", indication = indi)

#### C
comb <- most_rep_comb(merge_reac = T, merge_indi = T, top_reac = 2, top_indi = 2, roles = c("C"))
reac <- comb$reac$reac[2]
indi <- comb$indi$indi[2]
drugs_df <- make_df(comb, n_rows = 10)

o_war_anemia_c <- komega(drugs = drugs_df, reactions = reac, title_reac = "Nausea (C)", indication = indi)

## Anticoagulant therapy
# ...

# Other drug
comb <- most_rep_comb(top_drug = 2)

#######

