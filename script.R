# Authors: Jacopo Palombarini, Angela Boccia
# Date: 13/10/2025

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
pids_inter_susp <- unique(Drug[role_cod %in% c("SS", "I")]$primaryid) # 5.789.396

## Count number of drugs per report within SS or I, with at least one I -----
drug_counts_all <- Drug[primaryid %in% pids_inter & role_cod %in% c("SS", "I"), 
                        .N, by = primaryid]
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
bars[, role_cod := factor(role_cod, levels = c("I", "SS"),
                          labels = c("I", "SS"))]
bars <- bars[, total := sum(count), by = total_drugs
][, perc := round(100 * count / total, 1)][]

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
           
## Count number of drugs per report within PS, SS, C or I, with at least one I -----
drug_counts_all_2 <- Drug[primaryid %in% pids_inter & !is.na(role_cod),.N, by = primaryid]
setnames(drug_counts_all_2, "N", "total_drugs")

# table of reports by total_drugs (for 95% cutoff on reports)
tab_2 <- drug_counts_all_2[, .N, by = total_drugs][order(total_drugs)]
tab_2[, cumprop := cumsum(N) / sum(N)]
cutoffN_2 <- tab_2[cumprop >= 0.95, min(total_drugs)]

# join total_drugs to each drug row, then count drugs by total_drugs and role_cod
bars_2 <- Drug[primaryid %in% pids_inter & !is.na(role_cod)][
  drug_counts_all_2, on = "primaryid"
][, .(count = .N), by = .(total_drugs, role_cod)]

# tidy labels & percentages within each bar
bars_2[, role_cod := factor(role_cod, levels = c("PS", "SS", "I", "C"), 
                            labels = c("PS", "SS", "I", "C"))]
bars_2 <- bars_2[, total := sum(count), by = total_drugs
][, perc := round(100 * count / total, 1)][]

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

p2 <- ggplot(bars_2, aes(x = total_drugs, y = count, fill = role_cod)) +
  geom_col(position = "stack") +
  geom_segment(aes(x = cutoffN_2, xend = cutoffN_2, y = 0, yend = y_top_2 - 1),
               linetype = "dashed", color = "black", size = 0.7, alpha = 0.6) +
  geom_text(aes(x = cutoffN_2 + 0.2, y = y_top_2 - 50, 
                label = "95% of reports"), 
            hjust = 0, size = 3, color = "black") +
  labs(x = "Number of drugs (PS or SS or I or C)",
       y = "Number of reports (pids)", fill = "Drug role") +
  scale_x_continuous(breaks = 0:25, limits = c(0, 25),
                     expand = expansion(mult = c(0, 0.02))) +
  scale_y_continuous(breaks = seq(0, y_top_2, by = step_guess),
                     limits = c(0, y_top_2 + step_guess * 0.25)) +
  scale_fill_manual(values = c("I" = "steelblue", "SS" = "grey70",
                               "PS" = "lightgreen", "C" = "orange")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  # segment percentages (per role)
  geom_text(
    aes(label = paste0(perc, "%")),
    position = position_stack(vjust = 0.5),
    size = 3,
    color = "white"
  ) +
  # TOTAL percentage per bar (top label) — no inherited aes
  geom_text(
    data = totals_2,
    mapping = aes(x = total_drugs, y = count_total, label = paste0(perc_total, "%")),
    vjust = -0.3,
    size = 3.5,
    fontface = "bold",
    color = "black",
    inherit.aes = FALSE
  ) + 
  # geom text on the top right corner with total number of reports with at least one I and total number of reports
  annotate("text", x = 19.5, y = y_top_2*0.75, 
           label = paste0("Within total reports with at least one I: ", formatC(sum(totals_2$count_total), format = "d", big.mark = ","), 
                          "\n- At least one PS: ", 
                          formatC(sum(bars_2$count[bars_2$role_cod == "PS"]), format = "d", big.mark = ","),
                          "\n- At least one SS: ", 
                          formatC(sum(bars_2$count[bars_2$role_cod == "SS"]), format = "d", big.mark = ","),
                          "\n- At least one C: ",
                          formatC(sum(bars_2$count[bars_2$role_cod == "C"]), format = "d", big.mark = ","),
                          "\n- No PS, SS or C: ",
                          formatC(sum(bars_2$count[bars_2$role_cod == "I"]), format = "d", big.mark = ",")),
           hjust = 0, size = 3.5, color = "black")

p1 / p2 + plot_layout(heights = c(1, 1.2)) + 
  plot_annotation(title = "Reports by number of reported drugs", 
                  subtitle = "(1986-2024)")

ggplot(bars_2, aes(x = total_drugs, y = perc, 
                   group = role_cod, color = role_cod)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 3) +
  scale_x_continuous(breaks = 2:35, limits = c(2, 35),
                     expand = expansion(mult = c(0, 0.02))) +
  scale_y_continuous(breaks = seq(0, 60, by = 5),
                     limits = c(0, 60)) +
  labs(title = "Percentages by drug role per number of reported drugs",
       x = "Number of drugs (PS or SS or I or C)",
       y = "% by drug role", color = "Drug role") +
  scale_color_manual(values = c("I" = "steelblue", 
                                "SS" = "grey70",
                                "PS" = "lightgreen", 
                                "C" = "orange")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# Descriptive - all data -----
descriptive(pids_inter, file_name = "results/Descriptive_all_data.xlsx")
  
# Most rreported interacting drugs -----
Drug_inter_count <- Drug[role_cod == "I" & !is.na(substance), .N, by= .(substance)][order(-N)]
Drug_inter_count <- Drug_inter_count[1:10, ]
write.xlsx(Drug_inter_count, file = "results/Most_reported_interacting_drugs.xlsx")

# Descriptive for the most reported interacting drugs --------------------------
descriptive(pids_cases = pids_inter, drug = "warfarin", file_name = "results/Descriptive_warfarin.xlsx")
descriptive(pids_cases = pids_inter, drug = "tacrolimus", file_name = "results/Descriptive_tacrolimus.xlsx")
descriptive(pids_cases = pids_inter, drug = "acetylsalicylic acid", file_name = "results/Descriptive_acetylsalicylic_acid.xlsx")
descriptive(pids_cases = pids_inter, drug = "quetiapine", file_name = "results/Descriptive_quetiapine.xlsx")
descriptive(pids_cases = pids_inter, drug = "valproic acid", file_name = "results/Descriptive_valproic_acid.xlsx")
descriptive(pids_cases = pids_inter, drug = "rivaroxaban", file_name = "results/Descriptive_rivaroxaban.xlsx")
descriptive(pids_cases = pids_inter, drug = "paracetamol", file_name = "results/Descriptive_paracetamol.xlsx")
descriptive(pids_cases = pids_inter, drug = "ritonavir", file_name = "results/Descriptive_ritonavir.xlsx")
descriptive(pids_cases = pids_inter, drug = "clozapine", file_name = "results/Descriptive_clozapine.xlsx")
descriptive(pids_cases = pids_inter, drug = "furosemide", file_name = "results/Descriptive_furosemide_acid.xlsx")

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


## ADR -----
# Count ADR repetitions for interacting drugs ---------------------
Reac_inter_count <- Reac[primaryid %in% pids_inter & !is.na(pt), .N, by= .(pt)][order(-N)]
Reac_inter_count <- Reac_inter_count[1:10,]


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
  colnames(react_count[[i]]) <- c(as.character(Drug_inter_count$substance[i]), "N") 
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

# Plot for interacting drugs by classes -----
atc_inter <- ATC[substance %in% Drug_inter_count$substance]
atc_inter_drug <- left_join(Drug_inter_count, atc_inter,  by = "substance")

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


#####

# Omega -----

#####

#######

