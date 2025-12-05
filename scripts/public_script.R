# Authors: Jacopo Palombarini, Angela Boccia
# Last update: 5/12/2025

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
pids_inter <- unique(Drug[role_cod == "I"]$primaryid) # 95.947
# Descriptive, full dataset -----
descriptive(pids_inter, file_name = "results/Descriptive_all_data_inter.xlsx")
# Most reported interacting drugs -----
Drug_inter_count <- Drug[role_cod == "I" & !is.na(substance), .N, by= .(substance)][order(-N)]
Drug_inter_count$perc <- round(Drug_inter_count$N / sum(Drug_inter_count$N), 3) * 100
# Ratio for most reported interacting drugs over non interacting reports --------------------------
subs <- as.character(Drug_inter_count$substance[1:10])
get_ratio <- function(s) {
  pids_inter     <- unique(Drug[role_cod == "I"  & substance == s]$primaryid)
  pids_non_inter <- unique(Drug[role_cod != "I" & substance == s]$primaryid)
  length(pids_inter) / length(pids_non_inter)
}
ratios <- sapply(subs, get_ratio)
ratios_summary <- data.frame(
  substance = subs,
  ratio_int_vs_nonint = round(ratios * 100, 3)
)
most_reported <- merge(Drug_inter_count, ratios_summary, by = "substance")
write.xlsx(most_reported, file = "results/Most_reported_interacting_drugs.xlsx")
# Count ADR repetitions for interacting drugs ---------------------
Reac_inter_count <- Reac[primaryid %in% pids_inter & !is.na(pt), .N, by= .(pt)][order(-N)]
Reac_inter_count <- Reac_inter_count[1:10,]
write.xlsx(Reac_inter_count, file = "results/Most_reported_ADR_for_inter_drugs.xlsx")
# Most reported ADRs for most reported interacting drugs --------------------------
pids_list <- lapply(subs, function(s) {
  unique(Drug[role_cod == "I" & substance == s]$primaryid)
})
names(pids_list) <- subs
get_react_count <- function(pids) {
  as.data.frame(
    Reac[primaryid %in% pids & !is.na(pt), .N, by = pt][order(-N)][1:10]
  )
}
react_count <- lapply(pids_list, get_react_count)
react_count <- lapply(react_count, function(x) {
  colnames(x) <- c("substance", "N")
  x
})
react_count_unique <- dplyr::bind_rows(react_count)
write.xlsx(react_count_unique, file = "results/Most_reported_ADR_for_most_reported_inter_drugs.xlsx")






#####

# Plots -----
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

p1 <- ggplot(bars_2, aes(x = total_drugs, y = count, fill = role_cod)) +
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
        plot.title  = element_text(size = 10, face = "bold"))

## Plot only roles trend per numbers of reported drugs
p2 <- ggplot(bars_2, aes(x = total_drugs, y = perc, 
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

p1 / p2 + plot_layout(heights = c(1, 1.2)) + 
  plot_annotation(title = "Reports by number of reported drugs and roles", 
                  subtitle = "(1986-2024, Reports with at least one interacting drug (I))")

ggsave("results/plots/plot_count.tiff", plot = p2 / p3, width = 14, height = 10, units = "in", dpi = 1000, compression = "lzw")

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

pastel_palette0 <- c( "#A6CEF4", "#529FE4", "#9BD9F0", "#A4E5EB", "#B3DAB5", "#6DB671", "#ECEFF1",
                     "#FFE3A3", "#FFC899", "#FFAF70", "#FCE4EC",  "#EC98BC", "#D7B8EB", "#D6DCE0")

p0 <- ggplot(atc_inter_drug_c1, aes(x = reorder(Class1, -N), y = N)) +
  geom_col(linewidth = 0.15, fill = pastel_palette0, colour = "black") +
  labs(title = "Interacting Drugs by ATC Class 1",
       x = "ATC Class",
       y = "Count of Interacting Drugs") +
  theme_void() +
  theme(
    axis.title.y  = element_text(size = 10, face = "bold", angle = 90),
    axis.title.x  = element_text(size = 10, face = "bold"),
    axis.text     = element_text(size = 10, face = "bold"),
    axis.text.x   = element_text(size = 10, face = "bold", angle = 45, hjust = 1, vjust = 1),
    axis.text.y   = element_text(size = 10, face = "bold", angle = 45, hjust = 1),
    plot.title    = element_text(size = 10, face = "bold"))


ggsave("results/plots/plot0.tiff", plot = p0, width = 12, height = 8, units = "in", dpi = 1000, compression = "lzw")

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
    size = 3.3,
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
    size = 3.3,
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
    size = 3.3,
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

ggsave("results/plots/plot5.tiff", plot = p1, width = 12, height = 8, units = "in", dpi = 1000, compression = "lzw")
ggsave("results/plots/plot6.tiff", plot = p2, width = 12, height = 8, units = "in", dpi = 1000, compression = "lzw")
ggsave("results/plots/plot7.tiff", plot = p3, width = 12, height = 8, units = "in", dpi = 1000, compression = "lzw")

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

# Statistical analyses -----
# Chisq. test ----------------------
descriptive(pids_non_inter, file_name = "results/Descriptive_all_data_non_inter.xlsx")

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
tab_pvals <- data.frame(Variable = c("Sex", "Age", "Outcome", "Continent", ),
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
    geom_tile(, colour = "black") +
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
      axis.title.y  = element_text(size = 13, face = "bold", angle = 90),
      axis.title.x  = element_text(size = 13, face = "bold"),
      axis.text   = element_text(size = 13, face = "bold"),
      axis.text.x = element_text(size = 13, face = "bold", angle = 45, hjust = 1, vjust = 1),
      axis.text.y = element_text(size = 13, face = "bold", angle = 45, hjust = 1),
      plot.title  = element_text(size = 13, face = "bold"))
  
}
p_chisq <- plot_residuals(test_sex$residuals, "Sex") + 
  plot_residuals(test_age$residuals, "Age") +
  plot_residuals(test_out$residuals, "Outcome") +
  plot_residuals(test_con$residuals, "Continent") +
  plot_layout(ncol = 2)

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


