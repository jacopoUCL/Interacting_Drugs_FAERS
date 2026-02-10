# Authors: Jacopo Palombarini, Angela Boccia
# Last update: 21/1/2025
# RUN THE "helpers.R" SCRPT FIRST

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
## Selection of reports ----------------------
pids_inter <- unique(Drug[role_cod == "I"]$primaryid) # 95.947
pids_susp <- unique(Drug[role_cod %in% c("SS", "PS")]$primaryid) # 14.738.725
pids_conc <- unique(Drug[role_cod == "C"]$primaryid) # 5.122.564
pids_non_inter <- unique(Drug[role_cod != "I"]$primaryid) # 14.739.289
subs_inter <- unique(Drug[role_cod == "I"]$substance) # 6359
subs_non_inter <- unique(Drug[role_cod != "I"]$substance) # 6359

## NAs analysis ------
# int
Demo_int <- copy(Demo)[primaryid %in% pids_inter]
Drug_int <- copy(Drug)[primaryid %in% pids_inter]
Reac_int <- copy(Reac)[primaryid %in% pids_inter]
Indi_int <- copy(Indi)[primaryid %in% pids_inter]
Outc_int <- copy(Outc)[primaryid %in% pids_inter]
Ther_int <- copy(Ther)[primaryid %in% pids_inter]
datasets_int <- list(Demo_int, Drug_int, Indi_int, Reac_int, Outc_int, Ther_int)
nas_counts_list_int <- data.frame("variable" = character(), "NAs_int" = integer())
for (i in seq_along(datasets_int)) {
  dataset <- datasets_int[[i]]
  for (j in 1:ncol(dataset)) {
    na_count <- sum(is.na(dataset[[j]]))
    nas_counts_list_int <- rbind(nas_counts_list_int, 
                                 data.frame("variable" = colnames(dataset)[j], 
                                            "NAs_int" = na_count))
  }
}

# conc
Demo_conc <- copy(Demo)[primaryid %in% pids_conc]
Drug_conc <- copy(Drug)[primaryid %in% pids_conc]
Reac_conc <- copy(Reac)[primaryid %in% pids_conc]
Indi_conc <- copy(Indi)[primaryid %in% pids_conc]
Outc_conc <- copy(Outc)[primaryid %in% pids_conc]
Ther_conc <- copy(Ther)[primaryid %in% pids_conc]
datasets_conc <- list(Demo_conc, Drug_conc, Indi_conc, Reac_conc, Outc_conc, Ther_conc)
nas_counts_list_conc <- data.frame("variable" = character(), "NAs_conc" = integer())
for (i in seq_along(datasets_conc)) {
  dataset <- datasets_conc[[i]]
  for (j in 1:ncol(dataset)) {
    na_count <- sum(is.na(dataset[[j]]))
    nas_counts_list_conc <- rbind(nas_counts_list_conc, 
                                  data.frame("variable" = colnames(dataset)[j], 
                                             "NAs_conc" = na_count))
  }
}

# susp
Demo_susp <- copy(Demo)[primaryid %in% pids_susp]
Drug_susp <- copy(Drug)[primaryid %in% pids_susp]
Reac_susp <- copy(Reac)[primaryid %in% pids_susp]
Indi_susp <- copy(Indi)[primaryid %in% pids_susp]
Outc_susp <- copy(Outc)[primaryid %in% pids_susp]
Ther_susp <- copy(Ther)[primaryid %in% pids_susp]
datasets_susp <- list(Demo_susp, Drug_susp, Indi_susp, Reac_susp, Outc_susp, Ther_susp)
nas_counts_list_susp <- data.frame("variable" = character(), "NAs_susp" = integer())
for (i in seq_along(datasets_susp)) {
  dataset <- datasets_susp[[i]]
  for (j in 1:ncol(dataset)) {
    na_count <- sum(is.na(dataset[[j]]))
    nas_counts_list_susp <- rbind(nas_counts_list_susp, 
                                  data.frame("variable" = colnames(dataset)[j], 
                                             "NAs_susp" = na_count))
  }
}

setDT(nas_counts_list_int)
setDT(nas_counts_list_conc)
setDT(nas_counts_list_susp)

# remove specific elements
nas_counts_list_int <- nas_counts_list_int[!variable %in% c("primaryid", "wt_in_kgs", "reporter_country", "init_fda_dt", 
                                                            "event_dt", "premarketing", "literature", "RB_duplicates", 
                                                            "RB_duplicates_only_susp", "drug_seq", "drug_rec_act", 
                                                            "start_dt", "dur_in_days","end_dt", "rept_cod", "indi_pt", 
                                                            "role_cod", "pt", "outc_cod", "fda_dt")]
nas_counts_list_conc <- nas_counts_list_conc[!variable %in% c("primaryid", "wt_in_kgs", "reporter_country", "init_fda_dt", 
                                                              "event_dt", "premarketing", "literature", "RB_duplicates", 
                                                              "RB_duplicates_only_susp", "drug_seq", "drug_rec_act", 
                                                              "start_dt", "dur_in_days","end_dt", "rept_cod", "indi_pt", 
                                                              "role_cod", "pt", "outc_cod", "fda_dt")]
nas_counts_list_susp <- nas_counts_list_susp[!variable %in% c("primaryid", "wt_in_kgs", "reporter_country", "init_fda_dt", 
                                                              "event_dt", "premarketing", "literature", "RB_duplicates", 
                                                              "RB_duplicates_only_susp", "drug_seq", "drug_rec_act", 
                                                              "start_dt", "dur_in_days","end_dt", "rept_cod", "indi_pt", 
                                                              "role_cod", "pt", "outc_cod", "fda_dt")]

tot_obs <- data.frame("variable" = nas_counts_list_int$variable, 
                      "N_var" = c(rep(nrow(Demo), 4), nrow(Drug), nrow(Ther)))

df_nas <- merge(
  nas_counts_list_int,
  merge(nas_counts_list_conc,
        merge(nas_counts_list_susp, tot_obs, by = "variable"),
        by = "variable"),
  by = "variable"
)

df_nas$variable <- c("Age", "Reporter", "Country", "Sex", "Substance", "TTO")

pids_I <- unique(Drug[role_cod == "I"]$primaryid)
pids_C <- unique(Drug[role_cod == "C"]$primaryid)
pids_S <- unique(Drug[role_cod %in% c("SS","PS")]$primaryid)

N_role_I <- length(pids_I)
N_role_C <- length(pids_C)
N_role_S <- length(pids_S)

N_role_var_I_demo <- nrow(Demo[primaryid %in% pids_I])
N_role_var_C_demo <- nrow(Demo[primaryid %in% pids_C])
N_role_var_S_demo <- nrow(Demo[primaryid %in% pids_S])

N_role_var_I_drug <- nrow(Drug[role_cod == "I"])
N_role_var_C_drug <- nrow(Drug[role_cod == "C"])
N_role_var_S_drug <- nrow(Drug[role_cod %in% c("SS","PS")])

N_role_var_I_ther <- nrow(Ther[primaryid %in% pids_I])
N_role_var_C_ther <- nrow(Ther[primaryid %in% pids_C])
N_role_var_S_ther <- nrow(Ther[primaryid %in% pids_S])

N_variable_vec <- c(rep(nrow(Demo), 4), nrow(Drug), nrow(Ther))
names(N_variable_vec) <- c("Age","Reporter","Country","Sex","Substance","TTO")

df_plot <- df_nas %>%
  mutate(
    N_variable = N_variable_vec[variable],
    NAs_variable = NAs_int + NAs_conc + NAs_susp
  ) %>%
  tidyr::pivot_longer(
    cols = c(NAs_int, NAs_conc, NAs_susp),
    names_to = "role",
    values_to = "NAs_role_variable"
  ) %>%
  mutate(
    role = dplyr::recode(role,
                         NAs_int = "I",
                         NAs_conc = "C",
                         NAs_susp = "S"),
    N_role = dplyr::case_when(
      role == "I" ~ N_role_I,
      role == "C" ~ N_role_C,
      role == "S" ~ N_role_S,
      TRUE ~ NA_real_
    ),
    N_role_variable = dplyr::case_when(
      variable %in% c("Age","Reporter","Country","Sex") & role == "I" ~ N_role_var_I_demo,
      variable %in% c("Age","Reporter","Country","Sex") & role == "C" ~ N_role_var_C_demo,
      variable %in% c("Age","Reporter","Country","Sex") & role == "S" ~ N_role_var_S_demo,
      
      variable %in% c("Substance") & role == "I" ~ N_role_var_I_drug,
      variable %in% c("Substance") & role == "C" ~ N_role_var_C_drug,
      variable %in% c("Substance") & role == "S" ~ N_role_var_S_drug,
      
      variable %in% c("TTO") & role == "I" ~ N_role_var_I_ther,
      variable %in% c("TTO") & role == "C" ~ N_role_var_C_ther,
      variable %in% c("TTO") & role == "S" ~ N_role_var_S_ther,
      
      TRUE ~ NA_real_
    )
  ) %>%
  dplyr::group_by(role) %>%
  dplyr::ungroup() %>%
  mutate(
    NAs_relative_variable = NAs_variable / N_variable,
    NAs_relative_role_variable = NAs_role_variable / N_role_variable
  ) %>%
  dplyr::select(
    variable,
    role,
    N_variable,
    N_role_variable,
    NAs_variable,
    NAs_role_variable,
    NAs_relative_variable,
    NAs_relative_role_variable
  )

write.xlsx(df_plot, file = "results/tables/nas_by_role.xlsx")

## Descriptive, full dataset -----
descriptive(pids_inter, file_name = "results/tables/Descriptive_all_data_inter.xlsx")
## Most reported interacting drugs -----
Drug_inter_count <- Drug[role_cod == "I" & !is.na(substance), .N, by= .(substance)][order(-N)]
Drug_inter_count$perc <- round(Drug_inter_count$N / sum(Drug_inter_count$N), 3) * 100
## Ratio for most reported interacting drugs over non interacting reports --------------------------
subs <- as.character(Drug_inter_count$substance[1:10])
ratios <- sapply(subs, get_ratio)
ratios_summary <- data.frame(
  substance = subs,
  ratio_int_vs_nonint = round(ratios * 100, 3)
)
most_reported <- merge(Drug_inter_count, ratios_summary, by = "substance")
write.xlsx(most_reported, file = "results/tables/Most_reported_interacting_drugs.xlsx")
## Count ADR repetitions for interacting drugs ---------------------
Reac_inter_count <- Reac[primaryid %in% pids_inter & !is.na(pt), .N, by= .(pt)][order(-N)]
Reac_inter_count <- Reac_inter_count[1:10,]
write.xlsx(Reac_inter_count, file = "results/tables/Most_reported_ADR_for_inter_drugs.xlsx")
## Most reported ADRs for most reported interacting drugs --------------------------
pids_list <- lapply(subs, function(s) {
  unique(Drug[role_cod == "I" & substance == s]$primaryid)
})
names(pids_list) <- subs
react_count <- lapply(pids_list, get_react_count)
react_count <- lapply(react_count, function(x) {
  colnames(x) <- c("substance", "N")
  x
})
react_count_unique <- dplyr::bind_rows(react_count)
reacs <- vector()
for (i in 1:length(names(pids_list))) {
  reacs <- c(reacs, rep(names(pids_list)[i], 10))
}
react_count_unique$mrid <- reacs
write.xlsx(react_count_unique, file = "results/tables/Most_reported_ADR_for_most_reported_inter_drugs.xlsx")

#####

# Plots -----
## NAs -----
pd <- position_dodge(width = 0.75)

df_mark <- df_plot %>%
  filter(role == last(levels(factor(role))) & NAs_variable > 0) %>%
  distinct(variable, NAs_variable, NAs_relative_variable)

df_legend <- tibble::tibble(
  variable = df_plot$variable[1],
  NAs_variable = df_plot$NAs_variable[1],
  NAs_relative_role_variable = 0,
  role = df_plot$role[1],
  leg_line = "NAs per variable / N variable",
  leg_square = "NAs per variable|role / N per variable|role"
)

plot_nas <- ggplot(df_plot, aes(x = reorder(variable, NAs_variable), y = NAs_relative_role_variable, fill = role)) +
  geom_col(width = 0.75, position = pd, colour = "black") +
  geom_segment(data = df_mark, inherit.aes = FALSE, 
               aes(x = reorder(variable, NAs_variable), xend = reorder(variable, NAs_variable), 
                   y = 0, yend = NAs_relative_variable), linewidth = 0.6, color = "black") +
  geom_errorbar(data = df_mark, inherit.aes = FALSE, 
                aes(x = reorder(variable, NAs_variable), ymin = NAs_relative_variable, 
                    ymax = NAs_relative_variable), width = 0.1, linewidth = 0.9, color = "black") +
  geom_segment(
    data = df_legend,
    inherit.aes = FALSE,
    aes(x = 1, xend = 1.6, y = 0, yend = 0, linetype = leg_line),
    color = "black",
    linewidth = 0.8) +
  geom_point(
    data = df_legend,
    inherit.aes = FALSE,
    aes(x = 1, y = 0, shape = leg_square),
    size = 4,
    stroke = 1,
    color = "black",
    fill = "white") +
  scale_fill_manual(values = c("I" = "#BBDAF8", "S" = "#90C893", "C" = "#FFC899")) +
  scale_linetype_manual(
    name = NULL,
    values = c("NAs per variable / N variable" = "solid")) +
  scale_shape_manual(
    name = NULL,
    values = c("NAs per variable|role / N per variable|role" = 22)) +
  scale_y_continuous(
    labels = scales::percent_format(accuracy = 1),
    breaks = seq(0, 1, by = 0.1),
    expand = expansion(mult = c(0, 0.15))) +
  labs(x = "Variables", y = "Relative NAs", fill = "Role") +
  coord_flip() +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "top",
    panel.grid.major.y = element_blank(),
    axis.title = element_text(face = "bold"),
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(face = "bold"),
    axis.text.y = element_text(
      angle = 45, hjust = 1, vjust = -0.5,
      size = 12, face = "bold"),
    axis.text.x = element_text(
      angle = 45, hjust = 1, vjust = 1,
      size = 12, face = "bold")) +
  guides(
    fill = guide_legend(order = 1),
    linetype = guide_legend(order = 2, override.aes = list(color = "black", linewidth = 0.8)),
    shape = guide_legend(order = 3, override.aes = list(fill = "white", color = "black", size = 4)))

ggsave("results/plots/nas_plot.tiff", plot = plot_nas, width = 14, 
       height = 10, units = "in", dpi = 1000, compression = "lzw")

## Count number of drugs per report within PS, SS, C or I, with at least one I -----
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

p1 <- ggplot(bars_2, aes(x = total_drugs, y = count, fill = role_cod)) +
  geom_col(position = "stack") +
  geom_segment(aes(x = cutoffN_2, xend = cutoffN_2, y = 0, yend = y_top_2 - 1),
               linetype = "dashed", color = "black", size = 0.7, alpha = 0.6) +
  geom_text(aes(x = cutoffN_2 + 0.2, y = y_top_2 - 50, 
                label = "95% of reports"), 
            hjust = 0, size = 3, color = "black") +
  labs(x = "Number of drugs (S (PS or SS) or I or C)",
       y = "Number of reports (pids)", fill = "Drug role") +
  scale_x_continuous(breaks = 2:23, limits = c(1, 23),
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
        axis.text   = element_text(size = 10, face = "bold"),
        axis.text.y = element_text(size = 10, face = "bold", angle = 45, hjust = 1.5),
        plot.title  = element_text(size = 10, face = "bold"), 
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(), 
        axis.title.x = element_blank())

## Plot only roles trend per numbers of reported drugs
p2 <- ggplot(bars_2, aes(x = total_drugs, y = perc, 
                         group = role_cod, color = role_cod)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 3) +
  scale_x_continuous(breaks = 2:23, limits = c(1, 23),
                     expand = expansion(mult = c(0, 0.02))) +
  scale_y_continuous(breaks = seq(0, 60, by = 5),
                     limits = c(0, 60)) +
  labs(x = "Number of coreported drugs",
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

ggsave("results/plots/plot_count.tiff", plot = p1 / p2, width = 14, height = 10, units = "in", dpi = 1000, compression = "lzw")

## Plot for interacting drugs by classes -----
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

## Plot trend for most reported interacting drugs -----
int_data <- read_excel("results/tables/rankings/Most_reported_interacting_drugs.xlsx")
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
## Chisq. test ----------------------
descriptive(pids_non_inter, file_name = "results/Descriptive_all_data_non_inter.xlsx")

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
p_chisq <- plot_residuals(test_sex$residuals, "Sex") + 
  plot_residuals(test_age$residuals, "Age") +
  plot_residuals(test_out$residuals, "Outcome") +
  plot_residuals(test_con$residuals, "Continent") +
  plot_layout(ncol = 2)

ggsave("results/plots/plot_chisq.tiff", plot = p_chisq, width = 17, height = 13, units = "in", dpi = 1000, compression = "lzw")


# Cramer's V of imbalance
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
add_chi_sheet(wb, "Sex", test_sex, t(tab_sex))
add_chi_sheet(wb, "Age", test_age, t(tab_age))
add_chi_sheet(wb, "Outcome", test_out, t(tab_out))
add_chi_sheet(wb, "Continent", test_con, t(tab_con))
saveWorkbook(wb, "results/chi_tables.xlsx", overwrite = TRUE)

#####
