# Run this entire script before running the main one

get_ratio <- function(s) {
  pids_inter     <- unique(Drug[role_cod == "I"  & substance == s]$primaryid)
  pids_non_inter <- unique(Drug[role_cod != "I" & substance == s]$primaryid)
  length(pids_inter) / length(pids_non_inter)
}

get_react_count <- function(pids) {
  as.data.frame(
    Reac[primaryid %in% pids & !is.na(pt), .N, by = pt][order(-N)][1:10]
  )
}

getN   <- function(role) freq_roles[role_cod == role, N]

getPct <- function(role) freq_roles[role_cod == role, perc]

pastel_palette <- c("#E3F2FD", "#D0E6FB", "#BBDAF8", "#A6CEF4", "#90C2F0", "#7BB7EC", "#66ABE8", "#529FE4",
                    "#D7F1FA", "#C3E9F7", "#AFE1F4", "#9BD9F0", "#87D1ED", "#73C9E9", "#E0F7FA", "#CCF1F5", 
                    "#B8EBF0", "#A4E5EB", "#90DFE6", "#7CD9E1", "#E8F5E9", "#D6ECD7", "#C5E3C6", "#B3DAB5", 
                    "#A2D1A4", "#90C893", "#7FBF82", "#6DB671", "#FFF8E1", "#FFF1CC", "#FFEAB7", "#FFE3A3", 
                    "#FFDC8E", "#FFD57A", "#FFEBD6", "#FFE0C2", "#FFD5AE", "#FFC899", "#FFBC85", "#FFAF70", 
                    "#FCE4EC", "#F8D1E0", "#F4BED4", "#F0ABC8", "#EC98BC", "#F3E5F5", "#E5CEF0", "#D7B8EB", 
                    "#ECEFF1", "#D6DCE0")

pastel_palette0 <- c( "#A6CEF4", "#529FE4", "#9BD9F0", "#A4E5EB", "#B3DAB5", "#6DB671", "#ECEFF1",
                      "#FFE3A3", "#FFC899", "#FFAF70", "#FCE4EC",  "#EC98BC", "#D7B8EB", "#D6DCE0")

years_ <- function(path = NA) {
  d <- read_excel(path)
  d <- t(d)
  colnames(d) <- d[1,]
  d <- as.data.frame(d)[-c(1,3), (ncol(d)-15):ncol(d)]
  d <- as.data.frame(t(d))
  colnames(d) <- "N"
  return(d)
}

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

plot_residuals <- function(resid_mat, title) {
  df_resid <- as.data.frame(resid_mat)
  df_resid$Group <- rownames(df_resid)
  
  df_long <- df_resid |>
    dplyr::mutate(dplyr::across(where(is.factor), as.character)) |>
    tidyr::pivot_longer(
      cols      = where(is.numeric),   # <-- only numeric columns become Residual
      names_to  = "Category",
      values_to = "Residual"
    )
  
  ggplot2::ggplot(df_long, ggplot2::aes(x = Category, y = Group, fill = Residual)) +
    ggplot2::geom_tile(colour = "black") +
    ggplot2::geom_text(ggplot2::aes(label = round(Residual, 1))) +
    ggplot2::scale_fill_gradient2(
      midpoint = 0,
      low = "#BBDAF8",
      mid = "#FFFFFF",
      high = "#FFC899",
      limits = range(df_long$Residual, na.rm = TRUE)
    ) +
    ggplot2::labs(
      title = title,
      x = NULL,
      y = NULL,
      fill = "Std. residual"
    ) +
    ggplot2::theme_void() +
    ggplot2::theme(
      axis.title.y  = ggplot2::element_text(size = 13, face = "bold", angle = 90),
      axis.title.x  = ggplot2::element_text(size = 13, face = "bold"),
      axis.text     = ggplot2::element_text(size = 13, face = "bold"),
      axis.text.x   = ggplot2::element_text(size = 13, face = "bold", angle = 45, hjust = 1, vjust = 1),
      axis.text.y   = ggplot2::element_text(size = 13, face = "bold", angle = 45, hjust = 1),
      plot.title    = ggplot2::element_text(size = 13, face = "bold")
    )
}

cramers_v <- function(tab) {
  chi       <- suppressWarnings(chisq.test(tab, correct = FALSE))
  N         <- sum(tab)
  k         <- min(nrow(tab), ncol(tab))
  V         <- sqrt(chi$statistic / (N * (k - 1)))
  as.numeric(V)
}

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