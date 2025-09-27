set.seed(42)
##########preprocessing ############
missing_summary <- function(df) {
  total_missing_pct <- sum(is.na(df)) / (nrow(df) * ncol(df)) * 100
  
  col_missing <- sapply(df, function(x) mean(is.na(x)) * 100)
  col_missing <- col_missing[col_missing > 0]
  col_missing_df <- data.frame(
    column = names(col_missing),
    missing_pct = col_missing
  )
  
  # Print results
  cat("Total missing percentage:", round(total_missing_pct, 2), "%\n\n")
  if (nrow(col_missing_df) > 0) {
    cat("Columns with missing values:\n")
    print(col_missing_df[order(-col_missing_df$missing_pct), ])
  } else {
    cat("No columns have missing values.\n")
  }
}

summarize_data <- function(df) {
  print(length(unique(df$id)))
  max_tab <- df %>%
    dplyr::group_by(id) %>%
    dplyr::summarise(max_meds = max(Total_Meds, na.rm = TRUE), .groups = "drop") %>%
    dplyr::count(max_meds)
  print("Distribution of max Total_Meds per participant:")
  print(max_tab)
  visit_count_summary(df, id_col = "id")
}

clean_by_measure <- function(df, measure = c("MMSE", "CDR"), id_col = "id") {
  measure <- match.arg(measure)
  
  if (measure == "MMSE") {
    # remove rows with missing MMSE
    df <- df[!is.na(df[["MMSE"]]), ]
    # drop the CDR column if present
    if ("CDR" %in% names(df)) {
      df$CDR <- NULL
    }
  } else if (measure == "CDR") {
    # remove rows with missing CDR
    df <- df[!is.na(df[["CDR"]]), ]
    # drop the MMSE column if present
    if ("MMSE" %in% names(df)) {
      df$MMSE <- NULL
    }
  }
  
  # keep only ids with at least 2 records
  id_counts <- table(df[[id_col]])
  valid_ids <- names(id_counts[id_counts >= 2])
  df <- df[df[[id_col]] %in% valid_ids, ]
  
  return(df)
}


missranger_impute <- function(df, columns_ignores) {
  if (!is.data.frame(df)) {
    stop("The input must be a dataframe.")
  }
  if (!all(columns_ignores %in% colnames(df))) {
    stop("Some ID columns specified are not present in the dataframe.")
  }
  
  
  ignore_cols <- df[, columns_ignores, drop = FALSE]
  data_to_impute <- df[, !colnames(df) %in% columns_ignores, drop = FALSE]
  set.seed(42)
  imputed_data <- missRanger(data_to_impute, pmm.k = 10)
  result <- cbind(ignore_cols, imputed_data)
  return(result)
}

visit_count_summary <- function(df, id_col = "id", time_col = "months_since_baseline") {
  visit_counts <- table(df[[id_col]])
  cat("Number of visits per participant:\n")
  for (n in sort(unique(visit_counts))) {
    cat(n, "visit(s):", sum(visit_counts == n), "participants\n")
  }
  
  # duration per participant = max months_since_baseline
  dur_by_id <- tapply(df[[time_col]], df[[id_col]], max, na.rm = TRUE)
  # handle ids with all-NA times (tapply -> -Inf); set to NA and drop
  dur_by_id[is.infinite(dur_by_id)] <- NA_real_
  dur_by_id <- dur_by_id[!is.na(dur_by_id)]
  if (length(dur_by_id)) {
    cat("\nFollow-up duration (months_since_baseline):\n")
    cat("Average:", round(mean(dur_by_id), 2), "\n")
    cat("Minimum:", round(min(dur_by_id), 2), "\n")
    cat("Maximum:", round(max(dur_by_id), 2), "\n")
  } else {
    cat("\nFollow-up duration (", time_col, ") not available (all missing).\n", sep = "")
  }
}

####################### preparing for different data structure #################
split_by_baseline <- function(df, id_col = "id", visit_col = "visit", status_col = "status") {
  # Find IDs with baseline AD
  baseline_ids <- df[df[[visit_col]] == 1 & df[[status_col]] == "AD", id_col]
  
  # Remove all records for those IDs
  df_filtered <- df[!df[[id_col]] %in% baseline_ids, ]
  
  # Identify baseline HC and MCI IDs
  hc_ids <- df_filtered[df_filtered[[visit_col]] == 1 & df_filtered[[status_col]] == "HC", id_col]
  mci_ids <- df_filtered[df_filtered[[visit_col]] == 1 & df_filtered[[status_col]] == "MCI", id_col]
  ad_ids <- df_filtered[df_filtered[[visit_col]] == 1 & df_filtered[[status_col]] == "AD", id_col]
  
  # Subset longitudinal data for those IDs
  hc_df <- df_filtered[df_filtered[[id_col]] %in% hc_ids, ]
  mci_df <- df_filtered[df_filtered[[id_col]] %in% mci_ids, ]
  ad_df <- df_filtered[df_filtered[[id_col]] %in% ad_ids, ]
  
  list(HC = hc_df, MCI = mci_df, AD=ad_df)
}


# ordinal
mmse_stage <- function(mmse) {
  labs <- c("CU", "Questionable", "Mild", "Moderate", "Severe")
  cut(
    mmse,
    breaks = c(-Inf, 10, 20, 25, 29, 30),        # note ordering low→high
    labels = labs,                          # align labels high→low cut
    right = TRUE,
    ordered_result = TRUE
  )
}

cdr_stage <- function(cdr_sb) {
  labs <- c("CU", "Questionable", "Very mild", "Mild", "Moderate", "Severe")
  cut(
    cdr_sb,
    breaks = c(-Inf, 0, 2.5, 4.0, 9.0, 15.5, Inf),
    labels = rev(labs),
    right = TRUE,
    ordered_result = TRUE
  )
}
