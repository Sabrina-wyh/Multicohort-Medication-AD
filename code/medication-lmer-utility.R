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


# missranger_impute <- function(df, columns_ignores) {
#   if (!is.data.frame(df)) {
#     stop("The input must be a dataframe.")
#   }
#   if (!all(columns_ignores %in% colnames(df))) {
#     stop("Some ID columns specified are not present in the dataframe.")
#   }
#   
#   
#   ignore_cols <- df[, columns_ignores, drop = FALSE]
#   data_to_impute <- df[, !colnames(df) %in% columns_ignores, drop = FALSE]
#   set.seed(42)
#   imputed_data <- missRanger(data_to_impute, pmm.k = 10)
#   result <- cbind(ignore_cols, imputed_data)
#   return(result)
# }

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

# enforce_static_and_reindex <- function(df, id_col="id", visit_col="visit_no", static_cols, prefer_non_na=TRUE){
#   stopifnot(all(c(id_col, visit_col, static_cols) %in% names(df)))
#   suppressPackageStartupMessages(require(dplyr))
#   df2 <- df %>%
#     group_by(.data[[id_col]]) %>%
#     mutate("{visit_col}" := if (any(.data[[visit_col]] == 1, na.rm=TRUE)) .data[[visit_col]] 
#            else .data[[visit_col]] - min(.data[[visit_col]], na.rm=TRUE) + 1L) %>%
#     ungroup()
#   base <- df2 %>% group_by(.data[[id_col]]) %>%
#     { if (prefer_non_na) {
#       summarise(., across(all_of(static_cols), ~{x <- .; x[match(TRUE, !is.na(x))] %||% NA}, .names="{.col}"))
#     } else {
#       filter(., .data[[visit_col]]==1) %>% slice_head(n=1) %>% select(all_of(static_cols))
#     }
#     } %>% ungroup() %>% bind_cols(df2 %>% distinct(.data[[id_col]]))
#   names(base)[ncol(base)] <- id_col
#   df2 %>% select(-all_of(static_cols)) %>% left_join(base, by=id_col)
# }

enforce_static_and_reindex <- function(df, id_col="id", visit_col="visit_no", static_cols, prefer_non_na=TRUE){
  stopifnot(all(c(id_col, visit_col, static_cols) %in% names(df))); library(dplyr)
  df2 <- df %>% group_by(.data[[id_col]]) %>%
    mutate("{visit_col}" := if (any(.data[[visit_col]]==1, na.rm=TRUE)) .data[[visit_col]] else .data[[visit_col]]-min(.data[[visit_col]], na.rm=TRUE)+1L) %>% ungroup()
  base <- df2 %>% group_by(.data[[id_col]]) %>%
    { if (prefer_non_na) summarise(., across(all_of(static_cols), ~{x<-.; x[match(TRUE, !is.na(x))]}, .names="{.col}"), .groups="drop")
      else filter(., .data[[visit_col]]==1) %>% slice_head(n=1) %>% ungroup() %>% select(all_of(c(id_col, static_cols))) }
  df2 %>% select(-all_of(static_cols)) %>% left_join(base, by=id_col)
}

# missranger_impute <- function(df, columns_ignores, id_col = "id", static_cols = NULL, seed = 42, ...) {
#   stopifnot(is.data.frame(df))
#   stopifnot(all(columns_ignores %in% names(df)))
#   stopifnot(id_col %in% names(df))
#   if (!is.null(static_cols)) stopifnot(all(static_cols %in% names(df)))
#   
#   # 1) Build per-ID static values (first non-NA); warn if conflicting values exist
#   if (!is.null(static_cols) && length(static_cols)) {
#     warn_conflict <- function(v) length(unique(na.omit(v))) > 1
#     by_id <- lapply(split(df[, c(id_col, static_cols), drop = FALSE], df[[id_col]]), function(d) {
#       vals <- vapply(static_cols, function(col) { x <- d[[col]]; x[which(!is.na(x))[1]] %||% NA }, FUN.VALUE = NA)
#       if (any(vapply(static_cols, function(col) warn_conflict(d[[col]]), logical(1)))) {
#         warning(sprintf("Static columns have conflicting values within some %s groups; using first non-NA.", id_col), call. = FALSE)
#       }
#       c(setNames(list(d[[id_col]][1]), id_col), as.list(vals))
#     })
#     static_map <- do.call(rbind, lapply(by_id, as.data.frame.list, stringsAsFactors = FALSE))
#     rownames(static_map) <- NULL
#     for (col in static_cols) if (is.factor(df[[col]])) static_map[[col]] <- factor(static_map[[col]], levels = levels(df[[col]]))
#   }
#   
#   # 2) Prepare data for missRanger (exclude ID/ignore cols)
#   ignore_cols <- df[, columns_ignores, drop = FALSE]
#   data_to_impute <- df[, setdiff(names(df), columns_ignores), drop = FALSE]
#   
#   # 3) Impute
#   set.seed(seed)
#   imputed_data <- missRanger(data_to_impute, pmm.k = 10, ...)
#   
#   # 4) Re-enforce static columns as constant within id AFTER imputation
#   if (!is.null(static_cols) && length(static_cols)) {
#     imputed_data[[id_col]] <- df[[id_col]]  # ensure join key present
#     imputed_data <- merge(imputed_data, static_map, by = id_col, all.x = TRUE, suffixes = c("", ".static"))
#     for (col in static_cols) {
#       sc <- paste0(col)
#       imputed_data[[sc]] <- imputed_data[[col]] <- imputed_data[[col]] %||% imputed_data[[paste0(col)]]
#       # overwrite with static (even if non-NA) to guarantee constancy
#       imputed_data[[col]] <- imputed_data[[col %>% paste0(".static")]]
#       imputed_data[[paste0(col, ".static")]] <- NULL
#     }
#     imputed_data[[id_col]] <- NULL  # remove duplicate id if it wasn’t originally in data_to_impute
#   }
#   
#   # 5) Restore original column order
#   result <- cbind(ignore_cols, imputed_data)
#   result <- result[, names(df), drop = FALSE]
#   result
# }


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
# split_by_baseline <- function(df, id_col = "id", visit_col = "visit", status_col = "status") {
#   # Find IDs with baseline AD
#   baseline_ids <- df[df[[visit_col]] == 1 & df[[status_col]] == "AD", id_col]
#   
#   # Remove all records for those IDs
#   df_filtered <- df[!df[[id_col]] %in% baseline_ids, ]
#   
#   # Identify baseline HC and MCI IDs
#   hc_ids <- df_filtered[df_filtered[[visit_col]] == 1 & df_filtered[[status_col]] == "HC", id_col]
#   mci_ids <- df_filtered[df_filtered[[visit_col]] == 1 & df_filtered[[status_col]] == "MCI", id_col]
#   ad_ids <- df_filtered[df_filtered[[visit_col]] == 1 & df_filtered[[status_col]] == "AD", id_col]
#   
#   # Subset longitudinal data for those IDs
#   hc_df <- df_filtered[df_filtered[[id_col]] %in% hc_ids, ]
#   mci_df <- df_filtered[df_filtered[[id_col]] %in% mci_ids, ]
#   ad_df <- df_filtered[df_filtered[[id_col]] %in% ad_ids, ]
#   
#   list(HC = hc_df, MCI = mci_df, AD=ad_df)
# }

split_by_baseline <- function(df, id_col = "id", visit_col = "visit", status_col = "status") {
  # Identify baseline HC and MCI IDs
  hc_ids <- df[df[[visit_col]] == 1 & df[[status_col]] == "HC", id_col]
  mci_ids <- df[df[[visit_col]] == 1 & df[[status_col]] == "MCI", id_col]
  ad_ids <- df[df[[visit_col]] == 1 & df[[status_col]] == "AD", id_col]
  
  # Subset longitudinal data for those IDs
  hc_df <- df[df[[id_col]] %in% hc_ids, ]
  mci_df <- df[df[[id_col]] %in% mci_ids, ]
  ad_df <- df[df[[id_col]] %in% ad_ids, ]
  
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
