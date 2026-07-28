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

enforce_static_and_reindex <- function(df, id_col="id", visit_col="visit_no", static_cols, prefer_non_na=TRUE){
  stopifnot(all(c(id_col, visit_col, static_cols) %in% names(df))); library(dplyr)
  df2 <- df %>% group_by(.data[[id_col]]) %>%
    mutate("{visit_col}" := if (any(.data[[visit_col]]==1, na.rm=TRUE)) .data[[visit_col]] else .data[[visit_col]]-min(.data[[visit_col]], na.rm=TRUE)+1L) %>% ungroup()
  base <- df2 %>% group_by(.data[[id_col]]) %>%
    { if (prefer_non_na) summarise(., across(all_of(static_cols), ~{x<-.; x[match(TRUE, !is.na(x))]}, .names="{.col}"), .groups="drop")
      else filter(., .data[[visit_col]]==1) %>% slice_head(n=1) %>% ungroup() %>% select(all_of(c(id_col, static_cols))) }
  df2 %>% select(-all_of(static_cols)) %>% left_join(base, by=id_col)
}


propagate_baseline_comorbidities <- function(
    df,
    id_col           = "id",
    visit_col        = "visit_no",
    comorbidity_cols = c("CVD", "Endocrine", "Psychiatric"),
    output_type      = c("logical", "binary"),
    verbose          = TRUE
) {
  
  library(dplyr)
  
  output_type <- match.arg(output_type)
  
  existing_cols <- intersect(comorbidity_cols, names(df))
  
  if (length(existing_cols) == 0) {
    warning("None of the specified comorbidity columns found in the data.")
    return(df)
  }
  
  if (!visit_col %in% names(df)) {
    stop("`", visit_col, "` column not found. This function requires longitudinal data.")
  }
  
  # ── Step 1: Normalise to 0/1 for processing ──────────────────────────────
  to_binary <- function(x) {
    if (is.character(x)) return(as.integer(toupper(trimws(x)) == "TRUE"))
    if (is.logical(x))   return(as.integer(x))
    return(as.integer(x))
  }
  
  df_work <- df %>%
    mutate(across(all_of(existing_cols), to_binary))
  
  # ── Step 2: Extract baseline (visit 1) values only ───────────────────────
  baseline_vals <- df_work %>%
    group_by(.data[[id_col]]) %>%
    arrange(.data[[visit_col]]) %>%
    slice(1) %>%                              # ← first visit = baseline
    ungroup() %>%
    select(all_of(c(id_col, existing_cols)))  # id + comorbidity cols only
  
  if (verbose) {
    n_participants <- n_distinct(df[[id_col]])
    n_baseline_found <- nrow(baseline_vals)
    cat("✔ [Baseline] Extracted from", n_baseline_found, "of", n_participants, "participants\n")
    cat("  Columns propagated:", paste(existing_cols, collapse = ", "), "\n")
  }
  
  # ── Step 3: Replace comorbidity cols with baseline values across all visits 
  df_out <- df_work %>%
    select(-all_of(existing_cols)) %>%        # drop original comorbidity cols
    left_join(baseline_vals, by = id_col)     # rejoin baseline values to ALL rows
  
  # ── Step 4: Restore to requested output type ─────────────────────────────
  if (output_type == "logical") {
    df_out <- df_out %>%
      mutate(across(all_of(existing_cols), as.logical))
  } else {
    df_out <- df_out %>%
      mutate(across(all_of(existing_cols), as.integer))
  }
  
  # ── Step 5: Restore original column order ────────────────────────────────
  df_out <- df_out %>%
    select(all_of(names(df)))
  
  # ── Verbose summary ───────────────────────────────────────────────────────
  if (verbose) {
    
    # Check how many rows actually changed
    original_binary <- df %>%
      mutate(across(all_of(existing_cols), to_binary))
    
    changed_counts <- sapply(existing_cols, function(col) {
      sum(original_binary[[col]] != as.integer(as.logical(df_out[[col]])),
          na.rm = TRUE)
    })
    
    cat("\n── Summary ──────────────────────────────────────────\n")
    cat("  Participants      :", n_distinct(df[[id_col]]),  "\n")
    cat("  Total rows        :", nrow(df),                  "\n")
    cat("  Output type       :", output_type,               "\n")
    cat("  Rows changed per column:\n")
    for (col in existing_cols) {
      cat(sprintf("    %-20s: %d rows corrected\n", col, changed_counts[col]))
    }
    cat("─────────────────────────────────────────────────────\n")
  }
  
  return(df_out)
}

##   Cleans and standardises longitudinal clinical data by:
##   1. Fixing comorbidity columns to reflect baseline visit values only
##      (a participant's comorbidity status is determined at their first visit)
##   2. Fixing medication columns so that if a participant was ever on a
##      medication across any visit, all their visits are marked as 1/TRUE
##   3. Recalculating Total_Meds, Total_Disease, Med_Control, Disease_Control
##      to stay consistent with the updated values
##   4. Restoring each column group back to its original input type
##      (binary 0/1 or logical TRUE/FALSE) independently per group
clean_longitudinal_data <- function(df, 
                                    id_col             = "id",
                                    visit_col          = "visit_no",
                                    comorbidity_cols   = c("CVD", "Endocrine", "Psychiatric", "Arthritis"),
                                    medication_cols    = c("ACEi", "ARB", "BetaBlk", "CCB", 
                                                           "Diuretic", "Statin", "Metformin"),
                                    control_cols       = c("Med_Control", "Disease_Control"),
                                    comorbidity_type   = c("binary", "logical"),  # ← type per group
                                    medication_type    = c("binary", "logical"),  # ← type per group
                                    control_type       = c("binary", "logical"),  # ← type per group
                                    fix_comorbidities  = TRUE,
                                    fix_medications    = TRUE,
                                    verbose            = TRUE) {
  
  comorbidity_type <- match.arg(comorbidity_type)
  medication_type  <- match.arg(medication_type)
  control_type     <- match.arg(control_type)
  
  df_clean <- df
  
  existing_comorbidity_cols <- intersect(comorbidity_cols, names(df_clean))
  existing_medication_cols  <- intersect(medication_cols,  names(df_clean))
  existing_control_cols     <- intersect(control_cols,     names(df_clean))
  
  
  # ── Step 1: Normalise each group to 0/1 internally for processing ────────────
  to_binary <- function(x) {
    if (is.character(x)) return(as.integer(toupper(trimws(x)) == "TRUE"))
    if (is.logical(x))   return(as.integer(x))
    return(as.integer(x))
  }
  
  df_clean <- df_clean %>%
    mutate(across(all_of(c(existing_comorbidity_cols,
                           existing_medication_cols,
                           existing_control_cols)), to_binary))
  
  if (verbose) cat("✔ [Normalise] All target columns converted to 0/1 internally\n")
  
  
  # ── Step 2: Comorbidities — use BASELINE visit only ──────────────────────────
  if (fix_comorbidities && length(existing_comorbidity_cols) > 0) {
    
    baseline_df <- df_clean %>%
      filter(.data[[visit_col]] == 1) %>%
      select(all_of(c(id_col, existing_comorbidity_cols)))
    
    df_clean <- df_clean %>%
      select(-all_of(existing_comorbidity_cols)) %>%
      left_join(baseline_df, by = id_col)
    
    if (verbose) cat("✔ [Comorbidities] Set to baseline values for:",
                     paste(existing_comorbidity_cols, collapse = ", "), "\n")
    
  } else if (verbose) {
    cat("⏭  [Comorbidities] Skipped (fix_comorbidities = FALSE)\n")
  }
  
  
  # ── Step 3: Medications — if EVER 1, set ALL visits to 1 ─────────────────────
  if (fix_medications && length(existing_medication_cols) > 0) {
    
    ever_df <- df_clean %>%
      group_by(.data[[id_col]]) %>%
      summarise(across(all_of(existing_medication_cols),
                       ~ as.integer(any(.x == 1, na.rm = TRUE))),
                .groups = "drop")
    
    df_clean <- df_clean %>%
      select(-all_of(existing_medication_cols)) %>%
      left_join(ever_df, by = id_col)
    
    if (verbose) cat("✔ [Medications] Set to 1 if ever prescribed for:",
                     paste(existing_medication_cols, collapse = ", "), "\n")
    
  } else if (verbose) {
    cat("⏭  [Medications] Skipped (fix_medications = FALSE)\n")
  }
  
  
  # ── Step 4: Recalculate Total_Meds and Total_Disease ─────────────────────────
  if ("Total_Meds" %in% names(df_clean) && length(existing_medication_cols) > 0) {
    df_clean <- df_clean %>%
      mutate(Total_Meds = rowSums(pick(all_of(existing_medication_cols)), na.rm = TRUE))
    if (verbose) cat("✔ [Totals] Total_Meds recalculated\n")
  }
  
  if ("Total_Disease" %in% names(df_clean) && length(existing_comorbidity_cols) > 0) {
    df_clean <- df_clean %>%
      mutate(Total_Disease = rowSums(pick(all_of(existing_comorbidity_cols)), na.rm = TRUE))
    if (verbose) cat("✔ [Totals] Total_Disease recalculated\n")
  }
  
  
  # ── Step 5: Recalculate Med_Control / Disease_Control ────────────────────────
  if ("Med_Control" %in% existing_control_cols && length(existing_medication_cols) > 0) {
    df_clean <- df_clean %>%
      mutate(Med_Control = as.integer(
        rowSums(pick(all_of(existing_medication_cols)), na.rm = TRUE) > 0))
    if (verbose) cat("✔ [Controls] Med_Control recalculated\n")
  }
  
  if ("Disease_Control" %in% existing_control_cols && length(existing_comorbidity_cols) > 0) {
    df_clean <- df_clean %>%
      mutate(Disease_Control = as.integer(
        rowSums(pick(all_of(existing_comorbidity_cols)), na.rm = TRUE) > 0))
    if (verbose) cat("✔ [Controls] Disease_Control recalculated\n")
  }
  
  
  # ── Step 6: Restore each group back to its original input type ───────────────
  # Helper converters
  restore_type <- function(df, cols, type) {
    cols <- intersect(cols, names(df))
    if (length(cols) == 0) return(df)
    if (type == "logical") {
      df <- df %>% mutate(across(all_of(cols), as.logical))
    } else {
      df <- df %>% mutate(across(all_of(cols), as.integer))
    }
    return(df)
  }
  
  df_clean <- df_clean %>%
    restore_type(existing_comorbidity_cols, comorbidity_type) %>%
    restore_type(existing_medication_cols,  medication_type)  %>%
    restore_type(existing_control_cols,     control_type)
  
  if (verbose) {
    cat("✔ [Output] Comorbidities restored to :", comorbidity_type, "\n")
    cat("✔ [Output] Medications restored to   :", medication_type,  "\n")
    cat("✔ [Output] Controls restored to      :", control_type,     "\n")
  }
  
  
  if (verbose) {
    cat("\n── Summary ──────────────────────────────────────────\n")
    cat("  Input rows        :", nrow(df),           "\n")
    cat("  Output rows       :", nrow(df_clean),      "\n")
    cat("  fix_comorbidities :", fix_comorbidities,   "\n")
    cat("  fix_medications   :", fix_medications,     "\n")
    cat("  Comorbidity type  :", comorbidity_type,    "\n")
    cat("  Medication type   :", medication_type,     "\n")
    cat("  Control type      :", control_type,        "\n")
    cat("─────────────────────────────────────────────────────\n")
  }
  
  return(df_clean)
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
# 
# split_by_progression <- function(df,
#                                  id_col     = "id",
#                                  status_col = "status",
#                                  visit_col  = "visit_no",
#                                  cu_labels  = c("HC", "CU"),
#                                  mci_labels = c("MCI"),
#                                  ad_labels  = c("AD", "Dementia"),
#                                  verbose    = TRUE) {
#   
#   # ── Step 1: Summarise each participant's trajectory ──────────────────────────
#   trajectory <- df %>%
#     arrange(.data[[id_col]], .data[[visit_col]]) %>%
#     group_by(.data[[id_col]]) %>%
#     summarise(
#       baseline_status = first(.data[[status_col]]),
#       final_status    = last(.data[[status_col]]),
#       all_statuses    = list(unique(.data[[status_col]])),
#       n_visits        = n(),
#       .groups         = "drop"
#     )
#   
#   
#   # ── Step 2: Classify each participant into a progression group ───────────────
#   classify_progression <- function(baseline, final, all_statuses) {
#     
#     is_cu  <- baseline %in% cu_labels
#     is_mci <- baseline %in% mci_labels
#     is_ad  <- baseline %in% ad_labels
#     
#     progressed_to_mci <- any(all_statuses %in% mci_labels)
#     progressed_to_ad  <- any(all_statuses %in% ad_labels)
#     
#     # CU stable: starts CU, never reaches MCI or AD
#     if (is_cu && !progressed_to_mci && !progressed_to_ad) return("CU_stable")
#     
#     # MCI stable: starts MCI, never reaches AD
#     if (is_mci && !progressed_to_ad) return("MCI_stable")
#     
#     # AD stable: starts AD, stays AD
#     if (is_ad && final %in% ad_labels) return("AD_stable")
#     
#     # Progressors: CU→MCI, CU→AD, MCI→AD
#     if ((is_cu  && (progressed_to_mci || progressed_to_ad)) ||
#         (is_mci && progressed_to_ad)) return("CU_MCI_prog")
#     
#     # Fallback
#     return(NA_character_)
#   }
#   
#   trajectory <- trajectory %>%
#     rowwise() %>%
#     mutate(progression_group = classify_progression(
#       baseline_status, final_status, all_statuses)) %>%
#     ungroup() %>%
#     select(.data[[id_col]], progression_group)
#   
#   
#   # ── Step 3: Join group labels back to full dataset ───────────────────────────
#   df_labelled <- df %>%
#     left_join(trajectory, by = id_col)
#   
#   
#   # ── Step 4: Split into named list ────────────────────────────────────────────
#   result <- list(
#     CU_stable   = df_labelled %>% filter(progression_group == "CU_stable")   %>% select(-progression_group),
#     MCI_stable  = df_labelled %>% filter(progression_group == "MCI_stable")  %>% select(-progression_group),
#     CU_MCI_prog = df_labelled %>% filter(progression_group == "CU_MCI_prog") %>% select(-progression_group),
#     AD_stable   = df_labelled %>% filter(progression_group == "AD_stable")   %>% select(-progression_group)
#   )
#   
#   
#   # ── Step 5: Verbose summary ───────────────────────────────────────────────────
#   if (verbose) {
#     
#     n_ids <- function(x, id) length(unique(x[[id]]))
#     
#     cat("\n── Progression Group Summary ────────────────────────────\n")
#     cat(sprintf("  %-20s : %4d participants  |  %5d rows\n", 
#                 "CU_stable",   n_ids(result$CU_stable,   id_col), nrow(result$CU_stable)))
#     cat(sprintf("  %-20s : %4d participants  |  %5d rows\n", 
#                 "MCI_stable",  n_ids(result$MCI_stable,  id_col), nrow(result$MCI_stable)))
#     cat(sprintf("  %-20s : %4d participants  |  %5d rows\n", 
#                 "CU_MCI_prog", n_ids(result$CU_MCI_prog, id_col), nrow(result$CU_MCI_prog)))
#     cat(sprintf("  %-20s : %4d participants  |  %5d rows\n", 
#                 "AD_stable",   n_ids(result$AD_stable,   id_col), nrow(result$AD_stable)))
#     cat(sprintf("  %-20s : %4d participants  |  %5d rows\n",
#                 "TOTAL",
#                 n_ids(df_labelled, id_col),
#                 nrow(df)))
#     
#     # Warn about any unclassified participants
#     unclassified <- df_labelled %>% filter(is.na(progression_group))
#     if (nrow(unclassified) > 0) {
#       cat(sprintf("\n  ⚠ WARNING: %d rows could not be classified (NA progression_group)\n",
#                   nrow(unclassified)))
#     }
#     cat("─────────────────────────────────────────────────────────\n\n")
#   }
#   
#   return(result)
# }


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







split_by_progression <- function(
    df,
    id_col = "id",
    status_col = "status",
    visit_col = "visit_no",
    cu_labels = c("HC", "CU"),
    mci_labels = c("MCI"),
    ad_labels = c("AD", "Dementia"),
    verbose = TRUE
) {
  
  library(dplyr)
  
  # ---- trajectories ----
  traj <- df %>%
    arrange(.data[[id_col]], .data[[visit_col]]) %>%
    group_by(.data[[id_col]]) %>%
    summarise(
      baseline_status = first(.data[[status_col]]),
      final_status = last(.data[[status_col]]),
      all_statuses = list(unique(.data[[status_col]])),
      n_visits = n(),
      .groups = "drop"
    )
  
  # ---- classify ----
  classify <- function(baseline, statuses, final) {
    
    is_cu  <- baseline %in% cu_labels
    is_mci <- baseline %in% mci_labels
    is_ad  <- baseline %in% ad_labels
    
    has_cu  <- any(statuses %in% cu_labels)
    has_mci <- any(statuses %in% mci_labels)
    has_ad  <- any(statuses %in% ad_labels)
    
    # Priority 1: Anyone who ever reaches AD goes to AD groups
    if (has_ad) {
      # If started as AD, classify as AD_stable
      if (is_ad) {
        return("AD_stable")
      }
      # If started as CU or MCI and reached AD, classify as progression to AD
      if (is_cu || is_mci) {
        return("CU_MCI_AD_progression")
      }
    }
    
    # Priority 2: CU to MCI progression (never reached AD)
    if (is_cu && has_mci && !has_ad) {
      return("CU_MCI_progression")
    }
    
    # Priority 3: Non-progression (stable or reverters)
    # This includes: stable CU, stable MCI, and MCI reverters (MCI->CU)
    if (!has_ad) {
      return("Non_progression")
    }
    
    # Fallback (should rarely happen)
    return("Non_progression")
  }
  
  traj <- traj %>%
    rowwise() %>%
    mutate(
      progression_group = classify(
        baseline_status,
        all_statuses,
        final_status
      )
    ) %>%
    ungroup()
  
  # ---- merge back ----
  df_labelled <- df %>%
    left_join(
      traj %>% select(all_of(id_col), progression_group),
      by = id_col
    )
  
  # ---- split ----
  result <- split(df_labelled, df_labelled$progression_group)
  
  # ---- summary ----
  if (verbose) {
    
    cat("\n── Progression Group Summary ─────────────────────\n")
    
    for (g in names(result)) {
      
      x <- result[[g]]
      
      cat(
        sprintf(
          "%-22s: %4d participants | %5d rows\n",
          g,
          length(unique(x[[id_col]])),
          nrow(x)
        )
      )
    }
    
    cat("──────────────────────────────────────────────────\n\n")
  }
  
  return(list(
    groups = result,
    labelled_df = df_labelled
  ))
}



first_age_by_med <- function(
    df,
    med_cols,
    age_col = "age",
    id_col = "id",
    visit_col = "visit_no"
) {
  
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(broom)
  
  # ---- reshape long ----
  long <- df %>%
    select(all_of(c(id_col, visit_col, age_col, med_cols))) %>%
    pivot_longer(
      cols = all_of(med_cols),
      names_to = "drug",
      values_to = "user"
    ) %>%
    mutate(user = as.logical(user))
  
  # ---- users: first medication age ----
  users <- long %>%
    filter(user == TRUE) %>%
    arrange(.data[[id_col]], drug, .data[[visit_col]]) %>%
    group_by(.data[[id_col]], drug) %>%
    slice(1) %>%
    ungroup() %>%
    mutate(group = "user")
  
  # ---- non-users: baseline age ----
  non_users <- long %>%
    group_by(.data[[id_col]], drug) %>%
    filter(!any(user == TRUE, na.rm = TRUE)) %>%
    arrange(.data[[visit_col]]) %>%
    slice(1) %>%
    ungroup() %>%
    mutate(group = "non-user")
  
  # ---- combine ----
  dat <- bind_rows(users, non_users)
  
  # ---- summary ----
  out <- dat %>%
    group_by(drug, group) %>%
    summarise(
      n = n(),
      mean = mean(.data[[age_col]], na.rm = TRUE),
      sd = sd(.data[[age_col]], na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      mean = round(mean, 2),
      sd = round(sd, 2),
      mean_sd = paste0(mean, " (", sd, ")")
    )
  
  # ---- Wald test ----
  wald_df <- dat %>%
    group_by(drug) %>%
    group_modify(~{
      
      g <- .x
      
      if (length(unique(g$group)) < 2) {
        
        return(
          tibble(
            wald_chi2 = NA_real_,
            wald_p = NA_real_
          )
        )
      }
      
      g$user_bin <- ifelse(g$group == "user", 1, 0)
      
      fit <- lm(
        as.formula(
          paste(age_col, "~ user_bin")
        ),
        data = g
      )
      
      an <- car::Anova(fit, type = 3)
      
      tibble(
        wald_chi2 = unname(an$`F value`[1]),
        wald_p = unname(an$`Pr(>F)`[1])
      )
    }) %>%
    ungroup()
  
  out <- out %>%
    left_join(wald_df, by = "drug")
  
  return(out)
}

baseline_summary <- function(
    df,
    med_cols = c(
      "ACEi", "ARB", "BetaBlk",
      "CCB", "Diuretic",
      "Statin", "Metformin"
    ),
    disease_cols = c(
      "CVD", "Endocrine", "Psychiatric"
    )
) {
  
  library(dplyr)
  library(tidyr)
  
  df <- df %>% mutate(across(everything(), identity))
  
  # ---- baseline ----
  if ("visit_no" %in% names(df)) {
    
    baseline <- df %>%
      filter(visit_no == 1)
    
  } else {
    
    baseline <- df %>%
      filter(months_since_baseline == 0)
  }
  
  baseline <- baseline %>%
    arrange(id, months_since_baseline) %>%
    distinct(id, .keep_all = TRUE)
  
  n <- n_distinct(baseline$id)
  
  # ---- helpers ----
  fmt_mean_sd <- function(x) {
    
    x <- x[!is.na(x)]
    
    if (length(x) == 0) return(NA)
    
    sprintf(
      "%.1f ± %.1f",
      mean(x),
      sd(x)
    )
  }
  
  fmt_count_pct <- function(x) {
    
    sprintf(
      "%d (%.1f%%)",
      x,
      x / n * 100
    )
  }
  
  fmt_minmax <- function(x) {
    
    x <- x[!is.na(x)]
    
    if (length(x) == 0) return(NA)
    
    sprintf(
      "%.1f – %.1f",
      min(x),
      max(x)
    )
  }
  
  # ---- demographics ----
  age <- fmt_mean_sd(baseline$age)
  
  female <- fmt_count_pct(
    sum(
      toupper(as.character(baseline$sex)) %in% c("F", "FEMALE") |
        baseline$sex == 1,
      na.rm = TRUE
    )
  )
  
  edu <- fmt_mean_sd(baseline$edu)
  
  apoe <- fmt_count_pct(
    sum(
      toupper(as.character(baseline$APOE4)) %in% c("1", "YES", "TRUE") |
        baseline$APOE4 == 1,
      na.rm = TRUE
    )
  )
  
  visits <- fmt_minmax(
    df %>%
      group_by(id) %>%
      summarise(n = n()) %>%
      pull(n)
  )
  
  followup <- fmt_mean_sd(
    df %>%
      group_by(id) %>%
      summarise(
        fu = max(months_since_baseline, na.rm = TRUE) -
          min(months_since_baseline, na.rm = TRUE)
      ) %>%
      pull(fu)
  )
  
  meds <- fmt_mean_sd(baseline$Total_Meds)
  
  # ---- progression groups ----
  prog_counts <- list()
  
  if ("progression_group" %in% names(baseline)) {
    
    prog_name_map <- c(
      "CU_stable" = "Stable CU",
      "MCI_stable" = "Stable MCI",
      "CU_MCI_progression" = "CU/MCI progression",
      "MCI_reverter" = "MCI reversion",
      "AD_stable" = "AD dementia",
      "Mixed" = "Mixed trajectory"
    )
    
    for (k in names(prog_name_map)) {
      
      prog_counts[[prog_name_map[[k]]]] <-
        fmt_count_pct(
          sum(
            baseline$progression_group == k,
            na.rm = TRUE
          )
        )
    }
  }
  
  # ---- medication prevalence ----
  med_prev <- list()
  
  for (m in med_cols) {
    
    if (m %in% names(baseline)) {
      
      med_prev[[paste0(m, " prevalence")]] <-
        fmt_count_pct(
          sum(
            as.logical(baseline[[m]]),
            na.rm = TRUE
          )
        )
    }
  }
  
  # ---- disease prevalence ----
  disease_prev <- list()
  
  for (d in disease_cols) {
    
    if (d %in% names(baseline)) {
      
      disease_prev[[paste0(d, " prevalence")]] <-
        fmt_count_pct(
          sum(
            as.logical(baseline[[d]]),
            na.rm = TRUE
          )
        )
    }
  }
  
  # ---- final table ----
  summary_list <- c(
    list(
      "Age at baseline (year)" = age,
      "Gender (Female)" = female,
      "Education (year)" = edu,
      "APOE4 (YES)" = apoe,
      "Visits (Record)" = visits,
      "Follow up intervals (Month)" = followup,
      "Average medication taken at baseline" = meds
    ),
    prog_counts,
    med_prev,
    disease_prev
  )
  
  tibble(
    Measure = names(summary_list),
    Value = unlist(summary_list)
  )
}


# ---- 2) First medication age by cohort ----
p_format <- function(p) {
  case_when(
    is.na(p) ~ NA_character_,
    p < 0.001 ~ "<0.001***",
    p < 0.01  ~ paste0(sprintf("%.3f", p), "**"),
    p < 0.05  ~ paste0(sprintf("%.3f", p), "*"),
    TRUE      ~ sprintf("%.3f", p)
  )
}

first_age_by_med_table <- function(df, med_cols, cohort_name) {
  
  first_age_by_med(df, med_cols) %>%
    mutate(
      value = paste0(mean, " ± ", sd),
      p_value = p_format(wald_p),
      cohort = cohort_name
    ) %>%
    select(cohort, drug, group, value, p_value) %>%
    pivot_wider(
      names_from = group,
      values_from = value
    ) %>%
    select(cohort, drug, `non-user`, user, p_value)
}


# ============================================================================
# Drug-Drug Interaction Utilities
# ============================================================================
#
# OVERVIEW:
# These functions help identify and model drug-drug interactions in longitudinal data.
# 
# KEY CONSIDERATIONS:
# 1. With 7 medications, there are 21 possible pairwise interactions (7 choose 2)
# 2. Adding all interactions may lead to:
#    - Model convergence issues
#    - Overfitting
#    - Multiple testing concerns
#    - Interpretability challenges
#
# RECOMMENDED APPROACH:
# Option A: Data-driven selection
#   - Use identify_drug_combinations() to find most common pairs
#   - Include only interactions with prevalence > 5%
#   - Typically results in 3-5 key interactions
#
# Option B: Clinically-guided selection
#   - Use get_clinical_drug_pairs() for evidence-based pairs
#   - Focus on known synergistic/antagonistic combinations
#
# Option C: Exploratory analysis
#   - Run models with all interactions separately
#   - Use model comparison (AIC/BIC) to select most important
#
# EXAMPLE USAGE:
# ```r
# # Identify common combinations
# combos <- identify_drug_combinations(NACC_full_df_cdr, min_prevalence = 0.05)
# 
# # Build formula with top 3 interactions
# base_formula <- "CDR ~ year_since_baseline + age + sex + edu + APOE4 + 
#                  CVD + Endocrine + Psychiatric +
#                  ACEi + ARB + BetaBlk + CCB + Diuretic + Metformin + Statin"
# 
# extended_formula <- add_drug_interactions(
#   base_formula, 
#   combos, 
#   n_interactions = 3,
#   include_time_interaction = TRUE
# )
# ```
#


# identify_drug_combinations <- function(
#     df,
#     med_cols = c("ACEi", "ARB", "BetaBlk", "CCB", "Diuretic", "Metformin", "Statin"),
#     min_prevalence = 0.05
# ) {
#   
#   library(dplyr)
#   library(tidyr)
#   
#   # Ensure binary
#   df_binary <- df %>%
#     mutate(across(all_of(med_cols), ~ as.integer(as.logical(.))))
#   
#   # Get baseline or collapse to participant level
#   if ("visit_no" %in% names(df_binary)) {
#     df_baseline <- df_binary %>%
#       group_by(id) %>%
#       summarise(across(all_of(med_cols), ~ max(., na.rm = TRUE)), .groups = "drop")
#   } else {
#     df_baseline <- df_binary
#   }
#   
#   n_total <- nrow(df_baseline)
#   
#   # Calculate all pairwise combinations
#   combinations <- expand.grid(
#     drug1 = med_cols,
#     drug2 = med_cols,
#     stringsAsFactors = FALSE
#   ) %>%
#     filter(drug1 < drug2) # avoid duplicates
#   
#   results <- combinations %>%
#     rowwise() %>%
#     mutate(
#       both = sum(df_baseline[[drug1]] == 1 & df_baseline[[drug2]] == 1, na.rm = TRUE),
#       drug1_only = sum(df_baseline[[drug1]] == 1 & df_baseline[[drug2]] == 0, na.rm = TRUE),
#       drug2_only = sum(df_baseline[[drug1]] == 0 & df_baseline[[drug2]] == 1, na.rm = TRUE),
#       neither = sum(df_baseline[[drug1]] == 0 & df_baseline[[drug2]] == 0, na.rm = TRUE),
#       prevalence = both / n_total,
#       # Conditional probability: P(drug2|drug1)
#       cond_prob_2_given_1 = both / (both + drug1_only),
#       # Conditional probability: P(drug1|drug2)
#       cond_prob_1_given_2 = both / (both + drug2_only)
#     ) %>%
#     ungroup() %>%
#     filter(prevalence >= min_prevalence) %>%
#     arrange(desc(prevalence))
#   
#   return(results)
# }
get_common_drug_pairs <- function(
    cohort_list,
    med_cols       = c("ACEi", "ARB", "BetaBlk", "CCB", "Diuretic", "Metformin", "Statin"),
    min_prevalence = 0.00,
    verbose        = TRUE
) {
  
  library(dplyr)
  library(purrr)
  
  if (is.null(names(cohort_list))) {
    stop("cohort_list must be a named list, e.g. list(NACC = df1, AIBL = df2, HABS = df3)")
  }
  
  # ── Step 1: Run identify_drug_combinations on each cohort ────────────────
  all_results <- cohort_list %>%
    imap_dfr(~ identify_drug_combinations(
      df             = .x,
      med_cols       = med_cols,
      min_prevalence = min_prevalence,
      level          = "record"
    ) %>% mutate(cohort = .y))
  
  if (verbose) {
    cat("── Per-cohort pairs found (above min_prevalence = ", min_prevalence, ") ──\n")
    all_results %>%
      count(cohort) %>%
      { cat(sprintf("  %-25s: %d pairs\n", .$cohort, .$n)); invisible(.) }
  }
  
  # ── Step 2: Keep only pairs present in ALL cohorts ───────────────────────
  n_cohorts <- length(cohort_list)
  
  common_pairs <- all_results %>%
    group_by(drug1, drug2) %>%
    summarise(
      n_cohorts_present   = n_distinct(cohort),
      cohorts_present     = paste(sort(unique(cohort)), collapse = ", "),
      mean_prevalence     = mean(prevalence,  na.rm = TRUE),
      min_prevalence_obs  = min(prevalence,   na.rm = TRUE),
      max_prevalence_obs  = max(prevalence,   na.rm = TRUE),
      mean_lift           = mean(lift,         na.rm = TRUE),
      mean_cond_2_given_1 = mean(cond_prob_2_given_1, na.rm = TRUE),
      mean_cond_1_given_2 = mean(cond_prob_1_given_2, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    filter(n_cohorts_present == n_cohorts) %>%
    arrange(desc(mean_prevalence))
  
  if (verbose) {
    cat("\n── Pairs present in ALL", n_cohorts, "cohorts ──────────────────────────\n")
    cat("  Total common pairs:", nrow(common_pairs), "\n\n")
    print(common_pairs %>% select(drug1, drug2, mean_prevalence,
                                  min_prevalence_obs, max_prevalence_obs, mean_lift))
    cat("────────────────────────────────────────────────────────────\n")
  }
  
  # ── Step 3: Attach per-cohort breakdown as attribute ─────────────────────
  per_cohort_detail <- all_results %>%
    semi_join(common_pairs, by = c("drug1", "drug2")) %>%
    select(drug1, drug2, cohort, n_total, both, prevalence,
           cond_prob_2_given_1, cond_prob_1_given_2, lift) %>%
    arrange(drug1, drug2, cohort)
  
  attr(common_pairs, "per_cohort_detail") <- per_cohort_detail
  attr(common_pairs, "all_results")       <- all_results
  
  return(common_pairs)
}



identify_drug_combinations <- function(
    df,
    med_cols = c("ACEi", "ARB", "BetaBlk", "CCB", "Diuretic", "Metformin", "Statin"),
    min_prevalence = 0.00,
    level = c("record", "participant_baseline", "participant_ever", "participant_majority")
) {
  
  library(dplyr)
  library(tidyr)
  
  level <- match.arg(level)
  
  # ── Step 1: Ensure binary ─────────────────────────────────────────────────
  df_binary <- df %>%
    mutate(across(all_of(med_cols), ~ as.integer(as.logical(.))))
  
  # ── Step 2: Collapse based on level ──────────────────────────────────────
  has_visits <- "visit_no" %in% names(df_binary)
  
  df_analysis <- switch(level,
                        
                        "record" = {
                          if (verbose_flag <- FALSE) cat("✔ Using record-level (all rows)\n")
                          df_binary                                          # ← no collapse, use every row
                        },
                        
                        "participant_baseline" = {
                          if (!has_visits) stop("`visit_no` column required for participant_baseline level")
                          df_binary %>%
                            group_by(id) %>%
                            arrange(visit_no) %>%
                            slice(1) %>%
                            ungroup()
                        },
                        
                        "participant_ever" = {
                          if (!has_visits) stop("`visit_no` column required for participant_ever level")
                          df_binary %>%
                            group_by(id) %>%
                            summarise(across(all_of(med_cols),
                                             ~ as.integer(any(.x == 1, na.rm = TRUE))),
                                      .groups = "drop")
                        },
                        
                        "participant_majority" = {
                          if (!has_visits) stop("`visit_no` column required for participant_majority level")
                          df_binary %>%
                            group_by(id) %>%
                            summarise(across(all_of(med_cols),
                                             ~ as.integer(mean(.x, na.rm = TRUE) >= 0.5)),
                                      .groups = "drop")
                        }
  )
  
  n_total <- nrow(df_analysis)
  
  # ── Step 3: All pairwise combinations ────────────────────────────────────
  combinations <- expand.grid(
    drug1 = med_cols,
    drug2 = med_cols,
    stringsAsFactors = FALSE
  ) %>%
    filter(drug1 < drug2)
  
  results <- combinations %>%
    rowwise() %>%
    mutate(
      level              = level,
      n_total            = n_total,
      both               = sum(df_analysis[[drug1]] == 1 & df_analysis[[drug2]] == 1, na.rm = TRUE),
      drug1_only         = sum(df_analysis[[drug1]] == 1 & df_analysis[[drug2]] == 0, na.rm = TRUE),
      drug2_only         = sum(df_analysis[[drug1]] == 0 & df_analysis[[drug2]] == 1, na.rm = TRUE),
      neither            = sum(df_analysis[[drug1]] == 0 & df_analysis[[drug2]] == 0, na.rm = TRUE),
      prevalence         = both / n_total,
      cond_prob_2_given_1 = both / (both + drug1_only),
      cond_prob_1_given_2 = both / (both + drug2_only),
      prev_drug1         = (both + drug1_only) / n_total,
      prev_drug2         = (both + drug2_only) / n_total,
      lift               = prevalence / (prev_drug1 * prev_drug2)
    ) %>%
    ungroup() %>%
    filter(prevalence > min_prevalence) %>%
    arrange(desc(prevalence))
  
  return(results)
}

#' Add drug-drug interaction terms to model formula
#' 
#' @param base_formula Character string of base formula
#' @param drug_pairs Data frame from identify_drug_combinations (top N rows will be used)
#' @param n_interactions Number of top interactions to include
#' @param include_time_interaction Include three-way interaction with time?
#' @return Updated formula as character string
add_drug_interactions <- function(
    base_formula,
    drug_pairs,
    n_interactions = 5,
    include_time_interaction = TRUE
) {
  
  # Select top N combinations
  top_pairs <- drug_pairs %>%
    slice_head(n = n_interactions)
  
  # Create interaction terms
  interaction_terms <- paste0(
    top_pairs$drug1, ":", top_pairs$drug2
  )
  
  # Add to formula
  new_formula <- paste(base_formula, "+", paste(interaction_terms, collapse = " + "))
  
  # Add three-way interactions with time if requested
  if (include_time_interaction) {
    time_interactions <- paste0(
      "year_since_baseline:", top_pairs$drug1, ":", top_pairs$drug2
    )
    new_formula <- paste(new_formula, "+", paste(time_interactions, collapse = " + "))
  }
  
  return(new_formula)
}


#' Generate recommended drug-drug interaction pairs based on clinical relevance
#' 
#' @return Data frame with clinically relevant drug pairs and rationale
get_clinical_drug_pairs <- function() {
  
  tibble::tribble(
    ~drug1,      ~drug2,        ~rationale,
    "ACEi",      "ARB",         "Often combined for synergistic BP control; contraindicated in some cases",
    "ACEi",      "Diuretic",    "Common combination for hypertension and heart failure",
    "ARB",       "Diuretic",    "Common combination for hypertension",
    "BetaBlk",   "CCB",         "Combined for angina and arrhythmias",
    "Statin",    "Metformin",   "Common in metabolic syndrome/diabetes",
    "Statin",    "BetaBlk",     "Common post-MI or in CVD",
    "Diuretic",  "Metformin",   "Diabetes with hypertension"
  )
}


###############################################################################
# STRATIFIED FOREST PLOT UTILITIES
###############################################################################

#' Extract Stratified Comparison Data from dt_
#'
#' Helper function to extract estimates, CIs, and p-values for stratified comparisons
#'
#' @param dt_ Data frame from make_dt_from_fits
#' @param feature Feature name to extract
#' @param group1 First group name (e.g., "Female")
#' @param group2 Second group name (e.g., "Male")
#' @param cohort Cohort name (e.g., "NACC") or "META" for meta-analysis
#' @param is_meta Logical, whether this is a meta-analysis row
#'
#' @return Named list with est1, low1, hi1, p1, est2, low2, hi2, p2
extract_stratified_estimates <- function(dt_, feature, group1, group2, cohort, is_meta = FALSE) {
  
  # Find feature index
  feat_idx <- which(trimws(gsub("^\\s+", "", dt_$Measurements)) == feature)
  
  if (length(feat_idx) == 0) {
    return(list(
      est1 = NA, low1 = NA, hi1 = NA, p1 = NA,
      est2 = NA, low2 = NA, hi2 = NA, p2 = NA
    ))
  }
  
  if (is_meta) {
    # Meta-analysis estimates
    est1 <- dt_[[paste0(group1, "_META_est")]][feat_idx]
    low1 <- dt_[[paste0(group1, "_META_low")]][feat_idx]
    hi1 <- dt_[[paste0(group1, "_META_hi")]][feat_idx]
    p1 <- dt_[[paste0(group1, "_META_p")]][feat_idx]
    
    est2 <- dt_[[paste0(group2, "_META_est")]][feat_idx]
    low2 <- dt_[[paste0(group2, "_META_low")]][feat_idx]
    hi2 <- dt_[[paste0(group2, "_META_hi")]][feat_idx]
    p2 <- dt_[[paste0(group2, "_META_p")]][feat_idx]
  } else {
    # Cohort-specific estimates
    est1 <- dt_[[paste0(group1, "_", cohort, "_est")]][feat_idx]
    low1 <- dt_[[paste0(group1, "_", cohort, "_low")]][feat_idx]
    hi1 <- dt_[[paste0(group1, "_", cohort, "_hi")]][feat_idx]
    p1 <- dt_[[paste0(group1, "_", cohort, "_p")]][feat_idx]
    
    est2 <- dt_[[paste0(group2, "_", cohort, "_est")]][feat_idx]
    low2 <- dt_[[paste0(group2, "_", cohort, "_low")]][feat_idx]
    hi2 <- dt_[[paste0(group2, "_", cohort, "_hi")]][feat_idx]
    p2 <- dt_[[paste0(group2, "_", cohort, "_p")]][feat_idx]
  }
  
  list(
    est1 = est1, low1 = low1, hi1 = hi1, p1 = p1,
    est2 = est2, low2 = low2, hi2 = hi2, p2 = p2
  )
}


#' Format Estimate with CI and Significance Stars
#'
#' @param est Point estimate
#' @param low Lower CI bound
#' @param hi Upper CI bound
#' @param p P-value
#' @param digits Number of decimal places (default 3)
#'
#' @return Formatted string with estimate [CI] and stars
format_est_ci_stars <- function(est, low, hi, p, digits = 3) {
  stars <- ifelse(is.na(p), "",
                  ifelse(p < 0.001, "***",
                         ifelse(p < 0.01, "**",
                                ifelse(p < 0.05, "*", ""))))
  
  ifelse(is.na(est) | is.na(low) | is.na(hi), 
         "",
         sprintf(paste0("%+.", digits, "f [%.", digits, "f, %.", digits, "f]%s"), 
                 est, low, hi, stars))
}
