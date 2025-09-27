library(dplyr)
library(purrr)
library(forestploter)
library(grid)
library(forestploter)
library(grid)

############## general description ##############


###### Upset ######
# helper: coerce medication columns to logical
bool_to_binary <- function(df, cols) {
  df[cols] <- lapply(df[cols], function(x) as.integer(as.logical(x)))
  df
}

select_baseline <- function(df, cols) {
  df %>%
    arrange(visit_no) %>%
    distinct(id, .keep_all = TRUE) %>%
    select(all_of(cols)) %>%
    mutate(across(everything(), as.vector))  # Ensure all are vectors
}

# select_baseline1 <- function(df, cols) {
#   df[df$months_since_baseline == 0, cols, drop = FALSE]
# }


collapse_longitudinal <- function(df, med_cols = NULL) {
  if (is.data.frame(med_cols)) med_cols <- names(med_cols)
  if (is.null(med_cols)) med_cols <- setdiff(names(df), c("id", "dataset"))
  missing_cols <- setdiff(c("id","dataset", med_cols), names(df))
  if (length(missing_cols)) stop("Missing columns: ", paste(missing_cols, collapse = ", "))
  
  df %>%
    select(id, dataset, all_of(med_cols)) %>%
    # coerce meds to 0/1 scalars (handles logicals, numerics with NA, etc.)
    mutate(across(all_of(med_cols), ~ as.integer(replace(.x, is.na(.x), 0) > 0))) %>%
    group_by(id) %>%
    summarise(
      dataset = first(dataset),
      across(all_of(med_cols), ~ max(.x, na.rm = TRUE)),
      .groups = "drop"
    ) %>%
    mutate(dataset = as.factor(dataset)) %>%   # fine for ComplexUpset queries
    as.data.frame()
}


make_ids_unique <- function(df, cohort_name, id_col = "id") {
  df %>%
    mutate(!!id_col := paste0(cohort_name, "-", .data[[id_col]]))
}

###### interacton heatmap ######
# ============================================================
# make_participant_summary()
# - Computes annualized change via baseline_last for a chosen outcome (MMSE/CDR/etc.)
# - Produces one row per id with either:
#     * "baseline" values for med_cols + demo_cols, or
#     * "overall" (ever/first non-missing) collapsed values
# - med_cols can be a character vector or a data.frame (names() taken).
# ============================================================

make_participant_summary <- function(
    df,
    id_col        = "id",
    time_col      = "months_since_baseline",
    visit_col     = NULL,                # optional; if NULL uses time_col to find baseline
    dataset_col   = NULL,                # OPTIONAL now
    value_col     = "MMSE",              # e.g., "MMSE" or "CDR"
    out_col       = NULL,                # e.g., "dMMSE_per_year" / "dCDR_per_year"
    med_cols      = NULL,                # char vec OR data.frame (names taken)
    demo_cols     = NULL,                # char vec of demographic columns
    collapse_mode = c("baseline","overall"),  # how to summarize med+demo across visits
    na_to_zero_meds = TRUE               # NA meds→0 before collapsing when "overall"
) {
  collapse_mode <- match.arg(collapse_mode)
  if (is.null(out_col)) out_col <- paste0("d", value_col, "_per_year")

  # Normalize columns input
  if (is.data.frame(med_cols))  med_cols  <- names(med_cols)
  if (is.null(med_cols))         med_cols  <- character(0)
  if (is.null(demo_cols))        demo_cols <- character(0)

  # Do we have a dataset column to keep?
  use_dataset <- !is.null(dataset_col) && dataset_col %in% names(df)

  # Check required columns
  required <- c(id_col, time_col, value_col, med_cols, demo_cols)
  if (use_dataset) required <- c(required, dataset_col)
  miss <- setdiff(unique(required), names(df))
  if (length(miss)) stop("Missing columns: ", paste(miss, collapse = ", "))

  # Keep only needed columns
  keep_cols <- unique(c(id_col, time_col, value_col, med_cols, demo_cols,
                        if (use_dataset) dataset_col,
                        if (!is.null(visit_col)) visit_col))
  dat <- df[, keep_cols, drop = FALSE]

  # Basic typing
  dat[[id_col]]    <- as.character(dat[[id_col]])
  dat[[time_col]]  <- as.numeric(dat[[time_col]])
  dat[[value_col]] <- as.numeric(dat[[value_col]])
  if (!is.null(visit_col)) dat[[visit_col]] <- as.numeric(dat[[visit_col]])

  # ---------- 1) Annualized change (baseline_last) ----------
  # Define baseline/last by visit if provided, else by time
  if (is.null(visit_col)) {
    dat <- dplyr::arrange(dat, .data[[id_col]], .data[[time_col]])
  } else {
    dat <- dplyr::arrange(dat, .data[[id_col]], .data[[visit_col]])
  }

  slopes <- dat |>
    dplyr::group_by(.data[[id_col]]) |>
    dplyr::summarise(
      n_visits   = dplyr::n(),
      t0         = dplyr::first(.data[[time_col]])/12,
      t1         = dplyr::last(.data[[time_col]])/12,
      y0         = dplyr::first(.data[[value_col]]),
      y1         = dplyr::last(.data[[value_col]]),
      years_span = pmax(0, t1 - t0),
      !!out_col  := dplyr::if_else(years_span > 0, (y1 - y0)/years_span, NA_real_),
      .groups = "drop"
    )

  # ---------- 2) Summarise med_cols + demo_cols (baseline OR overall) ----------
  # helper: pick baseline row per id
  select_baseline_row <- function(d) {
    if (!is.null(visit_col)) {
      d |>
        dplyr::arrange(.data[[visit_col]]) |>
        dplyr::slice(1)
    } else {
      d |>
        dplyr::arrange(.data[[time_col]]) |>
        dplyr::slice(1)
    }
  }

  # helper: first non-missing
  first_non_missing <- function(x) {
    idx <- which(!is.na(x))[1]
    if (length(idx) == 0) NA else x[idx]
  }

  summary_cols <- c(id_col, med_cols, demo_cols, if (use_dataset) dataset_col)

  if (collapse_mode == "baseline") {
    summary_tbl <- dat |>
      dplyr::group_by(.data[[id_col]]) |>
      dplyr::group_modify(~select_baseline_row(.x)) |>
      dplyr::ungroup() |>
      dplyr::select(dplyr::all_of(summary_cols))

    # Coerce meds to 0/1
    if (length(med_cols)) {
      summary_tbl <- summary_tbl |>
        dplyr::mutate(dplyr::across(dplyr::all_of(med_cols), ~ as.integer(.x > 0)))
    }

  } else { # "overall"
    tmp <- dat |>
      dplyr::select(dplyr::all_of(summary_cols))

    # Coerce med cols to 0/1, and optionally set NA->0 before collapsing
    if (length(med_cols)) {
      if (na_to_zero_meds) {
        tmp <- tmp |>
          dplyr::mutate(dplyr::across(dplyr::all_of(med_cols),
                                      ~ as.integer(replace(.x, is.na(.x), 0) > 0)))
      } else {
        tmp <- tmp |>
          dplyr::mutate(dplyr::across(dplyr::all_of(med_cols), ~ as.integer(.x > 0)))
      }
    }

    summary_tbl <- tmp |>
      dplyr::group_by(.data[[id_col]]) |>
      dplyr::summarise(
        !!!(if (use_dataset) rlang::list2(!!dataset_col := dplyr::first(.data[[dataset_col]])) else rlang::list2()),
        dplyr::across(dplyr::all_of(med_cols), ~ max(.x, na.rm = TRUE)),
        dplyr::across(dplyr::all_of(demo_cols), ~ first_non_missing(.x)),
        .groups = "drop"
      )
  }

  # Factorize dataset only if present
  if (use_dataset && dataset_col %in% names(summary_tbl)) {
    summary_tbl[[dataset_col]] <- as.factor(summary_tbl[[dataset_col]])
  }

  # ---------- 3) Merge & return ----------
  # Note: 'slopes' keeps the original id column name (id_col)
  out <- dplyr::left_join(summary_tbl, slopes, by = id_col) |>
    as.data.frame()

  out
}



.effect_one_group <- function(data, meds, outcome = "dMMSE_per_year",
                              covars = c("age"), group_label = "All") {
  stopifnot(outcome %in% names(data))
  map_dfr(meds, function(m) {
    if (!m %in% names(data) || dplyr::n_distinct(data[[m]]) < 2) {
      tibble(med = m, beta = NA_real_, p = NA_real_)
    } else {
      rhs <- paste(c(sprintf("`%s`", m), covars), collapse = " + ")
      fml <- as.formula(paste0("`", outcome, "` ~ ", rhs))
      fit <- lm(fml, data = data)
      s   <- summary(fit)$coefficients
      rn  <- rownames(s)
      target <- which(rn == sprintf("`%s`", m))
      if (length(target) == 0) target <- 2L
      tibble(med = m, beta = s[target,1], p = s[target,4])
    }
  }) |>
    mutate(subgroup = group_label)
}


# Build effects for an arbitrary list of subgroup rows
# groups_df must have columns: label, filter_expr (string)
build_effects <- function(data, meds, outcome,
                          groups_df,
                          covars = c("age"),
                          min_n = 30) {
  stopifnot(all(c("label","filter_expr") %in% names(groups_df)))
  res <- purrr::map_dfr(seq_len(nrow(groups_df)), function(i) {
    lab <- groups_df$label[i]
    expr_str <- groups_df$filter_expr[i]
    expr <- rlang::parse_expr(expr_str)
    d_sub <- dplyr::filter(data, rlang::eval_tidy(expr, data = data))
    if (nrow(d_sub) < min_n) return(NULL)
    .effect_one_group(d_sub, meds, outcome, covars, group_label = lab)
  })
  res |>
    mutate(star = .p_to_star(p))
}

# Plot heatmap
add_pvalue_categories <- function(effects_tbl) {
  effects_tbl |>
    mutate(
      p_category = case_when(
        is.na(p) ~ "No data",
        p < 0.001 ~ "p < 0.001",
        p < 0.01 ~ "p < 0.01", 
        p < 0.05 ~ "p < 0.05",
        p < 0.1 ~ "p < 0.1",
        TRUE ~ "p ≥ 0.1"
      ),
      p_category = factor(p_category, levels = c("p < 0.001", "p < 0.01", "p < 0.05", "p < 0.1", "p ≥ 0.1", "No data"))
    )
}

# NEW: P-value only heatmap function
plot_pvalue_continuous_heatmap <- function(effects_tbl, meds_order = NULL,
                                           title = "P-value Significance Heatmap",
                                           subtitle = "Cell color represents statistical significance (continuous)",
                                           row_order = NULL,
                                           p_limit = 0.1) {  # Cap p-values for better visualization
  df <- effects_tbl
  
  if (is.null(meds_order)) meds_order <- unique(df$med)
  
  # if row_order not given, fall back to median-based sort
  if (is.null(row_order)) {
    row_order <- df |>
      dplyr::group_by(subgroup) |>
      dplyr::summarize(med_beta_median = median(beta, na.rm = TRUE)) |>
      dplyr::arrange(dplyr::desc(med_beta_median)) |>
      dplyr::pull(subgroup)
  }
  
  df$subgroup <- factor(df$subgroup, levels = row_order)
  df$med <- factor(df$med, levels = meds_order)
  
  # Cap p-values for better color scaling (optional)
  df <- df |>
    mutate(p_capped = pmin(p, p_limit, na.rm = TRUE))
  
  ggplot(df, aes(x = med, y = subgroup, fill = p_capped)) +
    geom_tile(color = "white") +
    # Continuous color scale - low p-values (significant) = dark colors
    scale_fill_gradient(
      low = "#E28187",      # Dark red for p ≈ 0 (highly significant)
      high = "#ffffff",     # Light blue for p ≈ 0.1 (not significant)
      limits = c(0, p_limit),
      name = "P-value",
      na.value = "grey90",
      breaks = c(0, 0.01, 0.05, 0.1),
      labels = c("0.00000000001", "0.01", "0.05", "≥0.1")
    ) +
    labs(title = title, subtitle = subtitle, x = " ", y = NULL) +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.title = element_text(size = 7),
      legend.text = element_text(size = 5) 
    )
}


################## Forest plot #################
# make_dt_from_fits <- function(fits, terms, term_labels = NULL) {
#   cohorts <- c("NACC", "AIBL", "HABS")
#   groups <- names(fits)
#   
#   if (is.null(term_labels)) term_labels <- setNames(terms, terms)
#   stopifnot(all(terms %in% names(term_labels)))
#   
#   clean_labels <- trimws(unname(term_labels[terms]))
#   out <- data.frame(Measurements = clean_labels, stringsAsFactors = FALSE)
#   
#   get_sig_stars <- function(p) dplyr::case_when(
#     is.na(p) ~ "", p < 0.001 ~ "***", p < 0.01 ~ "**",
#     p < 0.05 ~ "*", p < 0.1 ~ "^", TRUE ~ ""
#   )
#   
#   for (grp in groups) {
#     grp_fits <- fits[[grp]]
#     stopifnot(length(grp_fits) == length(cohorts))
#     
#     for (i in seq_along(cohorts)) {
#       co <- summary(grp_fits[[i]])$coefficients
#       est <- rep(NA_real_, length(terms))
#       se  <- rep(NA_real_, length(terms))
#       names(est) <- terms; names(se) <- terms
#       
#       available_terms <- intersect(terms, rownames(co))
#       if (length(available_terms) > 0) {
#         est[available_terms] <- co[available_terms, "Estimate"]
#         se[available_terms]  <- co[available_terms, "Std. Error"]
#       }
#       
#       nm <- paste0(grp, "_", cohorts[i])
#       out[[paste0(nm, "_est")]] <- est
#       out[[paste0(nm, "_low")]] <- est - 1.96 * se
#       out[[paste0(nm, "_hi")]]  <- est + 1.96 * se
#     }
#     
#     est_mat  <- matrix(NA_real_, nrow = length(terms), ncol = length(cohorts), dimnames = list(terms, cohorts))
#     se_mat   <- est_mat
#     star_mat <- matrix("", nrow = length(terms), ncol = length(cohorts), dimnames = list(terms, cohorts))
#     
#     for (i in seq_along(cohorts)) {
#       co <- summary(grp_fits[[i]])$coefficients
#       available_terms <- intersect(terms, rownames(co))
#       if (length(available_terms) > 0) {
#         est_mat[available_terms, i] <- co[available_terms, "Estimate"]
#         se_mat[available_terms, i]  <- co[available_terms, "Std. Error"]
#         
#         pvcoli <- intersect(c("Pr(>|t|)", "Pr(>|z|)", "Pr(>|W|)", "Pr(>Chisq)"), colnames(co))[1]
#         if (!is.na(pvcoli)) {
#           pvals <- co[available_terms, pvcoli]
#           star_mat[available_terms, i] <- get_sig_stars(pvals)
#         }
#       }
#     }
#     
#     txt <- vapply(seq_along(terms), function(j) {
#       has_data <- !is.na(est_mat[j, ])
#       if (!any(has_data)) return(paste(rep("", length(cohorts)), collapse = "\n"))
#       lines <- character(length(cohorts))
#       for (k in seq_along(cohorts)) {
#         lines[k] <- if (has_data[k]) {
#           sprintf("%+.3f (%.3f)%s", est_mat[j, k], se_mat[j, k], star_mat[j, k])
#         } else {
#           " "  # cleaner than "NA" or em-dash
#         }
#       }
#       paste(lines, collapse = "\n")
#     }, character(1))
#     
#     out[[paste0(grp, "_Est_SE")]] <- txt
#     out[[grp]] <- paste(rep(" ", 20), collapse = " ")  # empty CI column
#   }
#   
#   out
# }

make_dt_from_fits <- function(fits, terms, term_labels = NULL) {
  cohorts <- c("NACC","AIBL","HABS"); groups <- names(fits)
  if (is.null(term_labels)) term_labels <- setNames(terms, terms)
  stopifnot(all(terms %in% names(term_labels)))
  # clean_labels <- trimws(unname(term_labels[terms]))
  clean_labels <- paste0("     ", trimws(unname(term_labels[terms])))
  out <- data.frame(Measurements = clean_labels, stringsAsFactors = FALSE)
  
  get_sig_stars <- function(p) dplyr::case_when(
    is.na(p) ~ "", p < 0.001 ~ "***", p < 0.01 ~ "**",
    p < 0.05 ~ "*", p < 0.1 ~ "^", TRUE ~ ""
  )
  
  for (grp in groups) {
    grp_fits <- fits[[grp]]; stopifnot(length(grp_fits) == length(cohorts))
    for (i in seq_along(cohorts)) {
      co <- summary(grp_fits[[i]])$coefficients
      est <- se <- rep(NA_real_, length(terms)); names(est) <- names(se) <- terms
      available <- intersect(terms, rownames(co))
      if (length(available) > 0) {
        est[available] <- co[available, "Estimate"]
        se [available] <- co[available, "Std. Error"]
      }
      nm <- paste0(grp,"_",cohorts[i])
      low <- hi <- rep(NA_real_, length(terms)); ok <- !is.na(est) & !is.na(se)
      low[ok] <- est[ok] - 1.96*se[ok]; hi[ok] <- est[ok] + 1.96*se[ok]
      out[[paste0(nm,"_est")]] <- est
      out[[paste0(nm,"_low")]] <- low
      out[[paste0(nm,"_hi")]]  <- hi
    }
    
    est_mat  <- matrix(NA_real_, length(terms), length(cohorts), dimnames=list(terms, cohorts))
    se_mat   <- est_mat; star_mat <- matrix("", length(terms), length(cohorts), dimnames=list(terms, cohorts))
    for (i in seq_along(cohorts)) {
      co <- summary(grp_fits[[i]])$coefficients
      available <- intersect(terms, rownames(co))
      if (length(available) > 0) {
        est_mat[available, i] <- co[available, "Estimate"]
        se_mat [available, i] <- co[available, "Std. Error"]
        pv <- intersect(c("Pr(>|t|)","Pr(>|z|)","Pr(>|W|)","Pr(>Chisq)"), colnames(co))[1]
        if (!is.na(pv)) star_mat[available, i] <- get_sig_stars(co[available, pv])
      }
    }
    
    txt <- vapply(seq_along(terms), function(j){
      has <- !is.na(est_mat[j, ])
      lines <- character(length(cohorts))
      for (k in seq_along(cohorts)) {
        lines[k] <- if (has[k]) sprintf("%+.3f (%.3f)%s", est_mat[j,k], se_mat[j,k], star_mat[j,k]) else "NA"
      }
      paste(lines, collapse = "\n")
    }, character(1))
    
    out[[paste0(grp,"_Est_SE")]] <- txt
    out[[grp]] <- paste(rep(" ", 20), collapse = " ")
  }
  out
}

# create_forest_data <- function(dt_, groups) {
#   # forest_dt <- data.frame(Measurements = trimws(dt_$Measurements), stringsAsFactors = FALSE)
#   forest_dt <- data.frame(Measurements = dt_$Measurements, stringsAsFactors = FALSE)
#   
#   for (grp in groups) {
#     forest_dt[[paste0(grp, "_Est_SE")]] <- trimws(dt_[[paste0(grp, "_Est_SE")]])
#     forest_dt[[grp]] <- paste(rep(" ", 20), collapse = " ")
#   }
#   hdr <- c("Measurements")
#   for (g in groups) hdr <- c(hdr, paste(g, "Estimate (SE)", sep = "\n"), "")
#   colnames(forest_dt) <- hdr
#   forest_dt
# }

create_forest_data <- function(dt_, groups, add_gutter = TRUE, gutter_label = " ") {
  forest_dt <- data.frame(Measurements = dt_$Measurements, stringsAsFactors = FALSE)
  for (grp in groups) {
    forest_dt[[paste0(grp, "_Est_SE")]] <- dt_[[paste0(grp, "_Est_SE")]]
    forest_dt[[grp]] <- paste(rep(" ", 20), collapse = " ")
    if (add_gutter) forest_dt[[paste0(grp, "_sp")]] <- gutter_label
  }
  hdr <- c("Measurements")
  for (g in groups) {
    hdr <- c(hdr, paste(g, "Estimate (SE)", sep = "\n"), "")                 # Est(SE), CI
    if (add_gutter) hdr <- c(hdr, "")                                        # spacer
  }
  colnames(forest_dt) <- hdr
  forest_dt
}


# create_multi_group_forest <- function(
#     dt_,
#     groups  = c("Full","Female","Male","APOE4pos","APOE4neg"),
#     cohorts = c("NACC","AIBL","HABS"),
#     ref_line = 0,
#     title = "",
#     cohort_colors = c("#DAA87C","#92A5D1","#C9DCC4"),
#     cohort_shapes = c(15,16,18),
#     filename = NULL, width = NULL, height = NULL,
#     row_padding_mm = 10, x_pad = 0.15, cohort_nudge = 0.35,
#     point_size = 0.6, line_width = 1.6,
#     xlim_fixed = c(-0.5, 0.5)   # << fixed axis as requested
# ){
#   library(forestploter); library(grid)
#   
#   forest_dt <- create_forest_data(dt_, groups)
#   
#   # Build lists in COHORT-MAJOR order: NACC(all groups), AIBL(all), HABS(all)
#   est_list <- lower_list <- upper_list <- list()
#   for (cohort in cohorts) {
#     for (grp in groups) {
#       est <- dt_[[paste0(grp,"_",cohort,"_est")]]
#       low <- dt_[[paste0(grp,"_",cohort,"_low")]]
#       hi  <- dt_[[paste0(grp,"_",cohort,"_hi")]]
#       if (is.null(est) || length(est)!=nrow(dt_)) est <- rep(NA_real_, nrow(dt_))
#       if (is.null(low) || length(low)!=nrow(dt_)) low <- rep(NA_real_, nrow(dt_))
#       if (is.null(hi)  || length(hi) !=nrow(dt_)) hi  <- rep(NA_real_, nrow(dt_))
#       est_list   <- c(est_list, list(est))
#       lower_list <- c(lower_list, list(low))
#       upper_list <- c(upper_list, list(hi))
#     }
#   }
#   
#   # One CI column per group (not repeated)
#   ci_columns <- seq(3, by = 2, length.out = length(groups))
#   
#   # Use fixed xlim
#   xlim_use <- xlim_fixed
#   
#   # Column widths & theme tweaks
#   subgroup_w <- 0.38
#   group_w    <- 0.62 / length(groups)
#   col_w <- c(subgroup_w, rep(c(group_w*0.72, group_w*0.28), length(groups)))
#   
#   tm <- forest_theme(
#     base_size = 10,
#     refline_gp = gpar(lty = "solid", col = "grey70", lwd = 1.1),  # softer grey ref line
#     ci_pch = cohort_shapes,
#     ci_col = cohort_colors,
#     ci_lwd = line_width,                                          # stronger CI line
#     legend_name = "Cohort",
#     legend_value = cohorts,
#     core = list(
#       fg_params = list(hjust = rep(0, ncol(forest_dt)), x = rep(0.01, ncol(forest_dt)), lineheight = 1.05),
#       padding = unit(c(row_padding_mm, row_padding_mm), "mm")
#     ),
#     colwidth = col_w
#   )
#   
#   p <- forest(
#     forest_dt,
#     est = est_list, lower = lower_list, upper = upper_list,
#     ci_column = ci_columns,
#     ref_line = ref_line, xlim = xlim_use,
#     nudge_y = cohort_nudge, sizes = point_size,
#     theme = tm, title = title,
#     arrow_lab = c("Decrease","Increase"),
#     indent = rep(0, nrow(forest_dt))
#   )
#   
#   p <- edit_plot(p, row = which(forest_dt$Measurements != ""), col = 1, gp = gpar(fontface = "bold"))
#   
#   if (!is.null(filename)) { grDevices::pdf(filename, width = width, height = height); plot(p); grDevices::dev.off() }
#   plot(p); invisible(p)
# }

create_multi_group_forest <- function(
    dt_, groups, cohorts,
    ref_line = 0, title = "",
    cohort_colors = c("#DAA87C","#92A5D1","#C9DCC4"),
    cohort_shapes = c(15,16,18),
    filename = NULL, width = NULL, height = NULL,
    row_padding_mm = 10, xlim_fixed = c(-0.5, 0.5),
    cohort_nudge = 0.35, point_size = 0.6, line_width = 1.6,
    gutter = TRUE, gutter_width = 0.04,         # <<< control the inter-group gap here
    first_col_width = 0.40, est_frac = 0.22, ci_frac = 0.16
){
  library(forestploter); library(grid)
  
  forest_dt <- create_forest_data(dt_, groups, add_gutter = gutter)
  
  ## Build est/low/hi in COHORT-MAJOR order (NACC all groups, AIBL all, HABS all)
  est_list <- lower_list <- upper_list <- list()
  for (cohort in cohorts) for (grp in groups) {
    est <- dt_[[paste0(grp,"_",cohort,"_est")]]
    low <- dt_[[paste0(grp,"_",cohort,"_low")]]
    hi  <- dt_[[paste0(grp,"_",cohort,"_hi")]]
    if (is.null(est) || length(est)!=nrow(dt_)) est <- rep(NA_real_, nrow(dt_))
    if (is.null(low) || length(low)!=nrow(dt_)) low <- rep(NA_real_, nrow(dt_))
    if (is.null(hi)  || length(hi) !=nrow(dt_)) hi  <- rep(NA_real_, nrow(dt_))
    est_list <- c(est_list, list(est)); lower_list <- c(lower_list, list(low)); upper_list <- c(upper_list, list(hi))
  }
  
  ## Map CI columns correctly whether gutter is present or not
  chunk <- if (gutter) 3 else 2                       # per group: Est, CI, (spacer?)
  first_ci_col <- 1 + 2                               # Measurements=1; Est is +1; CI is +2
  ci_columns <- first_ci_col + (0:(length(groups)-1)) * chunk
  
  ## X-axis (fixed)
  xlim_use <- xlim_fixed
  
  ## Column widths: [Measurements] + per-group (Est, CI, [spacer])
  per_group <- if (gutter) c(est_frac, ci_frac, gutter_width) else c(est_frac, ci_frac)
  col_w <- c(first_col_width, rep(per_group, length(groups)))
  
  tm <- forest_theme(
    base_size = 10,
    refline_gp = gpar(lty = "dotted", col = "grey70", lwd = 1.1),
    ci_pch = cohort_shapes, ci_col = cohort_colors, ci_lwd = line_width,
    legend_name = "Cohort", legend_value = cohorts,
    core = list(
      fg_params = list(hjust = rep(0, ncol(forest_dt)), x = rep(0.01, ncol(forest_dt)), lineheight = 1.05),
      padding = unit(c(row_padding_mm, row_padding_mm), "mm")
    ),
    colwidth = col_w
  )
  
  p <- forest(
    forest_dt,
    est = est_list, lower = lower_list, upper = upper_list,
    ci_column = ci_columns,
    ref_line = ref_line, xlim = xlim_use,
    nudge_y = cohort_nudge, sizes = point_size,
    theme = tm, title = title, arrow_lab = c("Decrease","Increase"),
    indent = rep(0, nrow(forest_dt))
  )
  
  p <- edit_plot(p, row = which(forest_dt$Measurements != ""), col = 1, gp = gpar(fontface = "bold"))
  if (!is.null(filename)) { grDevices::pdf(filename, width = width, height = height); plot(p); grDevices::dev.off() }
  plot(p); invisible(p)
}



# plot_forest_subset <- function(
#     fits, terms,
#     groups = c("Full", "Female", "Male"),
#     cohorts = c("NACC", "AIBL", "HABS"),
#     term_labels = NULL,
#     ref_line = 0,
#     cohort_colors = c("#377eb8", "#e41a1c", "#4daf4a"),
#     cohort_shapes = c(15, 16, 18),
#     title = NULL,
#     filename = NULL,
#     width = NULL,
#     height = NULL,
#     row_padding_mm = 10,
#     x_pad = 0.15,
#     cohort_nudge = 0.35,
#     point_size = 0.35,
#     line_width = 1
# ) {
#   fits_subset <- fits[names(fits) %in% groups]
#   fits_subset <- fits_subset[groups]
#   
#   if (is.null(title)) title <- paste("Forest Plot:", paste(groups, collapse = ", "))
#   
#   dt_ <- make_dt_from_fits(fits_subset, terms, term_labels = term_labels)
#   
#   create_multi_group_forest(
#     dt_, groups = groups, cohorts = cohorts,
#     ref_line = ref_line, title = title,
#     cohort_colors = cohort_colors, cohort_shapes = cohort_shapes,
#     filename = filename, width = width, height = height,
#     row_padding_mm = row_padding_mm, x_pad = x_pad,
#     cohort_nudge = cohort_nudge, point_size = point_size, line_width = line_width
#   )
# }


plot_forest_subset <- function(
    fits, terms,
    groups  = c("Full","Female","Male","APOE4pos","APOE4neg"),
    cohorts = c("NACC","AIBL","HABS"),
    term_labels = NULL,
    title = "",
    filename = NULL, width = NULL, height = NULL,
    # layout/style passthrough
    ref_line = 0,
    cohort_colors = c("#DAA87C","#92A5D1","#C9DCC4"),
    cohort_shapes = c(15,16,18),
    row_padding_mm = 10,
    xlim_fixed = c(-0.5, 0.5),
    cohort_nudge = 0.35,
    point_size = 0.6,
    line_width = 1.6,
    gutter = TRUE,
    gutter_width = 0.04,
    first_col_width = 0.40,
    est_frac = 0.22,
    ci_frac = 0.16
){
  fits_subset <- fits[names(fits) %in% groups]
  fits_subset <- fits_subset[groups]
  dt_ <- make_dt_from_fits(fits_subset, terms, term_labels = term_labels)
  
  create_multi_group_forest(
    dt_, groups = groups, cohorts = cohorts,
    ref_line = ref_line, title = title,
    cohort_colors = cohort_colors, cohort_shapes = cohort_shapes,
    filename = filename, width = width, height = height,
    row_padding_mm = row_padding_mm, xlim_fixed = xlim_fixed,
    cohort_nudge = cohort_nudge, point_size = point_size, line_width = line_width,
    gutter = gutter, gutter_width = gutter_width,
    first_col_width = first_col_width, est_frac = est_frac, ci_frac = ci_frac
  )
}








