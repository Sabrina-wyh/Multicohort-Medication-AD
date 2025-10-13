library(dplyr)
library(purrr)
library(forestploter)
library(grid)
library(forestploter)
library(tidyr)
library(ggplot2)
library(forcats)
library(stringr)
library(scales)

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


###### interacton heatmap 
# ============================================================
# make_participant_summary()
# - Computes annualized change via baseline_last for a chosen outcome (MMSE/CDR/etc.)
# - Produces one row per id with either:
#     * "baseline" values for med_cols + demo_cols, or
#     * "overall" (ever/first non-missing) collapsed values
# - med_cols can be a character vector or a data.frame (names() taken).
# ============================================================

# make_participant_summary <- function(
#     df,
#     id_col        = "id",
#     time_col      = "months_since_baseline",
#     visit_col     = NULL,                # optional; if NULL uses time_col to find baseline
#     dataset_col   = NULL,                # OPTIONAL now
#     value_col     = "MMSE",              # e.g., "MMSE" or "CDR"
#     out_col       = NULL,                # e.g., "dMMSE_per_year" / "dCDR_per_year"
#     med_cols      = NULL,                # char vec OR data.frame (names taken)
#     demo_cols     = NULL,                # char vec of demographic columns
#     collapse_mode = c("baseline","overall"),  # how to summarize med+demo across visits
#     na_to_zero_meds = TRUE               # NA meds→0 before collapsing when "overall"
# ) {
#   collapse_mode <- match.arg(collapse_mode)
#   if (is.null(out_col)) out_col <- paste0("d", value_col, "_per_year")
# 
#   # Normalize columns input
#   if (is.data.frame(med_cols))  med_cols  <- names(med_cols)
#   if (is.null(med_cols))         med_cols  <- character(0)
#   if (is.null(demo_cols))        demo_cols <- character(0)
# 
#   # Do we have a dataset column to keep?
#   use_dataset <- !is.null(dataset_col) && dataset_col %in% names(df)
# 
#   # Check required columns
#   required <- c(id_col, time_col, value_col, med_cols, demo_cols)
#   if (use_dataset) required <- c(required, dataset_col)
#   miss <- setdiff(unique(required), names(df))
#   if (length(miss)) stop("Missing columns: ", paste(miss, collapse = ", "))
# 
#   # Keep only needed columns
#   keep_cols <- unique(c(id_col, time_col, value_col, med_cols, demo_cols,
#                         if (use_dataset) dataset_col,
#                         if (!is.null(visit_col)) visit_col))
#   dat <- df[, keep_cols, drop = FALSE]
# 
#   # Basic typing
#   dat[[id_col]]    <- as.character(dat[[id_col]])
#   dat[[time_col]]  <- as.numeric(dat[[time_col]])
#   dat[[value_col]] <- as.numeric(dat[[value_col]])
#   if (!is.null(visit_col)) dat[[visit_col]] <- as.numeric(dat[[visit_col]])
# 
#   # ---------- 1) Annualized change (baseline_last) ----------
#   # Define baseline/last by visit if provided, else by time
#   if (is.null(visit_col)) {
#     dat <- dplyr::arrange(dat, .data[[id_col]], .data[[time_col]])
#   } else {
#     dat <- dplyr::arrange(dat, .data[[id_col]], .data[[visit_col]])
#   }
# 
#   slopes <- dat |>
#     dplyr::group_by(.data[[id_col]]) |>
#     dplyr::summarise(
#       n_visits   = dplyr::n(),
#       t0         = dplyr::first(.data[[time_col]])/12,
#       t1         = dplyr::last(.data[[time_col]])/12,
#       y0         = dplyr::first(.data[[value_col]]),
#       y1         = dplyr::last(.data[[value_col]]),
#       years_span = pmax(0, t1 - t0),
#       !!out_col  := dplyr::if_else(years_span > 0, (y1 - y0)/years_span, NA_real_),
#       .groups = "drop"
#     )
# 
#   # ---------- 2) Summarise med_cols + demo_cols (baseline OR overall) ----------
#   # helper: pick baseline row per id
#   select_baseline_row <- function(d) {
#     if (!is.null(visit_col)) {
#       d |>
#         dplyr::arrange(.data[[visit_col]]) |>
#         dplyr::slice(1)
#     } else {
#       d |>
#         dplyr::arrange(.data[[time_col]]) |>
#         dplyr::slice(1)
#     }
#   }
# 
#   # helper: first non-missing
#   first_non_missing <- function(x) {
#     idx <- which(!is.na(x))[1]
#     if (length(idx) == 0) NA else x[idx]
#   }
# 
#   summary_cols <- c(id_col, med_cols, demo_cols, if (use_dataset) dataset_col)
# 
#   if (collapse_mode == "baseline") {
#     summary_tbl <- dat |>
#       dplyr::group_by(.data[[id_col]]) |>
#       dplyr::group_modify(~select_baseline_row(.x)) |>
#       dplyr::ungroup() |>
#       dplyr::select(dplyr::all_of(summary_cols))
# 
#     # Coerce meds to 0/1
#     if (length(med_cols)) {
#       summary_tbl <- summary_tbl |>
#         dplyr::mutate(dplyr::across(dplyr::all_of(med_cols), ~ as.integer(.x > 0)))
#     }
# 
#   } else { # "overall"
#     tmp <- dat |>
#       dplyr::select(dplyr::all_of(summary_cols))
# 
#     # Coerce med cols to 0/1, and optionally set NA->0 before collapsing
#     if (length(med_cols)) {
#       if (na_to_zero_meds) {
#         tmp <- tmp |>
#           dplyr::mutate(dplyr::across(dplyr::all_of(med_cols),
#                                       ~ as.integer(replace(.x, is.na(.x), 0) > 0)))
#       } else {
#         tmp <- tmp |>
#           dplyr::mutate(dplyr::across(dplyr::all_of(med_cols), ~ as.integer(.x > 0)))
#       }
#     }
# 
#     summary_tbl <- tmp |>
#       dplyr::group_by(.data[[id_col]]) |>
#       dplyr::summarise(
#         !!!(if (use_dataset) rlang::list2(!!dataset_col := dplyr::first(.data[[dataset_col]])) else rlang::list2()),
#         dplyr::across(dplyr::all_of(med_cols), ~ max(.x, na.rm = TRUE)),
#         dplyr::across(dplyr::all_of(demo_cols), ~ first_non_missing(.x)),
#         .groups = "drop"
#       )
#   }
# 
#   # Factorize dataset only if present
#   if (use_dataset && dataset_col %in% names(summary_tbl)) {
#     summary_tbl[[dataset_col]] <- as.factor(summary_tbl[[dataset_col]])
#   }
# 
#   # ---------- 3) Merge & return ----------
#   # Note: 'slopes' keeps the original id column name (id_col)
#   out <- dplyr::left_join(summary_tbl, slopes, by = id_col) |>
#     as.data.frame()
# 
#   out
# }
# 
# 
# 
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
# plot_pvalue_continuous_heatmap <- function(effects_tbl, meds_order = NULL,
#                                            title = "P-value Significance Heatmap",
#                                            subtitle = "Cell color represents statistical significance (continuous)",
#                                            row_order = NULL,
#                                            p_limit = 0.1) {  # Cap p-values for better visualization
#   df <- effects_tbl
# 
#   if (is.null(meds_order)) meds_order <- unique(df$med)
# 
#   # if row_order not given, fall back to median-based sort
#   if (is.null(row_order)) {
#     row_order <- df |>
#       dplyr::group_by(subgroup) |>
#       dplyr::summarize(med_beta_median = median(beta, na.rm = TRUE)) |>
#       dplyr::arrange(dplyr::desc(med_beta_median)) |>
#       dplyr::pull(subgroup)
#   }
# 
#   df$subgroup <- factor(df$subgroup, levels = row_order)
#   df$med <- factor(df$med, levels = meds_order)
# 
#   # Cap p-values for better color scaling (optional)
#   df <- df |>
#     mutate(p_capped = pmin(p, p_limit, na.rm = TRUE))
# 
#   ggplot(df, aes(x = med, y = subgroup, fill = p_capped)) +
#     geom_tile(color = "white") +
#     # Continuous color scale - low p-values (significant) = dark colors
#     scale_fill_gradient(
#       low = "#E28187",
#       high = "#ffffff",
#       limits = c(0, p_limit),
#       name = "P-value",
#       na.value = "grey90",
#       breaks = c(0.00001, 0.01, 0.05, 0.1),
#       labels = c("0.00001", "0.01", "0.05", "≥0.1")
#     ) +
#     labs(title = title, subtitle = subtitle, x = " ", y = NULL) +
#     theme_minimal(base_size = 12) +
#     theme(
#       panel.grid = element_blank(),
#       axis.text.x = element_text(angle = 45, hjust = 1),
#       legend.title = element_text(size = 7),
#       legend.text = element_text(size = 5)
#     )
# }

plot_pvalue_continuous_heatmap <- function(
    effects_tbl,
    meds_order   = NULL,
    title        = "P-value Significance Heatmap",
    subtitle     = "Cell color represents statistical significance (continuous)",
    row_order    = NULL,
    p_limit      = 0.1,                # cap p-values for better visualization
    dataset_col  = "dataset",          # column holding cohort labels
    dataset_order= c("NACC","AIBL","HABS"),
    nrow_facets  = 1                   # 1 => side-by-side
) {
  df <- effects_tbl
  
  # Require dataset column for faceting
  if (!dataset_col %in% names(df)) {
    stop(sprintf("Column '%s' not found in effects_tbl", dataset_col))
  }
  
  # medication order
  if (is.null(meds_order)) meds_order <- unique(df$med)
  
  # row (subgroup) order: median beta descending if not provided
  if (is.null(row_order)) {
    row_order <- df |>
      dplyr::group_by(subgroup) |>
      dplyr::summarize(med_beta_median = median(beta, na.rm = TRUE), .groups = "drop") |>
      dplyr::arrange(dplyr::desc(med_beta_median)) |>
      dplyr::pull(subgroup)
  }
  
  # apply factor orders
  df$subgroup <- factor(df$subgroup, levels = row_order)
  df$med      <- factor(df$med,      levels = meds_order)
  df[[dataset_col]] <- factor(df[[dataset_col]], levels = dataset_order)
  
  # cap p-values for color scaling
  df <- df |>
    dplyr::mutate(p_capped = pmin(p, p_limit, na.rm = TRUE))
  
  # build heatmap and facet by dataset (ordered NACC, AIBL, HABS)
  ggplot(df, aes(x = med, y = subgroup, fill = p_capped)) +
    geom_tile(color = "white") +
    scale_fill_gradient(
      low    = "#E28187",
      high   = "#FFFFFF",
      limits = c(0, p_limit),
      name   = "P-value",
      na.value = "grey90",
      breaks = c(0.00001, 0.01, 0.05, 0.1),
      labels = c("0.00001", "0.01", "0.05", "≥0.1")
    ) +
    facet_wrap(stats::as.formula(paste("~", dataset_col)), nrow = nrow_facets) +
    labs(title = title, subtitle = subtitle, x = NULL, y = NULL) +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid  = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1),
      strip.text  = element_text(face = "bold"),
      legend.title= element_text(size = 7),
      legend.text = element_text(size = 5)
    )
}

################################################

make_participant_summary <- function(
    df,
    id_col        = "id",
    time_col      = "months_since_baseline",
    visit_col     = NULL,
    dataset_col   = NULL,
    value_col     = "MMSE",
    out_col       = NULL,
    med_cols      = NULL,
    demo_cols     = NULL,
    collapse_mode = c("baseline","overall"),
    na_to_zero_meds = TRUE
) {
  collapse_mode <- match.arg(collapse_mode)
  if (is.null(out_col)) out_col <- paste0("d", value_col, "_per_year")
  
  if (is.data.frame(med_cols))  med_cols  <- names(med_cols)
  if (is.null(med_cols))         med_cols  <- character(0)
  if (is.null(demo_cols))        demo_cols <- character(0)
  
  use_dataset <- !is.null(dataset_col) && dataset_col %in% names(df)
  
  required <- c(id_col, time_col, value_col, med_cols, demo_cols)
  if (use_dataset) required <- c(required, dataset_col)
  miss <- setdiff(unique(required), names(df))
  if (length(miss)) stop("Missing columns: ", paste(miss, collapse = ", "))
  
  keep_cols <- unique(c(id_col, time_col, value_col, med_cols, demo_cols,
                        if (use_dataset) dataset_col,
                        if (!is.null(visit_col)) visit_col))
  dat <- df[, keep_cols, drop = FALSE]
  
  dat[[id_col]]    <- as.character(dat[[id_col]])
  dat[[time_col]]  <- as.numeric(dat[[time_col]])
  dat[[value_col]] <- as.numeric(dat[[value_col]])
  if (!is.null(visit_col)) dat[[visit_col]] <- as.numeric(dat[[visit_col]])
  
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
  
  select_baseline_row <- function(d) {
    if (!is.null(visit_col)) {
      d |> dplyr::arrange(.data[[visit_col]]) |> dplyr::slice(1)
    } else {
      d |> dplyr::arrange(.data[[time_col]])  |> dplyr::slice(1)
    }
  }
  first_non_missing <- function(x) { idx <- which(!is.na(x))[1]; if (length(idx)==0) NA else x[idx] }
  
  summary_cols <- c(id_col, med_cols, demo_cols, if (use_dataset) dataset_col)
  
  if (collapse_mode == "baseline") {
    summary_tbl <- dat |>
      dplyr::group_by(.data[[id_col]]) |>
      dplyr::group_modify(~select_baseline_row(.x)) |>
      dplyr::ungroup() |>
      dplyr::select(dplyr::all_of(summary_cols))
    
    if (length(med_cols)) {
      summary_tbl <- summary_tbl |>
        dplyr::mutate(dplyr::across(dplyr::all_of(med_cols), ~ as.integer(.x > 0)))
    }
  } else {
    tmp <- dat |> dplyr::select(dplyr::all_of(summary_cols))
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
  
  if (use_dataset && dataset_col %in% names(summary_tbl)) {
    summary_tbl[[dataset_col]] <- as.factor(summary_tbl[[dataset_col]])
  }
  
  out <- dplyr::left_join(summary_tbl, slopes, by = id_col) |> as.data.frame()
  
  # store med order for plotting helpers (not used by this function otherwise)
  attr(out, "med_order") <- med_cols
  out
}


###### forest plot ######
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


######### drug proportion
available_meds <- function(df, med_cols) {
  keep <- intersect(med_cols, names(df))
  if (!length(keep)) stop("None of the specified med_cols were found.")
  keep
}

normalize_meds <- function(df, meds) {
  df %>% mutate(across(all_of(meds), ~ as.numeric(replace(., is.na(.), 0))))
}

tag_dataset <- function(df, label, meds) {
  if (is.null(df) || nrow(df) == 0) return(tibble())
  meds_in_data <- available_meds(df, meds)
  df %>%
    select(any_of(c("id")), all_of(meds_in_data)) %>%
    mutate(dataset = label) %>%
    relocate(dataset, id) %>%
    normalize_meds(meds_in_data)
}


###### drug proportion ######



available_meds <- function(df, med_cols) {
  keep <- intersect(med_cols, names(df))
  if (!length(keep)) stop("None of the specified med_cols were found.")
  keep
}

normalize_meds <- function(df, meds) {
  df %>% dplyr::mutate(dplyr::across(tidyselect::all_of(meds), ~ as.numeric(replace(., is.na(.), 0))))
}

tag_dataset <- function(df, label, meds) {
  if (is.null(df) || !nrow(df)) return(tibble::tibble())
  meds_in_data <- available_meds(df, meds)
  df %>%
    dplyr::select(tidyselect::any_of("id"), dplyr::all_of(meds_in_data)) %>%
    dplyr::mutate(dataset = label) %>%
    dplyr::relocate(dataset, id) %>%
    normalize_meds(meds_in_data)
}

# Long table with per-med and "Any medication" proportions at baseline
make_prop_table_baseline <- function(baseline_df, meds) {
  med_long <- baseline_df %>%
    tidyr::pivot_longer(tidyselect::all_of(meds), names_to = "Medication", values_to = "has_med") %>%
    dplyr::group_by(dataset, Medication) %>%
    dplyr::summarise(n = dplyr::n_distinct(id), prop_yes = mean(has_med == 1), .groups = "drop") %>%
    dplyr::mutate(prop_no = 1 - prop_yes) %>%
    tidyr::pivot_longer(c(prop_yes, prop_no), names_to = "Status", values_to = "Proportion") %>%
    dplyr::mutate(Status = dplyr::if_else(Status == "prop_yes", "Yes", "No"),
                  Scope  = "Baseline")
  
  any_df <- baseline_df %>%
    dplyr::mutate(any_med = as.integer(rowSums(as.matrix(dplyr::across(tidyselect::all_of(meds)))) > 0)) %>%
    dplyr::group_by(dataset) %>%
    dplyr::summarise(n = dplyr::n_distinct(id), prop_yes = mean(any_med == 1), .groups = "drop") %>%
    dplyr::mutate(prop_no = 1 - prop_yes, Medication = "Any medication") %>%
    tidyr::pivot_longer(c(prop_yes, prop_no), names_to = "Status", values_to = "Proportion") %>%
    dplyr::mutate(Status = dplyr::if_else(Status == "prop_yes", "Yes", "No"),
                  Scope  = "Baseline")
  
  dplyr::bind_rows(med_long, any_df) %>%
    dplyr::select(dataset, Scope, Medication, Status, Proportion, n)
}

# Plot (keeps med_cols order top->bottom; "Any medication" is last)
plot_baseline_props <- function(prop_df,
                                meds_order,
                                dataset_order = c("AIBL","HABS","NACC"),
                                show_only_any = FALSE,
                                fill_values = c("Yes" = "#2C7BB6", "No" = "#BDBDBD"),
                                add_borders = FALSE,
                                text_size = 9) {
  
  med_levels <- c(meds_order, "Any medication")
  
  df <- prop_df %>%
    dplyr::mutate(
      dataset    = factor(dataset, levels = dataset_order),
      Medication = factor(Medication, levels = med_levels),
      Status     = factor(Status, levels = c("Yes","No"))
    )
  
  if (show_only_any) df <- dplyr::filter(df, Medication == "Any medication")
  
  ggplot2::ggplot(df, ggplot2::aes(x = Medication, y = Proportion, fill = Status)) +
    ggplot2::geom_col(position = "fill",
                      color = if (add_borders) "white" else NA,
                      linewidth = 0.2) +
    ggplot2::scale_fill_manual(values = fill_values, drop = FALSE) +
    ggplot2::scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    ggplot2::labs(x = NULL, y = "Proportion of participants", fill = "Has medication",
                  title = " ", subtitle = " ") +
    ggplot2::facet_wrap(~ dataset, nrow = 1) +
    ggplot2::coord_flip() +
    # Key: reverse x (which becomes y after flip) so top→bottom matches med_cols; "Any" at bottom
    ggplot2::scale_x_discrete(limits = rev(med_levels), drop = FALSE) +
    ggplot2::theme_minimal(base_size = text_size) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_text(hjust = 1, size = text_size - 1),
      axis.text.x = ggplot2::element_text(size = text_size - 1),
      legend.text = ggplot2::element_text(size = text_size - 1),
      legend.title = ggplot2::element_text(size = text_size),
      plot.title = ggplot2::element_text(size = text_size + 1, face = "bold"),
      plot.subtitle = ggplot2::element_text(size = text_size - 1),
      legend.position = "bottom",
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "bold", size = text_size)
    )
}

# PDF wrapper
baseline_props_to_pdf <- function(
    NACC_df, AIBL_df, HABS_df,
    out_pdf = "../plots/heatmap/MMSE_baselineMed_pvalue_heatmap_groups.pdf",
    width = 10, height = 5,
    med_cols = c("ACEi","ARB","BetaBlk","CCB","Diuretic","Metformin","Statin"),
    dataset_order = c("NACC","AIBL","HABS"),
    show_only_any = FALSE,
    fill_values = c("Yes" = "#2C7BB6", "No" = "#BDBDBD"),
    add_borders = FALSE,
    text_size = 9) {
  
  aibl_b <- tag_dataset(AIBL_df, "AIBL", med_cols)
  habs_b <- tag_dataset(HABS_df, "HABS", med_cols)
  nacc_b <- tag_dataset(NACC_df, "NACC", med_cols)
  baseline_all <- dplyr::bind_rows(aibl_b, habs_b, nacc_b)
  
  meds_in <- available_meds(baseline_all, med_cols)
  tbl_baseline <- make_prop_table_baseline(baseline_all, meds_in)
  
  p <- plot_baseline_props(
    prop_df = tbl_baseline,
    meds_order = meds_in,
    dataset_order = dataset_order,
    show_only_any = show_only_any,
    fill_values = fill_values,
    add_borders = add_borders,
    text_size = text_size
  )
  
  out_dir <- dirname(out_pdf)
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  
  grDevices::pdf(out_pdf, width = width, height = height)
  on.exit(grDevices::dev.off(), add = TRUE)
  print(p)
  
  invisible(list(plot = p, table = tbl_baseline, meds_in = meds_in,
                 path = out_pdf, width = width, height = height))
}

# Neat table (with "Any medication" last within each dataset)
baseline_props_table <- function(
    NACC_df, AIBL_df, HABS_df,
    med_cols = c("ACEi","ARB","BetaBlk","CCB","Diuretic","Statin","Metformin"),
    dataset_order = c("AIBL","HABS","NACC"),
    include_individual_meds = TRUE,
    include_any_med = TRUE,
    digits = 1,
    out_csv = NULL) {
  
  stopifnot(include_individual_meds || include_any_med)
  
  aibl_b <- tag_dataset(AIBL_df, "AIBL", med_cols)
  habs_b <- tag_dataset(HABS_df, "HABS", med_cols)
  nacc_b <- tag_dataset(NACC_df, "NACC", med_cols)
  baseline_all <- dplyr::bind_rows(aibl_b, habs_b, nacc_b)
  
  meds_in <- available_meds(baseline_all, med_cols)
  
  med_rows <- NULL
  if (include_individual_meds) {
    med_rows <- baseline_all %>%
      tidyr::pivot_longer(tidyselect::all_of(meds_in),
                          names_to = "Medication", values_to = "has_med") %>%
      dplyr::group_by(dataset, Medication) %>%
      dplyr::summarise(
        n_total = dplyr::n_distinct(id),
        n_yes   = sum(has_med == 1, na.rm = TRUE),
        n_no    = n_total - n_yes,
        prop_yes = n_yes / n_total,
        prop_no  = n_no / n_total,
        .groups = "drop"
      )
  }
  
  any_rows <- NULL
  if (include_any_med) {
    any_rows <- baseline_all %>%
      dplyr::mutate(any_med = as.integer(rowSums(as.matrix(dplyr::across(tidyselect::all_of(meds_in)))) > 0)) %>%
      dplyr::group_by(dataset) %>%
      dplyr::summarise(
        Medication = "Any medication",
        n_total = dplyr::n_distinct(id),
        n_yes   = sum(any_med == 1, na.rm = TRUE),
        n_no    = n_total - n_yes,
        prop_yes = n_yes / n_total,
        prop_no  = n_no / n_total,
        .groups = "drop"
      )
  }
  
  tbl <- dplyr::bind_rows(med_rows, any_rows) %>%
    dplyr::mutate(
      dataset   = factor(dataset, levels = dataset_order),
      Medication = factor(Medication, levels = c(meds_in, "Any medication"))
    ) %>%
    dplyr::arrange(dataset, Medication) %>%
    dplyr::mutate(
      prop_yes = round(100 * prop_yes, digits),
      prop_no  = round(100 * prop_no,  digits),
      Yes  = sprintf("%d (%.1f%%)", n_yes, prop_yes),
      No   = sprintf("%d (%.1f%%)", n_no,  prop_no)
    ) %>%
    dplyr::select(dataset, Medication, n_total, n_yes, n_no, prop_yes, prop_no, Yes, No)
  
  if (!is.null(out_csv)) {
    dir.create(dirname(out_csv), showWarnings = FALSE, recursive = TRUE)
    readr::write_csv(tbl, out_csv)
  }
  tbl
}



############# medication trajectories with age
# ---------- helpers ----------
available_meds <- function(df, med_cols) intersect(med_cols, names(df))

normalize_meds <- function(df, meds) {
  dplyr::mutate(df, dplyr::across(dplyr::all_of(meds), ~ as.numeric(replace(., is.na(.), 0))))
}

tag_dataset <- function(df, label, meds) {
  stopifnot("id" %in% names(df))
  meds_in <- available_meds(df, meds)
  df |>
    dplyr::select(dplyr::any_of(c("id","age")), dplyr::all_of(meds_in)) |>
    dplyr::mutate(dataset = label) |>
    dplyr::relocate(dataset, id) |>
    normalize_meds(meds_in)
}

tag_dataset_long <- function(df, label, med_cols) {
  stopifnot(all(c("id","age") %in% names(df)))
  meds_in <- available_meds(df, med_cols)
  df %>%
    dplyr::select(dplyr::any_of(c("id","age")), dplyr::all_of(meds_in)) %>%
    dplyr::mutate(dataset = label) %>%
    dplyr::relocate(dataset, id, age) %>%
    normalize_meds(meds_in)
}

# ---------- 1) Baseline-age smoother (wide -> long inside) ----------
plot_med_prevalence_by_age_smooth_multi <- function(
    NACC_df, AIBL_df, HABS_df,
    med_cols       = c("ACEi","ARB","BetaBlk","CCB","Diuretic","Statin","Metformin"),
    dataset_order  = c("NACC","AIBL","HABS"),
    include_any_med = TRUE,
    show_ci         = FALSE,
    text_size       = 8,
    line_size_any   = 1.1,
    line_size_med   = 0.8,
    palette         = NULL,      # named vector c(ACEi="#...", ..., "Any medication"="#...")
    add_panel_frame = TRUE,
    frame_color     = "black",
    frame_size      = 0.4,
    free_x          = TRUE       # << NEW: use free x-scale per cohort
){
  # combine labeled baselines
  aibl_b <- tag_dataset(AIBL_df, "AIBL", med_cols)
  habs_b <- tag_dataset(HABS_df, "HABS", med_cols)
  nacc_b <- tag_dataset(NACC_df, "NACC", med_cols)
  df <- dplyr::bind_rows(aibl_b, habs_b, nacc_b)
  if (!"age" %in% names(df)) stop("Age column not found. Please include an 'age' column in each df.")
  
  meds_in <- available_meds(df, med_cols)
  
  # long (7 meds)
  long_meds <- tidyr::pivot_longer(
    df, tidyselect::all_of(meds_in),
    names_to = "Medication", values_to = "has_med"
  )
  
  # append Any medication if requested
  if (isTRUE(include_any_med)) {
    any_df <- df |>
      dplyr::mutate(
        has_med = as.integer(rowSums(dplyr::across(tidyselect::all_of(meds_in))) > 0),
        Medication = "Any medication"
      )
    long_all <- dplyr::bind_rows(
      long_meds[, c("dataset","id","age","Medication","has_med")],
      any_df[,   c("dataset","id","age","Medication","has_med")]
    )
  } else {
    long_all <- long_meds
  }
  
  # factor levels
  med_levels <- unique(c(meds_in, if (include_any_med) "Any medication"))
  long_all <- long_all |>
    dplyr::filter(!is.na(age)) |>
    dplyr::mutate(
      dataset    = factor(dataset, levels = dataset_order),
      Medication = factor(Medication, levels = med_levels)
    )
  
  # colors
  if (is.null(palette)) {
    base_cols <- RColorBrewer::brewer.pal(8, "Dark2")
    # recycle if needed
    palette <- setNames(rep(base_cols, length.out = length(med_levels)), med_levels)
  } else {
    stopifnot(!is.null(names(palette)))
    palette <- palette[med_levels]
  }
  
  # linewidth mapping (avoids ggplot linewidth error)
  lw_values <- setNames(rep(line_size_med, length(med_levels)), med_levels)
  if ("Any medication" %in% names(lw_values)) lw_values["Any medication"] <- line_size_any
  
  p <- ggplot2::ggplot(
    long_all,
    ggplot2::aes(x = age, y = has_med, color = Medication, linewidth = Medication, group = Medication)
  ) +
    ggplot2::geom_smooth(
      method = "glm", method.args = list(family = "binomial"),
      se = show_ci
    ) +
    ggplot2::scale_color_manual(values = palette, drop = FALSE) +
    ggplot2::scale_linewidth_manual(values = lw_values, guide = "none") +
    ggplot2::scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0,1)) +
    ggplot2::labs(
      x = "Age (years)", y = "Prevalence",
      title = " "
    ) +
    (if (free_x)
      ggplot2::facet_grid(. ~ dataset, scales = "free_x", space = "fixed")
     else
       ggplot2::facet_wrap(~ dataset, nrow = 1)
    ) +
    ggplot2::theme_minimal(base_size = text_size) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "bold", size = text_size),
      legend.position = "bottom",
      legend.title = ggplot2::element_blank(),
      panel.border = if (add_panel_frame)
        ggplot2::element_rect(color = frame_color, fill = NA, linewidth = frame_size)
      else ggplot2::element_blank(),
      panel.spacing = grid::unit(6, "pt")
    )
  
  return(p)
}

# ---------- 2) Attained-age (visit-level) smoother ----------
plot_attained_age_prevalence_multi <- function(
    NACC_long, AIBL_long, HABS_long,
    med_cols       = c("ACEi","ARB","BetaBlk","CCB","Diuretic","Statin","Metformin"),
    dataset_order  = c("AIBL","HABS","NACC"),
    include_any_med = TRUE,
    show_ci         = FALSE,
    text_size       = 8,
    line_size_any   = 1.1,
    line_size_med   = 0.8,
    palette         = NULL,
    add_panel_frame = TRUE,
    frame_color     = "black",
    frame_size      = 0.4,
    free_x          = TRUE
){
  # stack visit-level with tags
  aibl_v <- tag_dataset_long(AIBL_long, "AIBL", med_cols)
  habs_v <- tag_dataset_long(HABS_long, "HABS", med_cols)
  nacc_v <- tag_dataset_long(NACC_long, "NACC", med_cols)
  df <- dplyr::bind_rows(aibl_v, habs_v, nacc_v)
  
  meds_in <- available_meds(df, med_cols)
  
  # long form for 7 meds
  long_meds <- tidyr::pivot_longer(
    df, tidyselect::all_of(meds_in),
    names_to = "Medication", values_to = "has_med"
  )
  
  # append Any if requested
  if (isTRUE(include_any_med)) {
    any_df <- df %>%
      dplyr::mutate(
        has_med = as.integer(rowSums(dplyr::across(dplyr::all_of(meds_in))) > 0),
        Medication = "Any medication"
      )
    long_all <- dplyr::bind_rows(
      long_meds[, c("dataset","id","age","Medication","has_med")],
      any_df[,   c("dataset","id","age","Medication","has_med")]
    )
  } else {
    long_all <- long_meds
  }
  
  # factors
  med_levels <- unique(c(meds_in, if (include_any_med) "Any medication"))
  long_all <- long_all %>%
    dplyr::filter(!is.na(age)) %>%
    dplyr::mutate(
      dataset    = factor(dataset, levels = dataset_order),
      Medication = factor(Medication, levels = med_levels)
    )
  
  # colors
  if (is.null(palette)) {
    base_cols <- RColorBrewer::brewer.pal(8, "Dark2")
    palette <- setNames(rep(base_cols, length.out = length(med_levels)), med_levels)
  } else {
    stopifnot(!is.null(names(palette)))
    palette <- palette[med_levels]
  }
  
  # linewidth mapping
  lw_values <- setNames(rep(line_size_med, length(med_levels)), med_levels)
  if ("Any medication" %in% names(lw_values)) lw_values["Any medication"] <- line_size_any
  
  p <- ggplot2::ggplot(
    long_all,
    ggplot2::aes(x = age, y = has_med, color = Medication, linewidth = Medication, group = Medication)
  ) +
    ggplot2::geom_smooth(
      method = "glm", method.args = list(family = "binomial"),
      se = show_ci
    ) +
    ggplot2::scale_color_manual(values = palette, drop = FALSE) +
    ggplot2::scale_linewidth_manual(values = lw_values, guide = "none") +
    ggplot2::scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0,1)) +
    ggplot2::labs(
      x = "Age (years)", y = "Prevalence",
      title = " "
    ) +
    (if (free_x)
      ggplot2::facet_grid(. ~ dataset, scales = "free_x", space = "fixed")
     else
       ggplot2::facet_wrap(~ dataset, nrow = 1)
    ) +
    ggplot2::theme_minimal(base_size = text_size) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "bold", size = text_size),
      legend.position = "bottom",
      legend.title = ggplot2::element_blank(),
      panel.border = if (add_panel_frame)
        ggplot2::element_rect(color = frame_color, fill = NA, linewidth = frame_size)
      else ggplot2::element_blank(),
      panel.spacing = grid::unit(6, "pt")
    )
  
  return(p)
}

############# medication correlation plots
# =========================
# Correlation heatmaps by cohort (7 meds) — brace-safe version
# =========================
# suppressPackageStartupMessages({
#   library(dplyr)
#   library(tidyr)
#   library(ggplot2)
#   library(purrr)
# })
# 
# # (optional) tetrachoric helper
# .safe_tetrachoric <- function(mat) {
#   if (!requireNamespace("psych", quietly = TRUE)) return(NULL)
#   out <- try(psych::tetrachoric(mat), silent = TRUE)
#   if (inherits(out, "try-error")) return(NULL)
#   out$rho
# }
# 
# # Pearson (φ on 0/1) or tetrachoric, plus p-values for Pearson
# .cor_with_p <- function(df, method = c("pearson","tetrachoric"),
#                         use = "pairwise.complete.obs") {
#   method <- match.arg(method)
#   mat <- as.matrix(df)
#   
#   if (method == "tetrachoric") {
#     rho <- .safe_tetrachoric(mat)
#     if (!is.null(rho)) {
#       pmat <- matrix(NA_real_, nrow = ncol(mat), ncol = ncol(mat),
#                      dimnames = list(colnames(mat), colnames(mat)))
#       return(list(r = rho, p = pmat))
#     }
#     method <- "pearson"  # fallback
#   }
#   
#   r <- suppressWarnings(stats::cor(mat, use = use, method = "pearson"))
#   p <- matrix(NA_real_, nrow = ncol(mat), ncol = ncol(mat),
#               dimnames = list(colnames(mat), colnames(mat)))
#   for (i in seq_len(ncol(mat))) {
#     for (j in seq_len(ncol(mat))) {
#       if (i == j) {
#         p[i, j] <- 0
#       } else if (i < j) {
#         ct <- suppressWarnings(stats::cor.test(mat[, i], mat[, j]))
#         p[i, j] <- ct$p.value
#         p[j, i] <- ct$p.value
#       }
#     }
#   }
#   list(r = r, p = p)
# }
# 
# # Tidy correlation matrix
# .tidy_corr <- function(r_mat, p_mat = NULL, med_order) {
#   r_long <- as.data.frame(r_mat) |>
#     tibble::rownames_to_column("Var1") |>
#     tidyr::pivot_longer(-Var1, names_to = "Var2", values_to = "r") |>
#     mutate(Var1 = factor(Var1, levels = med_order),
#            Var2 = factor(Var2, levels = med_order))
#   
#   if (is.null(p_mat)) return(r_long)
#   
#   p_long <- as.data.frame(p_mat) |>
#     tibble::rownames_to_column("Var1") |>
#     tidyr::pivot_longer(-Var1, names_to = "Var2", values_to = "p")
#   
#   left_join(r_long, p_long, by = c("Var1","Var2"))
# }
# 
# # Main plotting function
# plot_med_corr_by_cohort <- function(
#     NACC_df, AIBL_df, HABS_df,
#     med_cols = c("ACEi","ARB","BetaBlk","CCB","Diuretic","Metformin","Statin"),
#     dataset_order = c("NACC","AIBL","HABS"),
#     method = c("pearson","tetrachoric"),
#     mask_upper_triangle = FALSE,
#     annotate_values = TRUE,
#     digits = 2,
#     breaks = seq(-1, 1, by = 0.2),
#     limits = c(-1, 1),
#     out_pdf = NULL, width = 10, height = 4.2, dpi = 300
# ) {
#   method <- match.arg(method)
#   
#   # Ensure meds exist and are 0/1
#   make_binary_subset <- function(df, nm) {
#     keep <- intersect(nm, names(df))
#     if (length(keep) < length(nm)) {
#       missing <- setdiff(nm, keep)
#       warning("Missing medication columns dropped: ", paste(missing, collapse = ", "))
#     }
#     out <- df[, keep, drop = FALSE]
#     out |>
#       mutate(across(everything(), ~ {
#         v <- .
#         if (is.logical(v)) as.integer(v)
#         else if (is.numeric(v)) as.integer(v != 0)
#         else as.integer(tolower(as.character(v)) %in% c("1","yes","y","true","t"))
#       }))
#   }
#   
#   cohorts <- list(
#     "NACC" = make_binary_subset(NACC_df, med_cols),
#     "AIBL" = make_binary_subset(AIBL_df, med_cols),
#     "HABS" = make_binary_subset(HABS_df, med_cols)
#   )
#   
#   res_list <- purrr::imap(cohorts, function(dat, nm) {
#     cp <- .cor_with_p(dat, method = method)
#     .tidy_corr(cp$r, cp$p, med_order = med_cols) |>
#       mutate(dataset = nm)
#   })
#   
#   plot_df <- bind_rows(res_list) |>
#     mutate(dataset = factor(dataset, levels = dataset_order))
#   
#   if (mask_upper_triangle) {
#     plot_df <- plot_df |>
#       group_by(dataset) |>
#       mutate(i = as.integer(Var1), j = as.integer(Var2)) |>
#       ungroup() |>
#       filter(i >= j) |>
#       select(-i, -j)
#   }
#   
#   plot_df <- plot_df |>
#     mutate(
#       stars = dplyr::case_when(
#         is.na(p) ~ "",
#         p < 0.001 ~ "***",
#         p < 0.01  ~ "**",
#         p < 0.05  ~ "*",
#         p < 0.1   ~ "†",
#         TRUE      ~ ""
#       ),
#       r_lab = sprintf(paste0("%.", digits, "f"), r),
#       lab = ifelse(annotate_values, ifelse(is.na(r), "", r_lab), "")
#     )
#   
#   # Build base heatmap (no braces in the chain)
#   plt <- ggplot(plot_df, aes(x = Var2, y = Var1, fill = r)) +
#     geom_tile(color = "white", linewidth = 0.3) +
#     scale_fill_gradientn(
#       colours = c("#313695","#4575B4","#74ADD1","#ABD9E9","#E0F3F8",
#                   "#FFFFFF",
#                   "#FEE090","#FDAE61","#F46D43","#D73027","#A50026"),
#       limits = limits, breaks = breaks, oob = scales::squish
#     ) +
#     coord_fixed() +
#     facet_wrap(~ dataset, nrow = 1) +
#     labs(x = NULL, y = NULL,
#          fill = ifelse(method == "tetrachoric", "ρ (tetrachoric)", "r (φ)"),
#          title = "Medication co-usage correlations by cohort",
#          subtitle = ifelse(method == "tetrachoric",
#                            "Tetrachoric correlation (fallback to φ if psych unavailable)",
#                            "Pearson correlation on 0/1 indicators (φ)")) +
#     theme_minimal(base_size = 10) +
#     theme(
#       panel.grid = element_blank(),
#       axis.text.x = element_text(angle = 45, hjust = 1),
#       strip.text = element_text(face = "bold"),
#       plot.title = element_text(face = "bold")
#     )
#   
#   # Conditionally add text layer (no `{}` inside pipeline)
#   if (annotate_values) {
#     plt <- plt + geom_text(aes(label = lab), size = 2.8)
#   }
#   
#   if (!is.null(out_pdf)) {
#     ggsave(out_pdf, plt, width = width, height = height, dpi = dpi)  # device inferred by extension
#     message("Saved to: ", out_pdf)
#   }
#   plt
# }


# =========================
# Correlation heatmaps by cohort (7 meds)
# - No text annotations on tiles
# - Uses scale_fill_gradient2() for colours
# =========================
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(purrr)
})

# (optional) tetrachoric helper
.safe_tetrachoric <- function(mat) {
  if (!requireNamespace("psych", quietly = TRUE)) return(NULL)
  out <- try(psych::tetrachoric(mat), silent = TRUE)
  if (inherits(out, "try-error")) return(NULL)
  out$rho
}

# Pearson (φ on 0/1) or tetrachoric; p-values only for Pearson
.cor_with_p <- function(df, method = c("pearson","tetrachoric"),
                        use = "pairwise.complete.obs") {
  method <- match.arg(method)
  mat <- as.matrix(df)
  
  if (method == "tetrachoric") {
    rho <- .safe_tetrachoric(mat)
    if (!is.null(rho)) {
      pmat <- matrix(NA_real_, nrow = ncol(mat), ncol = ncol(mat),
                     dimnames = list(colnames(mat), colnames(mat)))
      return(list(r = rho, p = pmat))
    }
    method <- "pearson"  # fallback
  }
  
  r <- suppressWarnings(stats::cor(mat, use = use, method = "pearson"))
  p <- matrix(NA_real_, nrow = ncol(mat), ncol = ncol(mat),
              dimnames = list(colnames(mat), colnames(mat)))
  for (i in seq_len(ncol(mat))) {
    for (j in seq_len(ncol(mat))) {
      if (i == j) {
        p[i, j] <- 0
      } else if (i < j) {
        ct <- suppressWarnings(stats::cor.test(mat[, i], mat[, j]))
        p[i, j] <- ct$p.value
        p[j, i] <- ct$p.value
      }
    }
  }
  list(r = r, p = p)
}

# Tidy a square matrix to long format (keeps a fixed med order)
.tidy_corr <- function(r_mat, p_mat = NULL, med_order) {
  r_long <- as.data.frame(r_mat) |>
    tibble::rownames_to_column("Var1") |>
    tidyr::pivot_longer(-Var1, names_to = "Var2", values_to = "r") |>
    mutate(Var1 = factor(Var1, levels = med_order),
           Var2 = factor(Var2, levels = med_order))
  
  if (is.null(p_mat)) return(r_long)
  
  p_long <- as.data.frame(p_mat) |>
    tibble::rownames_to_column("Var1") |>
    tidyr::pivot_longer(-Var1, names_to = "Var2", values_to = "p")
  
  dplyr::left_join(r_long, p_long, by = c("Var1","Var2"))
}

# -------------------------
# Main plotting function
# -------------------------
plot_med_corr_by_cohort <- function(
    NACC_df, AIBL_df, HABS_df,
    med_cols = c("ACEi","ARB","BetaBlk","CCB","Diuretic","Metformin","Statin"),
    dataset_order = c("NACC","AIBL","HABS"),
    method = c("pearson","tetrachoric"),
    mask_upper_triangle = FALSE,
    # colour control (scale_fill_gradient2)
    low = "#4575B4", mid = "white", high = "#D73027", midpoint = 0,
    limits = c(-1, 1), breaks = seq(-1, 1, by = 0.2),
    out_pdf = NULL, width = 10, height = 4.2, dpi = 300
) {
  method <- match.arg(method)
  
  # ensure meds exist and are 0/1
  make_binary_subset <- function(df, nm) {
    keep <- intersect(nm, names(df))
    if (length(keep) < length(nm)) {
      missing <- setdiff(nm, keep)
      warning("Missing medication columns dropped: ", paste(missing, collapse = ", "))
    }
    out <- df[, keep, drop = FALSE]
    out |>
      mutate(across(everything(), ~ {
        v <- .
        if (is.logical(v)) as.integer(v)
        else if (is.numeric(v)) as.integer(v != 0)
        else as.integer(tolower(as.character(v)) %in% c("1","yes","y","true","t"))
      }))
  }
  
  cohorts <- list(
    "NACC" = make_binary_subset(NACC_df, med_cols),
    "AIBL" = make_binary_subset(AIBL_df, med_cols),
    "HABS" = make_binary_subset(HABS_df, med_cols)
  )
  
  res_list <- purrr::imap(cohorts, function(dat, nm) {
    cp <- .cor_with_p(dat, method = method)
    .tidy_corr(cp$r, cp$p, med_order = med_cols) |>
      mutate(dataset = nm)
  })
  
  plot_df <- bind_rows(res_list) |>
    mutate(dataset = factor(dataset, levels = dataset_order))
  
  if (mask_upper_triangle) {
    plot_df <- plot_df |>
      group_by(dataset) |>
      mutate(i = as.integer(Var1), j = as.integer(Var2)) |>
      ungroup() |>
      filter(i >= j) |>
      select(-i, -j)
  }
  
  # Build heatmap (NO geom_text -> no values on tiles)
  plt <- ggplot(plot_df, aes(x = Var2, y = Var1, fill = r)) +
    geom_tile(color = "white", linewidth = 0.3) +
    scale_fill_gradient2(
      low = low, mid = mid, high = high, midpoint = midpoint,
      limits = limits, breaks = breaks
    ) +
    coord_fixed() +
    facet_wrap(~ dataset, nrow = 1) +
    labs(
      x = NULL, y = NULL,
      fill = ifelse(method == "tetrachoric", "ρ (tetrachoric)", "r"),
      title = " "
    ) +
    theme_minimal(base_size = 10) +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1),
      strip.text = element_text(face = "bold"),
      plot.title = element_text(face = "bold")
    )
  
  if (!is.null(out_pdf)) {
    ggsave(out_pdf, plt, width = width, height = height, dpi = dpi)  # device from extension
    message("Saved to: ", out_pdf)
  }
  plt
}
