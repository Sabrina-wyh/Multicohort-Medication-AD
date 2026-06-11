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
library(meta)
library(metafor)
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
    dataset_order= c("NACC","AIBL","HABS-HD"),
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
meta_one <- \(est, se, method="REML"){ 
  ok <- is.finite(est)&is.finite(se)
  if(sum(ok)<2) return(list(b=NA,lo=NA,hi=NA,w=rep(NA,length(est))))
  fit <- rma.uni(yi=est[ok], sei=se[ok], method=method)
  w <- rep(NA_real_, length(est))
  w[ok] <- weights(fit)/sum(weights(fit))*100
  list(b=as.numeric(fit$b), lo=fit$ci.lb, hi=fit$ci.ub, w=w, I2=fit$I2, tau2=fit$tau2) 
  }


# make_dt_from_fits <- function(
#     fits,
#     terms,
#     cohorts = c("NACC","AIBL","HABS"),
#     term_labels = NULL,
#     method = "REML",          # "REML" (random-effects) or "FE" (common-effect)
#     digits = 3,
#     add_pred_int = TRUE,
#     p_adjust = c("none","BH") # NEW: control p-value adjustment
# ){
#   p_adjust <- match.arg(p_adjust)
#   
#   get_sig_stars <- function(p){
#     s <- rep("", length(p))
#     s[is.finite(p) & p < 0.001] <- "***"
#     s[is.finite(p) & p >= 0.001 & p < 0.01] <- "**"
#     s[is.finite(p) & p >= 0.01  & p < 0.05] <- "*"
#     s[is.finite(p) & p >= 0.05  & p <= 0.10] <- "^"
#     s
#   }
#   
#   meta_one <- function(est, se, method="REML", add_pred_int=TRUE){
#     ok <- is.finite(est) & is.finite(se)
#     if (sum(ok) < 2L){
#       return(list(b=NA_real_, se=NA_real_, lo=NA_real_, hi=NA_real_, p=NA_real_,
#                   w=rep(NA_real_, length(est)), I2=NA_real_, tau2=NA_real_,
#                   pi_lo=NA_real_, pi_hi=NA_real_))
#     }
#     fit <- metafor::rma.uni(yi = est[ok], sei = se[ok], method = method)
#     ww  <- stats::weights(fit); w <- rep(NA_real_, length(est)); w[ok] <- ww / sum(ww) * 100
#     pi_lo <- pi_hi <- NA_real_
#     if (add_pred_int){
#       pr <- tryCatch(stats::predict(fit), error = function(e) NULL)
#       if (!is.null(pr) && all(c("pi.lb","pi.ub") %in% names(pr))) {
#         pi_lo <- pr$pi.lb; pi_hi <- pr$pi.ub
#       }
#     }
#     list(b=as.numeric(fit$b), se=as.numeric(fit$se), lo=fit$ci.lb, hi=fit$ci.ub, p=fit$pval,
#          w=w, I2=fit$I2, tau2=fit$tau2, pi_lo=pi_lo, pi_hi=pi_hi)
#   }
#   
#   groups <- names(fits)
#   stopifnot(length(groups) > 0)
#   stopifnot(all(vapply(fits, length, 1L) == length(cohorts)))
#   if (is.null(term_labels)) term_labels <- stats::setNames(terms, terms)
#   stopifnot(all(terms %in% names(term_labels)))
#   clean_labels <- paste0("     ", trimws(unname(term_labels[terms])))
#   
#   out <- data.frame(Measurements = clean_labels, stringsAsFactors = FALSE, check.names = FALSE)
#   
#   for (grp in groups){
#     grp_fits <- fits[[grp]]
#     
#     est_mat <- matrix(NA_real_, nrow = length(terms), ncol = length(cohorts),
#                       dimnames = list(terms, cohorts))
#     se_mat  <- est_mat
#     pv_mat  <- est_mat
#     
#     for (i in seq_along(cohorts)){
#       nm <- paste0(grp, "_", cohorts[i])
#       co <- tryCatch(summary(grp_fits[[i]])$coefficients, error=function(e) NULL)
#       
#       est <- se <- low <- hi <- rep(NA_real_, length(terms))
#       names(est) <- names(se) <- names(low) <- names(hi) <- terms
#       
#       if (!is.null(co) && is.matrix(co)){
#         available <- intersect(terms, rownames(co))
#         if (length(available) > 0){
#           est[available] <- co[available, "Estimate"]
#           se [available] <- co[available, "Std. Error"]
#           pv_col <- intersect(colnames(co), c("Pr(>|t|)","Pr(>|z|)","Pr(>|W|)","Pr(>Chisq)"))
#           if (length(pv_col)) pv_mat[available, cohorts[i]] <- co[available, pv_col[1L]]
#         }
#       }
#       ok <- is.finite(est) & is.finite(se)
#       low[ok] <- est[ok] - 1.96 * se[ok]
#       hi [ok] <- est[ok] + 1.96 * se[ok]
#       
#       out[[paste0(nm, "_est")]] <- est
#       out[[paste0(nm, "_low")]] <- low
#       out[[paste0(nm, "_hi") ]] <- hi
#       
#       est_mat[, cohorts[i]] <- est
#       se_mat [, cohorts[i]] <- se
#     }
#     
#     ## --- NEW: adjust cohort p-values column-wise with BH, if requested ---
#     if (p_adjust == "BH") {
#       for (k in seq_along(cohorts)){
#         colp <- pv_mat[, k]
#         if (any(is.finite(colp))) {
#           pv_mat[, k] <- stats::p.adjust(colp, method = "BH")
#         }
#       }
#     }
#     
#     star_mat <- matrix("", nrow = length(terms), ncol = length(cohorts),
#                        dimnames = list(terms, cohorts))
#     has_p <- is.finite(pv_mat)
#     if (any(has_p)) star_mat[has_p] <- get_sig_stars(pv_mat[has_p])
#     
#     estse_txt <- vapply(seq_along(terms), function(j){
#       paste(vapply(seq_along(cohorts), function(k){
#         e <- est_mat[j,k]; s <- se_mat[j,k]; st <- star_mat[j,k]
#         if (is.finite(e) && is.finite(s)) sprintf("%+.*f (%.3f)%s", digits, e, s, st) else "NA"
#       }, character(1)), collapse = "\n")
#     }, character(1))
#     out[[paste0(grp, "_Est_SE")]] <- estse_txt
#     out[[grp]] <- paste(rep(" ", 20), collapse = " ")
#     
#     meta_list <- lapply(seq_len(nrow(est_mat)), function(j)
#       meta_one(est_mat[j,], se_mat[j,], method = method, add_pred_int = add_pred_int))
#     
#     out[[paste0(grp, "_META_est")]] <- vapply(meta_list, `[[`, numeric(1), "b")
#     out[[paste0(grp, "_META_low")]] <- vapply(meta_list, `[[`, numeric(1), "lo")
#     out[[paste0(grp, "_META_hi") ]] <- vapply(meta_list, `[[`, numeric(1), "hi")
#     out[[paste0(grp, "_META_se") ]] <- vapply(meta_list, `[[`, numeric(1), "se")
#     
#     meta_p_vec <- vapply(meta_list, `[[`, numeric(1), "p")
#     
#     ## --- NEW: adjust META p-values across terms with BH, if requested ---
#     if (p_adjust == "BH" && any(is.finite(meta_p_vec))) {
#       meta_p_vec <- stats::p.adjust(meta_p_vec, method = "BH")
#     }
#     out[[paste0(grp, "_META_p")  ]] <- meta_p_vec
#     
#     out[[paste0(grp, "_META_I2") ]] <- vapply(meta_list, `[[`, numeric(1), "I2")
#     out[[paste0(grp, "_META_tau2")]] <- vapply(meta_list, `[[`, numeric(1), "tau2")
#     if (add_pred_int){
#       out[[paste0(grp, "_META_pi_low")]] <- vapply(meta_list, `[[`, numeric(1), "pi_lo")
#       out[[paste0(grp, "_META_pi_hi") ]] <- vapply(meta_list, `[[`, numeric(1), "pi_hi")
#     }
#     
#     meta_stars <- get_sig_stars(out[[paste0(grp, "_META_p")]])
#     out[[paste0(grp, "_META_Est_SE")]] <-
#       ifelse(is.finite(out[[paste0(grp, "_META_est")]]) & is.finite(out[[paste0(grp, "_META_se")]]),
#              sprintf("%+.*f (%.3f)%s",
#                      digits, out[[paste0(grp, "_META_est")]], out[[paste0(grp, "_META_se")]], meta_stars),
#              "NA")
#     
#     for (i in seq_along(cohorts)){
#       w_num <- vapply(meta_list, function(x) x$w[i], NA_real_)
#       out[[paste0(grp, "_", cohorts[i], "_w") ]]  <- w_num
#       out[[paste0(grp, "_", cohorts[i], "_W%")]] <- ifelse(is.finite(w_num), sprintf("%.1f%%", w_num), "NA")
#     }
#   }
#   
#   ## Optional: if you also want the adjusted cohort p's exported, uncomment below:
#   # for (k in seq_along(cohorts)){
#   #   out[[paste0("AdjP_", cohorts[k])]] <- pv_mat[, k]
#   # }
#   
#   out
# }


make_dt_from_fits <- function(
    fits,
    terms,
    cohorts = c("NACC","AIBL","HABS-HD"),
    term_labels = NULL,
    method = "REML",                # "REML" (random-effects) or "FE" (common-effect)
    digits = 3,
    add_pred_int = TRUE,
    p_adjust = c("none","BH")       # controls which p's drive stars and META_p
){
  p_adjust <- match.arg(p_adjust)
  
  get_sig_stars <- function(p){
    s <- rep("", length(p))
    s[is.finite(p) & p < 0.001] <- "***"
    s[is.finite(p) & p >= 0.001 & p < 0.01] <- "**"
    s[is.finite(p) & p >= 0.01  & p < 0.05] <- "*"
    s[is.finite(p) & p >= 0.05  & p <= 0.10] <- "^"
    s
  }
  
  meta_one <- function(est, se, method="REML", add_pred_int=TRUE){
    ok <- is.finite(est) & is.finite(se)
    if (sum(ok) < 2L){
      return(list(
        b=NA_real_, se=NA_real_, lo=NA_real_, hi=NA_real_, p=NA_real_,
        w=rep(NA_real_, length(est)), I2=NA_real_, tau2=NA_real_,
        pi_lo=NA_real_, pi_hi=NA_real_
      ))
    }
    fit <- metafor::rma.uni(yi = est[ok], sei = se[ok], method = method)
    ww  <- stats::weights(fit); w <- rep(NA_real_, length(est)); w[ok] <- ww / sum(ww) * 100
    pi_lo <- pi_hi <- NA_real_
    if (add_pred_int){
      pr <- tryCatch(stats::predict(fit), error=function(e) NULL)
      if (!is.null(pr) && all(c("pi.lb","pi.ub") %in% names(pr))) {
        pi_lo <- pr$pi.lb; pi_hi <- pr$pi.ub
      }
    }
    list(b=as.numeric(fit$b), se=as.numeric(fit$se), lo=fit$ci.lb, hi=fit$ci.ub, p=fit$pval,
         w=w, I2=fit$I2, tau2=fit$tau2, pi_lo=pi_lo, pi_hi=pi_hi)
  }
  
  groups <- names(fits)
  stopifnot(length(groups) > 0)
  stopifnot(all(vapply(fits, length, 1L) == length(cohorts)))
  if (is.null(term_labels)) term_labels <- stats::setNames(terms, terms)
  stopifnot(all(terms %in% names(term_labels)))
  clean_labels <- paste0("     ", trimws(unname(term_labels[terms])))
  
  out <- data.frame(Measurements = clean_labels, stringsAsFactors = FALSE, check.names = FALSE)
  
  for (grp in groups){
    grp_fits <- fits[[grp]]
    
    est_mat <- matrix(NA_real_, nrow = length(terms), ncol = length(cohorts),
                      dimnames = list(terms, cohorts))
    se_mat  <- est_mat
    pv_mat  <- est_mat
    
    for (i in seq_along(cohorts)){
      nm <- paste0(grp, "_", cohorts[i])
      co <- tryCatch(summary(grp_fits[[i]])$coefficients, error=function(e) NULL)
      
      est <- se <- low <- hi <- rep(NA_real_, length(terms))
      names(est) <- names(se) <- names(low) <- names(hi) <- terms
      
      if (!is.null(co) && is.matrix(co)){
        available <- intersect(terms, rownames(co))
        if (length(available) > 0){
          est[available] <- co[available, "Estimate"]
          se [available] <- co[available, "Std. Error"]
          pv_col <- intersect(colnames(co), c("Pr(>|t|)","Pr(>|z|)","Pr(>|W|)","Pr(>Chisq)"))
          if (length(pv_col)) pv_mat[available, cohorts[i]] <- co[available, pv_col[1L]]
        }
      }
      ok <- is.finite(est) & is.finite(se)
      low[ok] <- est[ok] - 1.96 * se[ok]
      hi [ok] <- est[ok] + 1.96 * se[ok]
      
      out[[paste0(nm, "_est")]] <- est
      out[[paste0(nm, "_low")]] <- low
      out[[paste0(nm, "_hi") ]] <- hi
      
      est_mat[, cohorts[i]] <- est
      se_mat [, cohorts[i]] <- se
    }
    
    ## raw & BH-adjusted cohort p's
    pv_mat_raw <- pv_mat
    pv_mat_adj <- pv_mat
    for (k in seq_along(cohorts)){
      colp <- pv_mat_raw[, k]
      pv_mat_adj[, k] <- if (any(is.finite(colp))) stats::p.adjust(colp, method = "BH") else colp
    }
    pv_use <- if (p_adjust == "BH") pv_mat_adj else pv_mat_raw
    
    ## stars from chosen p's
    star_mat <- matrix("", nrow = length(terms), ncol = length(cohorts),
                       dimnames = list(terms, cohorts))
    has_p <- is.finite(pv_use)
    if (any(has_p)) star_mat[has_p] <- get_sig_stars(pv_use[has_p])
    
    ## compact text "est (se){stars}" stacked by cohort
    estse_txt <- vapply(seq_along(terms), function(j){
      paste(vapply(seq_along(cohorts), function(k){
        e <- est_mat[j,k]; s <- se_mat[j,k]; st <- star_mat[j,k]
        if (is.finite(e) && is.finite(s)) sprintf("%+.*f (%.3f)%s", digits, e, s, st) else "NA"
      }, character(1)), collapse = "\n")
    }, character(1))
    out[[paste0(grp, "_Est_SE")]] <- estse_txt
    out[[grp]] <- paste(rep(" ", 20), collapse = " ")
    
    ## meta per-term
    meta_list <- lapply(seq_len(nrow(est_mat)), function(j)
      meta_one(est_mat[j,], se_mat[j,], method = method, add_pred_int = add_pred_int))
    
    out[[paste0(grp, "_META_est")]] <- vapply(meta_list, `[[`, numeric(1), "b")
    out[[paste0(grp, "_META_low")]] <- vapply(meta_list, `[[`, numeric(1), "lo")
    out[[paste0(grp, "_META_hi") ]] <- vapply(meta_list, `[[`, numeric(1), "hi")
    out[[paste0(grp, "_META_se") ]] <- vapply(meta_list, `[[`, numeric(1), "se")
    
    ## META p's (raw & BH-adjusted) + chosen p
    meta_p_raw <- vapply(meta_list, `[[`, numeric(1), "p")
    meta_p_adj <- if (any(is.finite(meta_p_raw))) stats::p.adjust(meta_p_raw, method = "BH") else meta_p_raw
    meta_p_vec <- if (p_adjust == "BH") meta_p_adj else meta_p_raw
    out[[paste0(grp, "_META_p_raw")]] <- meta_p_raw
    out[[paste0(grp, "_META_p_adj")]] <- meta_p_adj
    out[[paste0(grp, "_META_p")    ]] <- meta_p_vec
    
    out[[paste0(grp, "_META_I2") ]] <- vapply(meta_list, `[[`, numeric(1), "I2")
    out[[paste0(grp, "_META_tau2")]] <- vapply(meta_list, `[[`, numeric(1), "tau2")
    if (add_pred_int){
      out[[paste0(grp, "_META_pi_low")]] <- vapply(meta_list, `[[`, numeric(1), "pi_lo")
      out[[paste0(grp, "_META_pi_hi") ]] <- vapply(meta_list, `[[`, numeric(1), "pi_hi")
    }
    
    ## META stars from chosen p
    meta_stars <- get_sig_stars(out[[paste0(grp, "_META_p")]])
    out[[paste0(grp, "_META_Est_SE")]] <-
      ifelse(is.finite(out[[paste0(grp, "_META_est")]]) & is.finite(out[[paste0(grp, "_META_se")]]),
             sprintf("%+.*f (%.3f)%s",
                     digits,
                     out[[paste0(grp, "_META_est")]],
                     out[[paste0(grp, "_META_se")]],
                     meta_stars),
             "NA")
    
    ## weights per cohort
    for (i in seq_along(cohorts)){
      w_num <- vapply(meta_list, function(x) x$w[i], NA_real_)
      out[[paste0(grp, "_", cohorts[i], "_w") ]]  <- w_num
      out[[paste0(grp, "_", cohorts[i], "_W%")]] <- ifelse(is.finite(w_num), sprintf("%.1f%%", w_num), "NA")
    }
    
    ## export cohort p's (raw / adj / chosen) per cohort
    for (k in seq_along(cohorts)) {
      co <- cohorts[k]
      out[[paste0(grp, "_", co, "_p_raw")]] <- pv_mat_raw[, k]
      out[[paste0(grp, "_", co, "_p_adj")]] <- pv_mat_adj[, k]
      out[[paste0(grp, "_", co, "_p")    ]] <- pv_use[, k]
    }
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
    row_padding_mm = 10, 
    # xlim_fixed = c(-0.5, 0.5),
    xlim_per_group = NULL,  # NEW: list or named list with xlims per group
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
  
  # ## X-axis (fixed)
  # xlim_use <- xlim_fixed
  ## NEW: Handle xlim_per_group
  # Convert to list format if needed
  if (is.null(xlim_per_group)) {
    # Default: same range for all groups
    xlim_use <- rep(list(c(-0.5, 0.5)), length(groups))
  } else if (is.list(xlim_per_group) && !is.null(names(xlim_per_group))) {
    # Named list: extract in order of groups
    xlim_use <- lapply(groups, function(g) {
      if (g %in% names(xlim_per_group)) xlim_per_group[[g]] else c(-0.5, 0.5)
    })
  } else if (is.list(xlim_per_group) && length(xlim_per_group) == length(groups)) {
    # Unnamed list: use as-is (assumed to match groups order)
    xlim_use <- xlim_per_group
  } else {
    stop("xlim_per_group must be NULL, a named list with group names, or a list of length equal to number of groups")
  }
  
  ## Column widths: [Measurements] + per-group (Est, CI, [spacer])
  per_group <- if (gutter) c(est_frac, ci_frac, gutter_width) else c(est_frac, ci_frac)
  col_w <- c(first_col_width, rep(per_group, length(groups)))
  
  tm <- forest_theme(
    base_size = 11,
    refline_gp = gpar(lty = "dotted", col = "grey70", lwd = 1.1),
    ci_pch = cohort_shapes, ci_col = cohort_colors, ci_lwd = line_width,
    legend_name = "Cohort", legend_value = cohorts,
    core = list(
      fg_params = list(hjust = rep(0, ncol(forest_dt)), x = rep(0.01, ncol(forest_dt)), lineheight = 1.05),
      padding = unit(c(row_padding_mm, row_padding_mm), "mm")
    ),
    colwidth = col_w
  )
  
  p <- forestploter::forest(
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
    cohorts = c("NACC","AIBL","HABS-HD"),
    term_labels = NULL,
    title = "",
    filename = NULL, width = NULL, height = NULL,
    # layout/style passthrough
    ref_line = 0,
    cohort_colors = c("#DAA87C","#92A5D1","#C9DCC4"),
    cohort_shapes = c(15,16,18),
    row_padding_mm = 10,
    # xlim_fixed = c(-0.5, 0.5),
    xlim_per_group = NULL,  # NEW: replaces xlim_fixed
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
  dt_ <- make_dt_from_fits(fits_subset, terms, cohorts = cohorts, term_labels = term_labels)
  
  create_multi_group_forest(
    dt_, groups = groups, cohorts = cohorts,
    ref_line = ref_line, title = title,
    cohort_colors = cohort_colors, cohort_shapes = cohort_shapes,
    filename = filename, width = width, height = height,
    row_padding_mm = row_padding_mm, 
    # xlim_fixed = xlim_fixed,
    xlim_per_group = xlim_per_group,  # NEW
    cohort_nudge = cohort_nudge, point_size = point_size, line_width = line_width,
    gutter = gutter, gutter_width = gutter_width,
    first_col_width = first_col_width, est_frac = est_frac, ci_frac = ci_frac
  )
}

plot_forest_single_cohort <- function(
    fits, terms, cohort = "NACC",
    groups  = c("Full","Female","Male","APOE4pos","APOE4neg"),
    term_labels = NULL, title = "",
    filename = NULL, width = NULL, height = NULL,
    ref_line = 0, cohort_color = "#92A5D1", cohort_shape = 16,
    row_padding_mm = 10, 
    # xlim_fixed = c(-0.5, 0.5),
    xlim_per_group = NULL,  # NEW: replaces xlim_fixed
    cohort_nudge = 0.35, point_size = 0.6, line_width = 1.6,
    gutter = TRUE, gutter_width = 0.04,
    first_col_width = 0.40, est_frac = 0.22, ci_frac = 0.16
){
  # reuse your existing pipeline, but restrict cohorts to the single choice
  plot_forest_subset(
    fits = fits, terms = terms, groups = groups,
    cohorts = c(cohort), term_labels = term_labels, title = title,
    filename = filename, width = width, height = height,
    ref_line = ref_line,
    cohort_colors = c(cohort_color), cohort_shapes = c(cohort_shape),
    row_padding_mm = row_padding_mm, xlim_fixed = xlim_fixed,
    cohort_nudge = cohort_nudge, point_size = point_size, line_width = line_width,
    gutter = gutter, gutter_width = gutter_width,
    first_col_width = first_col_width, est_frac = est_frac, ci_frac = ci_frac
  )
}

###### add for meta analysis
# 1) Expand each Measurement into 4 sub-rows (NACC, AIBL, HABS, Meta)
#    NO gap rows - spacing handled by row_padding_mm
expand_dt_with_meta_rows <- function(dt_, groups,
                                     cohorts = c("NACC","AIBL","HABS-HD"),
                                     meta_label = "Meta",
                                     add_block_gap = FALSE){
  K <- 4L
  order4 <- c(cohorts, meta_label)
  out <- vector("list", nrow(dt_) * K); k <- 0L
  
  for (i in seq_len(nrow(dt_))){
    ## 4 cohort/meta rows
    for (j in seq_along(order4)){
      co <- order4[j]; k <- k + 1L
      # Measurement label ONCE (first sub-row of the block)
      lab <- if (j == 1L) dt_$Measurements[i] else " "
      row <- list(Measurements = lab, .__cohort__ = co)
      
      for (g in groups){
        if (co == meta_label){
          e  <- dt_[[paste0(g,"_META_est")]][i]
          se <- dt_[[paste0(g,"_META_se")]][i]
          stars <- sub("^[^)]*\\)", "", dt_[[paste0(g,"_META_Est_SE")]][i])
          row[[paste0(g,"__EST")]] <- if (is.finite(e) && is.finite(se)) sprintf("%+.3f (%.3f)%s", e, se, stars) else "NA"
          row[[paste0(g,"__",co,"_W")]]   <- if (is.finite(e)) "100%" else "NA"
          row[[paste0(g,"__",co,"_est")]] <- e
          row[[paste0(g,"__",co,"_low")]] <- dt_[[paste0(g,"_META_low")]][i]
          row[[paste0(g,"__",co,"_hi") ]] <- dt_[[paste0(g,"_META_hi") ]][i]
        } else {
          e <- dt_[[paste0(g,"_",co,"_est")]][i]
          lines <- strsplit(dt_[[paste0(g,"_Est_SE")]][i], "\n", fixed = TRUE)[[1]]
          idx <- match(co, cohorts)
          row[[paste0(g,"__EST")]] <- if (!is.na(idx) && idx <= length(lines)) lines[idx] else if (is.finite(e)) sprintf("%+.3f (NA)", e) else "NA"
          row[[paste0(g,"__",co,"_W")]]   <- dt_[[paste0(g,"_",co,"_W%")]][i]
          row[[paste0(g,"__",co,"_est")]] <- e
          row[[paste0(g,"__",co,"_low")]] <- dt_[[paste0(g,"_",co,"_low")]][i]
          row[[paste0(g,"__",co,"_hi") ]] <- dt_[[paste0(g,"_",co,"_hi") ]][i]
        }
        row[[paste0(g,"__CI")]] <- " "
      }
      out[[k]] <- row
    }
  }
  
  df <- dplyr::bind_rows(out)
  rownames(df) <- NULL
  df
}

# 2) Build forest table: per-group = Est | Weight | CI | (Spacer)
#    Group name ONLY in Estimate column header
create_forest_table_expanded <- function(exp_dt, groups, add_gutter = TRUE, ci_spaces = 60){
  forest_dt <- data.frame(Measurements = exp_dt$Measurements, check.names = FALSE, stringsAsFactors = FALSE)
  hdr <- c("Measurements")
  
  # a long run of spaces reserves horizontal room for the CI panel
  ci_pad <- strrep(" ", ci_spaces)
  
  for (g in groups){
    forest_dt[[paste0(g,"__ESTCOL")]] <- exp_dt[[paste0(g,"__EST")]]
    forest_dt[[paste0(g,"__WC")]]     <- vapply(seq_len(nrow(exp_dt)), function(i){
      co <- exp_dt$.__cohort__[i]
      v  <- exp_dt[[paste0(g,"__", co, "_W")]][i]
      if (is.null(v) || is.na(v) || !nzchar(v)) "NA" else as.character(v)
    }, character(1))
    
    # force a *wide* blank CI column using spaces
    forest_dt[[paste0(g,"__CIC")]] <- ci_pad
    
    if (add_gutter) forest_dt[[paste0(g,"__SP")]] <- " "
    
    # headers
    hdr <- c(hdr,
             paste(g, "Estimate (SE)", sep = "\n"),
             "Weight",
             "",                    # CI plot has no header
             if (add_gutter) "")
  }
  names(forest_dt) <- hdr
  forest_dt
}

# 3) Plot (expanded rows): fixed widths, narrower rows, per-block bg
# NO gap rows - cleaner alternating backgrounds
plot_forest_expanded <- function(
    dt_, groups, cohorts = c("NACC","AIBL","HABS-HD"),
    IF_META_ANALYSIS = TRUE, meta_label = "Meta",
    title = "", filename = NULL, width = NULL, height = NULL,
    ref_line = 0, 
    # xlim_fixed = c(-0.8, 0.8),
    xlim_per_group = NULL,  # NEW: replaces xlim_fixed
    # point styles
    cohort_colors = c("#DAA87C","#92A5D1","#C9DCC4"), cohort_shapes = c(15,16,18),
    meta_color = "#000000", meta_shape = 17,
    # tighter rows
    row_padding_mm = 4, lineheight_txt = 0.96,
    # symbol sizing
    point_size = 0.85, line_width = 2.2,
    # column widths - ABSOLUTE widths in inches
    add_gutter = TRUE, gutter_width = 0.5,
    first_col_width = 3.5,   # Measurements column in inches
    est_frac = 2.5,          # Estimate column in inches
    wgt_frac = 1.5,          # Weight column in inches
    ci_frac = 8.0,           # CI plot in inches (WIDE!)
    # background fill per 4-row block
    block_col1 = "#C9DCC4", block_col2 = "#FFFFFF"){
  
  stopifnot(requireNamespace("forestploter", quietly = TRUE))
  stopifnot(requireNamespace("grid", quietly = TRUE))
  
  cohorts2 <- if (IF_META_ANALYSIS) c(cohorts, meta_label) else cohorts
  cols2    <- if (IF_META_ANALYSIS) c(cohort_colors, meta_color) else cohort_colors
  pchs2    <- if (IF_META_ANALYSIS) c(cohort_shapes, meta_shape) else cohort_shapes
  
  # expand rows (4 per measurement, NO gap row) and build table
  exp_dt    <- expand_dt_with_meta_rows(dt_, groups, cohorts = cohorts, meta_label = meta_label, add_block_gap = FALSE)
  forest_dt <- create_forest_table_expanded(exp_dt, groups, add_gutter = add_gutter)
  
  # series: one per (cohort × group); values only on that cohort's rows
  est_list <- lower_list <- upper_list <- list()
  add_na <- function(x) if (is.null(x)) rep(NA_real_, nrow(exp_dt)) else x
  row_co <- exp_dt$.__cohort__
  
  for (co in cohorts2){
    sel <- row_co == co
    for (g in groups){
      e  <- add_na(exp_dt[[paste0(g,"__",co,"_est")]]); e[!sel]  <- NA_real_
      lo <- add_na(exp_dt[[paste0(g,"__",co,"_low")]]); lo[!sel] <- NA_real_
      hi <- add_na(exp_dt[[paste0(g,"__",co,"_hi") ]]); hi[!sel] <- NA_real_
      est_list   <- c(est_list,   list(e))
      lower_list <- c(lower_list, list(lo))
      upper_list <- c(upper_list, list(hi))
    }
  }
  
  # CI columns: per-group chunk = Est, Weight, CI, (Spacer)
  chunk <- if (add_gutter) 4L else 3L
  ci_columns <- (2 + (seq_along(groups)-1L)*chunk) + 2L
  
  ## NEW: Handle xlim_per_group
  if (is.null(xlim_per_group)) {
    # Default: same range for all groups
    xlim_use <- rep(list(c(-0.8, 0.8)), length(groups))
  } else if (is.list(xlim_per_group) && !is.null(names(xlim_per_group))) {
    # Named list: extract in order of groups
    xlim_use <- lapply(groups, function(g) {
      if (g %in% names(xlim_per_group)) xlim_per_group[[g]] else c(-0.8, 0.8)
    })
  } else if (is.list(xlim_per_group) && length(xlim_per_group) == length(groups)) {
    # Unnamed list: use as-is (assumed to match groups order)
    xlim_use <- xlim_per_group
  } else {
    stop("xlim_per_group must be NULL, a named list with group names, or a list of length equal to number of groups")
  }
  
  # Use grid::unit() for column widths
  per_group <- if (add_gutter) c(est_frac, wgt_frac, ci_frac, gutter_width) else c(est_frac, wgt_frac, ci_frac)
  col_w <- grid::unit(c(first_col_width, rep(per_group, length(groups))), "inches")
  
  tm <- forestploter::forest_theme(
    base_size  = 12,
    refline_gp = grid::gpar(lty = "dotted", col = "grey70", lwd = 1.1),
    ci_pch     = pchs2, ci_col = cols2, ci_lwd = line_width,
    legend_name  = " ", legend_value = cohorts2,
    core = list(
      fg_params = list(hjust = rep(0, ncol(forest_dt)),
                       x = rep(0.01, ncol(forest_dt)),
                       lineheight = lineheight_txt),
      padding   = grid::unit(c(row_padding_mm, row_padding_mm), "mm")
    ),
    colwidth = col_w
  )
  
  p <- forestploter::forest(
    forest_dt,
    est   = est_list, lower = lower_list, upper = upper_list,
    ci_column = ci_columns,
    ref_line  = ref_line, 
    xlim = xlim_use,  # NEW: pass list of xlims
    # xlim = xlim_fixed,
    nudge_y   = 0, sizes = point_size, theme = tm, title = title,
    arrow_lab = c("Decrease","Increase"),
    indent = rep(0, nrow(forest_dt))
  )
  
  
  ## Bold significant cells in the printed table
  for (g in groups){
    col_idx  <- which(names(forest_dt) == paste(g, "Estimate (SE)", sep = "\n"))
    sig_rows <- which(grepl("\\*|\\^", exp_dt[[paste0(g,"__EST")]]))
    if (length(sig_rows)) p <- forestploter::edit_plot(p, row = sig_rows, col = col_idx, which="text", gp=grid::gpar(fontface="bold"))
  }
  
  p <- forestploter::edit_plot(p, row=1, col=seq_len(ncol(forest_dt)), which="text", y=grid::unit(0.5,"npc"), vjust=0.5)
  
  # Bold the (single) Measurement labels (top row of each block)
  p <- forestploter::edit_plot(
    p, row = which(forest_dt$Measurements != " "),
    col = 1, gp = grid::gpar(fontface = "bold")
  )
  
  # ---- Per-block background (each 4-row block) ----
  step <- 4L                           # 4 cohort rows (no gap)
  block_starts <- seq(1, nrow(forest_dt), by = step)
  for (i in block_starts){
    fill_col <- if ((((i-1)/step) %% 2) == 0) block_col1 else block_col2
    p <- forestploter::edit_plot(
      p, row = i:(i+3), which = "background",
      gp = grid::gpar(fill = fill_col, col = NA)
    )
  }
  
  if (!is.null(filename)) {
    # grDevices::pdf(filename, width = width, height = height)
    grDevices::cairo_pdf(filename, width = width, height = height, family = "Arial")
    plot(p);
    grDevices::dev.off()
  }
  
  invisible(p)
}


plot_forest_subset_meta_expanded <- function(
    fits, terms,
    groups  = c("Full","Female","Male","APOE4 carriers","APOE4 non-carriers"),
    cohorts = c("NACC","AIBL","HABS-HD"),
    term_labels = NULL, title = "",
    filename = NULL, width = NULL, height = NULL,
    # layout / styling knobs
    ref_line = 0, 
    # xlim_fixed = c(-0.8, 0.8),
    xlim_per_group = NULL,  # NEW: replaces xlim_fixed
    cohort_colors = c("#DAA87C", "#92A5D1", "#C9DCC4"),
    cohort_shapes = c(15, 16, 18),
    IF_META_ANALYSIS = TRUE, meta_label = "Meta", meta_color = "#000000", meta_shape = 17,
    row_padding_mm = 4, lineheight_txt = 0.96,
    point_size = 0.85, line_width = 2.2,
    add_gutter = TRUE, gutter_width = 0.5,
    first_col_width = 3.5, est_frac = 2.5, wgt_frac = 1.5, ci_frac = 8.0,
    block_col1 = "#f8f1ff", block_col2 = "#FFFFFF",
    p_adjust = c("none","BH"),              # keep passing through to make_dt_from_fits
    feature_subset = NULL                   # NEW: plot only these features (terms)
){
  stopifnot(all(vapply(fits, length, 1L) == length(cohorts)))
  fits_subset <- fits[names(fits) %in% groups]; fits_subset <- fits_subset[groups]
  
  # ## --- NEW: choose which terms to plot, preserving original order ---
  # terms_use <- terms
  # if (!is.null(feature_subset)) {
  #   miss <- setdiff(feature_subset, terms)
  #   if (length(miss)) warning("These requested features are not in 'terms': ", paste(miss, collapse=", "))
  #   terms_use <- intersect(terms, feature_subset)
  #   if (!length(terms_use)) stop("No overlap between 'feature_subset' and 'terms'.")
  # }
  
  dt_full <- make_dt_from_fits(
    fits_subset, terms,
    cohorts = cohorts, term_labels = term_labels,
    method = "REML", add_pred_int = TRUE, p_adjust = p_adjust
  )
  
  ## now subset rows for plotting only (preserve original order in `terms`)
  if (!is.null(feature_subset)) {
    miss <- setdiff(feature_subset, terms)
    if (length(miss)) warning("These requested features are not in 'terms': ", paste(miss, collapse=", "))
    keep_idx <- which(terms %in% feature_subset)
    if (!length(keep_idx)) stop("No overlap between 'feature_subset' and 'terms'.")
    dt_ <- dt_full[keep_idx, , drop = FALSE]
  } else {
    dt_ <- dt_full
  }
  
  plot_forest_expanded(
    dt_, groups = groups, cohorts = cohorts,
    IF_META_ANALYSIS = IF_META_ANALYSIS, meta_label = meta_label,
    title = title, filename = filename, width = width, height = height,
    ref_line = ref_line, 
    xlim_per_group = xlim_per_group,  # NEW
    # xlim_fixed = xlim_fixed,
    cohort_colors = cohort_colors, cohort_shapes = cohort_shapes,
    meta_color = meta_color, meta_shape = meta_shape,
    row_padding_mm = row_padding_mm, lineheight_txt = lineheight_txt,
    point_size = point_size, line_width = line_width,
    add_gutter = add_gutter, gutter_width = gutter_width,
    first_col_width = first_col_width, est_frac = est_frac, wgt_frac = wgt_frac, ci_frac = ci_frac,
    block_col1 = block_col1, block_col2 = block_col2
  )
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
                                dataset_order = c("AIBL","HABS-HD","NACC"),
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
    # ggplot2::scale_fill_manual(values = fill_values, drop = FALSE) +
    ggplot2::scale_fill_manual(values = fill_values, labels = c(Yes = "User", No = "Non-user"), name = NULL, drop = FALSE)+
    ggplot2::scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    ggplot2::labs(x = NULL, y = "Proportion of participants", fill = " ",
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
    dataset_order = c("NACC","AIBL","HABS-HD"),
    show_only_any = FALSE,
    fill_values = c("Yes" = "#2C7BB6", "No" = "#BDBDBD"),
    add_borders = FALSE,
    text_size = 9) {
  
  aibl_b <- tag_dataset(AIBL_df, "AIBL", med_cols)
  habs_b <- tag_dataset(HABS_df, "HABS-HD", med_cols)
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
  p + theme(plot.subtitle = element_text(size = 30), legend.text = element_text(size=20), legend.key.size = grid::unit(4,"mm"))
  p <- p + theme(strip.text.x = element_text(size = 14, face = "bold"))
  # p + scale_fill_manual(values=fill_values, labels=c("Yes"="User","No"="Non-user"), name=NULL)
  
  out_dir <- dirname(out_pdf)
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  
  # grDevices::pdf(out_pdf, width = width, height = height)
  grDevices::cairo_pdf(out_pdf, width = width, height = height, family = "Arial")
  on.exit(grDevices::dev.off(), add = TRUE)
  print(p)
  
  invisible(list(plot = p, table = tbl_baseline, meds_in = meds_in,
                 path = out_pdf, width = width, height = height))
}

# Neat table (with "Any medication" last within each dataset)
baseline_props_table <- function(
    NACC_df, AIBL_df, HABS_df,
    med_cols = c("ACEi","ARB","BetaBlk","CCB","Diuretic","Statin","Metformin"),
    dataset_order = c("AIBL","HABS-HD","NACC"),
    include_individual_meds = TRUE,
    include_any_med = TRUE,
    digits = 1,
    out_csv = NULL) {
  
  stopifnot(include_individual_meds || include_any_med)
  
  aibl_b <- tag_dataset(AIBL_df, "AIBL", med_cols)
  habs_b <- tag_dataset(HABS_df, "HABS-HD", med_cols)
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
    dataset_order  = c("NACC","AIBL","HABS-HD"),
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
  habs_b <- tag_dataset(HABS_df, "HABS-HD", med_cols)
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
  p <- p + theme(plot.subtitle = element_text(size = 12), legend.text = element_text(size = 12), text = element_text(size = 12))
  p <- p + theme(strip.text.x = element_text(size = 14, face = "bold"))
  return(p)
}

# ---------- 2) Attained-age (visit-level) smoother ----------
plot_attained_age_prevalence_multi <- function(
    NACC_long, AIBL_long, HABS_long,
    med_cols       = c("ACEi","ARB","BetaBlk","CCB","Diuretic","Statin","Metformin"),
    dataset_order  = c("AIBL","HABS-HD","NACC"),
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
  habs_v <- tag_dataset_long(HABS_long, "HABS-HD", med_cols)
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
  p <- p + theme(plot.subtitle = element_text(size = 12), legend.text = element_text(size = 12), text = element_text(size = 12))
  p <- p + theme(strip.text.x = element_text(size = 14, face = "bold"))
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
    dataset_order = c("NACC","AIBL","HABS-HD"),
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
    "HABS-HD" = make_binary_subset(HABS_df, med_cols)
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
