library(dplyr)
library(purrr)
library(rlang)
library(R6)
library(dplyr)
library(tidyr)
library(purrr)
library(tibble)
library(matrixStats)
library(glmmTMB)
library(parameters)
library(stringr)

min_except_zero_per_col <- function(data) {
  data <- as.data.frame(data)
  if (!requireNamespace("dplyr", quietly = TRUE))
    stop("Package 'dplyr' is required but not installed.", call. = FALSE)
  
  dplyr::summarise(
    data,
    dplyr::across(
      where(is.numeric),
      ~ {
        non_zero_vals <- .x[.x != 0]
        if (length(non_zero_vals) == 0) NA_real_ else min(non_zero_vals, na.rm = TRUE)
      }
    )
  )
}


ClonalAnalysis <- R6Class(
  "ClonalAnalysis",
  
  public = list(
    
    metadata       = NULL,
    counts         = NULL,
    params         = list(max_fract_thr = 2e-5,
                          min_dt        = 1e-6,
                          use_weights   = TRUE),
    last_meta      = NULL,
    analysis_info  = NULL,
    counts_filt    = NULL,
    dR_long        = NULL,
    modelling_data = NULL,
    model          = NULL,
    
    #####################. CONSTRUCTOR
    initialize = function(metadata, counts) {
      self$metadata <- private$.validate_and_augment_metadata(metadata)
      
      if (is.data.frame(counts)) counts <- as.matrix(counts)
      if (!is.matrix(counts)) stop("`counts` must be a matrix / data.frame.")
      if (!is.numeric(counts)) storage.mode(counts) <- "numeric"
      if (is.null(rownames(counts)) || any(rownames(counts) == ""))
        stop("Row names (lineage IDs) are required in `counts`.")
      
      bad <- setdiff(self$metadata$sample_name, colnames(counts))
      if (length(bad))
        stop("Counts matrix is missing columns for: ", paste(bad, collapse = ", "))
      
      self$counts <- counts
    },
    
    ### PUBLIC
    fit_model = function(max_fract_thr = 2e-5, min_dt = 1e-6,
                         treatments   = NULL,
                         control0     = NULL,
                         control1     = NULL,
                         use_weights  = TRUE) {
      
      self$params$max_fract_thr <- max_fract_thr
      self$params$min_dt        <- min_dt
      self$params$use_weights   <- use_weights
      
      local_meta   <- private$.build_local_metadata(self, treatments,
                                                    control0, control1)
      local_counts <- self$counts[, local_meta$sample_name, drop = FALSE]
      
      orig_meta <- self$metadata
      on.exit(self$metadata <- orig_meta, add = TRUE)
      self$metadata <- local_meta
      
      ## FRACTION matrix (computed once)
      self$counts_filt <- private$.prepare_counts(
        local_counts, local_meta, self$params$max_fract_thr)
      
      private$.prep_model(self)
      private$.fit(self)
      self$last_meta <- local_meta
      
      ctrl0_map <- local_meta %>%
        dplyr::filter(treatment == "control0") %>%
        dplyr::select(rep, ctrl0_sample = sample_name)
      
      ctrl1_map <- local_meta %>%
        dplyr::filter(treatment == "control1") %>%
        dplyr::select(rep, ctrl1_sample = sample_name)
      
      self$analysis_info <- local_meta %>%
        dplyr::select(sample_name, treatment, rep, time) %>%
        dplyr::left_join(ctrl0_map, by = "rep") %>%
        dplyr::left_join(ctrl1_map, by = "rep")
      
      invisible(local_meta)
    },
    
    summary = function(...) {
      if (is.null(self$model)) stop("fit_model() has not been run.")
      summary(self$model, ...)
    },
    
    get_model_rates = function(aggregate = TRUE) {
      if (is.null(self$model)) stop("fit_model() has not been run.")
      df <- dplyr::filter(self$modelling_data, sample_type != "control0")
      
      if (aggregate) {
        df %>%
          dplyr::group_by(lineage, condition, treatment) %>%
          dplyr::summarise(dR = mean(predicted_dR),
                           R  = mean(predicted_R, na.rm = TRUE),
                           .groups = "drop")
      } else {
        dplyr::select(df, lineage, rep, condition, treatment,
                      dR = predicted_dR, R = predicted_R)
      }
    },
    
## ── PUBLIC results() ────────────────────────────────────────────────
results = function(mode                    = "model",
                   delta                   = NULL,
                   alpha                   = NULL,
                   aggregate               = TRUE,
                   format                  = c("long", "list"),
                   enable.dev              = FALSE,
                   intercept_zero          = TRUE,
                   include_raw_versions    = FALSE) {
  if (is.null(self$model)) stop("fit_model() has not been run.")
  mode   <- match.arg(mode, c("model", "all", "raw"))
  format <- match.arg(format)

  if (!isTRUE(enable.dev) && mode != "model")
    stop('When enable.dev = FALSE, only mode = "model" is allowed.')

  internal_mode <- if (isTRUE(include_raw_versions)) "all" else mode

  long_tbl <- private$..results_internal(
    internal_mode, delta, alpha, aggregate, intercept_zero = intercept_zero
  )

  if (isTRUE(enable.dev)) {
    if (format == "long") return(long_tbl)

    unique_map <- long_tbl %>%
      dplyr::distinct(condition, treatment) %>%
      dplyr::add_count(treatment, name = "n_by_trt") %>%
      dplyr::filter(n_by_trt == 1) %>%
      dplyr::pull(condition, name = treatment)

    long_tbl2 <- long_tbl %>%
      dplyr::mutate(cond_label = ifelse(condition %in% unique_map, treatment, condition))

    metric_cols <- setdiff(
      names(long_tbl2),
      c("lineage", "rep", "condition", "cond_label", "treatment")
    )

    out <- purrr::map(metric_cols, \(m)
      long_tbl2 %>%
        dplyr::select(lineage, cond_label, !!rlang::sym(m)) %>%
        tidyr::pivot_wider(
          names_from  = cond_label,
          values_from = !!rlang::sym(m)
        ) %>%
        tibble::column_to_rownames("lineage") %>%
        as.data.frame()
    )
    names(out) <- metric_cols
    return(out)
  }

  if (isTRUE(include_raw_versions)) {
    metric_map <- c(
      ## model (pred) — shown under the usual names
      dR_glmnet_resid_pred = "cGR",
      R_pred               = "growth_rate",
      dR_pred              = "centered_growth_rate",
      R_ctrl_pred          = "growth_rate_control",
      dR_ctrl_pred         = "centered_growth_rate_control",
      fr_pred              = "fraction",
      fr_ctrl_pred         = "fraction_control",
      fr_t0_pred           = "fraction_t0",
      dR_glmnet_fit_pred   = "dR_glmnet_fit",
      ## raw companions
      R_raw                = "growth_rate_raw",
      dR_raw               = "centered_growth_rate_raw",
      R_ctrl_raw           = "growth_rate_control_raw",
      dR_ctrl_raw          = "centered_growth_rate_control_raw",
      fr_raw               = "fraction_raw",
      fr_ctrl_raw          = "fraction_control_raw",
      fr_t0_raw            = "fraction_t0_raw",
      ## keep these if they appear with suffixes in "all" mode
      growth_difference_pred = "growth_difference",
      FC_pred                = "FC",
      SR_delta_pred          = "SR_delta",
      SLR_alpha_pred         = "SLR_alpha",
      proportional_pred      = "proportional",
      growth_difference_raw  = "growth_difference_raw",
      FC_raw                 = "FC_raw",
      SR_delta_raw           = "SR_delta_raw",
      SLR_alpha_raw          = "SLR_alpha_raw",
      proportional_raw       = "proportional_raw"
    )
  } else {
    metric_map <- c(
      dR_glmnet_resid   = "cGR",
      R                 = "growth_rate",
      dR                = "centered_growth_rate",
      growth_difference = "growth_difference",
      FC                = "FC",
      R_ctrl            = "growth_rate_control",
      dR_ctrl           = "centered_growth_rate_control",
      p_value_pair      = "p_value",
      fr                = "fraction",
      fr_ctrl           = "fraction_control",
      fr_t0             = "fraction_t0",
      dR_glmnet_fit     = "dR_glmnet_fit"
    )
  }

  if (format == "long") {
    core_cols    <- intersect(
      c("lineage", "condition", "treatment", "rep"),
      names(long_tbl)
    )
    keep_metrics <- intersect(names(metric_map), names(long_tbl))

    result_main <- long_tbl %>%
      dplyr::select(dplyr::all_of(c(core_cols, keep_metrics))) %>%
      dplyr::rename_with(~ metric_map[.x], .cols = keep_metrics)

  } else {  # format == "list"
    unique_map <- long_tbl %>%
      dplyr::distinct(condition, treatment) %>%
      dplyr::add_count(treatment, name = "n_by_trt") %>%
      dplyr::filter(n_by_trt == 1) %>%
      dplyr::pull(condition, name = treatment)

    long_tbl2 <- long_tbl %>%
      dplyr::mutate(cond_label = ifelse(condition %in% unique_map, treatment, condition))

    metric_cols <- setdiff(
      names(long_tbl2),
      c("lineage", "rep", "condition", "cond_label", "treatment")
    )

    tmp <- purrr::map(metric_cols, \(m)
      long_tbl2 %>%
        dplyr::select(lineage, cond_label, !!rlang::sym(m)) %>%
        tidyr::pivot_wider(
          names_from  = cond_label,
          values_from = !!rlang::sym(m)
        ) %>%
        tibble::column_to_rownames("lineage") %>%
        as.data.frame()
    )
    names(tmp) <- metric_cols

    keep_metrics <- intersect(names(metric_map), names(tmp))
    result_main  <- tmp[keep_metrics]
    names(result_main) <- metric_map[keep_metrics]
  }

  list(
    result        = result_main,
    analysis_info = self$analysis_info,
    last_meta     = self$last_meta
  )
}
  ),
  
  
  
  
  ############################################################## 
  private = list(
    
    
    .build_local_metadata = function(self, treatments,
                                     control0, control1) {
      
      local_meta <- self$metadata
      
      ## 1. validate `treatments`
      if (!is.null(treatments)) {
        if (is.null(names(treatments)) || any(names(treatments) == ""))
          stop("`treatments` must be a *named* list; names become new labels.")
        clash <- intersect(names(treatments), unique(local_meta$treatment))
        if (length(clash))
          stop("Names in `treatments` collide with existing labels: ",
               paste(clash, collapse = ", "))
      }
      
      ## 2. build <sample_name → new label> map
      mapping <- character(0)
      add_map <- function(samples, lab) {
        if (is.null(samples)) return()
        dup <- intersect(names(mapping), samples)
        if (length(dup))
          stop("Sample(s) assigned twice: ", paste(dup, collapse = ", "))
        mapping[samples] <<- lab
      }
      add_map(control0, "control0")
      add_map(control1, "control1")
      if (!is.null(treatments)) purrr::imap(treatments, add_map)
      
      bad <- setdiff(names(mapping), local_meta$sample_name)
      if (length(bad))
        stop("Unknown sample_name(s): ", paste(bad, collapse = ", "))
      
      ## 3. relabel treatment column
      idx <- match(names(mapping), local_meta$sample_name)
      local_meta$treatment[idx] <- mapping
      
      ## 4. sub-set rows
      keep_row <- with(local_meta, ifelse(
        treatment == "control0",
        is.null(control0)  | sample_name %in% control0,
        ifelse(
          treatment == "control1",
          is.null(control1) | sample_name %in% control1,
          {
            if (is.null(treatments)) TRUE
            else treatment %in% names(treatments)
          }
        )
      ))
      local_meta <- local_meta[keep_row, , drop = FALSE]
      
      ## 4a. keep one control1 per replicate (closest in time)
      if (any(local_meta$treatment == "control1")) {
        
        ctrl_tbl <- local_meta %>%
          dplyr::filter(treatment == "control1") %>%
          dplyr::select(ctrl_sample = sample_name, rep, ctrl_time = time)
        
        trt_mean <- local_meta %>%
          dplyr::filter(!treatment %in% c("control0", "control1")) %>%
          dplyr::group_by(rep) %>%
          dplyr::summarise(t_bar = mean(time), .groups = "drop")
        
        match_tbl <- trt_mean %>%
          dplyr::left_join(ctrl_tbl, by = "rep") %>%
          dplyr::mutate(abs_dt = abs(ctrl_time - t_bar)) %>%
          dplyr::group_by(rep) %>%
          dplyr::slice_min(order_by = abs_dt, with_ties = TRUE) %>%
          dplyr::slice_min(order_by = ctrl_time) %>%
          dplyr::ungroup()
        
        keep_ctrl <- unique(match_tbl$ctrl_sample)
        
        local_meta <- local_meta %>%
          dplyr::filter(treatment != "control1" | sample_name %in% keep_ctrl)
        
        if (any(duplicated(
          dplyr::filter(local_meta, treatment == "control1")$rep)))
          stop("Control1 pruning failed: replicate still has > 1 vehicle.")
      }
      
      ## 5. sanity checks
      if (!any(local_meta$treatment == "control0"))
        stop("After sub-setting, no control0 samples remain.")
      if (!any(local_meta$treatment == "control1"))
        stop("After sub-setting, no control1 samples remain.")
      
      local_meta %>%
        dplyr::mutate(
          condition = paste(treatment, time, sep = "_"),
          rep       = factor(rep)
        )
    },
    
    ## ---------- misc utilities ------------------------------------------
    .need_cols = function(df, cols, label) {
      miss <- setdiff(cols, names(df))
      if (length(miss))
        stop(label, " missing columns: ", paste(miss, collapse = ", "))
    },
    
    .find_min_except_zero = function(x)
      min(x[x > 0], na.rm = TRUE),
    
    ## ---------- metadata validation / augmentation ----------------------
    .validate_and_augment_metadata = function(metadata) {
      
      core <- c("treatment", "time", "rep", "sample_name")
      private$.need_cols(metadata, core, "metadata")
      
      ## derive growth_rate if absent
      if (!"growth_rate" %in% names(metadata)) {
        
        if (!"fold_expansion" %in% names(metadata))
          stop("metadata needs `fold_expansion` to compute growth_rate.")
        
        metadata <- metadata %>%
          dplyr::group_by(rep) %>%
          dplyr::mutate(t0 = min(time[treatment == "control0"], na.rm = TRUE)) %>%
          dplyr::ungroup()
        
        if (any(is.infinite(metadata$t0)))
          stop("Each replicate needs ≥ 1 control0 to define baseline time.")
        
        metadata <- metadata %>%
          dplyr::mutate(delta_t = time - t0)
        
        bad_dt <- metadata %>%
          dplyr::filter(treatment != "control0" & delta_t <= 0)
        if (nrow(bad_dt))
          stop("Non-control0 sample(s) have delta_t ≤ 0: ",
               paste(bad_dt$sample_name, collapse = ", "))
        
        bad_fe <- metadata %>%
          dplyr::filter(treatment != "control0" &
                          (!is.finite(fold_expansion) | fold_expansion <= 0))
        if (nrow(bad_fe))
          stop("fold_expansion must be positive for: ",
               paste(bad_fe$sample_name, collapse = ", "))
        
        metadata <- metadata %>%
          dplyr::mutate(growth_rate = dplyr::case_when(
            treatment == "control0" ~ NA_real_,
            TRUE                    ~ log(fold_expansion) / delta_t
          ))
      }
      
      if (!all(c("control0", "control1") %in% metadata$treatment))
        stop("metadata must contain at least one control0 and one control1 row.")
      
      metadata %>%
        dplyr::mutate(
          condition = paste(treatment, time, sep = "_"),
          rep       = factor(rep)
        )
    },
    
    ## ---------- lineage filter
    .prepare_counts = function(counts, metadata, thr) {
      
      total_counts <- colSums(counts[, metadata$sample_name])
      cnt_pc  <- counts[, metadata$sample_name]
      col_specific_pseudo <- min_except_zero_per_col(cnt_pc)
      cnt_pc[cnt_pc == 0] <- rep(as.vector(unlist(col_specific_pseudo)),
                                 each = nrow(cnt_pc))[cnt_pc == 0]
      
      fr <- sweep(cnt_pc, 2, total_counts, "/")
      fr[matrixStats::rowMaxs(fr) > thr, , drop = FALSE]
    },
    
    ## ---------- delta R computation 
    .calc_dR = function(fr, metadata, min_dt) {
      if (any(abs(colSums(fr) - 1) > 1e-8))
        fr <- sweep(fr, 2, colSums(fr), "/")
      
      time_lut <- setNames(metadata$time, metadata$sample_name)
      
      purrr::map_dfr(seq_len(nrow(metadata)), function(i) {
        
        smp   <- metadata$sample_name[i]
        rep_i <- metadata$rep[i]
        
        c0 <- metadata %>%
          dplyr::filter(rep == rep_i, treatment == "control0") %>%
          dplyr::arrange(time) %>%
          dplyr::slice(1) %>%
          dplyr::pull(sample_name)
        
        dt <- max(abs(time_lut[smp] - time_lut[c0]), min_dt)
        
        tibble::tibble(
          sample_name  = smp,
          lineage      = rownames(fr),
          dR           = log(fr[, smp] / fr[, c0]) / dt,
          dt_to_ctrl0  = dt
        )
      }) %>%
        dplyr::left_join(
          metadata %>%
            dplyr::select(sample_name, treatment, time, rep,
                          condition, growth_rate),
          by = "sample_name"
        ) %>%
        dplyr::mutate(sample_type = dplyr::case_when(
          treatment == "control0" ~ "control0",
          treatment == "control1" ~ "control1",
          TRUE                    ~ "treatment"
        ))
    },
    
    ## ---------- weight-spline helper  ------------------------------------
    .fit_w = function(cnt_sub, dR_mat) {
      lm  <- log10(rowMeans(cnt_sub))
      lw  <- log10(1 / matrixStats::rowVars(dR_mat))
      idx <- is.finite(lm) & is.finite(lw)
      
      if (sum(idx) < 4)
        return(setNames(rep(1, length(lm)), names(lm)))
      
      sp  <- smooth.spline(lm[idx], lw[idx], df = 3)
      ap  <- approx(sp$x, sp$y, xout = lm, rule = 2)
      setNames(10**ap$y, names(lm))
    },
    
    ## ----------  modelling_data  ---------------------------------
    .prep_model = function(self) {
      dR_long <- private$.calc_dR(
        self$counts_filt, self$metadata, self$params$min_dt)
      
      cond_keep <- unique(dplyr::filter(
        dR_long, sample_type != "control0")$condition)
      w_vec     <- rep(NA_real_, nrow(dR_long))
      
      for (cond in cond_keep) {
        idx <- which(dR_long$condition == cond)
        dfc <- dR_long[idx, ]
        
        dR_mat <- dfc %>%
          dplyr::select(lineage, rep, dR) %>%
          tidyr::pivot_wider(names_from  = rep,
                             values_from = dR,
                             values_fn   = mean,
                             values_fill = NA) %>%
          tibble::column_to_rownames("lineage") %>%
          as.matrix()
        
        cnt_sub <- self$counts_filt[, unique(dfc$sample_name), drop = FALSE]
        w       <- private$.fit_w(cnt_sub, dR_mat)
        w_vec[idx] <- w[dfc$lineage] / mean(w, na.rm = TRUE)
      }
      
      dR_long$weight <- w_vec
      self$dR_long   <- dR_long
      
      self$modelling_data <- dplyr::filter(dR_long, sample_type != "control0") %>%
        dplyr::mutate(condition = factor(condition),
                      rep       = factor(rep))
    },
    
    ## ---------- mixed-effects fit 
    .run_mixed = function(df, use_weights = TRUE) {
      args <- list(
        formula = dR ~ 0 + condition + (0 + condition || lineage),
        data    = df,
        family  = gaussian(),
        control = glmmTMBControl(
          optimizer = stats::optim,
          optArgs   = list(method = "BFGS", maxit = 2e5),
          profile   = TRUE)
      )
      if (use_weights) args$weights <- df$weight
      do.call(glmmTMB, args)
    },
    
    ## ---------- fit wrappe
    .fit = function(self) {
      self$model <- private$.run_mixed(self$modelling_data,
                                       self$params$use_weights)
      
      ## predictions
      self$modelling_data <- self$modelling_data %>%
        dplyr::mutate(predicted_dR = predict(self$model, re.form = NULL),
                      predicted_R  = predicted_dR + growth_rate,
                      raw_R        = dR + growth_rate)
      
      ## lineage   condition p-values heuristic
      p_tbl <- broom.mixed::tidy(self$model, effects = "ran_vals", conf.int = FALSE) %>%
        dplyr::filter(component == "cond", group == "lineage") %>%
        dplyr::mutate(
          condition = sub("^condition", "", term),
          lineage   = level
        ) %>%
        dplyr::transmute(
          lineage,
          condition,
          p_value_pair = 2 * stats::pnorm(abs(estimate / std.error), lower.tail = FALSE)
        )
      
      ## merge everything
      self$modelling_data <- self$modelling_data %>%
        dplyr::left_join(p_tbl,  by = c("lineage", "condition"))
    },
    
    
    ..results_internal = function(mode, delta, alpha, aggregate, intercept_zero = TRUE) {
      # Bring modelling data and attach counts/fractions (columns added to `md`)
      md <- self$modelling_data
      md$counts <- self$counts[ cbind(md$lineage, md$sample_name) ]
      md$fr     <- self$counts_filt[ cbind(md$lineage, md$sample_name) ]
      
      # Choose columns based on mode
      use_pred <- mode != "raw"
      dRcol    <- if (use_pred) "predicted_dR" else "dR"
      Rcol     <- if (use_pred) "predicted_R"  else "raw_R"
      
      # Delta / alpha defaults + warning state
      if (use_pred && all(is.na(md[[Rcol]]))) {
        delta  <- delta %||% 1e-2
        alpha  <- alpha %||% 1
        warn_R <- TRUE
      } else {
        absR   <- abs(md[[Rcol]])
        delta  <- delta %||% max(1e-2, min(absR, na.rm = TRUE) / 2)
        alpha  <- alpha %||% stats::median(absR, na.rm = TRUE)
        warn_R <- FALSE
      }
      
      # control1 (vehicle-at-time) lookups from md
      ctrl <- dplyr::filter(md, sample_type == "control1") %>%
        dplyr::select(
          lineage, rep,
          R_ctrl  = dplyr::all_of(Rcol),
          dR_ctrl = dplyr::all_of(dRcol),
          fr_ctrl = fr
        )
      
      # control0 (baseline t0) fractions per replicate
      ctrl0_map <- self$metadata %>%
        dplyr::filter(treatment == "control0") %>%
        dplyr::arrange(rep, time) %>%
        dplyr::group_by(rep) %>%
        dplyr::slice(1) %>%
        dplyr::ungroup() %>%
        dplyr::select(rep, ctrl0_sample = sample_name)
      
      fr_t0_long <- purrr::map_dfr(seq_len(nrow(ctrl0_map)), function(i) {
        s <- ctrl0_map$ctrl0_sample[i]
        r <- ctrl0_map$rep[i]
        tibble::tibble(
          lineage = rownames(self$counts_filt),
          rep     = r,
          fr_t0   = self$counts_filt[, s]
        )
      })
      
      # Build per-replicate table (NO aggregation yet)
      trt <- dplyr::filter(md, sample_type != "control1") %>%
        dplyr::left_join(ctrl,       by = c("lineage", "rep")) %>%
        dplyr::left_join(fr_t0_long, by = c("lineage", "rep")) %>%
        dplyr::mutate(
          current_dR        = .data[[dRcol]],
          current_R         = .data[[Rcol]],
          growth_difference = (R_ctrl - current_R),
          FC                = dR_ctrl - current_dR,
          SR_delta          = growth_difference /
            (abs(current_R) + abs(R_ctrl) + 2 * delta),
          SLR_alpha         = asinh(current_R / alpha) - asinh(R_ctrl / alpha),
          proportional               = growth_difference / (abs(R_ctrl) + delta),
          R_ctrl            = R_ctrl,
          dR_ctrl           = dR_ctrl,
          fr_ctrl           = fr_ctrl
        )
      
      ## -------- GLMNET residualization on dR BEFORE aggregation (per replicate) ----
      ## Constrained lasso: slope(dR_ctrl) >= 0; dummy locked at 0.
      ## Optional zero-intercept via `intercept_zero` (default TRUE).
      trt <- trt %>%
        dplyr::group_by(treatment, condition, rep) %>%
        dplyr::group_modify(~{
          df   <- .x
          xc   <- as.numeric(df$dR_ctrl)       # predictor dR_ctrl
          yy   <- as.numeric(df$current_dR)    # response dR_ctrl
          keep <- stats::complete.cases(xc, yy)
          
          df$dR_glmnet_fit   <- NA_real_
          df$dR_glmnet_resid <- NA_real_
          
          # Need at least 2 points and variation in predictor
          if (sum(keep) >= 2L && isTRUE(stats::sd(xc[keep]) > 0)) {
            X <- cbind(dR_ctrl = xc[keep], dummy = 0 * xc[keep])
            y <- yy[keep]
            
            intercept_flag <- if (isTRUE(intercept_zero)) FALSE else TRUE
            
            cvfit <- glmnet::cv.glmnet(
              x = X, y = y, family = "gaussian",
              alpha = 1,
              standardize = FALSE,
              intercept = intercept_flag,       # <- controls intercept (0 if FALSE)
              lower.limits   = c(0,   0),       # slope >= 0 ; dummy fixed at 0
              upper.limits   = c(Inf, 0),
              penalty.factor = c(1,   0)
            )
            
            X_all <- cbind(dR_ctrl = xc, dummy = 0 * xc)
            yhat  <- as.numeric(stats::predict(cvfit, newx = X_all, s = "lambda.min"))
            df$dR_glmnet_fit   <- yhat
            df$dR_glmnet_resid <- yhat - yy
          }
          df
        }) %>%
        dplyr::ungroup()
      ## ------------------------------- end GLMNET block ----------------------------
      
      # Aggregate AFTER fitting (averaging residuals/fits across replicates)
      if (aggregate) {
        trt <- trt %>%
          dplyr::group_by(lineage, condition, treatment) %>%
          dplyr::summarise(
            dR                = mean(current_dR),
            R                 = mean(current_R, na.rm = TRUE),
            growth_difference = mean(growth_difference, na.rm = TRUE),
            FC                = mean(FC),
            SR_delta          = mean(SR_delta,  na.rm = TRUE),
            SLR_alpha         = mean(SLR_alpha, na.rm = TRUE),
            proportional               = mean(proportional,       na.rm = TRUE),
            counts            = mean(counts),
            fr                = mean(fr),
            fr_ctrl           = mean(fr_ctrl,   na.rm = TRUE),
            fr_t0             = mean(fr_t0,     na.rm = TRUE),
            dR_ctrl           = mean(dR_ctrl,   na.rm = TRUE),
            R_ctrl            = mean(R_ctrl,    na.rm = TRUE),
            p_value_pair      = mean(p_value_pair, na.rm = TRUE),
            # NEW: averaged glmnet outputs on dR
            dR_glmnet_fit     = mean(dR_glmnet_fit,   na.rm = TRUE),
            dR_glmnet_resid   = mean(dR_glmnet_resid, na.rm = TRUE),
            .groups = "drop"
          )
      } else {
        trt <- dplyr::select(
          trt,
          lineage, rep, condition, treatment,
          dR  = current_dR,
          R   = current_R,
          growth_difference, FC, SR_delta, SLR_alpha, proportional,
          counts, fr, fr_ctrl, fr_t0, dR_ctrl, R_ctrl, p_value_pair,
          dR_glmnet_fit, dR_glmnet_resid
        )
      }
      
      if (warn_R && mode != "raw")
        warning("R columns are NA because growth_rate data were absent.")
      
      # Developer "all" mode combines model + raw
      if (mode == "all") {
        pred_tbl <- private$..results_internal("model", delta, alpha, aggregate,
                                               intercept_zero = intercept_zero) %>%
          dplyr::select(-dplyr::any_of(c("counts")))
        raw_tbl  <- private$..results_internal("raw",   delta, alpha, aggregate,
                                               intercept_zero = intercept_zero)
        key_cols <- intersect(
          c("lineage", "condition", "rep", "treatment"),
          intersect(names(pred_tbl), names(raw_tbl))
        )
        return(
          dplyr::left_join(pred_tbl, raw_tbl, by = key_cols,
                           suffix = c("_pred", "_raw"))
        )
      }
      
      trt
    }
    
    
    
  )
)



############### MERGES 2 results() OUTPUTS ###############
merge_clonal_results <- function(...,
                                 id                 = NULL,
                                 intersect_barcodes = TRUE) {
  
  objects <- list(...)
  if (length(objects) == 1L &&
      all(c("result", "analysis_info", "last_meta") %in% names(objects[[1]])))
    objects <- list(objects[[1]])
  
  if (!length(objects))
    stop("No objects supplied.")
  
  needed_slots <- c("result", "analysis_info", "last_meta")
  bad <- purrr::map_lgl(objects, \(x) !all(needed_slots %in% names(x)))
  if (any(bad))
    stop("Every input must be the 3-slot list returned by ClonalAnalysis$results().")
  
  n_sets <- length(objects)
  if (is.null(id))
    id <- paste0("set", seq_len(n_sets))
  if (length(id) != n_sets)
    stop("`id` must be NULL or have one entry per input object.")
  
  shape <- if (is.data.frame(objects[[1]]$result)) "long" else "list"
  if (!all(purrr::map_lgl(objects, \(x)
                          { if (shape == "long") is.data.frame(x$result) else is.list(x$result) })))
    stop("All `result` slots must share the same shape (\"long\" or \"list\").")
  
  bc_sets <- purrr::map(objects, \(obj) {
    if (shape == "long")
      unique(obj$result$lineage)
    else
      unique(unlist(purrr::map(obj$result, rownames)))
  })
  
  bc_union      <- Reduce(union, bc_sets)
  bc_intersect  <- Reduce(intersect, bc_sets)
  barcodes_kept <- if (isTRUE(intersect_barcodes)) bc_intersect else bc_union
  
  lost_prop <- 1 - length(barcodes_kept) / length(bc_union)
  if (lost_prop > 0.30)
    warning(sprintf(
      "%.1f%% of barcodes were lost during merging (%d of %d retained).",
      100 * lost_prop, length(barcodes_kept), length(bc_union)
    ))
  
  if (shape == "long") {
    
    df_list <- purrr::map2(objects, id, \(obj, tag) {
      obj$result |>
        dplyr::filter(.data$lineage %in% barcodes_kept) |>
        dplyr::mutate(.dataset = tag, .before = 1)
    })
    
    merged_result <- dplyr::bind_rows(df_list)
    
    if (any(duplicated(names(merged_result)))) {
      names(merged_result) <- make.unique(names(merged_result), sep = "::")
    }
    
  } else {
    
    metrics <- unique(unlist(purrr::map(objects, \(x) names(x$result))))
    merged_result <- setNames(vector("list", length(metrics)), metrics)
    
    rename_dup_cols <- function(mat, obj, tag, taken) {
      
      cn <- colnames(mat)
      dup_now <- cn %in% taken
      if (any(dup_now)) {
        
        lm <- obj$last_meta
        for (j in which(dup_now)) {
          cur <- cn[j]
          
          cond_candidates <- unique(
            lm$condition[lm$treatment == cur | lm$condition == cur]
          )
          
          if (length(cond_candidates) == 1L &&
              !cond_candidates %in% c(taken, cn)) {
            
            cn[j] <- cond_candidates
            
          } else {
            cn[j] <- paste0(tag, "::", cur)
          }
        }
      }
      colnames(mat) <- cn
      cn
    }
    
    for (metric in metrics) {
      
      existing_names <- character(0)
      
      mats <- purrr::map2(objects, id, \(obj, tag) {
        
        m <- obj$result[[metric]]
        if (is.null(m))
          return(NULL)
        
        keep <- intersect(barcodes_kept, rownames(m))
        m    <- m[keep, , drop = FALSE]
        
        existing_names <<- c(existing_names,
                             rename_dup_cols(m, obj, tag, existing_names))
        
        m
      }) |>
        purrr::compact()
      
      if (length(mats))
        merged_result[[metric]] <- do.call(cbind, mats)
    }
  }
  
  analysis_info <- purrr::map2(objects, id, \(obj, tag)
                               obj$analysis_info |>
                                 dplyr::mutate(.dataset = tag, .before = 1)
  ) |>
    dplyr::bind_rows() |>
    dplyr::distinct()
  
  last_meta <- purrr::map2(objects, id, \(obj, tag)
                           obj$last_meta |>
                             dplyr::mutate(.dataset = tag, .before = 1)
  ) |>
    dplyr::bind_rows() |>
    dplyr::distinct()
  
  list(
    result        = merged_result,
    analysis_info = analysis_info,
    last_meta     = last_meta
  )
}
