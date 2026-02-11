# ============================================================================
# DPMirt Internal Utilities
# ============================================================================

# --------------------------------------------------------------------------
# Session Signature
# --------------------------------------------------------------------------

#' Generate a session signature hash for compiled model reuse validation
#'
#' Creates a hash based on model structure + data dimensions + missingness
#' pattern. Used to validate whether a compiled model can be reused with
#' new data.
#'
#' @param spec A \code{dpmirt_spec} object.
#' @return A character string hash.
#' @noRd
.session_signature <- function(spec) {
  key_parts <- paste0(
    spec$config$model, "_",
    spec$config$prior, "_",
    spec$config$parameterization, "_",
    spec$config$identification, "_",
    spec$constants$N, "_",
    spec$constants$I, "_",
    if (!is.null(spec$constants$M)) spec$constants$M else "NA", "_",
    sum(is.na(spec$data$y))   # missingness fingerprint

  )
  digest_simple(key_parts)
}


#' Simple string hashing (no external dependency)
#'
#' @param x A character string.
#' @return A hex digest string.
#' @noRd
digest_simple <- function(x) {
  # Use a simple but effective hash based on R's internal serialization
  raw_bytes <- serialize(x, connection = NULL)
  paste0(format(as.hexmode(as.integer(raw_bytes)), width = 2),
         collapse = "")
}


# --------------------------------------------------------------------------
# Input Validation Helpers
# --------------------------------------------------------------------------

#' Validate model argument
#' @noRd
.validate_model <- function(model) {
  valid <- c("rasch", "2pl", "3pl")
  model <- tolower(model)
  model <- match.arg(model, valid)
  model
}


#' Validate prior argument
#' @noRd
.validate_prior <- function(prior) {
  valid <- c("dpm", "normal")
  prior <- tolower(prior)
  prior <- match.arg(prior, valid)

  prior
}


#' Validate parameterization argument
#' @noRd
.validate_parameterization <- function(parameterization, model) {
  valid <- c("irt", "si")
  parameterization <- tolower(parameterization)
  parameterization <- match.arg(parameterization, valid)

  # SI is not meaningful for Rasch

  if (model == "rasch" && parameterization == "si") {
    stop("Slope-Intercept (SI) parameterization is not meaningful for ",
         "the Rasch model (discrimination is fixed at 1). ",
         "Use parameterization = 'irt'.",
         call. = FALSE)
  }

  parameterization
}


#' Resolve identification default based on model
#' @noRd
.resolve_identification <- function(identification, model, prior) {
  if (is.null(identification)) {
    # Model-specific defaults (Blueprint D3)
    identification <- switch(model,
      "rasch" = "constrained_item",
      "2pl"   = "unconstrained",
      "3pl"   = "unconstrained"
    )
  }

  valid <- c("constrained_item", "constrained_ability", "unconstrained")
  identification <- match.arg(identification, valid)

  # Enforce invalid combination rules (Blueprint Section 3.8)
  if (identification == "constrained_ability" && prior == "dpm") {
    stop("Cannot use 'constrained_ability' with DPM prior. ",
         "Constraining theta ~ N(0,1) defeats the purpose of the DPM prior.",
         call. = FALSE)
  }

  if (identification == "constrained_item" && model == "3pl") {
    stop("'constrained_item' identification is not implemented for 3PL. ",
         "Use 'unconstrained' (recommended) or 'constrained_ability'.",
         call. = FALSE)
  }

  identification
}


#' Validate model x prior x identification combination
#' @noRd
.validate_model_combination <- function(model, prior, parameterization,
                                        identification) {
  # All checks already done in individual validators, but this provides

  # a central enforcement point
  invisible(TRUE)
}


# --------------------------------------------------------------------------
# Data Validation and Transformation
# --------------------------------------------------------------------------

#' Validate and prepare input data
#'
#' @param data A matrix or data.frame of binary responses, or a long-format
#'   data.frame with columns (person, item, response).
#' @param data_format One of "auto", "matrix", "long".
#' @return A list with \code{y} (N x I matrix), \code{N}, \code{I},
#'   and \code{data_format}.
#' @noRd
.validate_data <- function(data, data_format = "auto") {
  data_format <- match.arg(data_format, c("auto", "matrix", "long"))

  if (data_format == "auto") {
    data_format <- .detect_data_format(data)
  }

  if (data_format == "matrix") {
    y <- .prepare_matrix_data(data)
  } else {
    y <- .prepare_long_data(data)
  }

  list(y = y, N = nrow(y), I = ncol(y), data_format = data_format)
}


#' Detect data format automatically
#' @noRd
.detect_data_format <- function(data) {
  if (is.matrix(data)) {
    return("matrix")
  }
  if (is.data.frame(data)) {
    # Long format: typically has 3 columns (person, item, response)
    if (ncol(data) <= 4 && nrow(data) > max(dim(data)) * 2) {
      return("long")
    }
    return("matrix")
  }
  stop("Data must be a matrix or data.frame.", call. = FALSE)
}


#' Prepare matrix-format data
#' @noRd
.prepare_matrix_data <- function(data) {
  y <- as.matrix(data)
  storage.mode(y) <- "double"

  # Check binary
  unique_vals <- unique(as.vector(y[!is.na(y)]))
  if (!all(unique_vals %in% c(0, 1))) {
    stop("Response data must be binary (0/1). Found values: ",
         paste(sort(unique_vals), collapse = ", "),
         call. = FALSE)
  }

  # Check dimensions
  if (nrow(y) < 2) {
    stop("Response matrix must have at least 2 persons (rows).",
         call. = FALSE)
  }
  if (ncol(y) < 2) {
    stop("Response matrix must have at least 2 items (columns).",
         call. = FALSE)
  }

  # Check for rows/cols that are all NA
  row_all_na <- apply(y, 1, function(x) all(is.na(x)))
  col_all_na <- apply(y, 2, function(x) all(is.na(x)))
  if (any(row_all_na)) {
    warning(sum(row_all_na), " person(s) have all-NA responses and will ",
            "produce uninformative posteriors.", call. = FALSE)
  }
  if (any(col_all_na)) {
    stop(sum(col_all_na), " item(s) have all-NA responses. ",
         "Remove these columns before fitting.", call. = FALSE)
  }

  # Check for items with zero variance (all 0s or all 1s)
  col_means <- colMeans(y, na.rm = TRUE)
  degenerate <- col_means == 0 | col_means == 1
  if (any(degenerate)) {
    warning(sum(degenerate), " item(s) have zero variance (all 0 or all 1). ",
            "These may cause estimation problems.", call. = FALSE)
  }

  y
}


#' Prepare long-format data (convert to matrix)
#' @noRd
.prepare_long_data <- function(data) {
  # Expects columns: person, item, response (in that order or by name)
  if (is.data.frame(data)) {
    # Try to identify columns by name
    names_lower <- tolower(names(data))
    person_col <- which(names_lower %in% c("person", "person_id", "respondent",
                                            "subject", "id"))
    item_col <- which(names_lower %in% c("item", "item_id", "question"))
    resp_col <- which(names_lower %in% c("response", "resp", "y", "score"))

    if (length(person_col) >= 1 && length(item_col) >= 1 &&
        length(resp_col) >= 1) {
      person_col <- person_col[1]
      item_col <- item_col[1]
      resp_col <- resp_col[1]
    } else {
      # Fall back to positional
      if (ncol(data) < 3) {
        stop("Long-format data must have at least 3 columns ",
             "(person, item, response).", call. = FALSE)
      }
      person_col <- 1
      item_col <- 2
      resp_col <- 3
    }

    # Convert to wide matrix
    persons <- sort(unique(data[[person_col]]))
    items   <- sort(unique(data[[item_col]]))
    y <- matrix(NA_real_, nrow = length(persons), ncol = length(items))

    person_idx <- match(data[[person_col]], persons)
    item_idx   <- match(data[[item_col]], items)
    y[cbind(person_idx, item_idx)] <- as.double(data[[resp_col]])

    colnames(y) <- as.character(items)
    rownames(y) <- as.character(persons)

    return(.prepare_matrix_data(y))
  }

  stop("Long-format data must be a data.frame.", call. = FALSE)
}


# --------------------------------------------------------------------------
# Seed Management
# --------------------------------------------------------------------------

#' Set seed reproducibly for NIMBLE
#' @noRd
.set_seed <- function(seed) {
  if (!is.null(seed)) {
    set.seed(seed)
  }
}


# --------------------------------------------------------------------------
# Timer Utilities
# --------------------------------------------------------------------------

#' Start a timer
#' @noRd
.start_timer <- function() {
  proc.time()
}


#' Stop a timer and return elapsed time in seconds
#' @noRd
.elapsed_time <- function(start_time) {
  as.numeric((proc.time() - start_time)["elapsed"])
}


# --------------------------------------------------------------------------
# Message Utilities
# --------------------------------------------------------------------------

#' Print a verbose message if verbose is TRUE
#' @noRd
.vmsg <- function(..., verbose = TRUE) {
  if (verbose) {
    message(paste0(...))
  }
}


#' Format time in human-readable form
#' @noRd
.format_time <- function(seconds) {
  if (seconds < 60) {
    sprintf("%.1f sec", seconds)
  } else if (seconds < 3600) {
    sprintf("%.1f min", seconds / 60)
  } else {
    sprintf("%.1f hr", seconds / 3600)
  }
}
