###############################################################################
# File:           utils.R
# Purpose:        Utility/helper functions for common tasks across Nextflow R scripts
#
# Description:
#   This file defines general-purpose functions used by multiple pipeline steps,
#   such as logging, safe file reading, error handling, etc.
#
# Usage:
#   Include this file in  R scripts via:
#       source(file.path(Sys.getenv("SCRIPT_DIR", "."), "utils.R"))
#
# Notes:
#   - Do not run this R file directly.
#   - Keep functions general and side-effect-free if possible.
#   - If adding new functions, document inline.
###############################################################################

#' @export
#' Log a message with a timestamp and level
#'
#' Prints a formatted log message to the console including a timestamp and a log level.
#'
#' @param msg A character string containing the message to log.
#' @param level A character string specifying the log level (e.g., "INFO", "ERROR", "DEBUG"). Default is "INFO".
#'
#' @return NULL. Called for side effects (prints to console).
#'
#' @examples
#' log_message("Starting analysis")
#' log_message("File not found", level = "ERROR")
log_message = function(msg, level = "INFO") {
  timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  cat(sprintf("[%s] [%s] %s\n", timestamp, level, msg))
  flush.console()
}

#' @export
#' Log a message with a timestamp and level
#'
#' Prints a formatted log message to the console including a timestamp and a log level.
#'
#' @param msg A character string containing the message to log.
#' @param level A character string specifying the log level (e.g., "INFO", "ERROR", "DEBUG"). Default is "INFO".
#' @param con A connection (e.g., a log file) specifying where the log should be written. Default: stdout().
#'
#' @return NULL. Called for side effects (prints to console).
#'
#' @examples
#' log_msg("Starting analysis")
#' log_msg("File not found", level = "ERROR")
# log_msg = function(msg, level = "INFO") {
#   timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S")
#   cat(sprintf("[%s] [%s] %s\n", timestamp, level, msg))
# }
log_msg = function(msg, level = "INFO", con = stdout()) {

  timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  
  # ANSI color codes
  colors = c(
    "START" = "\033[33m",
    "INFO" = "\033[97m",
    "PLOT" = "\033[38;5;240m",
    "WARN" = "\033[38;5;208m",
    "ERROR" = "\033[31m",
    "DEBUG" = "\033[34m",
    "MSG" = "\033[77m",
    "DONE" = "\033[32m"
  )

  # reset color
  reset = "\033[0m"

  level = toupper(level)
  level_fmt = sprintf("%5s", level)
  
  # Pick color for the level, default to white
  color = colors[[level]]
  if (is.null(color)) color = "\033[37m"

  # Convert tables/data frames to text
  if (is.data.frame(msg) || is.matrix(msg)) {
    msg = paste(capture.output(print(msg, row.names = FALSE)), collapse = "\n")
  }

  prefix = sprintf("%s[%s] [%s] ", color, timestamp, level_fmt)
  msg = gsub("\n", paste0(reset, "\n", prefix), msg, fixed = TRUE)
  
  # line = sprintf(
  #   "%s[%s] [%s] %s%s\n",
  #   color, timestamp, level_fmt, msg, reset
  # )
  line = sprintf(
    "%s%s%s\n",
    prefix, msg, reset
  )
  
  # Print colored message
  if (!is.null(con)) {
    cat(line, file = con, append = TRUE)
    flush(con)
  }
}


#' @export
#' Plot a PDF and PNG of the same image
#'
#' Attempts to plot a PDF and PNG. If it fails, logs the error and stops.
#'
#' @param plot A plot object
#' @param filename The name (without file extension) to be used for the output files
#' @param height The height of the plot in inches (cm)
#' @param width The width of the plot in inches (cm)
#' @return NULL (side effect: saves files)
#' @importFrom readr read_csv
plot_pdf_and_png = function(plot, filename, height, width, res=600) {
  tryCatch({
    # Write PDF
    pdf_file = paste0(filename, '.pdf')
    grDevices::pdf(pdf_file, height = height, width = width)
    print(plot)
    dev.off()

    # Write PNG
    png_file = paste0(filename, '.png')
    grDevices::png(png_file, height = height, width = width, res = res, units = 'in')
    print(plot)
    dev.off()

    log_message(sprintf("Saved plot to %s and %s", pdf_file, png_file), level = "INFO")
  },
    error = function(e) {
      log_message(sprintf("Could not write plot: %s", filename), level = "ERROR")
      log_message(e$message, level = "ERROR")
      stop(e)
    }
  )
}

#' @export
#' Replace one or more whitespace characters with the replacement
#'
cfn = function(
  x,
  replacement = "_"
) {  
  gsub(, replacement, trimws(x))
}


#' @export
#' Read a CSV or TSV file automatically
#'
#' Attempts to read a file based on its extension. Supports CSV (.csv) and TSV (.tsv).
#' Allows specifying whether the file contains a header row.
#' Logs the process and handles errors gracefully.
#'
#' @param file_path Path to the input file
#' @param has_header Logical indicating whether the file contains a header row (default TRUE)
#' @return A tibble containing the contents of the file
#' @importFrom readr read_csv read_tsv
read_delim_auto = function(file_path, has_header = TRUE) {
  tryCatch({
    ext = tolower(tools::file_ext(file_path))
    
    # readr's col_names argument corresponds to has_header
    col_names = has_header
    
    if (ext == "csv") {
      data = readr::read_csv(file_path, col_names = col_names, show_col_types = FALSE)
    } else if (ext == "tsv") {
      data = readr::read_tsv(file_path, col_names = col_names, show_col_types = FALSE)
    } else {
      stop(sprintf("Unsupported file extension: '%s'. Only .csv and .tsv are supported.", ext))
    }
    
    log_message(sprintf("Successfully read file: %s", file_path), level = "INFO")
    return(data)
  },
  error = function(e) {
    log_message(sprintf("Could not read file: %s", file_path), level = "ERROR")
    log_message(e$message, level = "ERROR")
    stop(e)
  })
}

#' Combine Sample IDs and Position IDs from Two Tibbles
#'
#' Takes two tibbles of identical dimensions: one containing sample IDs and the other containing position IDs,
#' and combines them into a long-format tibble with two columns, `sample_id` and `position_id`.
#' The function preserves the original row and column ordering to correctly pair sample and position IDs.
#'
#' @param sample_tbl A tibble containing sample IDs arranged in a grid-like structure.
#' @param position_tbl A tibble containing position IDs arranged in the same grid-like structure as `sample_tbl`.
#' 
#' @return A tibble with two columns:
#'   \describe{
#'     \item{sample_id}{Sample identifier from `sample_tbl`}
#'     \item{position_id}{Position identifier from `position_tbl`}
#'   }
#'   The rows correspond to the matched sample-position pairs in the original tibbles.
#' 
#' @importFrom dplyr mutate select inner_join row_number
#' @importFrom tidyr pivot_longer
#' 
#' @examples
#' sample_tbl = tibble::tribble(
#'   ~A, ~B, ~C,
#'   "S1", "S2", "S3",
#'   "S4", "S5", "S6"
#' )
#' position_tbl = tibble::tribble(
#'   ~A,  ~B,  ~C,
#'   "P1", "P2", "P3",
#'   "P4", "P5", "P6"
#' )
#' combine_sample_position(sample_tbl, position_tbl)
#'
#' @export
combine_sample_position = function(sample_tbl, position_tbl) {
  sample_long = sample_tbl %>%
    dplyr::mutate(row_id = dplyr::row_number()) %>%
    tidyr::pivot_longer(
      cols = -row_id,
      names_to = "col",
      values_to = "sample_id"
    )
  
  position_long = position_tbl %>%
    dplyr::mutate(row_id = dplyr::row_number()) %>%
    tidyr::pivot_longer(
      cols = -row_id,
      names_to = "col",
      values_to = "position_id"
    )
  
  combined = dplyr::inner_join(sample_long, position_long, by = c("row_id", "col")) %>%
    dplyr::select(sample_id, position_id)
  
  return(combined)
}

#' @title Safely Run an Expression with Logging and Error Handling
#'
#' @description
#' A utility wrapper to safely evaluate an expression with informative logging,
#' success/failure messages, and optional suppression of return values.
#'
#' This function is useful for consistent logging and error handling around
#' I/O or data-processing steps that may fail. On error, the function logs
#' the error message and stops execution.
#'
#' @param expr An unevaluated expression (use \code{quote(...)}) to be executed.
#' @param start_msg A character string to log before running the expression.
#' @param success_msg A character string to log upon successful completion.
#' @param error_msg A character string to log if an error occurs.
#' @param silent Logical; if \code{TRUE}, suppresses returning the result
#'        (useful for write operations that would print to console).
#'
#' @return The result of evaluating \code{expr}, unless \code{silent = TRUE}, in which case \code{NULL} is invisibly returned.
#'
#' @examples
#' \dontrun{
#' safe_run(
#'   expr = quote(readr::write_tsv(data, "output.tsv")),
#'   start_msg = "Writing data...",
#'   success_msg = "Write complete.",
#'   error_msg = "Failed to write data.",
#'   silent = TRUE
#' )
#' }
#'
#' @export
safe_run = function(expr, start_msg, success_msg, error_msg, silent = FALSE) {
  tryCatch({
    log_message(start_msg, level = "INFO")
    result = eval(expr)
    log_message(success_msg, level = "INFO")
    if (!silent) return(result) else invisible(NULL)
  }, error = function(e) {
    log_message(error_msg, level = "ERROR")
    log_message(e$message, level = "ERROR")
    stop(e)
  })
}

# TODO: ADD DOCUMENTATION
# Optional: Validate file extensions
check_extension = function(path, valid_ext) {
  ext = tolower(tools::file_ext(path))
  tolower(ext) %in% tolower(valid_ext)
}


# TODO: add documentation
#'
#' bins = c(0, 60, 120, 365, 730, 1825, Inf)
#' labels = c("1-60","60-120","120-365","365-730","730-1825","1825+")
#' df_new = bin_ages_df(df, "Age", bins = bins, labels = labels)
#' df_new
bin_ages_df = function(
  df,
  age_col,
  bins,
  labels,
  new_days_col = "age_days",
  new_bin_col = "age_bin")
{
  library(dplyr)
  library(stringr)
  # Check labels
  if (length(labels) != length(bins) - 1) {
    stop("Number of labels must be one less than number of bins!")
  }
  
  age_to_days = function(x) {
    if (is.na(x) || str_trim(x) == "") return(NA_real_)
    
    # Clean input
    x = str_to_lower(x)
    x = str_replace_all(x, "\\?|\\s", "")
    
    years = 0
    months = 0
    days = 0
    
    # Years
    if (str_detect(x, "y")) {
      years = as.numeric(str_extract(x, "\\d+(?=y)"))
    }
    
    # Months immediately after years (e.g., 14y10m)
    months_after_y = str_extract(x, "(?<=y)\\d+(?=m)")
    if (!is.na(months_after_y)) {
      months = as.numeric(months_after_y)
    }
    # Months with "mo" (standalone months)
    else if (str_detect(x, "mo")) {
      months = as.numeric(str_extract(x, "\\d+(?=mo)"))
    }
    # Standalone months "m" not after y
    else if (str_detect(x, "m") & !str_detect(x, "y")) {
      months = as.numeric(str_extract(x, "\\d+(?=m)"))
    }
    
    # Days
    if (str_detect(x, "d")) {
      days = as.numeric(str_extract(x, "\\d+(?=d)"))
    }
    
    # Total days
    return(years*365 + months*30 + days)
  }
  
  df = df %>%
    mutate(
      !!new_days_col := sapply(.data[[age_col]], age_to_days),
      !!new_bin_col := ifelse(
        is.na(.data[[new_days_col]]),
        NA_character_,
        as.character(
          cut(.data[[new_days_col]],
              breaks = bins,
              labels = labels,
              right = TRUE,
              include.lowest = TRUE)
        )
      )
    )
  
  return(df)
}


#' Validate Command-Line Arguments
#'
#' Checks that required arguments are provided, files exist, and file extensions are valid.
#'
#' @param opt A list of command-line options, typically returned by `optparse::parse_args`.
#' @param required_args Character vector of argument names that must be provided.
#' @param file_args Character vector of argument names that are expected to be file paths.
#' @param ext_rules Named list where names are argument names and values are character vectors
#'                  of allowed file extensions (without dot). Example: `list(input="csv")`.
#' @param opt_parser Optional, an `optparse` parser object. If provided, `print_help` will be called when missing arguments are detected.
#'
#' @return Invisibly returns TRUE if all checks pass.
#' @export
#'
#' @examples
#' \dontrun{
#' required_args <- c("input_seurat")
#' file_args <- c("input_seurat", "meta_file")
#' ext_rules <- list(
#'   input_seurat = c("qs2"),
#'   meta_file    = c("csv", "tsv")
#' )
#' validate_args(opt, required_args, file_args, ext_rules, opt_parser)
#' }
validate_args <- function(opt, required_args = NULL, file_args = NULL, ext_rules = list(), opt_parser = NULL) {
  
  # 1. Check for missing required arguments
  if (!is.null(required_args)) {
    missing_args <- required_args[!sapply(required_args, function(arg) !is.null(opt[[arg]]))]
    if (length(missing_args) > 0) {
      if (!is.null(opt_parser)) print_help(opt_parser)
      stop(sprintf("Missing required argument(s): %s\n", paste(missing_args, collapse = ", ")), call. = FALSE)
    }
  }
  
  # 2. Check that specified files exist
  if (!is.null(file_args)) {
    for (file_arg in file_args) {
      file_path <- opt[[file_arg]]
      if (!file.exists(file_path)) {
        stop(sprintf("File does not exist: --%s = %s", file_arg, file_path), call. = FALSE)
      }
    }
  }
  
  # 3. Check file extensions
  check_extension <- function(path, valid_ext) {
    ext <- tolower(tools::file_ext(path))
    tolower(ext) %in% tolower(valid_ext)
  }
  
  if (length(ext_rules) > 0) {
    for (arg_name in names(ext_rules)) {
      valid_ext <- ext_rules[[arg_name]]
      file_path <- opt[[arg_name]]
      if (!check_extension(file_path, valid_ext)) {
        stop(sprintf("The %s file must have one of these extensions: %s", arg_name, paste(valid_ext, collapse = ", ")), call. = FALSE)
      }
    }
  }
  
  invisible(TRUE)
}


check_opts <- function(
  opt,
  opt_parser,
  required = character(0),
  validators = list()
) {
  errors <- character(0)

  # required / non-null checks
  missing <- required[vapply(required, function(x) is.null(opt[[x]]), logical(1))]
  if (length(missing) > 0) {
    errors <- c(errors, sprintf("Missing required argument(s): %s", paste0("--", missing, collapse = ", ")))
  }

  # custom validators, e.g. list(input_file = file.exists, resolutions = jsonlite::validate)
  for (name in names(validators)) {
    val <- opt[[name]]
    if (!is.null(val) && !isTRUE(validators[[name]](val))) {
      errors <- c(errors, sprintf("Invalid value for --%s: %s", name, val))
    }
  }

  if (length(errors) > 0) {
    print_help(opt_parser)
    stop(paste(errors, collapse = "\n"), call. = FALSE)
  }
}

apply_json_mappings = function(df, column_to_match, json_file) {

    mapping_list = jsonlite::fromJSON(json_file)

    for (new_col in names(mapping_list)) {

        map_def  = mapping_list[[new_col]]
        method   = map_def$method
        fallback = map_def$fallback %||% NULL
        levels   = map_def$levels %||% NULL

        # ---- APPLY METHOD ----
        if (method == "pattern") {

            map_df = as.data.frame(map_def$mappings)

            df = df %>%
                dplyr::mutate(
                    !!new_col := purrr::map_chr(.data[[column_to_match]], function(x) {
                        hit = map_df$value[stringr::str_detect(x, map_df$pattern)]
                        if (length(hit) > 0) hit[1] else NA_character_
                    })
                )

            # default level order from mappings if not provided
            if (is.null(levels)) {
                levels = map_df$value
            }

        } else if (method == "exact") {

            map_df = readr::read_tsv(map_def$file, col_types = readr::cols())

            id_col    = map_def$id_col
            value_col = map_def$value_col

            if (!value_col %in% colnames(map_df)) {
                stop(glue::glue("Column '{value_col}' not found in {map_def$file}"))
            }

            map_df = map_df %>%
                dplyr::rename(value = !!rlang::sym(value_col))

            df = df %>%
                dplyr::left_join(map_df, by = setNames(id_col, column_to_match))

            df[[new_col]] = df$value
            df$value = NULL

            # default level order from observed values if not provided
            if (is.null(levels)) {
                levels = unique(na.omit(df[[new_col]]))
            }

        } else if (method == "relevel_only") {

            if (is.null(fallback)) {
                stop(glue::glue("'relevel_only' requires a fallback column"))
            }

            df[[new_col]] = df[[fallback]]

            # levels MUST be provided here
            if (is.null(levels)) {
                stop(glue::glue("'levels' must be provided for relevel_only ({new_col})"))
            }

        } else if (method == "run_function") {

            fn_name = map_def$function_name
            params  = map_def$params %||% list()

            fn = get(fn_name, mode = "function")

            df[[new_col]] = do.call(fn, c(list(df[[column_to_match]]), params))

            if (is.null(levels)) {
                levels = unique(na.omit(df[[new_col]]))
            }

        } else {
            stop(glue::glue("Unknown method: {method}"))
        }

        # ---- APPLY FALLBACK ----
        if (!is.null(fallback) && fallback %in% colnames(df)) {
            df[[new_col]] = dplyr::coalesce(df[[new_col]], df[[fallback]])
        }

        # ---- APPLY LEVELS (UNIVERSAL) ----
        if (!is.null(levels)) {
            df[[new_col]] = factor(df[[new_col]], levels = levels)
        }

    }

    return(df)
}

convert_age = \(x, to="months") {
  res = str_match(x, "(?:(\\d+)y)?(?:(\\d+)m)?")[, -1] %>%
    apply(1, \(x) {
      years = as.numeric(x[1]); if(is.na(years)) years = 0
      months = as.numeric(x[2]); if(is.na(months)) months = 0
      if(to=="months") years*12 + months else years + months/12
    })
  return(res)
}
