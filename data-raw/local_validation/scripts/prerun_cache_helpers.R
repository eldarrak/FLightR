# Local-only prerun object cache helpers for validation/benchmark scripts.
# Cache files live under data-raw/local_validation/outputs/prerun_cache/ and must not be committed.

truthy_env <- function(name) {
  tolower(Sys.getenv(name, unset = "false")) %in% c("1", "true", "yes", "y")
}

get_prerun_cache_dir <- function(output_root = "D:/GitHub/FLightR/data-raw/local_validation/outputs") {
  cache_dir <- file.path(output_root, "prerun_cache")
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  normalizePath(cache_dir, winslash = "/", mustWork = TRUE)
}

hash_object <- function(x) {
  if (requireNamespace("digest", quietly = TRUE)) {
    return(digest::digest(x, algo = "sha256"))
  }
  tf <- tempfile(fileext = ".rds")
  on.exit(unlink(tf), add = TRUE)
  saveRDS(x, tf)
  unname(tools::md5sum(tf))
}

file_fingerprint <- function(path) {
  info <- file.info(path)
  data.frame(
    path = normalizePath(path, winslash = "/", mustWork = FALSE),
    size = if (nrow(info)) info$size else NA_real_,
    mtime = if (nrow(info)) format(info$mtime, "%Y-%m-%d %H:%M:%S %Z") else NA_character_,
    md5 = if (file.exists(path)) unname(tools::md5sum(path)) else NA_character_,
    stringsAsFactors = FALSE
  )
}

make_prerun_cache_key <- function(components) {
  hash <- hash_object(components)
  list(hash = hash, components = components)
}

write_prerun_metadata <- function(path, key, object_summary, cache_status) {
  metadata_path <- sub("\\.rds$", "_metadata.csv", path)
  scalarize <- function(x) {
    if (length(x) == 0) return(NA_character_)
    paste(as.character(x), collapse = ";")
  }
  components <- key$components
  rows <- data.frame(
    name = names(components),
    value = vapply(components, scalarize, character(1)),
    stringsAsFactors = FALSE
  )
  header <- data.frame(
    name = c("cache_key", "cache_file", "created", names(cache_status), names(object_summary)),
    value = c(key$hash, path, format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
              vapply(cache_status, scalarize, character(1)),
              vapply(object_summary, scalarize, character(1))),
    stringsAsFactors = FALSE
  )
  utils::write.csv(rbind(header, rows), metadata_path, row.names = FALSE)
  metadata_path
}

read_prerun_metadata <- function(path) {
  metadata_path <- sub("\\.rds$", "_metadata.csv", path)
  if (!file.exists(metadata_path)) return(NULL)
  utils::read.csv(metadata_path, stringsAsFactors = FALSE)
}

summarize_prerun_object <- function(x) {
  list(
    grid_size = if (!is.null(x$Spatial$Grid)) nrow(x$Spatial$Grid) else NA_integer_,
    main_index_rows = if (!is.null(x$Indices$Main.Index)) nrow(x$Indices$Main.Index) else NA_integer_,
    matrix_index_rows = if (!is.null(x$Indices$Matrix.Index.Table)) nrow(x$Indices$Matrix.Index.Table) else NA_integer_,
    start_point = if (!is.null(x$Spatial$start.point)) x$Spatial$start.point else NA_integer_,
    stop_point = if (!is.null(x$Spatial$stop.point)) x$Spatial$stop.point else NA_integer_
  )
}

validate_cached_prerun <- function(x, metadata, key, expected) {
  if (is.null(metadata)) return("metadata_missing")
  key_row <- metadata$value[metadata$name == "cache_key"]
  if (length(key_row) != 1 || !identical(key_row, key$hash)) return("key_mismatch")
  summary <- summarize_prerun_object(x)
  if (!is.null(expected$grid_size) && is.finite(expected$grid_size) && !identical(as.integer(summary$grid_size), as.integer(expected$grid_size))) {
    return("grid_size_mismatch")
  }
  if (!is.null(expected$start_point) && is.finite(expected$start_point) && !identical(as.integer(summary$start_point), as.integer(expected$start_point))) {
    return("start_point_mismatch")
  }
  if (!is.null(expected$stop_point) && is.finite(expected$stop_point) && !identical(as.integer(summary$stop_point), as.integer(expected$stop_point))) {
    return("stop_point_mismatch")
  }
  "ok"
}

load_or_build_prerun <- function(key, build_fun, expected = list(), output_root = "D:/GitHub/FLightR/data-raw/local_validation/outputs") {
  cache_dir <- get_prerun_cache_dir(output_root)
  cache_file <- file.path(cache_dir, paste0("prerun_", substr(key$hash, 1, 24), ".rds"))
  reuse <- truthy_env("FLIGHTR_REUSE_PRERUN")
  force_rebuild <- truthy_env("FLIGHTR_FORCE_REBUILD_PRERUN")
  save_after_build <- reuse || truthy_env("FLIGHTR_SAVE_PRERUN")

  load_seconds <- 0
  build_seconds <- 0
  save_seconds <- 0
  cache_hit <- FALSE
  warning_message <- NA_character_
  value <- NULL

  if (reuse && !force_rebuild && file.exists(cache_file)) {
    load_start <- proc.time()[["elapsed"]]
    candidate <- try(readRDS(cache_file), silent = TRUE)
    load_seconds <- proc.time()[["elapsed"]] - load_start
    metadata <- read_prerun_metadata(cache_file)
    if (!inherits(candidate, "try-error")) {
      valid <- validate_cached_prerun(candidate, metadata, key, expected)
      if (identical(valid, "ok")) {
        value <- candidate
        cache_hit <- TRUE
      } else {
        warning_message <- paste("Cached prerun object invalid; rebuilding:", valid)
      }
    } else {
      warning_message <- paste("Cached prerun object could not be read; rebuilding:", conditionMessage(attr(candidate, "condition")))
    }
  }

  if (is.null(value)) {
    build_start <- proc.time()[["elapsed"]]
    value <- build_fun()
    build_seconds <- proc.time()[["elapsed"]] - build_start
    if (save_after_build) {
      save_start <- proc.time()[["elapsed"]]
      saveRDS(value, cache_file)
      object_summary <- summarize_prerun_object(value)
      cache_status <- list(
        reuse_enabled = reuse,
        force_rebuild = force_rebuild,
        save_after_build = save_after_build,
        cache_hit = FALSE,
        warning = warning_message
      )
      metadata_path <- write_prerun_metadata(cache_file, key, object_summary, cache_status)
      save_seconds <- proc.time()[["elapsed"]] - save_start
    } else {
      metadata_path <- NA_character_
    }
  } else {
    metadata_path <- sub("\\.rds$", "_metadata.csv", cache_file)
  }

  list(
    value = value,
    cache = list(
      enabled = reuse,
      hit = cache_hit,
      file = cache_file,
      key = key$hash,
      metadata_file = metadata_path,
      build_seconds = build_seconds,
      load_seconds = load_seconds,
      save_seconds = save_seconds,
      effective_seconds = if (cache_hit) load_seconds else build_seconds + save_seconds,
      warning = warning_message
    )
  )
}
