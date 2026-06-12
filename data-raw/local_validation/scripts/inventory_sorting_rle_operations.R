# Inventory sorting/RLE/sample operations in public R code.
# Does not read private data.

repo_root <- "D:/GitHub/FLightR"
output_root <- "D:/GitHub/FLightR/data-raw/local_validation/outputs"
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)

patterns <- c(
  "sort.int", "sort\\(", "rle\\(", "inverse.rle", "tabulate",
  "sample\\(", "sample.int", "rowProds", "cbind", "rbind",
  "get.transition.rle", "Transitions.rle", "Points.rle"
)

r_files <- list.files(file.path(repo_root, "R"), pattern = "\\.R$", full.names = TRUE, recursive = TRUE)

find_function <- function(lines, line_no) {
  candidates <- seq_len(line_no)
  candidates <- candidates[grepl("^[[:space:]]*[A-Za-z0-9_.]+[[:space:]]*<-[[:space:]]*function\\b", lines[candidates])]
  if (length(candidates) == 0) return(NA_character_)
  line <- lines[tail(candidates, 1)]
  sub("^[[:space:]]*([A-Za-z0-9_.]+)[[:space:]]*<-.*$", "\\1", line)
}

classify_occurrence <- function(fun, expr) {
  text <- paste(fun, expr)
  if (grepl("generate.points.dirs|pf.run.parallel.SO.resample|Current.Points|Last.State|New.Points|Transitions", text)) {
    return(c("live propagation / particle-history-sensitive", "yes", "no", "Do not replace with counts unless one-to-one particle history alignment is proven preserved."))
  }
  if (grepl("get.transition.rle|Points.rle|Transitions.rle", text)) {
    return(c("stored RLE output only", "no for stored marginal output; yes if feeding later particle identity", "maybe", "Candidate for tabulate/counting only after confirming output ordering and downstream decode semantics."))
  }
  if (grepl("get.coordinates.PF|get_ZI_distances|dist.fun|estimate.movement.parameters", text)) {
    return(c("post-processing only", "no live propagation identity, but result semantics matter", "maybe", "Profile first; optimize only aggregate/post-filter summaries."))
  }
  c("unknown, needs manual inspection", "unknown", "unknown", "Inspect manually before changing.")
}

rows <- list()
for (file in r_files) {
  lines <- readLines(file, warn = FALSE)
  for (i in seq_along(lines)) {
    expr <- lines[i]
    hits <- patterns[vapply(patterns, function(p) grepl(p, expr), logical(1))]
    if (length(hits) == 0) next
    fun <- find_function(lines, i)
    cls <- classify_occurrence(fun, expr)
    rows[[length(rows) + 1]] <- data.frame(
      file = gsub("\\\\", "/", file),
      line_number = i,
      function_name = fun,
      matched_pattern = paste(hits, collapse = ";"),
      expression = trimws(expr),
      classification = cls[1],
      particle_identity_history_must_be_preserved = cls[2],
      tabulate_counting_optimization_safe = cls[3],
      recommended_action = cls[4],
      stringsAsFactors = FALSE
    )
  }
}

inventory <- if (length(rows) == 0) {
  data.frame()
} else {
  do.call(rbind, rows)
}

out_path <- file.path(output_root, "sorting_rle_inventory.csv")
utils::write.csv(inventory, out_path, row.names = FALSE)

summary <- if (nrow(inventory)) {
  aggregate(line_number ~ classification, inventory, length)
} else {
  data.frame(classification = character(), line_number = integer())
}
names(summary) <- c("classification", "n_occurrences")
print(summary)
message("Inventory saved to: ", out_path)
