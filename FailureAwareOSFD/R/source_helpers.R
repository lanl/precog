# inverse_model/R/load_helpers.R
source_helpers <- function(
  helpers_dir = file.path("R", "helpers"),
  pattern = "\\.R$"
) {
  files <- list.files(helpers_dir, pattern = pattern, full.names = TRUE)
  if (!length(files)) stop("No helper files found in: ", helpers_dir)

  # Map: basename -> full path (assumes unique basenames)
  basenames <- basename(files)
  if (any(duplicated(basenames))) {
    dups <- unique(basenames[duplicated(basenames)])
    stop("Duplicate helper filenames found (must be unique): ",
         paste(dups, collapse = ", "))
  }
  file_map <- setNames(files, basenames)

  # Parse dependency header:  # depends: a.R, b.R
  deps_of <- function(path) {
    lines <- readLines(path, warn = FALSE, n = 50) # header is near top
    hit <- grep("^\\s*#\\s*depends\\s*:", lines, value = TRUE)
    if (!length(hit)) return(character())
    dep_str <- sub("^\\s*#\\s*depends\\s*:\\s*", "", hit[1])
    deps <- trimws(strsplit(dep_str, ",")[[1]])
    deps <- deps[deps != ""]
    # normalize: allow "a" or "a.R"
    deps <- ifelse(grepl("\\.R$", deps), deps, paste0(deps, ".R"))
    deps
  }

  deps <- lapply(file_map, deps_of)

  # Check missing deps
  declared <- unique(unlist(deps))
  missing <- setdiff(declared, names(file_map))
  if (length(missing)) {
    stop("Missing helpers referenced in depends: ",
         paste(missing, collapse = ", "))
  }

  # Topological sort with cycle detection
  state <- setNames(rep(0L, length(file_map)), names(file_map)) # 0=unseen,1=visiting,2=done
  order <- character()

  visit <- function(f) {
    if (state[[f]] == 1L) stop("Circular dependency detected involving: ", f)
    if (state[[f]] == 2L) return()
    state[[f]] <<- 1L
    for (d in deps[[f]]) visit(d)
    state[[f]] <<- 2L
    order <<- c(order, f)
  }

  for (f in names(file_map)) visit(f)

  # Source in order, but skip RcppExports.R (C++ functions loaded via sourceCpp instead)
  # RcppExports.R creates namespace-prefixed .Call() references that cause
  # "namespace 'epiOutputExpansion' is not available" warnings when loading saved objects
  for (f in order) {
    if (f == "RcppExports.R") {
      if (verbose_skip <- getOption("source_helpers.verbose", TRUE)) {
        cat("Skipping RcppExports.R (C++ functions already loaded via sourceCpp)\n")
      }
      next
    }
    source(file_map[[f]], local = FALSE)
  }

  invisible(list(order = order, files = file_map, deps = deps))
}
