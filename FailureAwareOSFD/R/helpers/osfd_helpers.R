# depends: 
#' @export
osfd_fun <- function(name, fallback = NULL) {
  ns <- asNamespace("OSFD")
  if (exists(name, envir = ns, mode = "function")) {
    get(name, envir = ns)
  } else {
    fallback
  }
}

# List likely helper/internal function names in your installed OSFD
osfd_list_helpers <- function() {
  ns <- asNamespace("OSFD")
  ls(ns, pattern = "(spanfill|fill|dist|perturb|ball|EI|UCB|nn|cpp)")
}


get_osfd_fn <- function(name) {
  if (!exists(name, envir = asNamespace("OSFD"), mode = "function")) {
    stop(sprintf("OSFD function '%s' not found in namespace.", name))
  }
  getFromNamespace(name, "OSFD")
}
