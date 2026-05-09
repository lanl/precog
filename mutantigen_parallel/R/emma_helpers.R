# Helper functions from Emma posterior outfile analysis.

#--------------------------------------------------
# 1. Modify the Mutantigen output file so that R can read it in
#--------------------------------------------------

tree_to_nexus <- function(infile, outfile)
{
  cat("# nexus\nbegin smptrees;\n\ttree = (", file=outfile, append=F)
  system(paste0("cat ", infile, " | tr -d '\n' >> ", outfile))
  cat(");\nend;\n", file=outfile, append=T)
  
  # there's an annoying inconsistency about extra surrounding (...)
  phy <- read.simmap(outfile)
  res <- try(print(phy), silent=T)
  if (class(res) == "try-error")
  {
    cat("# nexus\nbegin smptrees;\n\ttree = ", file=outfile,
        append=F)
    system(paste0("cat ", infile, " | tr -d '\n' >> ", outfile))
    cat(";\nend;\n", file=outfile, append=T)
  }
}

#--------------------------------------------------
# 2. Count tip sampling over time
#--------------------------------------------------

tips_to_bins <- function(tipfile, treefile, casefile)
{
  # all reported tip sampling times
  tips <- read_csv(tipfile) %>%
    distinct()
  
  # the ones actually on the tree
  # (it seems like we should be able to rely on tipfile, but that is not
  #   always true)
  phy <- read.simmap(treefile)
  tiplab <- tibble(name = phy$tip.label)
  tips <- tips %>%
    right_join(tiplab) %>%
    dplyr::select(name, year) %>%
    arrange(year)
  
  # cases
  case <- read_tsv(casefile) %>%
    dplyr::select(date, totalCases)
  
  # define time bins
  breaks <- case$date # every 7 days
  #breaks <- case$date[c(F,T)] # every two weeks
  
  # apply the bins
  # (un-tidy because of mysterious bug when using mutate here)
  case$bin <- cut(case$date, breaks=breaks)
  case$time <- as.integer(case$bin)
  tips$bin <- cut(tips$year, breaks=breaks)
  tips$time <- as.integer(tips$bin)
  
  ans <- tips %>%
    group_by(time) %>%
    dplyr::summarize(num_seq = dplyr::n()) %>%
    full_join(case) %>%
    dplyr::select(Time = date, Cases = totalCases, NumSeq = num_seq) %>%
    arrange(Time)
  
  return(ans)
}