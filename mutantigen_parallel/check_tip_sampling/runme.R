library(phytools)
library(tidyverse)

runnum <- "5408"
outdir <- "outfiles"

mafile <- file.path(outdir, paste0("out_", runnum,
                                   ".out.simmapAntigenic"))
nxfile <- file.path(outdir, paste0("out_", runnum,
                                   ".simmapAntigenic.nex"))
tpfile <- file.path(outdir, paste0("out_", runnum, ".tips"))
csfile <- file.path(outdir, paste0("out_", runnum, ".timeseries"))

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

tree_to_nexus(mafile, nxfile)

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
            select(name, year) %>%
            arrange(year)

    # cases
    case <- read_tsv(casefile) %>%
            select(date, totalCases)

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
           summarize(num_seq = n()) %>%
           full_join(case) %>%
           select(Time = date, Cases = totalCases, NumSeq = num_seq) %>%
           arrange(Time)

    return(ans)
}

ans <- tips_to_bins(tpfile, nxfile, csfile)
write_csv(ans, paste0("counts", runnum, ".csv"))

#--------------------------------------------------
# 3. Report whether the results seems "good enough"
#--------------------------------------------------
# I just made up a 90% rule.  Might be too strict?

# drop the first rows, where tips are not yet sampled
ans1 <- ans
while(anyNA(ans1[1,]) == T)
    ans1 <- ans1[-1,]

nt <- nrow(ans1)

# how many times have at least 5 cases?
c5 <- ans1 %>% filter(Cases > 4) %>% nrow
c5 / nt > 0.9 # want this to be true

# how many times have at least 5 sequences?
s5 <- ans1 %>% filter(NumSeq > 4) %>% nrow
s5 / nt > 0.9 # want this to be true
