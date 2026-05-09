library(phytools)
library(tidyverse)
library(this.path)
setwd(paste0(this.path::here(), "/../"))
source("R/emma_helpers.R")
runnum <- "16"
outdir <- paste0(this.path::here(), '/../outfiles')

mafile <- file.path(outdir, paste0("out_", runnum,
                                   ".out.simmapAntigenic"))
nxfile <- file.path(outdir, paste0("out_", runnum,
                                   ".simmapAntigenic.nex"))
tpfile <- file.path(outdir, paste0("out_", runnum, ".tips"))
csfile <- file.path(outdir, paste0("out_", runnum, ".timeseries"))

tree_to_nexus(mafile, nxfile)
ans <- tips_to_bins(tpfile, nxfile, csfile)
write_csv(ans, paste0("data_post_processing/counts", runnum, ".csv"))

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
