require(phytools)
require(stringr)
require(dplyr)
require(readr)

#--------------------------------------------------
# Get all existing run numbers
#--------------------------------------------------
# relying on the "out_NNN.summary" files

outdir <- file.path("../", "outfiles")

tmp <- basename(Sys.glob(file.path(outdir, "out_*.summary")))
all_runnum <- sort(as.integer(str_extract(tmp, "(\\d)+")))
#str_extract(tmp, "(?<=out_)(\\d)+(?=.summary)") # fancy :)

# Um, turns out not all these runs produced complete output.
# Will deal with that below.

#--------------------------------------------------
# Tree file format conversion
#--------------------------------------------------
# mutantigen output tree files are in almost-newick format,
# phytools needs them in nexus format

tree_to_nexus <- function(infile, outfile)
{
    cat("# nexus\nbegin smptrees;\n\ttree = (", file=outfile, append=F)
    system(paste0("cat ", infile, " | tr -d '\n' >> ", outfile))
    cat(");\nend;\n", file=outfile, append=T)
}

for (runnum in all_runnum)
{
    # trait = fitness load
    tree <- file.path(outdir,
               paste0("out_", runnum, ".simmapLoad"))
    if (file.exists(tree))
    {
        nex <- paste0(tree, ".nex")
        if (!file.exists(nex))
            tree_to_nexus(tree, nex)

        # trait = antigenic type
        # (only bother doing this one if load output exists)
        tree <- file.path(outdir,
                   paste0("out_", runnum, ".out.simmapAntigenic"))
        nex <- paste0(tree, ".nex")
        if (!file.exists(nex))
            tree_to_nexus(tree, nex)

        message(paste("did", runnum))
    }
}

#--------------------------------------------------
# Summarize the values of the load and antigenicity traits
#--------------------------------------------------

get_info <- function(nex)
{
    phy <- read.simmap(nex)
    x <- as.integer(unlist(lapply(phy$maps, function(x) names(x))))
    return(c("min" = min(x), "max" = max(x)))
}

get_nex_info <- function(runnum)
{
    nex <- file.path(outdir,
             paste0("out_", runnum, ".simmapLoad.nex"))

    if (file.exists(nex))
    {
        ans_load <- get_info(nex)

        nex <- file.path(outdir,
                 paste0("out_", runnum, ".out.simmapAntigenic.nex"))
        ans_anti <- get_info(nex)

        ans <- data.frame(rbind(ans_load, ans_anti))
        rownames(ans) <- NULL
        ans$run <- runnum
        ans$type <- c("load", "anti")

        message(paste("did", runnum))
    } else {
        ans <- rep(NA, 4)
        names(ans) <- c("min", "max", "run", "type")
    }

    return(ans)
}

ans0 <- lapply(all_runnum, get_nex_info)

ans <- bind_rows(ans0) %>%
       filter(!is.na(run)) %>%
       select(run, type, everything())

write_csv(ans, file=file.path(outdir, "trait_summaries.csv"))
