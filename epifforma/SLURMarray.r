# This is code to submit array jobs to the SLURM cluster
# Code developed by Lars Fritsche, modified by Lauren Beesley

slurmarray <- function(cmdLines,sname="Array",stime="60",smem="4G",
         soutdir=paste0(getwd(),"/temp/"),soutput="%a_%N_%j.log",sparallel=500,
         serror="%a_%N_%j.err.txt",scpus=1,extraOption=character(0))  {
    if(!file.exists(soutdir)) dir.create(soutdir)
    logdir <- paste0(soutdir,"arraylog")
    if(!file.exists(logdir)) dir.create(logdir)

    XTRA <- character(0)
    if(length(extraOption)>0){
        XTRA <- c(XTRA,paste0("#SBATCH",extraOption))
    }   


	file.prefix <- paste0(soutdir,sname,"_",format(Sys.time(), format="%Y%m%d_%H%M%S"))

    cmdfile <- paste0(file.prefix,".cmd.txt")
    write(cmdLines,cmdfile)    
    
    
    if(sparallel>0){
    	arrays <- paste0("#SBATCH --array=1-",length(cmdLines),"%",as.integer(sparallel))
    } else {
    	arrays <- paste0("#SBATCH --array=1-",length(cmdLines))
    }
    
    arrayLines <- c("#!/bin/bash -l",
        paste0("#SBATCH --error=",gsub('~','.',logdir),"/",sname,"_",serror),
        paste0("#SBATCH --output=",gsub('~','.',logdir),"/",sname,"_",soutput),
        paste0("#SBATCH --job-name=",sname),
        paste0("#SBATCH --time=",stime),
        paste0("#SBATCH --mem=",smem),
        paste0("#SBATCH --cpus-per-task=",scpus),
        XTRA,
        arrays,
        paste0("module load R"),
        paste0("bash -c \"$(head -n $SLURM_ARRAY_TASK_ID ",cmdfile," | tail -n 1)\"")
        )
    
    afile <- paste0(file.prefix,".sh")
    write(arrayLines,afile)

	Sys.sleep(1)

    system(paste("chmod u+rwx",afile))
    cat("Run SLURM array with the following command:\n")
    print(paste("sbatch",afile))
    return(paste("sbatch",afile))
}
