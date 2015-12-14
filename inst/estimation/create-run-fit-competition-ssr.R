cores_req <- "3"
mem_req <- "1000"
time_req <- "4:00"
queue_req <- "short"
lsfoutfilename <- "lag1.out"

files_path <- file.path("/home", "ngr67a", "ssr-influenza-competition")

requestCmds <- "#!/bin/bash\n"

requestCmds <- paste0(requestCmds, "#BSUB -n ", cores_req, " # how many cores we want for our job\n",
	"#BSUB -R span[hosts=1] # ask for all the cores on a single machine\n",
	"#BSUB -R rusage[mem=", mem_req, "] # ask for memory\n",
	"#BSUB -o ", lsfoutfilename, " # log LSF output to a file\n",
	"#BSUB -W ", time_req, " # run time\n",
	"#BSUB -q ", queue_req, " # which queue we want to run in\n")

phs <- 1:52
datasets <- c("ili_national",
              paste0("ili_region", 1:10,))
for(data_set in datasets) {
    for(prediction_horizon_limit in phs) {
        filename <- paste0(files_path, "/submit-ph", prediction_horizon_limit, "-", data_set, ".sh")
        
        cat(requestCmds, file = filename)
        cat("module load R/3.2.1\n", file = filename, append = TRUE)
        cat(paste0("R CMD BATCH --vanilla \'--args ", data_set, " ", 
                   prediction_horizon_limit,
                   "\' ", files_path, 
                   "/inst/estimation/fit-competition-ssr.R output-competition-ssr-ph", 
                   prediction_horizon_limit, "-", data_set, ".Rout"),
            file = filename, append = TRUE)

        bsubCmd <- paste0("bsub < ", filename)

        system(bsubCmd)
    }
}
