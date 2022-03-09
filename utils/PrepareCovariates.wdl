version 1.0

task prepare_covariates {
	input {
		File covariate_file
		Array[String] covariate_list
		Array[String] individuals_list
		String output_file
	}

	File covariate_list_file = write_lines(covariate_list)
	File individual_list_file = write_lines(individuals_list)

	command <<<
	Rscript --vanilla -<<EOF "~{covariate_file}" "~{covariate_list_file}" "~{individual_list_file}" "~{output_file}"


	args <- commandArgs(trailingOnly = TRUE)

	covariate_file <- args[1]
	sprintf("got the following arguments:")
	sprintf("args[1] (covariate file): %s", covariate_file)

	covariate_list <- read.delim(args[2], header=FALSE)[[1]]
	sprintf("args[2] (covariate list): %s", paste(collapse=",", covariate_list))
	
	individual_list <- read.delim(args[3], header=FALSE)[[1]]
	
	sprintf("args[3] (participants list): %s", paste(collapse=",", individual_list))
	
	covariates <- read.delim(covariate_file, header = TRUE, sep = '\t')
	covariates <- subset(covariates, subset= ID %in% individual_list, select=c("ID", covariate_list))
	

	# make sure that all requested samples are present:

	if(nrow(covariates) != length(individual_list)) {
		stop(sprintf("Got different number of individuals than requested: %d vs. %d", nrow(covariates), length(individual_list))
	}

	#prepare for transpose
	rownames(covariates) <- covariates$BQCID 
	covariates <- select(covariates,-BQCID)
	rotated <- t(covariates)
	rotated <- as.data.frame(rotated)
	rotated$ID <- rownames(rotated)

	write.table(x = rotated,file="extracted_coveriates.tsv", sep="\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

	EOF

	>>>
	output {
		File prepared_covariates = output_file
	}

	runtime {
		docker: "r-base:latest"
		memory: "2GB"
		disks: "local-disk 20 HDD"
	}

}


workflow PrepareCovariates{
	call prepare_covariates{}

	output {
		File covariates = prepare_covariates.prepared_covariates
	}
}
