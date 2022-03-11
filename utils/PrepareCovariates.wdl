version 1.0

task prepare_covariates {
	input {
		File covariate_file
		Array[String] covariate_list
		Array[String] individuals_list
	}

	File covariate_list_file = write_lines(covariate_list)
	File individual_list_file = write_lines(individuals_list)

	command <<<
	# need to quote the heredoc word so that the '$' in the R code don't get misunderstood.

	Rscript --vanilla -<<"EOF" "~{covariate_file}" "~{covariate_list_file}" "~{individual_list_file}" 


	args <- commandArgs(trailingOnly = TRUE)

	covariate_file <- args[1]
	cat(sprintf("got the following arguments:\n"))
	cat(sprintf("args[1] (covariate file): %s\n", covariate_file))

	covariate_list <- read.delim(args[2], header=FALSE)[[1]]
	cat(sprintf("args[2] (covariate list): %s\n", paste(collapse=",", covariate_list)))
	
	individual_list <- read.delim(args[3], header=FALSE)[[1]]
	
	cat(sprintf("args[3] (participants list): %s\n", paste(collapse=",", individual_list)))
	
	covariates <- read.delim(covariate_file, header = TRUE, sep = '\t')
	
	cat(names(covariates))
	cat("\n")

	covariates <- subset(covariates, subset = BQCID %in% individual_list, select=c("BQCID", covariate_list))
	

	# make sure that all requested samples are present:

	if(nrow(covariates) != length(individual_list)) {
		stop(sprintf("Got different number of individuals than requested: %d vs. %d.\n Samples requested that were not returned are: %s%n", 
			nrow(covariates), length(individual_list), setdiff(individual_list, covariates$BQCID)))
	}

	#debug
	
	# prepare for transpose
	rownames(covariates) <- covariates$BQCID 
	covariates <- subset(covariates, select=-BQCID)
	rotated <- t(covariates)
	rotated <- as.data.frame(rotated)
	rotated$ID <- rownames(rotated)
	# move ID to the first column
	rotated <- rotated[,c(ncol(rotated), seq(ncol(rotated) - 1))]
	write.table(x = rotated, file="extracted_covariates.tsv", sep="\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

	EOF

	>>>
	output {
		File prepared_covariates = "extracted_covariates.tsv"
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
