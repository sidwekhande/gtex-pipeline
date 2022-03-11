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
			nrow(covariates), length(individual_list), setdiff(individual_list,covariates$BQCID)))
	}

	#debug
	print(covariates)
	cat(covariates$BQCID)
	which(is.na(covariates$BQCID))

	# prepare for transpose
	rownames(covariates) <- covariates$BQCID 
	print("1")
	covariates <- subset(covariates,select=-BQCID)
	print("2")
	rotated <- t(covariates)
	print("3")
	rotated <- as.data.frame(rotated)
	print("4")
	rotated$ID <- rownames(rotated)
	# move ID to the first column
	print("5")
	rotated[,c(ncol(rotated),seq(ncol(rotated)-1))]
	print("6")
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
