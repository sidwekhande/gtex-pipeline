task echo_task {
	String message

	command {
		echo ${message}
	}

	runtime { 
		docker: "broadgdac/firecloud-ubuntu:15.10"
	}

	output {
		String echoed=read_string(stdout())
	}
}

workflow echo {
	call echo_task
}