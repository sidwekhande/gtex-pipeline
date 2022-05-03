version 1.0
import "write_array.wdl" as wa

task RenameVcfSamples{
	input {
		File vcf
		File rename_by
	}
	command <<<
		set -euo pipefail 
		bcftools reheader -s ~{rename_by}  -o output.vcf.gz 
		bcftools index -t output.vcf.gz 

	>>>
	output {
		File vcf_out="output.vcf.gz"
		File vcf_index_out="output.vcf.gz.tbi"
	}

	runtime {
		docker: "bschiffthaler/bcftools:latest"
		preemptible: 0
		disks: "local-disk " + (2*ceil(size(vcf,"GiB"))+20) + " HDD"
		bootDiskSizeGb: "16"
		memory: 20 + " GB"
	}
}


workflow RenameVcfSamplesWF{
	input {
		File vcf_in
		Array[String] samples_in_vcf
		Array[String] new_sample_names
	}

	call wa.write_array_to_tsv as array{
		input:
			array=[samples_in_vcf,new_sample_names]
	}

	call RenameVcfSamples{
		input:
			vcf=vcf_in,
			rename_by=array.tsv
	}
}