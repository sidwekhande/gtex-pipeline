version 1.0

workflow Cutadapt{
  input {
    String adapter_sequence1
    String adapter_sequence2
    String sample_id
    File pretrim_fastq1
    File pretrim_fastq2
  }
  
  call cutadapt{
    input:
    adapter_sequence1 = adapter_sequence1,
    adapter_sequence2 = adapter_sequence2,
    sample_id = sample_id,
    pretrim_fastq1 = pretrim_fastq1,
    pretrim_fastq2 = pretrim_fastq2
  }
}

task cutadapt {
  input { 
    String adapter_sequence1
    String adapter_sequence2
    String sample_id
    File pretrim_fastq1
    File pretrim_fastq2
  }

  command <<<
    
    curl -L -O https://github.com/broadinstitute/palantir-workflows/raw/main/Scripts/monitoring/cromwell_monitoring_script.sh 
    chmod a+x cromwell_monitoring_script.sh 
    ./cromwell_monitoring_script.sh &
        
    cutadapt -a ~{adapter_sequence1} -A ~{adapter_sequence2} \
    -o ~{sample_id}.trimmed.R1.fastq.gz -p ~{sample_id}.trimmed.R2.fastq.gz \
    ~{pretrim_fastq1} ~{pretrim_fastq2} \
    -q 20 --paired --phred33
  >>>
  
  output{
    File posttrim1 = sample_id +".trimmed.R1.fastq.gz"
    File posttrim2 = sample_id +".trimmed.R2.fastq.gz"
  }
  runtime {
       docker: "mskaccess/trim_galore:0.6.3"
       memory: "8 GB"
       cpu: 4
       disks: "local-disk 200 HDD"
  }
}