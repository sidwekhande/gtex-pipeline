version 1.0

workflow Cutadapt{
  input {
    String adapter_sequence1
    String sample_id
    File pretrim_fastq1
  }
  
  call cutadapt{
    input:
    adapter_sequence1 = adapter_sequence1,
    sample_id = sample_id,
    pretrim_fastq1 = pretrim_fastq1,
  }
}

task cutadapt {
  input { 
    String adapter_sequence1
    String sample_id
    File pretrim_fastq1
  }

  command <<<
    
    wget https://github.com/broadinstitute/palantir-workflows/raw/main/Scripts/monitoring/cromwell_monitoring_script.sh 
    chmod a+x cromwell_monitoring_script.sh 
    ./cromwell_monitoring_script.sh &
        
    cutadapt -a ~{adapter_sequence1} -o ~{sample_id}.trimmed.R1.fastq.gz \
    ~{pretrim_fastq1} \
    -q 20 --paired --phred33
  >>>
  
  output{
    File posttrim1 = sample_id +".trimmed.R1.fastq.gz"
  }
  runtime {
       docker: "mskaccess/trim_galore:0.6.3"
       memory: "8 GB"
       cpu: 4
       disks: "local-disk 200 HDD"
  }
}
