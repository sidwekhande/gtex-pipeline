version 1.0

workflow Kallisto_quant_step2{
  input {
    File reference
    String name
  }
  
  call Kallisto_quant{
    input:
    reference = reference,
    name=name
  }
}

task Kallisto_quant {
  input { 
    File reference
    String name
  }

  command <<<
    wget https://github.com/broadinstitute/palantir-workflows/raw/main/Scripts/monitoring/cromwell_monitoring_script.sh 
    chmod a+x cromwell_monitoring_script.sh 
    ./cromwell_monitoring_script.sh &
    
    kallisto index -i ~{name}	~{reference}  
  >>>
  output{
    File index = name
  }
  runtime {
       docker: "jjkrc/kallisto:0.46.1"
       memory: "8 GB"
       cpu: 4
       disks: "local-disk 200 HDD"
  }
}
