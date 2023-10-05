version 1.0

workflow Kallisto_quant_step2{
  input {
    File TranscriptIdx
    File read1
    File read2
    File gtf
    String sampleID
    String BootstrapIter
    Int disk_size
    Int? memoryMaybe
  }
  
  call Kallisto_quant{
    input:
    TranscriptIdx = TranscriptIdx,
    read1 = read1,
    read2 = read2,
    name = sampleID,
    BootstrapIter = BootstrapIter,
    gtf = gtf,
    disk_size = disk_size,
    memoryMaybe = memoryMaybe

  }
  
   output {
    File info = Kallisto_quant.info
    File tsv = Kallisto_quant.tsv
    File h5 = Kallisto_quant.h5
    File pseudobam = Kallisto_quant.pseudobam
    
  }
}

task Kallisto_quant {
  input { 
    File TranscriptIdx
    File read1
    File read2
    String name
    String BootstrapIter
    File gtf
    Int disk_size = 60 
    Int? memoryMaybe
  }

  
  Int memoryDefault=1
  Int memoryJava=select_first([memoryMaybe,memoryDefault])
  Int memoryRam=memoryJava+2

  command <<<
    set -xeuo pipefail 
    
    wget https://github.com/broadinstitute/palantir-workflows/raw/main/Scripts/monitoring/cromwell_monitoring_script.sh 
    chmod a+x cromwell_monitoring_script.sh 
    ./cromwell_monitoring_script.sh &
    
    kallisto quant -i ~{TranscriptIdx} -o ~{name} -b ~{BootstrapIter} ~{read1} ~{read2} -t 4 -g ~{gtf} --genomebam
  >>>
  
  output{
    File info = name + "/run_info.json"
    File tsv = name + "/abundance.tsv"
    File h5 = name + "/abundance.h5"
    File pseudobam = name + "/pseudoalignments.bam"
  }
  runtime {
       docker: "jjkrc/kallisto:0.46.1"
       memory: memoryRam + " GB"
       cpu: 4
       disks: "local-disk " + disk_size + " HDD"
  }
}
