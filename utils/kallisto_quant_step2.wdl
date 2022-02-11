version 1.0

workflow Kallisto_quant_step2{
  input {
    File TranscriptIdx
    File read1
    File read2
    String sampleID
    String BootstrapIter
    
    Int? memoryMaybe
  }
  
  output {
  
    File Kallisto_quant.info
    File Kallisto_quant.tsv
    File Kallisto_quant.h5

  }
  
  call Kallisto_quant{
    input:
    TranscriptIdx = TranscriptIdx,
    read1 = read1,
    read2 = read2,
    name = sampleID,
    BootstrapIter = BootstrapIter,
    memoryMaybe = memoryMaybe

  }
}

task Kallisto_quant {
  input { 
    File TranscriptIdx
    File read1
    File read2
    String name
    String BootstrapIter
    
    Int? memoryMaybe
  }
  
  Int memoryDefault=1
  Int memoryJava=select_first([memoryMaybe,memoryDefault])
  Int memoryRam=memoryJava+2
  Int disk_size = 10 + ceil(size([read1], "GB")) + ceil(size([read2], "GB"))

  command <<<
    wget https://github.com/broadinstitute/palantir-workflows/raw/main/Scripts/monitoring/cromwell_monitoring_script.sh 
    chmod a+x cromwell_monitoring_script.sh 
    ./cromwell_monitoring_script.sh &
    
    kallisto quant -i ~{TranscriptIdx} -o ~{name} -b ~{BootstrapIter} ~{read1} ~{read2} 
  >>>
  output{
    File info = name + ".run_info.json"
    File tsv = name + ".abundance.tsv"
    File h5 = name + "abundance.h5"
  }
  runtime {
       docker: "jjkrc/kallisto:0.46.1"
       memory: memoryRam + " GB"
       cpu: 4
       disks: "local-disk " + disk_size + " HDD"
  }
}
