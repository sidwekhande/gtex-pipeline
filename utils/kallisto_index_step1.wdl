version 1.0

workflow Kallisto_index_step1{
  input {
    File reference
    String name
  }
  
  call Kallisto_index{
    input:
    reference = reference,
    name=name
  }
}

task Kallisto_index {
  input { 
    File reference
    String name
  }

  command <<<
    kallisto index -i ~{name}	~{reference}  
  >>>
  output{
    File index = name
  }
  runtime {
       docker: "jjkrc/kallisto:0.46.1"
       memory: "8 GB"
       cpu: 4
       disks: "local-disk 200 SSD"
  }
}
