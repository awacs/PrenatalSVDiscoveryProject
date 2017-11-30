workflow Main{
    File Input
    String BAFSCRIPT
    String BATCH_KEY
    String BAFFILE
    call split{input:FILE=Input}
    # call Scatter.work {input: FILES=split.FILES}
    scatter (Chunks in split.output_chunks){
        call BAF{input: FILE=Chunks,BAFscript=BAFSCRIPT,Batch_key=BATCH_KEY,Baffile=BAFFILE}
    }
    call merge{input:FILES=BAF.BAFout,FILE=Input}
    output{
        File bafresult=merge.out
    }
}
task split{
    File FILE
    String prefix = basename(FILE)
    command{
        mkdir splits
        sed '/^#/d' ${FILE} \
          | sort -R \
          | split -a 4 -l 1000 - splits/${prefix}
    }
    output{
        Array[File] output_chunks = glob("splits/*")
    }
  runtime {
    sla: "-sla miket_sc"
    queue: "short"
    memory: "4 GB"
    cpu: "1"
  }
}
task BAF{
    File FILE
    String BAFscript
    String Batch_key
    String Baffile
    command{
        python ${BAFscript} ${FILE} ${Baffile} --batch ${Batch_key}> test.out
    }
    output{
        File BAFout ="test.out"
    }
  runtime {
    sla: "-sla miket_sc"
    queue: "short"
    memory: "4 GB"
    cpu: "1"
  }
}
task merge{
    Array[File] FILES
    File FILE
    String Name = basename(FILE,".bed")
    command{
        cat ${sep=" " FILES} > ${Name}.metrics
    }
  runtime {
    sla: "-sla miket_sc"
    queue: "short"
    memory: "4 GB"
    cpu: "1"
  }
  output {
    File out="${Name}.metrics"
  }
}