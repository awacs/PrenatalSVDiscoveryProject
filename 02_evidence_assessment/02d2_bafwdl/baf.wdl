import "baf_indv.wdl" as Main
workflow Wrap{
    String DIR
    String OUTDIR
    String BAFscript
    String Batch_key
    String Baffile
    call scatt{input:dir=DIR}
    scatter(bed in scatt.output_chunks){
        call Main.Main{input:Input=bed,BAFSCRIPT=BAFscript,BATCH_KEY=Batch_key,BAFFILE=Baffile}
    }
    # call Main.Main{input:Input=DIR,BAFSCRIPT=BAFscript,BATCH_KEY=Batch_key,BAFFILE=Baffile}
    call merge{input: FILES=Main.bafresult,outdir=OUTDIR}
    output{
        Array[File] metrics=merge.metrics
    }
}

task scatt{
    String dir
    command{
        ln -s ${dir} Beds
    }
  runtime {
    sla: "-sla miket_sc"
    queue: "short"
    memory: "4 GB"
    cpu: "1"
  }
    output{
        Array[File] output_chunks = glob("Beds/*")
    }
}
task merge{
    Array[File] FILES
    String outdir
    command<<<
        mkdir bafmetrics
        cp {${sep="," FILES}} bafmetrics
        cp {${sep="," FILES}} ${outdir}
    >>>
  runtime {
    sla: "-sla miket_sc"
    queue: "short"
    memory: "4 GB"
    cpu: "1"
  }
    output{
        Array[File] metrics = glob("bafmetrics/*")
    }
}