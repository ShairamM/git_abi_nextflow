nextflow.enable.dsl=2

//params.accession = "SRR1777174"
params.accession = "${projectDir}/sra_accessions.txt"
params.with_fastqc = false  // use false to stop

// Fastp trimming parameters
//params.cut_window_size = 4
//params.cut_mean_quality = 20
//params.length_required = 50
//params.average_quality = 20

process prefetch_sra {
   
    container "https://depot.galaxyproject.org/singularity/sra-tools%3A3.2.1--h4304569_1"
    publishDir "${projectDir}/SRA_files", mode:'copy', overwrite:true

    input:
        val accession
    output:
        path("${accession}/${accession}.sra") //tuple val(accession), 
    """
    prefetch ${accession}
    """
}

process FASTQ_dump {
    container "https://depot.galaxyproject.org/singularity/sra-tools%3A3.2.1--h4304569_1"
    publishDir "${projectDir}/FASTQ_files", mode:'copy', overwrite:true

    input:
        path(inputfile)  //tuple val(accession), 
    
    output:
        path("*.fastq")

    """
    fastq-dump --split-3 ${inputfile}
    """
}

process stats {
    container "https://depot.galaxyproject.org/singularity/ngsutils%3A0.5.9--py27h9801fc8_5"
    publishDir "${projectDir}/FASTQ_files", mode:'copy', overwrite:true
    
    input:
        path inputfile
    
    output:
        path "*.txt"

    """
    fastqutils stats ${inputfile} > ${inputfile.getSimpleName()}_stats.txt
    """

}

process FASTQC {
    container "https://depot.galaxyproject.org/singularity/fastqc%3A0.12.1--hdfd78af_0"
    publishDir "${projectDir}/FASTQ_files/FASTQC_files", mode:'copy', overwrite:true

    input:
        path inputfile
    
    output:
        path "${inputfile.getSimpleName()}_fastqc.zip"        //path "*_fastqc.{zip,html}", emit: fastqc_results
        path "${inputfile.getSimpleName()}_fastqc.html" 

    """
    fastqc ${inputfile}
    """

}

process trimming {
    container "https://depot.galaxyproject.org/singularity/fastp%3A1.0.1--heae3180_0"
    publishDir "${projectDir}/FASTQ_files/trimmed_files", mode:'copy', overwrite:true

    input:
        path(inputfile)
    
    output: 
        
        //path "${inputfile.getSimpleName()}_trimmed.fastq"
        //path "${inputfile.getSimpleName()}_trimmed.html"   
        path "${inputfile.getSimpleName()}_trimmed.json"  

    """
    fastp \
          --in1 "${inputfile}" --out1 "${inputfile.getSimpleName()}_trimmed.fastq" --json "${inputfile.getSimpleName()}_trimmed.json"
             
    """
}

    process multiQC {
    container "https://depot.galaxyproject.org/singularity/multiqc%3A1.29--pyhdfd78af_0"
    publishDir "${projectDir}/FASTQ_files/multiQC", mode:'copy', overwrite:true

    input:
        path(inputfile)
    
    output: 
         
        path "multiqc_*"  

    """
    multiqc .
             
    """

}
/*--json "${accession}_trimmed.json" --html "${accession}_trimmed.html" \
            --cut_front --cut_tail \
          --cut_window_size ${params.cut_window_size} \
          --cut_mean_quality ${params.cut_mean_quality} \
          --length_required ${params.length_required} \
          --average_qual ${params.average_quality}
*/

workflow {  
    accession_ch = Channel.fromPath(params.accession).splitText().map { it.trim() }
    prefetch_sra_ch = prefetch_sra(accession_ch)
    FASTQ_dump_ch = FASTQ_dump(prefetch_sra_ch)
    //fastq_files_ch = FASTQ_dump_ch.flatten()
    //fastq_paths_only_ch = FASTQ_dump_ch.flatMap { accession, files -> files } 
    
    stats(FASTQ_dump_ch)
        // This process is now conditional based on the --with_fastqc parameter
    if (params.with_fastqc) {
        FASTQC(FASTQ_dump_ch)
    }

    trim_ch = trimming(FASTQ_dump_ch)
    multiQC(trim_ch.collect())

    
}    
 /* 
   // Pair R1 + R2 per accession
    paired_fastq_ch = FASTQ_dump_ch
    .map { accession, files ->
        def read1 = files.find { it.name =~ /_1\.fastq$/ }
        def read2 = files.find { it.name =~ /_2\.fastq$/ }
        tuple(accession, read1, read2)
    }
    trimming_ch = trimming(paired_fastq_ch)
    */