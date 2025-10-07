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
        path "${inputfile.getSimpleName()}_fastqc.zip", emit: fastqc_zip        //path "*_fastqc.{zip,html}", emit: fastqc_results
        path "${inputfile.getSimpleName()}_fastqc.html", emit: fastqc_html

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
        
        path "${inputfile.getSimpleName()}_trimmed.fastq", emit: trimmed_fastq   
        //path "${inputfile.getSimpleName()}_trimmed.html"   
        path "${inputfile.getSimpleName()}_trimmed.json", emit: json_report  

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

    process spades {
    container "https://depot.galaxyproject.org/singularity/spades%3A4.2.0--haf24da9_0"
    publishDir "${projectDir}/FASTQ_files/spades", mode:'copy', overwrite:true

    input:
        path inputfile  // trimmed fastq

    output:
        path "${inputfile.getSimpleName()}_assembly", emit: assembly_dir
        path "${inputfile.getSimpleName()}_assembly/contigs.fasta", emit: contigs
    """
    spades.py -s "${inputfile}" -o ${inputfile.getSimpleName()}_assembly --isolate
    """
    }

    process bandage {
    container "https://depot.galaxyproject.org/singularity/bandage%3A0.9.0--h9948957_0"
    publishDir "${projectDir}/FASTQ_files/spades/bandage", mode:'copy', overwrite:true

    input:
        path inputfile  // assembly dir contains fastg

    output:
        path "${inputfile.getSimpleName()}_assembly_graph.png"
    """
    Bandage image ${inputfile}/assembly_graph.fastg ${inputfile.getSimpleName()}_assembly_graph.png
    """
    }

    process QUAST {
    container "https://depot.galaxyproject.org/singularity/quast%3A5.3.0--py313pl5321h5ca1c30_2"
    publishDir "${projectDir}/FASTQ_files/spades/quast_QC", mode:'copy', overwrite:true

    input:
        path inputfile  // takes fasta file

    output:
         path "quast_output/report.html", emit: quast_html
         path "quast_output", emit: quast_dir
    """
    quast.py "${inputfile}" -o quast_output
    """
    }

workflow {  
    accession_ch = Channel.fromPath(params.accession).splitText().map { it.trim() }
    prefetch_sra_ch = prefetch_sra(accession_ch)
    FASTQ_dump_ch = FASTQ_dump(prefetch_sra_ch)
    //fastq_files_ch = FASTQ_dump_ch.flatten()
    //fastq_paths_only_ch = FASTQ_dump_ch.flatMap { accession, files -> files } 
    
    stats(FASTQ_dump_ch)
        // This process is now conditional based on the --with_fastqc parameter
    if (params.with_fastqc) {
        fastqc_ch = FASTQC(FASTQ_dump_ch)
    }
    fastqc_report_ch = fastqc_ch.fastqc_html 

    trim_ch = trimming(FASTQ_dump_ch)
    // from emit
    trimmed_fastq_ch = trim_ch.trimmed_fastq
    trimmed_json_ch = trim_ch.json_report

    

    spades_ch = spades(trimmed_fastq_ch)
    // from emit
    assembly_dir_ch = spades_ch.assembly_dir     // contains .fastg and more
    contigs_ch = spades_ch.contigs        // contains .fasta
    
    bandage(assembly_dir_ch)

    QUAST_ch = QUAST(contigs_ch)
    // 1. Initialize the combined channel with the first report channel (Fix #2)
    all_reports_ch = fastqc_report_ch
    all_reports_ch = all_reports_ch.mix(trimmed_json_ch, QUAST_ch)
    multiQC_ch = all_reports_ch.collect()
    multiQC(multiQC_ch)

    
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