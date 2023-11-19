import pandas as pd

configfile: "config/files_config.yaml"
configfile: "config/parameter_config.yaml"
    
# envvars: "TMPDIR"

# 2021_09 version:
# - [x] simplify specification of output (target) files
# - [ ] reduce overreliance on specifying input values via wildcards from config file? wilcard matching from the target output should be able to suffice for lots of things
# - [ ] for rules that may or not involve multiple input files (merging fastq or bams...), add a check to see if there are indeed multiple files before unnecessarily running a slow command like cat or merging bams?

## stuff to add: 
# - mutsigcv or some such

wildcard_constraints:
    library="L-\d{4}",
    libcon_method="[a-z0-9]+",

rule all:
    input:
        config['output'],

##########################
# NEB LIBRARIES (NO UMI) #
##########################
rule concatenate_fastq_files:
# only necessary for older run configs where there were multiple fastqs per flowcell*lane*sample, but doing it for everything out of laziness and for sake of standardization
    input:
        all_fastqs = lambda wildcards:config["fastq"]["nebultra2"][wildcards.fastq]
        #         all_fastqs = lambda wildcards:TESTDICT[wildcards.fastq]
    output:
        r1_fastq_output = temp("tmp/{fastq}.nebultra2.R1.fq.gz"),
        r2_fastq_output = temp("tmp/{fastq}.nebultra2.R2.fq.gz")
    params:
        library = lambda wildcards:wildcards.fastq.split('.')[0] # just for cluster_slurm to know library?
    log:
        "logs/preprocessing_and_alignment/concatenate_fastq_files/{fastq}.nebultra2.concatenate_fastq_files.log"
    run:
        r1_fqs,r2_fqs = [ sorted([x for x in input.all_fastqs if re.search(y,x)]) for y in ("_R1_","_R2_") ]
        cmd = "(cat {r1_fqs} > {output.r1_fastq_output}; cat {r2_fqs} > {output.r2_fastq_output}) 2>> {log}"
        shell("echo \""+cmd+"\" > {log}")
        shell(cmd)
        
rule align_to_reference:
    input:
        r1_fq="tmp/{library}_{origlib}_{flowcell}_{lane}.nebultra2.R1.fq.gz",
        r2_fq="tmp/{library}_{origlib}_{flowcell}_{lane}.nebultra2.R2.fq.gz"
    output:
        temp("tmp/{library}_{origlib}_{flowcell}_{lane}.nebultra2.unsorted.bam")
    params:
        reference_genome = config["reference_genome_for_bwa"],
        bwa=config["bwa"],
        samtools=config["samtools"],
        bwa_threads=config["bwa_threads"]
    log:
        "logs/preprocessing_and_alignment/align_to_reference/{library}_{origlib}_{flowcell}_{lane}.nebultra2.align_to_reference.log"
    shell:
        "({params.bwa} mem -t {params.bwa_threads} -M "+\
        "-R '@RG\\tID:{wildcards.library}_{wildcards.origlib}_{wildcards.flowcell}_{wildcards.lane}\\t"+\
        "LB:{wildcards.library}\\tPL:ILLUMINA\\tPU:no_unit\\tSM:{wildcards.library}' "+\
        "{params.reference_genome} {input.r1_fq} {input.r2_fq} | {params.samtools} view -b - > {output}) 2> {log}"

rule sort_bam:
    input:
        "tmp/{library}_{origlib}_{flowcell}_{lane}.nebultra2.unsorted.bam"
    output:
        bam=temp("tmp/{library}_{origlib}_{flowcell}_{lane}.nebultra2.sorted.bam"),
        bai=temp("tmp/{library}_{origlib}_{flowcell}_{lane}.nebultra2.sorted.bai"),
    params:
        samtools=config["samtools"]
    log:
        "logs/preprocessing_and_alignment/sort/{library}_{origlib}_{flowcell}_{lane}.nebultra2.sort.log"
    shell:
        "{params.samtools} sort -o {output.bam} {input} 2> {log} "+\
        "&& {params.samtools} index {output.bam} {output.bai} 2>> {log}"

rule merge_same_sample_bams:
    input:
        bams=lambda wildcards:expand("tmp/{b}.nebultra2.sorted.bam",b=config["bams_to_merge"]['nebultra2'][wildcards.bams_to_merge]),
        indices=lambda wildcards:expand("tmp/{b}.nebultra2.sorted.bai",b=config["bams_to_merge"]['nebultra2'][wildcards.bams_to_merge])
    output:
        temp("tmp/{bams_to_merge}.nebultra2.sorted_and_merged.bam")
    params:
        samtools=config["samtools"]
    log:
        "logs/preprocessing_and_alignment/merge/{bams_to_merge}.nebultra2.merge.log"
    shell:
        "{params.samtools} merge {output} {input.bams} 2> {log}"

        
########################
# IDT LIBRARIES (+UMI) #
########################
rule concatenate_fastq_files_idtxgen:
# only necessary for older run configs where there were multiple fastqs per flowcell*lane*sample, but doing it for everything out of laziness and for sake of standardization
    input:
        all_fastqs = lambda wildcards:config["fastq"]["idtxgen"][wildcards.fastq]
        #         all_fastqs = lambda wildcards:TESTDICT[wildcards.fastq]
    output:
        r1_fastq_output = temp("tmp/{fastq}.idtxgen.R1.fq.gz"),
        r2_fastq_output = temp("tmp/{fastq}.idtxgen.R2.fq.gz")
    log:
        "logs/preprocessing_and_alignment/concatenate_fastq_files/{fastq}.idtxgen.concatenate_fastq_files.log"
    run:
        r1_fqs,r2_fqs = [ sorted([x for x in input.all_fastqs if re.search(y,x)]) for y in ("_R1_","_R2_") ]
        cmd = "(cat {r1_fqs} > {output.r1_fastq_output}; cat {r2_fqs} > {output.r2_fastq_output}) 2>> {log}"
        shell("echo \""+cmd+"\" > {log}")
        shell(cmd)

rule make_unmapped_bam_idtxgen:
    input:
        r1_fq="tmp/{library}_{origlib}_{flowcell}_{lane}.idtxgen.R1.fq.gz",
        r2_fq="tmp/{library}_{origlib}_{flowcell}_{lane}.idtxgen.R2.fq.gz",
    output:
        temp("tmp/{library}_{origlib}_{flowcell}_{lane}.idtxgen.unmapped.bam")
    params:
        picard_jar=config["picard_jar"],
        mem = config["default_gatk_mem"],
    log:
        "logs/preprocessing_and_alignment/make_unmapped_bam/{library}_{origlib}_{flowcell}_{lane}.idtxgen.log",
    shell:
        "java -Xmx{params.mem} -XX:ParallelGCThreads=1 -Djava.io.tmpdir=${{TMPDIR}} -jar {params.picard_jar} "+\
        "FastqToSam F1={input.r1_fq} F2={input.r2_fq} O={output} SM={wildcards.library} "+\
        "RG={wildcards.library}_{wildcards.origlib}_{wildcards.flowcell}_{wildcards.lane} LB={wildcards.library} PL=ILLUMINA PU=no_unit 2> {log}"

rule extract_umi_idtxgen:
    # need to update to make read-structure encoded in files_config on a library by library basis, not hard coded as a string in parameter_config
    # for example, this is how i got tumor frxn in cnvkit command below: purity = lambda wildcards:config['tumor_frxn'][wildcards.library],
    input:
        "tmp/{library}_{origlib}_{flowcell}_{lane}.idtxgen.unmapped.bam"
    output:
        temp("tmp/{library}_{origlib}_{flowcell}_{lane}.idtxgen.unmapped.umi_extracted.bam")
    params:
        fgbio_jar = config["fgbio_jar"],
        mem = config["default_gatk_mem"],
        fgbio_read_structure = config["fgbio_read_structure"],
    log:
        "logs/preprocessing_and_alignment/extract_umi/{library}_{origlib}_{flowcell}_{lane}.idtxgen.log",
    shell:
        "java -Xmx{params.mem} -XX:ParallelGCThreads=1 -Djava.io.tmpdir=${{TMPDIR}} -jar {params.fgbio_jar} "+\
        "ExtractUmisFromBam --input={input} --output={output} --read-structure={params.fgbio_read_structure} "+\
        "--molecular-index-tags=ZA ZB --single-tag=RX 2> {log}"

rule convert_bam_to_fastq_idtxgen:
    input:
        "tmp/{library}_{origlib}_{flowcell}_{lane}.idtxgen.unmapped.umi_extracted.bam"
    output:
        temp("tmp/{library}_{origlib}_{flowcell}_{lane}.idtxgen.fq")
    params:
        picard_jar = config["picard_jar"],
        mem = config["default_gatk_mem"],
    log:
        "logs/preprocessing_and_alignment/convert_bam_to_fastq/{library}_{origlib}_{flowcell}_{lane}.idtxgen.log"
    shell:
        "java -Xmx{params.mem} -XX:ParallelGCThreads=1 -Djava.io.tmpdir=${{TMPDIR}} -jar {params.picard_jar} "+\
        "SamToFastq I={input} F={output} INTERLEAVE=true 2> {log}"

rule align_to_reference_idtxgen:
    input:
        "tmp/{library}_{origlib}_{flowcell}_{lane}.idtxgen.fq",
    output:
        temp("tmp/{library}_{origlib}_{flowcell}_{lane}.idtxgen.tmp.bam")
    params:
        bwa = config["bwa"],
        bwa_threads = config["bwa_threads"],
        bwa_ref = config["reference_genome_for_bwa"], 
    log:
        "logs/preprocessing_and_alignment/align_to_reference/{library}_{origlib}_{flowcell}_{lane}.idtxgen.log"
    shell:
        "{params.bwa} mem -p -t {params.bwa_threads} {params.bwa_ref} {input} > {output} 2> {log}"

rule merge_aln_and_unmapped_bams_idtxgen:
    input:
        unmapped_bam="tmp/{library}_{origlib}_{flowcell}_{lane}.idtxgen.unmapped.umi_extracted.bam",
        mapped_bam="tmp/{library}_{origlib}_{flowcell}_{lane}.idtxgen.tmp.bam"
    output:
        temp("tmp/{library}_{origlib}_{flowcell}_{lane}.idtxgen.aln.bam")
    params:
        picard_jar = config["picard_jar"],
        mem = config["default_gatk_mem"],
        bwa = config["bwa"],
        bwa_threads = config["bwa_threads"],
        bwa_ref = config["reference_genome_for_bwa"],  # made sequence dict for bwa genome so that can be used for this?
        #         gatk_ref = config["reference_genome_for_gatk"],
    log:
        "logs/preprocessing_and_alignment/merge_aln_and_unmapped_bams/{library}_{origlib}_{flowcell}_{lane}.idtxgen.log"
    shell:
        "java -Xmx{params.mem} -XX:ParallelGCThreads=1 -Djava.io.tmpdir=${{TMPDIR}} -jar {params.picard_jar} "+\
        "MergeBamAlignment UNMAPPED={input.unmapped_bam} ALIGNED={input.mapped_bam} O={output} R={params.bwa_ref} SO=coordinate "+\
        "ALIGNER_PROPER_PAIR_FLAGS=true MAX_GAPS=-1 ORIENTATIONS=FR VALIDATION_STRINGENCY=SILENT 2> {log}"
    
rule merge_same_sample_bams_idtxgen:
    input:
        bams=lambda wildcards:expand("tmp/{b}.idtxgen.aln.bam",b=config["bams_to_merge"]['idtxgen'][wildcards.library]),
    output:
        temp("tmp/{library}.idtxgen.sorted_and_merged.bam")
    params:
        samtools=config["samtools"]
    log:
        "logs/preprocessing_and_alignment/merge/{library}.idtxgen.merge.log"
    shell:
        "{params.samtools} merge {output} {input.bams} 2> {log}"
        
###########
# UNIFIED #
###########

# BEGIN REMOVE MOUSE #
# this workflow is modeled after the ConcatRef published strategy and pipeline sent to me by Anna-Lisa Doebley on 10/14/20 #
# Jo S-Y, Kim E, Kim S. Impact of mouse contamination in genomic profiling of patient-derived models and best practice for robust analysis. Genome Biol [Internet]. 2019 Dec 11;20(1):231 #
rule remove_unpaired_reads_from_bam: # perhaps these cause problems in later filtering?
    input:
        "tmp/{library}.{libcon_method}.sorted_and_merged.bam"
    output:
        bam = "results/before_mouse_removal/{library}/{library}.{libcon_method}.bam",
        bai = "results/before_mouse_removal/{library}/{library}.{libcon_method}.bai"
    params:
        samtools=config["samtools"]
    log:
        "logs/remove_mouse_reads/remove_unpaired_reads_from_bam/{library}.{libcon_method}.log"
    shell:
        "({params.samtools} view -f 1 -b -o {output.bam} {input} && {params.samtools} index {output.bam} {output.bai}) 2> {log}"

rule make_unmapped_bam_from_mapped_bam: # need this output to preserve the read metadata, especially for the UMI libraries?
    input:
        bam = "results/before_mouse_removal/{library}/{library}.{libcon_method}.bam",
        bai = "results/before_mouse_removal/{library}/{library}.{libcon_method}.bai"
    output:
        temp( "tmp/{library}.{libcon_method}.sorted_and_merged.no_unpaired.unmapped.bam" )
    params:
        picard_jar=config["picard_jar"],
    log:
        "logs/remove_mouse_reads/make_unmapped_bam_from_mapped_bam/{library}.{libcon_method}.log"
    shell:
        "java -Xmx4G -XX:ParallelGCThreads=1 -Djava.io.tmpdir=${{TMPDIR}} -jar {params.picard_jar} "+\
        "RevertSam I={input.bam} O={output} 2> {log}"

rule make_fastq_from_unmapped_bam:
    input:
        "tmp/{library}.{libcon_method}.sorted_and_merged.no_unpaired.unmapped.bam"
    output:
        temp( "tmp/{library}.{libcon_method}.sorted_and_merged.no_unpaired.unmapped.fq" )
    params:
        picard_jar=config["picard_jar"],
    log:
        "logs/remove_mouse_reads/make_fastq_from_unmapped_bam/{library}.{libcon_method}.log"
    shell:
        "java -Xmx8G -XX:ParallelGCThreads=2 -Djava.io.tmpdir=${{TMPDIR}} -jar {params.picard_jar} "+
        "SamToFastq I={input} F={output} INTERLEAVE=true 2> {log}"
        
rule align_to_concatref:
    input:
        "tmp/{library}.{libcon_method}.sorted_and_merged.no_unpaired.unmapped.fq",
    output:
        temp("tmp/{library}.{libcon_method}.concatref.bam")
    params:
        bwa = config["bwa"],
        bwa_threads = config["bwa_threads"],
        bwa_ref = config["reference_genome_for_concatref"], 
    log:
        "logs/remove_mouse_reads/align_to_concatref/{library}.{libcon_method}.log"
    shell:
        "{params.bwa} mem -p -t {params.bwa_threads} {params.bwa_ref} {input} > {output} 2> {log}"

rule sort_and_index_concatref_bam:
    input:
        "tmp/{library}.{libcon_method}.concatref.bam"
    output:
        bam = temp("tmp/{library}.{libcon_method}.concatref.sorted.bam"),
        bai = temp("tmp/{library}.{libcon_method}.concatref.sorted.bai"),
    params:
        picard_jar=config["picard_jar"],
    log:
        "logs/remove_mouse_reads/sort_and_index_concatref_bam/{library}.{libcon_method}.log"
    shell:
        "java -Xmx15G -XX:ParallelGCThreads=4 -Djava.io.tmpdir=${{TMPDIR}} -jar {params.picard_jar} "+
        "SortSam I={input} O={output.bam} SORT_ORDER=coordinate CREATE_INDEX=true 2> {log}" 

rule remove_mouse_mapped_reads_from_concatref_bam: # both self and mate
    input:
        bam = "tmp/{library}.{libcon_method}.concatref.sorted.bam",
        bai = "tmp/{library}.{libcon_method}.concatref.sorted.bai",
    output:
        temp("tmp/{library}.{libcon_method}.concatref.sorted.cleaned.bam")
    params:
        samtools=config["samtools"],
        concatref_human_chromosome_bed=config['concatref_human_chromosome_bed'],
    log:
        "logs/remove_mouse_reads/remove_mouse_mapped_reads_from_concatref_bam/{library}.{libcon_method}.log"
    shell:
        "( {params.samtools} view -h -L {params.concatref_human_chromosome_bed} {input.bam} | "+
        "awk '$1~/^@/||($3!~/_GRCm38/&&$7!~/_GRCm38/)' | {params.samtools} view -bS - > {output} ) 2> {log}"

rule make_fastq_from_cleaned_bam:
    input:
        "tmp/{library}.{libcon_method}.concatref.sorted.cleaned.bam"
    output:
        temp( "tmp/{library}.{libcon_method}.concatref.sorted.cleaned.fq" )
    params:
        picard_jar=config["picard_jar"],
    log:
        "logs/remove_mouse_reads/make_fastq_from_cleaned_bam/{library}.{libcon_method}.log"
    shell:
        "java -Xmx16G -XX:ParallelGCThreads=4 -Djava.io.tmpdir=${{TMPDIR}} -jar {params.picard_jar} "+
        "SamToFastq I={input} F={output} INTERLEAVE=true 2> {log}"

rule realign_cleaned_bam_to_human:
    input:
        "tmp/{library}.{libcon_method}.concatref.sorted.cleaned.fq",
    output:
        temp("tmp/{library}.{libcon_method}.concatref.sorted.cleaned.realigned.bam")
    params:
        bwa = config["bwa"],
        bwa_threads = config["bwa_threads"],
        bwa_ref = config["reference_genome_for_bwa"], 
    log:
        "logs/remove_mouse_reads/realign_cleaned_bam_to_human/{library}.{libcon_method}.log"
    shell:
        "{params.bwa} mem -p -t {params.bwa_threads} {params.bwa_ref} {input} > {output} 2> {log}"
        
rule merge_aln_and_unmapped_bams_after_mouse_removal: # bring back in all the metadata?
    input:
        unmapped_bam="tmp/{library}.{libcon_method}.sorted_and_merged.no_unpaired.unmapped.bam",
        mapped_bam="tmp/{library}.{libcon_method}.concatref.sorted.cleaned.realigned.bam"
    output:
        temp("tmp/{library}.{libcon_method}.concatref.sorted.cleaned.realigned.merged.bam")
    params:
        picard_jar = config["picard_jar"],
        mem = config["default_gatk_mem"],
        bwa_ref = config["reference_genome_for_bwa"],  
    log:
        "logs/remove_mouse_reads/merge_aln_and_unmapped_bams_after_mouse_removal/{library}.{libcon_method}.log"
    shell:
        "java -Xmx{params.mem} -XX:ParallelGCThreads=2 -Djava.io.tmpdir=${{TMPDIR}} -jar {params.picard_jar} "+\
        "MergeBamAlignment UNMAPPED={input.unmapped_bam} ALIGNED={input.mapped_bam} O={output} R={params.bwa_ref} SO=coordinate "+\
        "ALIGNER_PROPER_PAIR_FLAGS=true MAX_GAPS=-1 ORIENTATIONS=FR VALIDATION_STRINGENCY=SILENT 2> {log}"

# END REMOVE MOUSE #

# BEGIN MARK DUPLICATES (LIBCON METHOD SPECIFIC) #
rule mark_duplicates:
    input:
        "tmp/{library}.nebultra2.concatref.sorted.cleaned.realigned.merged.bam"
    output:
        marked_bam="results/mark_duplicates/{library}/{library}.nebultra2.bam",
        marked_bai="results/mark_duplicates/{library}/{library}.nebultra2.bai",
        metrics="results/mark_duplicates/{library}/{library}.nebultra2.metrics.txt",
    params:
        picard_jar=config["picard_jar"],
        #         markdup_threads=config["markdup_threads"]
    log:
        "logs/mark_duplicates/{library}.nebultra2.mark_duplicates.log"
    shell:
        "java -Xmx15G -XX:ParallelGCThreads=4 -Djava.io.tmpdir=${{TMPDIR}} -jar {params.picard_jar} "+\
        "MarkDuplicates I={input} O={output.marked_bam} M={output.metrics} CREATE_INDEX=TRUE 2> {log}"

rule mark_duplicates_idtxgen:
    # use UMIs for duplicate marking ("BARCODE_TAG=RX")
    input:
        "tmp/{library}.idtxgen.concatref.sorted.cleaned.realigned.merged.bam"
    output:
        marked_bam="results/mark_duplicates/{library}/{library}.idtxgen.bam",
        metrics="results/mark_duplicates/{library}/{library}.idtxgen.metrics.txt",
        bai=touch("results/mark_duplicates/{library}/{library}.idtxgen.bai"),
    params:
        picard_jar=config["picard_jar"],
        #         markdup_threads=config["markdup_threads"]
    log:
        "logs/mark_duplicates/{library}.idtxgen.mark_duplicates.log"
    shell:
        "java -Xmx15G -XX:ParallelGCThreads=4 -Djava.io.tmpdir=${{TMPDIR}} -jar {params.picard_jar} "+\
        "MarkDuplicates I={input} O={output.marked_bam} M={output.metrics} BARCODE_TAG=RX CREATE_INDEX=TRUE 2> {log}"
# END MARK DUPLICATES (LIBCON METHOD SPECIFIC) #
        
rule base_recalibrator:
    input:
        bam="results/mark_duplicates/{library}/{library}.{libcon_method}.bam",
        index="results/mark_duplicates/{library}/{library}.{libcon_method}.bai"
    output:
        "results/base_recalibration/{library}/{library}.{libcon_method}.recal_data.table"
    params:
        gatk=config["gatk"],
        mem=config["default_gatk_mem"],
        reference_genome=config["reference_genome_for_gatk"],
        known_sites_list=expand("--known-sites {sites}",sites=config["base_recal_known_sites"]),
        #         tmpdir=os.environ["TMPDIR"],
    log:
        "logs/base_recalibration/{library}.{libcon_method}.base_recalibration.log"
    shell:
        "{params.gatk} --java-options \"-Xmx{params.mem} -Djava.io.tmpdir=${{TMPDIR}}\" BaseRecalibrator "+\
        "-I {input.bam} -R {params.reference_genome} {params.known_sites_list} -O {output} 2> {log}"

rule apply_base_recal:
    input:
        bam="results/mark_duplicates/{library}/{library}.{libcon_method}.bam",
        recal_table="results/base_recalibration/{library}/{library}.{libcon_method}.recal_data.table",
    output:
        bam="results/base_recalibration/{library}/{library}.{libcon_method}.bam",
    params:
        gatk=config["gatk"],
        mem=config["default_gatk_mem"],
        reference_genome=config["reference_genome_for_gatk"],
        #         tmpdir=os.environ["TMPDIR"],
    log:
        "logs/base_recalibration/{library}.{libcon_method}.apply_base_recal.log"
    shell:
        "{params.gatk} --java-options \"-Xmx{params.mem} -Djava.io.tmpdir=${{TMPDIR}}\" ApplyBQSR --create-output-bam-index \"false\" "+\
        "-I {input.bam} -R {params.reference_genome} --bqsr-recal-file {input.recal_table} -O {output.bam} 2> {log}"
        
rule index_base_recal_bam:
    # index generation ApplyBQSR didnt seem to work? using gatk for better error reporting? altho the issue with samtools is that i was using wrong input bam!
    input:
        bam="results/base_recalibration/{library}/{library}.{libcon_method}.bam",
    output:
        bai="results/base_recalibration/{library}/{library}.{libcon_method}.bai",
    params:
        picard_jar=config["picard_jar"],
        mem = config["default_gatk_mem"],
    log:
        "logs/index_base_recal_bam/{library}.{libcon_method}.index.log"
    shell:
        "java -Xmx{params.mem} -Djava.io.tmpdir=${{TMPDIR}} -jar {params.picard_jar} "+\
        "BuildBamIndex I={input.bam} O={output.bai} 2> {log}"

# BEGIN METRICS #
rule collect_alnmt_metrics:
    input:
        bam=lambda wildcards:"results/"+config['metrics']['collect_alnmt_metrics'][wildcards.alnmt_step]['path_string']+"/{library}/{library}.{libcon_method}.bam",
        reference = lambda wildcards:config['metrics']['collect_alnmt_metrics'][wildcards.alnmt_step]['fasta'],
    output:
        "results/collect_alnmt_metrics/{alnmt_step}/{library}/{library}.{libcon_method}.metrics.txt"
    params:
        picard_jar=config["picard_jar"],
    log:
        "logs/collect_alnmt_metrics/{alnmt_step}/{library}.{libcon_method}.log"
    shell:
        "java -Xmx4G -XX:ParallelGCThreads=1 -Djava.io.tmpdir=${{TMPDIR}} -jar {params.picard_jar} "+
        "CollectAlignmentSummaryMetrics I={input.bam} O={output} R={input.reference} 2> {log}" 
        
rule collect_wgs_metrics:
    input: # reusing config from "collect_alnmt_metrics" step
        bam=lambda wildcards:"results/"+config['metrics']['collect_alnmt_metrics'][wildcards.alnmt_step]['path_string']+"/{library}/{library}.{libcon_method}.bam",
        reference = lambda wildcards:config['metrics']['collect_alnmt_metrics'][wildcards.alnmt_step]['fasta'],
    output:
        "results/collect_wgs_metrics/{alnmt_step}/{library}/{library}.{libcon_method}.metrics.txt"
    params:
        picard_jar=config["picard_jar"],
        mem=config["default_gatk_mem"],
    log:
        "logs/collect_wgs_metrics/{alnmt_step}/{library}.{libcon_method}.log"
    shell:
        "java -Xmx{params.mem} -XX:ParallelGCThreads=1 -Djava.io.tmpdir=${{TMPDIR}} -jar {params.picard_jar} "+
        "CollectWgsMetrics I={input.bam} O={output} R={input.reference} 2> {log}" 

rule collect_isize_metrics:
    input:
        bam="results/base_recalibration/{library}/{library}.{libcon_method}.bam",
    output:
        metrics="results/collect_isize_metrics/{library}/{library}.{libcon_method}.metrics.txt",
        histogram="results/collect_isize_metrics/{library}/{library}.{libcon_method}.histogram.pdf",
    params:
        picard_jar=config["picard_jar"],
        mem = config["default_gatk_mem"],
        reference_genome=config["reference_genome_for_gatk"]
    log:
        "logs/collect_isize_metrics/{library}.{libcon_method}.collect_isize_metrics.log"
    shell:
        "java -Xmx{params.mem} -Djava.io.tmpdir=${{TMPDIR}} -jar {params.picard_jar} "+\
        "CollectInsertSizeMetrics I={input.bam} O={output.metrics} H={output.histogram} 2> {log}"

rule collect_hs_metrics:
    input:
        bam="results/base_recalibration/{library}/{library}.{libcon_method}.bam",
        regions=lambda wildcards:config['metrics']['collect_hs_metrics'][wildcards.target]['regions'],
        probes=lambda wildcards:config['metrics']['collect_hs_metrics'][wildcards.target]['probes'],
    output:
        "results/collect_hs_metrics/{target}/{library}/{library}.{libcon_method}.metrics.txt",
    params:
        picard_jar=config["picard_jar"],
        mem=config["default_gatk_mem"],
        reference_genome=config["reference_genome_for_gatk"],
        target=lambda wildcards:config['metrics']['collect_hs_metrics'][wildcards.target]['target'],
    log:
        "logs/collect_hs_metrics/{library}.{libcon_method}.{target}.collect_hs_metrics.log"
    shell:
        "java -Xmx{params.mem} -Djava.io.tmpdir=${{TMPDIR}} -jar {params.picard_jar} "+\
        "CollectHsMetrics I={input.bam} R={params.reference_genome} BAIT_INTERVALS={input.probes} TARGET_INTERVALS={input.regions} "+\
        "BAIT_SET_NAME=\"{params.target}\" O={output} 2> {log}"
        
# END METRICS #
        
###########
#         #
# MUTECT2 #
#         #
###########

#############################
# GENERATE PANEL OF NORMALS #
#############################
rule mutect2_for_panel_of_normals:
    input:
        bam = lambda wildcards:config["benign_lib_for_mutect_pon"][wildcards.library][wildcards.libcon_method]+".bam",
        bai = lambda wildcards:config["benign_lib_for_mutect_pon"][wildcards.library][wildcards.libcon_method]+".bai",
    output:
        "results/pon_for_mutect2/{library}/{library}.{libcon_method}.vcf.gz"
    params:
        gatk=config["gatk"],
        mem=config["default_gatk_mem"],
        reference_genome=config["reference_genome_for_gatk"]
    log:
        "logs/pon_for_mutect2/{library}.{libcon_method}.pon_for_mutect2.log"
    shell:
        "{params.gatk} --java-options \"-Xmx{params.mem} -Djava.io.tmpdir=${{TMPDIR}}\" Mutect2 --max-mnp-distance 0 "+\
        "-I {input.bam} -R {params.reference_genome} -O {output} 2> {log}"

rule pon_genomics_db_import:
    input:
        vcfs = [ "results/pon_for_mutect2/{l}/{l}.{m}.vcf.gz".format(l=l,m=m) for l in config["benign_lib_for_mutect_pon"] for m in config["benign_lib_for_mutect_pon"][l] ]
    output:
        touch(temp("tmp/pon_genomics_db_import.finished.txt"))
    params:
        gatk=config["gatk"],
        mem=config["default_gatk_mem"],
        reference_genome=config["reference_genome_for_gatk"],
        chrs=config["pon_genomics_db_import_chrs"]
    log:
        "logs/m2_pon_genomics_db_import/pon_genomics_db_import.log"
    run:
        vcf_string = expand("-V {v}",v=input.vcfs)
        cmd = "{params.gatk} --java-options \"-Xmx{params.mem} -Djava.io.tmpdir=${{TMPDIR}}\" GenomicsDBImport "
        cmd += "{vcf_string} -R {params.reference_genome} {params.chrs} --genomicsdb-workspace-path results/m2_pon_db/ 2> {log}"
        shell(cmd)
        
rule create_somatic_pons:
    input:
        "tmp/pon_genomics_db_import.finished.txt" # dummy input
    output:
        "results/pon_for_mutect2/panel_of_normals_vcf/panel_of_normals_vcf.vcf.gz" 
    params:
        gatk=config["gatk"],
        mem=config["default_gatk_mem"],
        reference_genome=config["reference_genome_for_gatk"]
    log:
        "logs/create_somatic_pons/create_somatic_pons.log"
    shell:
        "{params.gatk} --java-options \"-Xmx{params.mem} -Djava.io.tmpdir=${{TMPDIR}}\" CreateSomaticPanelOfNormals "+\
        "-R {params.reference_genome} -V gendb://results/m2_pon_db/ -O {output} 2> {log}"

######################
# NO/IGNORE GERMLINE #
######################
rule mutect2_somatic_only:
    input:
        pon_vcf = "results/pon_for_mutect2/panel_of_normals_vcf/panel_of_normals_vcf.vcf.gz",
        bam = "results/base_recalibration/{library}/{library}.{libcon_method}.bam",
        bai = "results/base_recalibration/{library}/{library}.{libcon_method}.bai",
    output:
        "results/mutect2/{library}/somatic_only/{library}.{libcon_method}.vcf.gz"
    params:
        gatk=config["gatk"],
        mem=config["default_gatk_mem"],
        reference_genome=config["reference_genome_for_gatk"],
        premade_germline_resource=config["mutect2_germline_resource"]
    log:
        "logs/mutect2/{library}.{libcon_method}.mutect2_somatic_only.log"
    shell:
        "{params.gatk} --java-options \"-Xmx{params.mem} -Djava.io.tmpdir=${{TMPDIR}}\" Mutect2 "+\
        "--germline-resource {params.premade_germline_resource} --panel-of-normals {input.pon_vcf} -I {input.bam} -R {params.reference_genome} -O {output} 2> {log}"
        
################
# USE GERMLINE #
################
rule mutect2_germline_matched:
    input:
        pon_vcf = "results/pon_for_mutect2/panel_of_normals_vcf/panel_of_normals_vcf.vcf.gz",
        somatic_bam =  lambda wildcards: config["germline_counterpart"][wildcards.library][wildcards.libcon_method]['somatic']+'.bam',
        somatic_bai =  lambda wildcards: config["germline_counterpart"][wildcards.library][wildcards.libcon_method]['somatic']+'.bai',
        germline_bam = lambda wildcards: config["germline_counterpart"][wildcards.library][wildcards.libcon_method]['germline']+'.bam',
        germline_bai = lambda wildcards: config["germline_counterpart"][wildcards.library][wildcards.libcon_method]['germline']+'.bam',
    output:
        vcf="results/mutect2/{library}/germline_matched/{library}.{libcon_method}.vcf.gz",
    params:
        gatk=config["gatk"],
        mem=config["default_gatk_mem"],
        reference_genome=config["reference_genome_for_gatk"],
        premade_germline_resource=config["mutect2_germline_resource"],
        germline_lib = lambda wildcards:config["germline_counterpart"][wildcards.library][wildcards.libcon_method]['germline_lib'],
    log:
        "logs/mutect2/{library}.{libcon_method}.mutect2_germline_matched.log"
    shell:
        "{params.gatk} --java-options \"-Xmx{params.mem} -Djava.io.tmpdir=${{TMPDIR}}\" Mutect2 "+\
        "--germline-resource {params.premade_germline_resource} --panel-of-normals {input.pon_vcf} -I {input.somatic_bam} -I {input.germline_bam} -normal {params.germline_lib} -R {params.reference_genome} -O {output.vcf} 2> {log}"
        
###########
# UNIFIED #
###########
rule filt_mutect2:
    input:
        "results/mutect2/{library}/{mutect2_run_type}/{library}.{libcon_method}.vcf.gz"
    output:
        "results/filt_mutect2/{library}/{mutect2_run_type}/{library}.{libcon_method}.vcf.gz"
    params:
        gatk=config["gatk"],
        mem=config["default_gatk_mem"],
        reference_genome=config["reference_genome_for_gatk"],
    log:
        "logs/filt_mutect2/{library}.{libcon_method}.{mutect2_run_type}.filt_mutect2.log"
    shell:
        "{params.gatk} --java-options \"-Xmx{params.mem} -Djava.io.tmpdir=${{TMPDIR}}\" FilterMutectCalls "+\
        "-R {params.reference_genome} -V {input} -O {output} 2> {log}"

rule funcotator:
    input:
        "results/filt_mutect2/{library}/{mutect2_run_type}/{library}.{libcon_method}.vcf.gz"
    output:
        "results/funcotator/{library}/{mutect2_run_type}/{library}.{libcon_method}.vcf.gz"
    params:
        gatk=config["gatk"],
        reference_genome=config["reference_genome_for_gatk"],
        ref_version=config["funcotator_ref_version"],
        data_sources=config["funcotator_data_sources"]
    log:
        "logs/funcotator/{library}.{libcon_method}.{mutect2_run_type}.funcotator.log"
    shell:
        "{params.gatk} --java-options \"-Xmx15G -XX:ParallelGCThreads=4 -Djava.io.tmpdir=${{TMPDIR}}\" Funcotator "+\
        "-R {params.reference_genome} --ref-version {params.ref_version} --output-file-format VCF --data-sources-path {params.data_sources} "+\
        "-V {input} -O {output} 2> {log}"

rule parse_funcotator_vcf:
    input:
        "results/funcotator/{library}/{mutect2_run_type}/{library}.{libcon_method}.vcf.gz"
    output:
        "results/parse_funcotator_vcf/{library}/{mutect2_run_type}/{library}.{libcon_method}.csv"
    params:
        python=config["python"],
        script=config["parse_funcotator_vcf_script"]
    log:
        "logs/parse_funcotator_vcf/{library}.{libcon_method}.{mutect2_run_type}.parse_funcotator_vcf.log"
    shell:
        "{params.python} {params.script} {input} {output} 2> {log}"

###########
#         #
# STRELKA #
#         #
###########
rule manta:
    input:
        somatic_bam =  lambda wildcards: config["germline_counterpart"][wildcards.library][wildcards.libcon_method]['somatic']+'.bam',
        somatic_bai =  lambda wildcards: config["germline_counterpart"][wildcards.library][wildcards.libcon_method]['somatic']+'.bai',
        germline_bam = lambda wildcards: config["germline_counterpart"][wildcards.library][wildcards.libcon_method]['germline']+'.bam',
        germline_bai = lambda wildcards: config["germline_counterpart"][wildcards.library][wildcards.libcon_method]['germline']+'.bai',
    output:
        manta_cmd_file = "results/manta/{library}/germline_matched/{library}.{libcon_method}/runWorkflow.py",
        small_indel_vcf="results/manta/{library}/germline_matched/{library}.{libcon_method}/results/variants/candidateSmallIndels.vcf.gz",
    params:
        manta_install_path=config["manta_install_path"],
        mem=config["default_gatk_mem"],
        reference_genome=config["reference_genome_for_gatk"],
        threads=config['bwa_threads'],
        strelka_regions=config['strelka_regions'],
    log:
        "logs/manta/{library}.{libcon_method}.manta.log"
    shell:
        "{params.manta_install_path}/configManta.py --exome --callRegions {params.strelka_regions} --normalBam {input.germline_bam} --tumorBam {input.somatic_bam} "+\
        "--referenceFasta {params.reference_genome} --runDir results/manta/{wildcards.library}/germline_matched/{wildcards.library}.{wildcards.libcon_method}/ 2> {log} && "+\
        "{output.manta_cmd_file} -j {params.threads} 2>> {log}"

rule strelka2:
    input:
        somatic_bam =  lambda wildcards: config["germline_counterpart"][wildcards.library][wildcards.libcon_method]['somatic']+'.bam',
        somatic_bai =  lambda wildcards: config["germline_counterpart"][wildcards.library][wildcards.libcon_method]['somatic']+'.bai',
        germline_bam = lambda wildcards: config["germline_counterpart"][wildcards.library][wildcards.libcon_method]['germline']+'.bam',
        germline_bai = lambda wildcards: config["germline_counterpart"][wildcards.library][wildcards.libcon_method]['germline']+'.bam',
        manta_small_indel_vcf = lambda wildcards: config["germline_counterpart"][wildcards.library][wildcards.libcon_method]['manta_small_indel_vcf'],
    output:
        strelka_cmd_file = "results/strelka2/{library}/germline_matched/{library}.{libcon_method}/runWorkflow.py",
        snv_vcf="results/strelka2/{library}/germline_matched/{library}.{libcon_method}/results/variants/somatic.snvs.vcf.gz",
        indel_vcf="results/strelka2/{library}/germline_matched/{library}.{libcon_method}/results/variants/somatic.indels.vcf.gz",
    params:
        strelka_install_path=config["strelka_install_path"],
        mem=config["default_gatk_mem"],
        reference_genome=config["reference_genome_for_gatk"],
        threads=config['bwa_threads'],
        strelka_regions=config['strelka_regions'],
    log:
        "logs/strelka2/{library}.{libcon_method}.strelka2.log"
    shell:
        "{params.strelka_install_path}/configureStrelkaSomaticWorkflow.py --exome --callRegions {params.strelka_regions} --normalBam {input.germline_bam} --tumorBam {input.somatic_bam} --indelCandidates {input.manta_small_indel_vcf} "+\
        "--referenceFasta {params.reference_genome} --runDir results/strelka2/{wildcards.library}/germline_matched/{wildcards.library}.{wildcards.libcon_method}/ 2> {log} && "+\
        "{output.strelka_cmd_file} -m local -j {params.threads} 2>> {log}"
        
rule funcotator_on_strelka2:
    input:
        "results/strelka2/{library}/germline_matched/{library}.{libcon_method}/results/variants/somatic.{var_type}.vcf.gz",
    output:
        "results/funcotator_on_strelka/{library}/germline_matched/{library}.{libcon_method}.{var_type}.vcf.gz",
    params:
        gatk=config["gatk"],
        reference_genome=config["reference_genome_for_gatk"],
        ref_version=config["funcotator_ref_version"],
        data_sources=config["funcotator_data_sources"],
    log:
        "logs/funcotator_on_strelka/{library}.{libcon_method}.germline_matched.{var_type}.funcotator.log"
    shell:
        "{params.gatk} --java-options \"-Xmx15G -XX:ParallelGCThreads=4 -Djava.io.tmpdir=${{TMPDIR}}\" Funcotator "+\
        "-R {params.reference_genome} --ref-version {params.ref_version} --output-file-format VCF --data-sources-path {params.data_sources} "+\
        "-V {input} -O {output} 2> {log}"

rule parse_funcotator_on_strelka_vcf:
    input:
        "results/funcotator_on_strelka/{library}/germline_matched/{library}.{libcon_method}.{var_type}.vcf.gz",
    output:
        "results/parse_funcotator_on_strelka_vcf/{library}/germline_matched/{library}.{libcon_method}.{var_type}.csv"
    params:
        python=config["python"],
        script=config["parse_funcotator_on_strelka_vcf_script"]
    log:
        "logs/parse_funcotator_vcf/{library}.{libcon_method}.{var_type}.parse_funcotator_vcf.log"
    shell:
        "{params.python} {params.script} {input} {output} 2> {log}"
        
##########
#        #
# CNVKIT #
#        #
##########
rule cnvkit_batch_benign_normal: # needs R/3.6.2 to produce the .cns files???
    input:
        tumor_bams = config["cnvkit_bams"]["tumor"]['input'],
        normal_bams = config["cnvkit_bams"]["benign"]['input'],
    output:
#         finished_txt = touch("results/cnvkit/tumor_v_benign/finished.txt"),
        touch(expand("results/cnvkit/batch/tumor_v_benign/{lib_and_method}.{cnx}",lib_and_method=config["cnvkit_bams"]["tumor"]['lib_and_method'],cnx=['cns','cnr'])),
    params:
        cnvkit_install_path = config['cnvkit_install_path'],
        mem = config['cnvkit_default_mem'],
        reference_genome = config["reference_genome_for_gatk"],
        cnvkit_tfbs_g_probes = config['cnvkit_tfbs_g_probes'],
        cnvkit_refflat = config['cnvkit_refflat'],
        cnvkit_access_excludeTSS = config['cnvkit_access_excludeTSS'],
    log:
        "logs/cnvkit/batch/tumor_v_benign/cnvkit.tumor_v_benign.log"
    shell:
        "{params.cnvkit_install_path} batch {input.tumor_bams} --normal {input.normal_bams} --targets {params.cnvkit_tfbs_g_probes} "+\
        "--fasta {params.reference_genome} --access {params.cnvkit_access_excludeTSS} --annotate {params.cnvkit_refflat} -p 0 "+\
        "--output-reference results/cnvkit/batch/tumor_v_benign/reference.cnn --output-dir results/cnvkit/batch/tumor_v_benign/ --scatter --diagram "+\
        "2> {log}"
        
rule cnvkit_batch_flat_normal:
    input:
        tumor_bams = [ x for k in config["cnvkit_bams"] for x in config["cnvkit_bams"][k]['input'] ], # all files
    output:
#         finished_txt = touch("results/cnvkit/flat_normal/finished.txt"),
        touch(expand("results/cnvkit/batch/flat_normal/{lib_and_method}.{cnx}",lib_and_method=[ x for k in config["cnvkit_bams"] for x in config["cnvkit_bams"][k]['lib_and_method'] ],cnx=['cns','cnr'])),
    params:
        cnvkit_install_path = config['cnvkit_install_path'],
        mem = config['cnvkit_default_mem'],
        reference_genome = config["reference_genome_for_gatk"],
        cnvkit_tfbs_g_probes = config['cnvkit_tfbs_g_probes'],
        cnvkit_refflat = config['cnvkit_refflat'],
        cnvkit_access_excludeTSS = config['cnvkit_access_excludeTSS'],
    log:
        "logs/cnvkit/batch/flat_normal/cnvkit.flat_normal.log"
    shell:
        "{params.cnvkit_install_path} batch {input.tumor_bams} --normal --targets {params.cnvkit_tfbs_g_probes} "+\
        "--fasta {params.reference_genome} --access {params.cnvkit_access_excludeTSS} --annotate {params.cnvkit_refflat} -p 0 "+\
        "--output-reference results/cnvkit/batch/flat_normal/reference.cnn --output-dir results/cnvkit/batch/flat_normal/ --scatter --diagram "+\
        "2> {log}"

rule cnvkit_batch_pbmc_normal:
    input:
        tumor_bams = config["cnvkit_bams"]["tumor"]['input'],
        normal_bams = config["cnvkit_bams"]["pbmc"]['input'],
    output:
        #         finished_txt = touch("results/cnvkit/tumor_v_pbmc/finished.txt"),
        touch(expand("results/cnvkit/batch/tumor_v_pbmc/{lib_and_method}.{cnx}",lib_and_method=config["cnvkit_bams"]["tumor"]['lib_and_method'],cnx=['cns','cnr'])),
    params:
        cnvkit_install_path = config['cnvkit_install_path'],
        mem = config['cnvkit_default_mem'],
        reference_genome = config["reference_genome_for_gatk"],
        cnvkit_tfbs_g_probes = config['cnvkit_tfbs_g_probes'],
        cnvkit_refflat = config['cnvkit_refflat'],
        cnvkit_access_excludeTSS = config['cnvkit_access_excludeTSS'],
    log:
        "logs/cnvkit/batch/tumor_v_pbmc/cnvkit.tumor_v_pbmc.log"
    shell:
        "{params.cnvkit_install_path} batch {input.tumor_bams} --normal {input.normal_bams} --targets {params.cnvkit_tfbs_g_probes} "+\
        "--fasta {params.reference_genome} --access {params.cnvkit_access_excludeTSS} --annotate {params.cnvkit_refflat} -p 0 "+\
        "--output-reference results/cnvkit/batch/tumor_v_pbmc/reference.cnn --output-dir results/cnvkit/batch/tumor_v_pbmc/ --scatter --diagram "+\
        "2> {log}"

rule cnvkit_metrics:
    input:
        cns = "results/cnvkit/batch/{comparison}/{library}.{libcon_method}.cns",
        cnr = "results/cnvkit/batch/{comparison}/{library}.{libcon_method}.cnr",
    output:
        "results/cnvkit/metrics/{comparison}/{library}.{libcon_method}.metrics.txt",
    params:
        cnvkit_install_path = config['cnvkit_install_path'],
    log:
        "logs/cnvkit/metrics/{comparison}/{library}.{libcon_method}.log"
    shell:
        "{params.cnvkit_install_path} metrics {input.cnr} -s {input.cns} > {output} 2> {log}"

rule cnvkit_adj_call_w_tumor_frxn:
    input:
        cns = "results/cnvkit/batch/{comparison}/{library}.{libcon_method}.cns",
    output:
        "results/cnvkit/adj_w_tumor_frxn/{comparison}/{library}.{libcon_method}.cns",
    params:
        cnvkit_install_path = config['cnvkit_install_path'],
        purity = lambda wildcards:config['tumor_frxn'][wildcards.library],
    log:
        "logs/cnvkit/adj_w_tumor_frxn/{comparison}/{library}.{libcon_method}.log"
    shell:
        "{params.cnvkit_install_path} call {input.cns} --purity {params.purity} -o {output} 2> {log}"

# rule cnvkit_heatmap:
#     input:
#         expand("results/cnvkit/{call_method}/{comparison}/{library}.{libcon_method}.cns",
#     output:
#         "results/cnvkit/{call_method}/{comparison}/{library}.{libcon_method}.cns",
#     params:
#         cnvkit_install_path = config['cnvkit_install_path'],
#         purity = lambda wildcards:config['tumor_frxn'][wildcards.library],
#     log:
#         "logs/cnvkit/adjust_call_w_tumor_frxn/{comparison}/{library}.{libcon_method}.log"
#     shell:
#         "{params.cnvkit_install_path} call {input.cns} --purity {params.purity} -o {output} 2> {log}"

# rule cnvkit_normal_only:
#     input:
#         normal_bams = config["cnvkit_bams"]["normal"],
#     output:
#         #         output_reference = "results/cnvkit/my_reference.cnn",
#         finished_txt = touch("results/cnvkit/normal_only/cnvkit.finished.txt"),
#     params:
#         cnvkit_install_path = config['cnvkit_install_path'],
#         mem = config['cnvkit_default_mem'],
#         reference_genome = config["reference_genome_for_gatk"],
#         cnvkit_tfbs_g_probes = config['cnvkit_tfbs_g_probes'],
#         cnvkit_refflat = config['cnvkit_refflat'],
#         cnvkit_access_excludeTSS = config['cnvkit_access_excludeTSS'],
#     log:
#         "logs/cnvkit/cnvkit.normal_only.log"
#     shell:
#         "{params.cnvkit_install_path} batch {input.tumor_bams} --normal {input.normal_bams} --targets {params.cnvkit_tfbs_g_probes} "+\
#         "--fasta {params.reference_genome} --access {params.cnvkit_access_excludeTSS} --annotate {params.cnvkit_refflat} -p 0 "+\
#         "--output-reference {output.output_reference} --output-dir results/cnvkit/normal_only/ --diagram --scatter "+\
#         "2> {log}"
    
# OLD ALL RULE TARGETS
# rule all:
#     input:
# #     # debug
# #         [ "results/mark_duplicates/{l}/{l}.{m}.bam".format(l=l,m=m) for m in config["bams_to_merge"] for l in config["bams_to_merge"][m] ],
#     # QC
#         [ "results/{metric}/{t}/{l}/{l}.{m}.metrics.txt".format(metric=metric,t=t,l=l,m=m) for metric in config['metrics'] for t in config['metrics'][metric] for m in config["bams_to_merge"] for l in config["bams_to_merge"][m] ],
#         [ "results/collect_isize_metrics/{l}/{l}.{m}.metrics.txt".format(l=l,m=m) for m in config["bams_to_merge"] for l in config["bams_to_merge"][m] ],
#     # parse vcf's
#         [ "results/parse_funcotator_vcf/{l}/somatic_only/{l}.{m}.csv".format(l=l,m=m) for m in config["bams_to_merge"] for l in config["bams_to_merge"][m] ],
#         [ "results/parse_funcotator_vcf/{l}/germline_matched/{l}.{m}.csv".format(l=l,m=m) for l in config["germline_counterpart"] for m in config["germline_counterpart"][l] ],
#         [ "results/parse_funcotator_on_strelka_vcf/{l}/germline_matched/{l}.{m}.{v}.csv".format(l=l,m=m,v=v) for l in config["germline_counterpart"] for m in config["germline_counterpart"][l] for v in ['snvs','indels'] ],
#     # cnvkit
#         expand("results/cnvkit/batch/{comparison}/{output}.cns",comparison=['tumor_v_benign','tumor_v_pbmc','flat_normal'],output=config["cnvkit_bams"]["tumor"]['lib_and_method']),
#         expand("results/cnvkit/metrics/{comparison}/{lib_and_method}.metrics.txt",comparison=['tumor_v_benign','tumor_v_pbmc'],lib_and_method=config["cnvkit_bams"]["tumor"]['lib_and_method']),
#         expand("results/cnvkit/metrics/{comparison}/{lib_and_method}.metrics.txt",comparison=['flat_normal'],lib_and_method=[ x for k in config["cnvkit_bams"] for x in config["cnvkit_bams"][k]['lib_and_method'] ]),
#         expand("results/cnvkit/adj_w_tumor_frxn/{comparison}/{lib_and_method}.cns",comparison=['tumor_v_benign','tumor_v_pbmc'],lib_and_method=config["cnvkit_bams"]["tumor"]['lib_and_method']),

        #         expand("results/base_recalibration/{library}/{library}.bam",library=config["bams_to_merge"]),
        #         expand("results/base_recalibration/{library}/{library}.bam.bai",library=config["bams_to_merge"]),
        #         expand("results/collect_hs_metrics/coding/{library}/{library}.metrics.txt",library=config["bams_to_merge"]),
        #         expand("results/collect_hs_metrics/all/{library}/{library}.metrics.txt",library=config["bams_to_merge"]),
        #         #expand("results/mutect2/{library}/{library}.vcf.gz",library=config["bams_to_merge"]), # no longer needed if funcotator is activated?
        #         #expand("results/funcotator/{library}/{library}.vcf.gz",library=config["bams_to_merge"]),
        #         expand("results/parse_funcotator_vcf/{library}/{library}.csv",library=config["bams_to_merge"])
        #         expand("tmp/idtxgen/{idtxgen_fastq}.R{num}.fq.gz",idtxgen_fastq=config["idtxgen_fastq"],num=['1','2'])
        #         expand("tmp/idtxgen/{idtxgen_fastq}.aln.bam",idtxgen_fastq=config["idtxgen_fastq"]),
        #         expand("results/idtxgen/mark_duplicates/{library}/{library}.bam",library=config["bams_to_merge"]['idtxgen']),
        #         expand("results/base_recalibration/{library}/{library}.idtxgen.recal_data.table",library=config["bams_to_merge"]['idtxgen'])+\
        #         expand("results/base_recalibration/{library}/{library}.nebultra2.recal_data.table",library=config["bams_to_merge"]['nebultra2'])
        #         [ "results/funcotator/{l}/somatic_only/{l}.{m}.vcf.gz".format(l=l,m=m) for m in config["bams_to_merge"] for l in config["bams_to_merge"][m] ],
#         [ "results/funcotator/{l}/germline_matched/{l}.{m}.vcf.gz".format(l=l,m=m) for l in config["germline_counterpart"] for m in config["germline_counterpart"][l] ],
#         [ "results/strelka2/{l}/germline_matched/{l}.{m}/results/variants/somatic.snvs.vcf.gz".format(l=l,m=m) for l in config["germline_counterpart"] for m in config["germline_counterpart"][l] ],
#         [ "results/funcotator_on_strelka/{l}/germline_matched/{l}.{m}.{v}.vcf.gz".format(l=l,m=m,v=v) for l in config["germline_counterpart"] for m in config["germline_counterpart"][l] for v in ['snvs','indels'] ],
#         expand("results/cnvkit/{c}/finished.txt",c=['tumor_v_benign','flat_normal','tumor_v_pbmc'])
#         [ "results/collect_hs_metrics/{t}/{l}/{l}.{m}.metrics.txt".format(t=t,l=l,m=m) for t in ['tfbs_g','tss','coding','all'] for m in config["bams_to_merge"] for l in config["bams_to_merge"][m] ],
#         [ "results/collect_wgs_metrics/{t}/{l}/{l}.{m}.metrics.txt" .format(t=t,l=l,m=m) for t in config['metrics']['collect_wgs_metrics'] for m in config["bams_to_merge"] for l in config["bams_to_merge"][m] ],

# OLD RULES #

# rule index_sorted_bams:
#     input:
#         "tmp/{library}_{origlib}_{flowcell}_{lane}.nebultra2.sorted.bam"
#     output:
#         temp("tmp/{library}_{origlib}_{flowcell}_{lane}.nebultra2.sorted.bai")
#     params:
#         samtools=config["samtools"]
#     log:
#         "logs/index/{library}_{origlib}_{flowcell}_{lane}.nebultra2.sorted.index.log"
#     shell:
#         "{params.samtools} index {input} {output} 2> {log}"

# rule collect_hs_metrics__all:
#     input:
#         bam="results/base_recalibration/{library}/{library}.{libcon_method}.bam",
#         #         bai="results/base_recalibration/{library}/{library}.{libcon_method}.bai",
#     output:
#         "results/collect_hs_metrics/all/{library}/{library}.{libcon_method}.metrics.txt"
#     params:
#         panel_type="all",
#         picard_jar=config["picard_jar"],
#         mem = config["default_gatk_mem"],
#         regions=config["capture_panel_regions__all"],
#         probes=config["capture_panel_probes__all"],
#         reference_genome=config["reference_genome_for_gatk"]
#     log:
#         "logs/collect_hs_metrics/{library}.{libcon_method}.all.collect_hs_metrics.log"
#     shell:
#         "java -Xmx{params.mem} -Djava.io.tmpdir=${{TMPDIR}} -jar {params.picard_jar} "+\
#         "CollectHsMetrics I={input.bam} R={params.reference_genome} BAIT_INTERVALS={params.probes} TARGET_INTERVALS={params.regions} "+\
#         "BAIT_SET_NAME=\"{params.panel_type}\" O={output} 2> {log}"

# rule collect_hs_metrics__tfbs_g:
#     input:
#         bam="results/base_recalibration/{library}/{library}.{libcon_method}.bam",
#         #         bai="results/base_recalibration/{library}/{library}.{libcon_method}.bai",
#     output:
#         "results/collect_hs_metrics/tfbs_g/{library}/{library}.{libcon_method}.metrics.txt"
#     params:
#         panel_type="tfbs_g",
#         picard_jar=config["picard_jar"],
#         mem=config["default_gatk_mem"],
#         regions=config["capture_panel_regions__tfbs_g"],
#         probes=config["capture_panel_probes__tfbs_g"],
#         reference_genome=config["reference_genome_for_gatk"],
#     log:
#         "logs/collect_hs_metrics/{library}.{libcon_method}.tfbs_g.collect_hs_metrics.log"
#     shell:
#         "java -Xmx{params.mem} -Djava.io.tmpdir=${{TMPDIR}} -jar {params.picard_jar} "+\
#         "CollectHsMetrics I={input.bam} R={params.reference_genome} BAIT_INTERVALS={params.probes} TARGET_INTERVALS={params.regions} "+\
#         "BAIT_SET_NAME=\"{params.panel_type}\" O={output} 2> {log}"
        
# rule collect_hs_metrics__coding:
#     input:
#         bam="results/base_recalibration/{library}/{library}.{libcon_method}.bam",
#         #         bai="results/base_recalibration/{library}/{library}.{libcon_method}.bai",
#     output:
#         "results/collect_hs_metrics/coding/{library}/{library}.{libcon_method}.metrics.txt"
#     params:
#         panel_type="coding",
#         picard_jar=config["picard_jar"],
#         mem=config["default_gatk_mem"],
#         regions=config["capture_panel_regions__coding_only"],
#         probes=config["capture_panel_regions__coding_only"], # havent filtered probes down for this...should at some point?
#         reference_genome=config["reference_genome_for_gatk"]
#     log:
#         "logs/collect_hs_metrics/{library}.{libcon_method}.coding.collect_hs_metrics.log"
#     shell:
#         "java -Xmx{params.mem} -Djava.io.tmpdir=${{TMPDIR}} -jar {params.picard_jar} "+\
#         "CollectHsMetrics I={input.bam} R={params.reference_genome} BAIT_INTERVALS={params.probes} TARGET_INTERVALS={params.regions} "+\
#         "BAIT_SET_NAME=\"{params.panel_type}\" O={output} 2> {log}"
        
# rule collect_hs_metrics__tss:
#     input:
#         bam="results/base_recalibration/{library}/{library}.{libcon_method}.bam",
#         #         bai="results/base_recalibration/{library}/{library}.{libcon_method}.bai",
#     output:
#         "results/collect_hs_metrics/tss/{library}/{library}.{libcon_method}.metrics.txt"
#     params:
#         panel_type="tss",
#         picard_jar=config["picard_jar"],
#         mem=config["default_gatk_mem"],
#         regions=config["capture_panel_regions__tss"],
#         probes=config["capture_panel_probes__tss"],
#         reference_genome=config["reference_genome_for_gatk"],
#     log:
#         "logs/collect_hs_metrics/{library}.{libcon_method}.tss.collect_hs_metrics.log"
#     shell:
#         "java -Xmx{params.mem} -Djava.io.tmpdir=${{TMPDIR}} -jar {params.picard_jar} "+\
#         "CollectHsMetrics I={input.bam} R={params.reference_genome} BAIT_INTERVALS={params.probes} TARGET_INTERVALS={params.regions} "+\
#         "BAIT_SET_NAME=\"{params.panel_type}\" O={output} 2> {log}"

# rule collect_wgs_metrics__concatRef:
#     input:
#         bam = "tmp/{library}.{libcon_method}.concatref.sorted.bam",
#         bai = "tmp/{library}.{libcon_method}.concatref.sorted.bai",
#     output:
#         "results/collect_wgs_metrics/concatRef/{library}/{library}.{libcon_method}.metrics.txt"
#     params:
#         picard_jar=config["picard_jar"],
#         reference=config['reference_genome_for_concatref'], # made sequence dict so that ref can be used w picard
#     log:
#         "logs/collect_wgs_metrics/concatRef/{library}.{libcon_method}.log"
#     shell:
#         "java -Xmx4G -XX:ParallelGCThreads=1 -Djava.io.tmpdir=${{TMPDIR}} -jar {params.picard_jar} "+
#         "CollectWgsMetrics I={input.bam} O={output} R={params.reference} 2> {log}" 
        
# rule collect_wgs_metrics__human:
#     input:
#         bam="results/base_recalibration/{library}/{library}.{libcon_method}.bam",
#     output:
#         "results/collect_wgs_metrics/human/{library}/{library}.{libcon_method}.metrics.txt"
#     params:
#         picard_jar=config["picard_jar"],
#         reference=config['reference_genome_for_gatk'], # made sequence dict so that ref can be used w picard
#     log:
#         "logs/collect_wgs_metrics/human/{library}.{libcon_method}.log"
#     shell:
#         "java -Xmx4G -XX:ParallelGCThreads=1 -Djava.io.tmpdir=${{TMPDIR}} -jar {params.picard_jar} "+
#         "CollectWgsMetrics I={input.bam} O={output} R={params.reference} 2> {log}" 

# rule collect_wgs_metrics:
#     input:
#         bam="results/base_recalibration/{library}/{library}.{libcon_method}.bam",
#         reference = lambda wildcards:config['metrics']['collect_wgs_metrics'][wildcards.target],
#     output:
#         "results/collect_wgs_metrics/{target}/{library}/{library}.{libcon_method}.metrics.txt"
#     params:
#         picard_jar=config["picard_jar"],
#     log:
#         "logs/collect_wgs_metrics/{target}/{library}.{libcon_method}.log"
#     shell:
#         "java -Xmx4G -XX:ParallelGCThreads=1 -Djava.io.tmpdir=${{TMPDIR}} -jar {params.picard_jar} "+
#         "CollectWgsMetrics I={input.bam} O={output} R={input.reference} 2> {log}"
# rule collect_wgs_metrics:
#     input:
#         bam=lambda wildcards:"results/"+config['metrics']['collect_wgs_metrics'][wildcards.target]['path_string']+"/{library}/{library}.{libcon_method}.bam",
#         reference = lambda wildcards:config['metrics']['collect_wgs_metrics'][wildcards.target]['fasta'],
#     output:
#         "results/collect_wgs_metrics/{target}/{library}/{library}.{libcon_method}.metrics.txt"
#     params:
#         picard_jar=config["picard_jar"],
#     log:
#         "logs/collect_wgs_metrics/{target}/{library}.{libcon_method}.log"
#     shell:
#         "java -Xmx4G -XX:ParallelGCThreads=1 -Djava.io.tmpdir=${{TMPDIR}} -jar {params.picard_jar} "+
#         "CollectWgsMetrics I={input.bam} O={output} R={input.reference} 2> {log}" 

# rule flag_stat:
#     input:
#         bam=lambda wildcards:"results/"+config['metrics']['collect_wgs_metrics'][wildcards.target]['path_string']+"/{library}/{library}.{libcon_method}.bam",
#     output:
#         "results/flag_stat/{target}/{library}/{library}.{libcon_method}.metrics.txt",
#     params:
#         gatk=config["gatk"],
#         mem=config["default_gatk_mem"],
#     log:
#         "logs/flag_stat/{library}.{libcon_method}.{target}.collect_hs_metrics.log"
#     shell:
#         "{params.gatk} --java-options \"-Xmx{params.mem} -Djava.io.tmpdir=${{TMPDIR}}\" FlagStat "+\
#         "-I {input.bam} > {output} 2> {log}"