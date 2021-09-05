## chimp RNA-seq snakemake
configfile: "env/rnaseq_config.yaml"

## star alignment and then feature counts
rule all:
    input:
        expand("data/samples/combined/{cell_line}.fastq.gz", cell_line=config["cell_line"]),
        expand('output/Alingment/{cell_line}/',cell_line=config['cell_line']),
        expand('output/PostAlignment/{cell_line}_featureCounts.txt',cell_line=config['cell_line']),
        expand('output/PostAlignment/{cell_line}_featureCounts_final.txt',cell_line=config['cell_line']),
        'output/PostAlignment/ipsc_chimp_counts.txt'
        
rule star_index:
    input:
        fasta=config['genome'],
        gff=config['gff']
    output:
        directory("data/Index_genome/")
    params:
        limitRAM=45000000000
    shell:
        'STAR \
         --runMode genomeGenerate \
         --limitGenomeGenerateRAM {params.limitRAM} \
         --genomeDir {output} \
         --genomeFastaFiles {input.fasta} \
         --sjdbGTFfile {input.gff}'

rule poolFiles:
    input:
        lambda wildcards: expand('data/samples/{sample}_1.fastq.gz',sample=config["cell_line"][wildcards.cell_line])
    output:
        "data/samples/combined/{cell_line}.fastq.gz"
    shell:
        "cat {input} > {output}"

rule star_map:
    input:
        sample=rules.poolFiles.output,
        index_dir="data/Index_genome/"
    output:
        aligned_reads=directory('output/Alingment/{cell_line}/')
    threads: 10
    shell:
        'STAR \
        --runThreadN {threads} \
        --readFilesIn {input.sample} \
        --readFilesCommand zcat \
        --genomeDir {input.index_dir} \
        --outFilterMultimapNmax 20 \
        --outFilterMatchNmin 16 \
        --outFilterMismatchNoverLmax 0.03 \
        --alignIntronMax 1 \
        --alignSJDBoverhangMin 1 \
        --outSAMtype BAM SortedByCoordinate \
        --outFileNamePrefix {output.aligned_reads} \
        --outReadsUnmapped Fastx'
        
rule featureCounts:
    input:
        bams='output/Alingment/{cell_line}/Aligned.sortedByCoord.out.bam',
        orthologous_exon_file=config['orthoexon']
    output:
        'output/PostAlignment/{cell_line}_featureCounts.txt'
    params:
        threads=10
    shell:
        "featureCounts -T {params.threads} \
        -F 'SAF' -s 0 \
        -a {input.orthologous_exon_file} -o {output} {input.bams}"

rule polishFiles:
    input:
        rules.featureCounts.output
    output:
        tmpfile=temp('output/PostAlignment/{cell_line}_featureCounts.tmp'),
        finalfile='output/PostAlignment/{cell_line}_featureCounts_final.txt'
    run:
        shell( "cut -f1,7,8,9,10,11,12 {input} | tail -n +3 > {output.tmpfile}")
        shell("echo -e 'GeneID\tCounts' | cat - {output.tmpfile} > {output.finalfile}")

## combine all counts into a single 'species' file
rule combineFiles:
    input:
        directory('output/PostAlignment/')
    output:
        combined_counts = 'output/PostAlignment/ipsc_chimp_counts.txt'
    script:
        '../code/combine_counts.py'