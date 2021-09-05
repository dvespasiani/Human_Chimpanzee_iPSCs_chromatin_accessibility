## use this to first pool together SRA tech rep files into a single sample file
## and run common analyses, i.e. fastqc and trimmomatic

configfile: "env/config.yaml"

rule all:
  input:  
## FastQC + Trimmomatic
    expand("data/input_files/samples/{samples}_1.fastq.gz",samples=config["samples"]),
    expand("data/input_files/samples/{samples}_2.fastq.gz",samples=config["samples"]),
    expand("data/input_files/FastQC/{samples}_{read}_fastqc.zip", samples=config["samples"],read=['1','2']),
    expand("data/input_files/Trimmed_reads/{samples}_1.fastq.gz", samples=config["samples"]),
    expand("data/input_files/Trimmed_reads/{samples}_2.fastq.gz", samples=config["samples"]),
    expand("data/input_files/Trimmed_reads/{samples}_1_unpaired.fastq.gz", samples=config["samples"]),
    expand("data/input_files/Trimmed_reads/{samples}_2_unpaired.fastq.gz", samples=config["samples"]),

rule poolfastq:
  input:
    R1=lambda wildcards: expand("data/input_files/SRA/{file}_1.fastq.gz",file=config["samples"][wildcards.samples]),
    R2=lambda wildcards: expand("data/input_files/SRA/{file}_2.fastq.gz",file=config["samples"][wildcards.samples])
  output:
    R1="data/input_files/samples/{samples}_1.fastq.gz",
    R2="data/input_files/samples/{samples}_2.fastq.gz"
  run:
    if len(input) > 1:
      shell("cat {input.R1} > {output.R1} & cat {input.R2} > {output.R2}")
    else:
      shell("cp {input.R1} {output.R1} && touch -h {output} & cp {input.R2} {output.R2} && touch -h {output.R2}")

rule fastqc:
  input:
    ["data/input_files/samples/{samples}_1.fastq.gz", "data/input_files/samples/{samples}_2.fastq.gz"]
  output:
    "data/input_files/FastQC/{samples}_1_fastqc.zip",
    "data/input_files/FastQC/{samples}_2_fastqc.zip"
  log:
    "logs/FastQC/pre_trimming/{samples}_pre_trimming.log"
  shell:
    "fastqc {input} --outdir=data/input_files/FastQC/ 2> {log}"

rule trimmomatic:
	input:
    # r1=rules.poolfastq.output.R1,
		# r2=rules.poolfastq.output.R2
		r1="data/input_files/samples/{samples}_1.fastq.gz",
		r2="data/input_files/samples/{samples}_2.fastq.gz"
	output:
		r1="data/input_files/Trimmed_reads/{samples}_1.fastq.gz",
		r2="data/input_files/Trimmed_reads/{samples}_2.fastq.gz",
		# reads where trimming entirely removed the mate
		r1_unpaired="data/input_files/Trimmed_reads/{samples}_1_unpaired.fastq.gz",
		r2_unpaired="data/input_files/Trimmed_reads/{samples}_2_unpaired.fastq.gz"
	log:
		"logs/Trimmomatic/{samples}.log"
	shell:
		"java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.39.jar  PE \
    -phred33 \
    {input.r1} {input.r2} \
    {output.r1} {output.r1_unpaired} \
    {output.r2} {output.r2_unpaired} \
    ILLUMINACLIP:./data/Adapters/NexteraPE-PE.fa:2:30:10 \
    LEADING:20 TRAILING:20 MINLEN:20 -trimlog {log}"