##################
# Snakemake rules
##################

rule download_latest_ECDC:
    output:
        ecdc="data/PREPROCESSING/results/ecdc.csv",
    shell:
        """
        wget -c https://opendata.ecdc.europa.eu/covid19/nationalcasedeath/csv -O {output.ecdc}
        """

rule belgium_seqs_fasta:
    input:
        fasta=config['sequences'],
        meta=config['raw_metadata'],
    output:
        be_fasta="data/PREPROCESSING/results/belgium_only.fasta",
    shell:
        """
        python data/PREPROCESSING/scripts/get_be_only.py \
            --sequences {input.fasta} \
            --metadata {input.meta} \
            --output {output.be_fasta}
        """

checkpoint partition_sequences_bel:
    input:
        sequences = rules.belgium_seqs_fasta.output.be_fasta,
    output:
        split_sequences = directory("data/PREPROCESSING/results/split_sequences/")
    log:
        ".logs/partition_sequences.txt"
    params:
        sequences_per_group = 100
    shell:
        """
        python3 scripts/partition-sequences.py \
            --sequences {input.sequences} \
            --sequences-per-group {params.sequences_per_group} \
            --output-dir {output.split_sequences} 2>&1 | tee {log}
        """

rule align_bel:
    message:
        """
        Aligning sequences to {input.reference}
          - gaps relative to reference are considered real
        Cluster:  {wildcards.cluster}
        """
    input:
        sequences = "data/PREPROCESSING/results/split_sequences/{cluster}.fasta",
        reference = config["reference_preproc"],
    output:
        alignment = "data/PREPROCESSING/results/split_alignments/{cluster}.fasta"
    log:
        "logs/align_{cluster}.txt"
    threads: 2
    shell:
        """
        augur align \
            --sequences {input.sequences} \
            --reference-sequence {input.reference} \
            --output {output.alignment} \
            --nthreads {threads} \
            --remove-reference 2>&1 | tee {log}
        """

def _get_alignments_bel(wildcards):
    checkpoint_output = checkpoints.partition_sequences_bel.get(**wildcards).output[0]
    return expand("data/PREPROCESSING/results/split_alignments/{i}.fasta",
                  i=glob_wildcards(os.path.join(checkpoint_output, "{i}.fasta")).i)

rule aggregate_alignments_bel:
    message: "Collecting alignments"
    input:
        alignments = _get_alignments_bel
    output:
        alignment = "data/PREPROCESSING/results/be_aligned.fasta"
    log:
        "logs/aggregate_alignments.txt"
    shell:
        """
        cat {input.alignments} > {output.alignment} 2> {log}
        """

rule pangolin_belgium:
    input:
        be_fasta=rules.aggregate_alignments_bel.output.alignment,
    output:
        be_pangolin="data/PREPROCESSING/results/be_pangolin.csv",
    params:
        pango_threads=config["pangolin_threads"],
    conda:
        "pangolin.yaml"
    # having issues with conda so for run everything from the env w pangolin
    # conda:
    #     #"data/PREPROCESSING/config/pangolin.yaml"
    shell:
        """
        pangolin --threads {params.pango_threads} --outfile {output.be_pangolin} {input.be_fasta}
        """

rule update_nextmeta:
    input:
        be_aln=rules.aggregate_alignments_bel.output.alignment,
        be_pango=rules.pangolin_belgium.output.be_pangolin,
        cases=rules.download_latest_ECDC.output.ecdc,
        nextmeta=config['raw_metadata'],
    output:
        newmeta=config['metadata'],
    shell:
        """
        python data/PREPROCESSING/scripts/modify_metadata_for_subsampling.py \
            --metadata {input.nextmeta} \
            --pangolin {input.be_pango} \
            --cases {input.cases} \
            --seqs {input.be_aln} \
            --out {output.newmeta}
        """
