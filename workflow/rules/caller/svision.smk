rule svision:
    input:
        bam=f"{MAPPER}/{{sample}}/{{sample}}.sorted.bam",
        csi=f"{MAPPER}/{{sample}}/{{sample}}.sorted.bam.csi",
        fasta=config["fasta"],
    output:
        vcf=f"svision/{{sample}}/chroms/{{chrom}}/{{sample}}.{{chrom}}.svision.s{config['min_reads']}.graph.vcf",
        graphs=directory("svision/{sample}/chroms/{chrom}/graphs"),
        segments=temp(directory("svision/{sample}/chroms/{chrom}/segments")),
        predict_results=temp(directory("svision/{sample}/chroms/{chrom}/predict_results")),
    log:
        "logs/{sample}/svision/{chrom}.log",
    wildcard_constraints:
        chrom=r"|".join(re.escape(chrom) for chrom in CHROMS_SVISION),
    container:
        "docker://jiadongxjtu/svision:latest"
    threads: 1
    params:
        dir="svision/{sample}/chroms/{chrom}",
        model=config["model_svision"],
        min_reads=config["min_reads"],
        min_quality_mapping=config["min_quality_mapping"],
        min_size=config["min_size"],
        max_size=config["max_size"],
    shell:
        """
        {{ SVision \\
            -t {threads} \\
            -s {params.min_reads} \\
            --min_mapq {params.min_quality_mapping} \\
            --min_sv_size {params.min_size} \\
            --max_sv_size {params.max_size} \\
            --qname \\
            --graph \\
            --min_gt_depth {params.min_reads} \\
            -o {params.dir} \\
            -b {input.bam} \\
            -m {params.model} \\
            -g {input.fasta} \\
            -n {wildcards.sample}.{wildcards.chrom} \\
            -c {wildcards.chrom}

        # SVision exits 0 without writing a VCF when a chromosome has no calls,
        # and removes segments/ and predict_results/ on success.
        if [ ! -f {output.vcf} ]; then
            echo -e "[INFO] SVision wrote no VCF for {wildcards.chrom}, leaving an empty one."
            touch {output.vcf}
        fi
        mkdir -p {output.graphs} {output.segments} {output.predict_results}; }} \\
        1> {log} 2>&1
        """


rule format_svision:
    input:
        vcfs=expand(
            "svision/{{sample}}/chroms/{chrom}/{{sample}}.{chrom}.svision.s{min_reads}.graph.vcf",
            chrom=CHROMS_SVISION,
            min_reads=config["min_reads"],
        ),
    output:
        tab=temp("svision/{sample}/rename.tab"),
        vcf=protected("svision/{sample}/svision.vcf"),
    log:
        "logs/{sample}/format_svision.log",
    conda:
        "../../envs/bcftools.yaml"
    params:
        chroms=CHROMS_SVISION,
    shell:
        """
        {{ chroms=({params.chroms})
        vcfs=({input.vcfs})
        chrom_lead=""
        vcf_lead=""
        for i in "${{!vcfs[@]}}"; do
            if [ -s "${{vcfs[$i]}}" ]; then
                chrom_lead="${{chroms[$i]}}"
                vcf_lead="${{vcfs[$i]}}"
                break
            fi
        done
        if [ -z "${{vcf_lead}}" ]; then
            echo -e "[ERROR] No non-empty SVision VCF found for sample {wildcards.sample}." >&2
            exit 1
        fi

        echo -e "{wildcards.sample}.${{chrom_lead}}\\t{wildcards.sample}" > {output.tab}

        {{ grep '^#' "${{vcf_lead}}"; cat {input.vcfs} | {{ grep -v '^#' || true; }}; }} \\
            | awk '/^##INFO=<ID=GFA_L/ && !f {{print "##INFO=<ID=GFA_ID,Number=.,Type=String,Description=\\"GFA_ID\\">"; f=1}} 1' \\
            | awk 'BEGIN{{FS=OFS="\\t"}} /^#/ || $5 == "<CSV>" {{print; next}} {{split($8, a, ";"); for(i in a) {{if(a[i] ~ /^SVTYPE=/) {{split(a[i], b, "="); if(b[2] == "tDUP") b[2] = "DUP:TANDEM"; gsub("<SV>", "<"b[2]">", $5)}}}}}}1' \\
            | awk '/^##ALT/ && !f {{print "##ALT=<ID=INS,Description=\\"INS\\">\\n##ALT=<ID=INV,Description=\\"INV\\">\\n##ALT=<ID=DUP,Description=\\"DUP\\">\\n##ALT=<ID=DUP:TANDEM,Description=\\"DUP:TANDEM\\">\\n##ALT=<ID=DEL,Description=\\"DEL\\">"; f=1}} 1' \\
            | awk '/^##INFO=<ID=SUPPORT,/ {{sub("Type=String", "Type=Integer")}} 1' \\
            | bcftools reheader -s {output.tab} \\
            | bcftools annotate --set-id '%CHROM\\_%POS\\_%REF\\_%FIRST_ALT' \\
            > {output.vcf}; }} \\
        1> {log} 2>&1
        """
