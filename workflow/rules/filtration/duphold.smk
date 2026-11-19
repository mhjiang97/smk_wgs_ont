rule duphold:
    input:
        vcf="{caller}/{sample}/{caller}.vcf",
        bam=f"{MAPPER}/{{sample}}/{{sample}}.sorted.bam",
        csi=f"{MAPPER}/{{sample}}/{{sample}}.sorted.bam.csi",
        fasta=config["fasta"],
    output:
        bcf=temp("{caller}/{sample}/{sample}.sorted.bcf"),
        vcf="{caller}/{sample}/{caller}.duphold.vcf",
    log:
        "logs/{sample}/duphold.{caller}.log",
    conda:
        "../../envs/bcftools.yaml"
    threads: 1
    shell:
        """
        exec 1> {log} 2>&1

        bcftools view -e 'CHROM~"^GL" || ALT~"[][]GL" || POS=0 || ALT~":0[][]"' {input.vcf} \\
            | bcftools sort -Ou - >{output.bcf}

        export DUPHOLD_SAMPLE_NAME={wildcards.sample}
        duphold \\
            -t {threads} \\
            -v {output.bcf} \\
            -b {input.bam} \\
            -f {input.fasta} \\
            -o {output.vcf}
        """
