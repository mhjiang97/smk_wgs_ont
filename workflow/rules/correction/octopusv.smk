rule octopusv_correct_bnds:
    input:
        vcf="{caller}/{sample}/{caller}.BND.vcf",
    output:
        svcf="{caller}/{sample}/{caller}.BND.svcf",
        vcf="{caller}/{sample}/{caller}.BND.octopusv.vcf",
    log:
        "logs/{sample}/octopusv_correct_bnds.{caller}.log",
    conda:
        "../../envs/octopusv.yaml"
    shell:
        """
        exec 1>{log} 2>&1

        octopusv correct {input.vcf} {output.svcf}
        octopusv svcf2vcf -i {output.svcf} -o {output.vcf}
        """


checkpoint split_octopusv_calls:
    input:
        vcf="{caller}/{sample}/{caller}.BND.octopusv.vcf",
    output:
        vcf="{caller}/{sample}/{caller}.octopusv.{type_sv}.vcf",
    log:
        "logs/{sample}/split_octopusv_calls.{caller}.{type_sv}.log",
    conda:
        "../../envs/bcftools.yaml"
    params:
        formula=get_split_octopusv_calls_formula,
    shell:
        """
        bcftools filter -i '{params.formula}' {input.vcf} 1>{output.vcf} 2>{log}
        """


rule concat_corrected_calls:
    input:
        vcf="{caller}/{sample}/{caller}.octopusv.{type_sv}.vcf",
        vcf_former="{caller}/{sample}/{caller}.{type_sv}.vcf",
    output:
        tmp_id=temp(touch("{caller}/{sample}/{type_sv}.id")),
        tmp_vcf=temp(
            touch("{caller}/{sample}/{caller}.octopusv.{type_sv}.vcf.tmp.gz")
        ),
        tmp_vcf_former=temp(touch("{caller}/{sample}/{caller}.{type_sv}.vcf.tmp.gz")),
        tmp_tbi=temp(
            touch("{caller}/{sample}/{caller}.octopusv.{type_sv}.vcf.tmp.gz.tbi")
        ),
        tmp_tbi_former=temp(
            touch("{caller}/{sample}/{caller}.{type_sv}.vcf.tmp.gz.tbi")
        ),
        vcf="{caller}/{sample}/{caller}.{type_sv}.corrected.vcf",
    log:
        "logs/{sample}/concat_corrected_calls.{caller}.{type_sv}.log",
    conda:
        "../../envs/bcftools.yaml"
    params:
        to_concatenate=whether_to_concatenate,
    shell:
        """
        exec 1>{log} 2>&1

        if [ "{params.to_concatenate}" == "False" ]; then
            cp {input.vcf_former} {output.vcf}
            exit 0
        fi

        if [ "{wildcards.type_sv}" == "BND" ]; then
            {{ grep -v '^#' {input.vcf} || true; }} \\
                | cut -f 3 - \\
                >{output.tmp_id}
            bcftools view -i "ID=@{output.tmp_id}" {input.vcf_former} \\
                | bcftools sort -O v - \\
                >{output.vcf}

        else
            bcftools sort -O z {input.vcf} >{output.tmp_vcf}
            bcftools sort -O z {input.vcf_former} >{output.tmp_vcf_former}
            tabix -p vcf {output.tmp_vcf}
            tabix -p vcf {output.tmp_vcf_former}
            bcftools concat \\
                -a \\
                {output.tmp_vcf_former} \\
                {output.tmp_vcf} \\
                >{output.vcf}
        fi
        """
