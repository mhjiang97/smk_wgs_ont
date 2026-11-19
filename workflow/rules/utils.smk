def get_targets():
    targets = []

    targets += [
        f"survivor/{sample}/final/{sample}.{type_sv}.merged.vcf"
        for sample in SAMPLES
        for type_sv in TYPES_SV
    ]

    return targets


def get_annotsv_cache_outputs():
    if SPECIES in ["homo_sapiens"]:
        return {
            "dir1": f"{config['cache_annotsv']}/Annotations_Human",
            "dir2": f"{config['cache_annotsv']}/Annotations_Exomiser",
        }
    else:
        raise ValueError("Unsupported species")


def get_annotsv_cache_parameters():
    if SPECIES in ["homo_sapiens"]:
        return {
            "arg_install": "install-human-annotation",
            "dirs": [
                "share/AnnotSV/Annotations_Human",
                "share/AnnotSV/Annotations_Exomiser",
            ],
        }
    else:
        raise ValueError("Unsupported species")


def get_convert_snpeff_arguments(wildcards):
    caller = wildcards.caller

    fields = CALLER2FMTS.get(caller)
    if fields is None:
        raise ValueError("Unsupported caller")

    arg = " ".join(f"GEN[*].{field}" for field in fields)

    return arg


def get_split_octopusv_calls_formula(wildcards):
    type_sv = wildcards.type_sv

    formula = f'INFO/SVTYPE == "{type_sv}"'

    if type_sv == "BND":
        formula = 'INFO/SVTYPE == "BND" | INFO/SVTYPE == "TRA"'

    return formula


def whether_to_concatenate(wildcards):
    sample = wildcards.sample
    type_sv = wildcards.type_sv
    caller = wildcards.caller

    path_vcf_new = checkpoints.split_octopusv_calls.get(
        sample=sample, caller=caller, type_sv=type_sv
    ).output.vcf
    path_vcf_old = f"{caller}/{sample}/{caller}.{type_sv}.vcf"

    l_n = sum(1 for line in open(path_vcf_new) if not line.startswith("#"))
    l_o = sum(1 for line in open(path_vcf_old) if not line.startswith("#"))

    if type_sv == "BND":
        if l_n == l_o:
            return False
    else:
        if l_n == 0:
            return False

    return True
