"""Snakemake wrapper for Bismark methylation extractor tool: bismark_methylation_extractor."""
# https://github.com/FelixKrueger/Bismark/blob/master/bismark_methylation_extractor

__author__ = "Roman Chernyatchik"
__copyright__ = "Copyright (c) 2019 JetBrains"
__email__ = "roman.chernyatchik@jetbrains.com"
__license__ = "MIT"


import os
import shutil


from pathlib import Path
from tempfile import TemporaryDirectory
from snakemake.shell import shell

params_extra = snakemake.params.get("extra", "")
cmdline_args = ["bismark_methylation_extractor {params_extra}"]

# trimming options
trimming_options = [
    "ignore",  # meth_bias_r1_5end
    "ignore_3prime",  # meth_bias_r1_3end
    "ignore_r2",  # meth_bias_r2_5end
    "ignore_3prime_r2",  # meth_bias_r2_3end
]
for key in trimming_options:
    value = snakemake.params.get(key, None)
    if value:
        cmdline_args.append("--{} {}".format(key, value))


with TemporaryDirectory() as tempdir:
    tempdir = Path(tempdir)
    input_files = []

    # copy input files and add to cmdline
    for input_path in snakemake.input:   
        input_bam_in_temp = Path(tempdir, os.path.basename(input_path))
        shutil.copy(input_path, input_bam_in_temp)
        input_files.append(str(input_bam_in_temp))
        
    cmdline_args.append("{input_files}")
    cmdline_args.append("-o {tempdir} ")  # always use tempdir for initial output
    # log
    log = snakemake.log_fmt_shell(stdout=True, stderr=True)
    cmdline_args.append("{log}")

    shell(" ".join(cmdline_args))

    key2prefix_suffix = [
        ("mbias_report", ("", ".M-bias.txt")),
        ("mbias_r1", ("", ".M-bias_R1.png")),
        ("mbias_r2", ("", ".M-bias_R2.png")),
        ("splitting_report", ("", "_splitting_report.txt")),
        ("methylome_CpG_cov", ("", ".bismark.cov.gz")),
        ("methylome_CpG_mlevel_bedGraph", ("", ".bedGraph.gz"))
    ]
    
    # add additional files is --comprehensive is set
    if "--comprehensive" in params_extra:
        key2prefix_suffix.extend([
            ("read_base_meth_state_cpg", ("CpG_context_", ".txt.gz")),
            ("read_base_meth_state_chg", ("CHG_context_", ".txt.gz")),
            ("read_base_meth_state_chh", ("CHH_context_", ".txt.gz"))
            ])


    log_append = snakemake.log_fmt_shell(stdout=True, stderr=True, append=True)

    for key, (prefix, suffix) in key2prefix_suffix:
        exp_path = snakemake.output.get(key, None)
        if exp_path:
            if len(snakemake.input) != 1:
                raise ValueError(
                    "bismark/bismark_methylation_extractor: Error: only one BAM file is"
                    " expected in input, but was <{}>".format(snakemake.input)
                )
            bam_file = input_files[0]
            bam_name = os.path.basename(bam_file)
            bam_wo_ext = os.path.splitext(bam_name)[0]

            actual_path = os.path.join(tempdir, prefix + bam_wo_ext + suffix)
            if exp_path != actual_path:
                shell("mv {actual_path:q} {exp_path:q} {log_append}")
