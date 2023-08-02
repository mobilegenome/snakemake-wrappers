"""Snakemake wrapper for Bismark aligned reads deduplication using deduplicate_bismark."""
# https://github.com/FelixKrueger/Bismark/blob/master/deduplicate_bismark

__author__ = "Roman Chernyatchik"
__copyright__ = "Copyright (c) 2019 JetBrains"
__email__ = "roman.chernyatchik@jetbrains.com"
__license__ = "MIT"

import os
import shutil

from pathlib import Path
from snakemake.shell import shell
from tempfile import TemporaryDirectory

bam_path = snakemake.output.get("bam", None)
report_path = snakemake.output.get("report", None)
if not bam_path or not report_path:
    raise ValueError(
        "bismark/deduplicate_bismark: Please specify both 'bam=..' and 'report=..' paths in output section"
    )

params_extra = snakemake.params.get("extra", "")
log = snakemake.log_fmt_shell(stdout=True, stderr=True)
log_append = snakemake.log_fmt_shell(stdout=True, stderr=True, append=True)


with TemporaryDirectory() as tempdir:
    tempdir = Path(tempdir)

    input_files = []
    for input_path in snakemake.input:
    
        input_bam = Path(tempdir, os.path.basename(input_path))
        # copy input to tempdir
        shutil.copy(os.path.basename(input_path), input_bam)
        input_files.append(str(input_bam))

    input_files = " ".join(input_files)

    shell(
        "deduplicate_bismark {params_extra} --bam "
        " --output_dir {tempdir} {input_files} {log}"
    )

    # Move outputs into proper position.
    fst_input_filename = os.path.basename(snakemake.input[0])
    fst_input_basename = os.path.splitext(fst_input_filename)[0]
    prefix = os.path.join(tempdir, fst_input_basename)

    deduplicated_bam_actual_name = prefix + ".deduplicated.bam"

    expected_2_actual_paths = [
        (bam_path, deduplicated_bam_actual_name),
        (
            report_path,
            prefix + ".deduplication_report.txt",
        ),
    ]
    for exp_path, actual_path in expected_2_actual_paths:
        if exp_path and (exp_path != actual_path):
            shell("mv {actual_path:q} {exp_path:q} {log_append}")
