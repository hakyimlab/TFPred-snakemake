
rule bwa_align_single:
    input:
        #r0="data/raw_samples/trimmed_fastq/{sample}_trimmed.fq.gz"
        r0=rules.trim_fastq_single.output.rOut
    output:
        sam_file=os.path.join(DATA_DIR, "sam_files/{file_basename}.sam")
        # sorted_bam_file=os.path.join(DATA_DIR, "sorted_bam/{file_basename}.bam"),
        # sorted_bai_file=os.path.join(DATA_DIR, "sorted_bam/{file_basename}.bai")
    params:
        sample_name="{file_basename}",
        bwa_index_file=bwa_index_file,
        temporary_dir=TMP_DIR,
        jobname='{file_basename}'
    message: "ALIGNING - single sample for {params.sample_name}"
    threads: 8
    resources:
        mem_mb=20000
    log:
        os.path.join(LOG_DIR, "bwa_align/{file_basename}_alignment.log")
    shell:
        """
        (bwa mem -M -t 8 {params.bwa_index_file} {input.r0} > {output.sam_file}) 2> {log}
        """