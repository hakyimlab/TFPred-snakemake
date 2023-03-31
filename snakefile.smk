include: "workflow/rules/common.smk"

# === one rule to rule them all ===
rule all:
    input:
        # expand("{data_dir}/trimmed_fastq/{fastq_file}", fastq_file=fastq_files, data_dir=DATA_DIR),
        # expand("{data_dir}/sam_files/{sample}.sam", sample=SAMPLES, data_dir=DATA_DIR),
        # expand("{data_dir}/sorted_bam/{sample}.bam", sample=SAMPLES, data_dir=DATA_DIR)
        # expand("{data_dir}/merged_bam/{individual}.bam", individual=individuals, data_dir=DATA_DIR),
        # expand("{data_dir}/deduplicated_bam/{individual}.bam", individual=individuals, data_dir=DATA_DIR),
        # expand("{data_dir}/recalibrated_bam/{individual}_recal_data.table", individual=individuals, data_dir=DATA_DIR),
        # expand("{data_dir}/recalibrated_bam/{individual}.bam", individual=individuals, data_dir=DATA_DIR),
        # expand(os.path.join(GVCF_DIR, "{individual}.g.vcf.gz"), individual=individuals),
        # expand(os.path.join(REBLOCKED_GVCF_DIR, "{individual}.reblocked.g.vcf.gz"), individual=individuals),
        # expand(os.path.join(scratch_folder, "{chrom}_db"), chrom=chromosomes),
        expand(os.path.join(FINAL_VCFS_DIR, f"unphased_chromosomes/{project_name}_{{chrom}}_unphased_genotypes.vcf.gz"), chrom=chromosomes),
        expand(os.path.join(FINAL_VCFS_DIR, f"phased_chromosomes/{project_name}_{{chrom}}_phased_genotypes.vcf.gz"), chrom=chromosomes)
        # expand(os.path.join(FINAL_VCFS_DIR, f"phased/{project_name}_merged_phased_genotypes.vcf.gz"), project_name=project_name)