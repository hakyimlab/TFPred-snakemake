#!/bin/sh
# properties = {"type": "single", "rule": "find_homer_motifs", "local": false, "input": ["data/homer_files/FOXA1/foxa1.lncap.motif"], "output": ["data/homer_files/FOXA1/scanMotifsGenomeWide_foxa1.lncap.txt"], "wildcards": {"tf": "FOXA1", "motif_file": "foxa1.lncap"}, "params": {"homer_cmd": "/lus/grand/projects/covid-ct/imlab/users/temi/software/homer/bin/scanMotifsGenome.pl", "genome": "/lus/grand/projects/covid-ct/imlab/users/temi/software/homer/data/genomes/hg38", "jobname": "FOXA1"}, "log": [], "threads": 1, "resources": {"mem_mb": 10000, "mem_mib": 9537, "disk_mb": 1000, "disk_mib": 954, "tmpdir": "<TBD>", "queue": "single-gpu", "nodes": 1, "time": "01:00:00", "project": "covid-ct"}, "jobid": 1, "cluster": {}}
cd '/lus/theta-fs0/projects/covid-ct/imlab/users/temi/projects/TFPred-snakemake' && /lus/theta-fs0/projects/covid-ct/imlab/users/temi/software/conda_envs/TFXcan-tools/bin/python3.11 -m snakemake --snakefile '/lus/theta-fs0/projects/covid-ct/imlab/users/temi/projects/TFPred-snakemake/workflow/rules/common.smk' --target-jobs 'find_homer_motifs:tf=FOXA1,motif_file=foxa1.lncap' --allowed-rules 'find_homer_motifs' --cores 'all' --attempt 2 --force-use-threads  --resources 'mem_mb=10000' 'mem_mib=9537' 'disk_mb=1000' 'disk_mib=954' 'nodes=1' --wait-for-files '/lus/theta-fs0/projects/covid-ct/imlab/users/temi/projects/TFPred-snakemake/.snakemake/tmp.gsoimkhc' 'data/homer_files/FOXA1/foxa1.lncap.motif' --force --keep-target-files --keep-remote --max-inventory-time 0 --nocolor --notemp --no-hooks --nolock --ignore-incomplete --rerun-triggers 'mtime' 'params' 'code' 'software-env' 'input' --skip-script-cleanup  --conda-frontend 'mamba' --wrapper-prefix 'https://github.com/snakemake/snakemake-wrappers/raw/' --printshellcmds  --latency-wait 60 --scheduler 'greedy' --scheduler-solver-path '/lus/theta-fs0/projects/covid-ct/imlab/users/temi/software/conda_envs/TFXcan-tools/bin' --default-resources 'mem_mb=8000' 'disk_mb=max(2*input.size_mb, 1000)' 'tmpdir=system_tmpdir' 'queue=single-gpu' 'nodes=1' 'time="01:00:00"' 'project="covid-ct"' --mode 2 && touch '/lus/theta-fs0/projects/covid-ct/imlab/users/temi/projects/TFPred-snakemake/.snakemake/tmp.gsoimkhc/1.jobfinished' || (touch '/lus/theta-fs0/projects/covid-ct/imlab/users/temi/projects/TFPred-snakemake/.snakemake/tmp.gsoimkhc/1.jobfailed'; exit 1)

