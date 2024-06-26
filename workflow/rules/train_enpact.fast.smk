
def catch(func, *args, handle=lambda e : e, **kwargs):
    try:
        return func(*args, **kwargs)
    except Exception as e:
        return(None)
        #return handle(e)


def count_number_of_lines(wildcards):
    # tf, tissue, = os.path.join(PREDICTORS_DIR, "{wildcards.tf}_{wildcards.tissue}.predictors.txt")
    # chckd = rules.create_training_set.output.f1 #get(**wildcards)
    chckd = os.path.join(PREDICTORS_DIR, f"{wildcards.tf}_{wildcards.tissue}.predictors.txt")
    nlines = catch(int, subprocess.run(f"wc -l {chckd} | awk '{{print $1}}'", shell=True, capture_output = True).stdout.decode().strip())
    if nlines is None:
        return(0)
    return(nlines)

rule all:
    input:
        motif_outputs,
        expand(os.path.join(HOMERFILES_DIR, '{tf}', 'merged_motif_file.txt'), tf = set(TF_list)),
        expand(os.path.join(PREDICTORS_DIR, '{tf}_{tissue}.predictors.txt'), zip, tf=TF_list, tissue=tissue_list),
        expand(os.path.join(PREDICTORS_DIR, '{tf}_{tissue}.ground_truth.txt'), zip, tf=TF_list, tissue=tissue_list),
        expand(os.path.join(AGGREGATION_DIR, f'{runname}_{config["enformer"]["aggtype"]}_{{tf}}_{{tissue}}.csv.gz'), zip, tf = TF_list, tissue = tissue_list),
        expand(os.path.join(AGGREGATION_DIR, f'train_{run}_{config["enformer"]["aggtype"]}.{{tf}}_{{tissue}}.prepared.csv.gz'), zip, tf = TF_list, tissue = tissue_list),
        expand(os.path.join(AGGREGATION_DIR, f'test_{run}_{config["enformer"]["aggtype"]}.{{tf}}_{{tissue}}.prepared.csv.gz'), zip, tf = TF_list, tissue = tissue_list),
        expand(os.path.join(MODELS_DIR, "{tf}_{tissue}", f'{{tf}}_{{tissue}}_{config["date"]}.logistic.rds'), zip, tf = TF_list, tissue = tissue_list),
        expand(os.path.join(MODELS_DIR, "{tf}_{tissue}", f'{{tf}}_{{tissue}}_{config["date"]}.linear.rds'), zip, tf = TF_list, tissue = tissue_list),
        expand(os.path.join(MODELS_EVAL_DIR, f'{{tf}}_{{tissue}}_{config["date"]}.linear.train_eval.txt.gz'), zip, tf = TF_list, tissue = tissue_list),
        expand(os.path.join(MODELS_EVAL_DIR, f'{{tf}}_{{tissue}}_{config["date"]}.logistic.train_eval.txt.gz'), zip, tf = TF_list, tissue = tissue_list),
        expand(os.path.join(MODELS_EVAL_DIR, f'{{tf}}_{{tissue}}_{config["date"]}.linear.test_eval.txt.gz'), zip, tf = TF_list, tissue = tissue_list),
        expand(os.path.join(MODELS_EVAL_DIR, f'{{tf}}_{{tissue}}_{config["date"]}.logistic.test_eval.txt.gz'), zip, tf = TF_list, tissue = tissue_list)

rule find_homer_motifs:
    # input:
    #     os.path.join(config["homer"]['motifs_database'], '{motif_file}')
    output: 
        os.path.join(HOMERFILES_DIR, '{tf}', 'scanMotifsGenomeWide.{motif_file}.txt')
    params:
        run = run,
        homer_cmd = config['homer']['scanMotifsGenome'],
        genome = config['homer']['genome'],
        mfile = os.path.join(config['homer']['motifs_database'], '{motif_file}'),
        #ofile = lambda output: os.path.join(output, 'scanMotifsGenomeWide_{motif_file}.txt'),
        jobname = '{tf}'
    message: "working on {wildcards}" 
    resources:
        mem_mb = 10000
    shell:
        """
        perl {params.homer_cmd} {params.mfile} {params.genome} > {output}
        """

rule merge_homer_motifs:
    input:
        #gatherMotifFiles
        lambda wildcards: expand(os.path.join(HOMERFILES_DIR, wildcards.tf, 'scanMotifsGenomeWide.{motif_file}.txt'), tf = set(homer_motifs_dict.keys()), motif_file = homer_motifs_dict[wildcards.tf])
    output:
        os.path.join(HOMERFILES_DIR, '{tf}', 'merged_motif_file.txt')
    message: "working on {input}" 
    params:
        jobname = '{tf}',
        run = run,
        # ifiles = lambda wildcards: expand(os.path.join(HOMERFILES_DIR, '{tf}','scanMotifsGenomeWide.{motif_file}.txt'), tf = set(TF_list), motif_file = homer_motifs_dict[wildcards.tf]),
        #ifiles = gatherMotifFiles,

    resources:
        mem_mb = 10000
    run:
        with open(output[0], 'w') as outfile:
            for fname in input:
                with open(fname) as infile:
                    for i, line in enumerate(infile):
                        outfile.write(line)


rule create_training_set:
    input: rules.merge_homer_motifs.output
    #     lambda wildcards: os.path.join(HOMERFILES_DIR, wildcards.tf, 'merged_motif_file.txt')
    output:
        #directory(PREDICTORS_DIR)
        f1=os.path.join(PREDICTORS_DIR, '{tf}_{tissue}.predictors.txt'),
        f2=os.path.join(PREDICTORS_DIR, '{tf}_{tissue}.ground_truth.txt'),
        # f3=os.path.join(PREDICTORS_DIR, '{tf}_{tissue}.info.txt.gz'),
        # f4=os.path.join(PREDICTORS_DIR, '{tf}_{tissue}.summary.txt')
    params:
        run = run,
        rscript = config['rscript'],
        tf_tissue = '{tf}_{tissue}',
        bedfiles_dir = config['peaks']['directory'],
        sortedbeds_dir = os.path.join(SORTEDBEDS_DIR, '{tf}_{tissue}'),
        jobname = '{tf}_{tissue}',
        basename = os.path.join(PREDICTORS_DIR, '{tf}_{tissue}'),
        test_chromosomes = ','.join(config['train_by_chromosome']['test']),
        peaks_files = lambda wildcards: ','.join(model_config[wildcards.tf]["peakFiles"][wildcards.tissue]) if isinstance(model_config[wildcards.tf]["peakFiles"][wildcards.tissue], list) else model_config[wildcards.tf]["peakFiles"][wildcards.tissue],
        #f1=os.path.join(PREDICTORS_DIR, '{tf}_{tissue}.predictors.txt'),
        # f2=os.path.join(PREDICTORS_DIR, '{tf}_{tissue}.ground_truth.txt'),
        f3=os.path.join(PREDICTORS_DIR, '{tf}_{tissue}.info.txt.gz'),
        f4=os.path.join(PREDICTORS_DIR, '{tf}_{tissue}.summary.txt'),
        mfile = rules.merge_homer_motifs.output
    message: "working on {wildcards}"
    resources:
        partition = "caslake", #if params.nfiles > 200 else "caslake",
        #attempt: attempt * 100,
        mem_cpu = 4, #lambda wildcards, attempt: attempt * 8,
        nodes = 1,
        # mem_cpu = 4, #lambda wildcards, attempt: attempt * 8,
        # nodes = 1,
        #load= 50 #if resources.partition == 'bigmem' else 1
    threads: 8
    shell:
        """
        {params.rscript} workflow/src/create_training_sets_bedtools.R --transcription_factor {wildcards.tf} --tissue {wildcards.tissue} --predicted_motif_file {params.mfile} --sorted_bedfiles_directory {params.sortedbeds_dir} --bedlinks_directory {params.bedfiles_dir} --predictors_file {output.f1} --ground_truth_file {output.f2} --info_file {params.f3} --peaks_files {params.peaks_files} --test_chromosomes {params.test_chromosomes} --summary_file {params.f4}; sleep 5
        """

rule aggregate_epigenomes:
    input:
        #rules.create_training_set.output.f1
        #os.path.join(PREDICTORS_DIR, '{tf}_{tissue}.predictors.txt')
        rules.create_training_set.output.f1
        # ifNotRunEnformer = collect_valid_predictor_files()
    output:
        os.path.join(AGGREGATION_DIR, f'{runname}_{config["enformer"]["aggtype"]}_{{tf}}_{{tissue}}.csv.gz')
    message: 
        "working on {wildcards}"
    resources:
        partition="caslake",
        mem_cpu=8,
        cpu_task=8,
        mem_mb=24000
    params:
        run = run,
        jobname = '{tf}_{tissue}',
        enformer_epigenome_directory = config['enformer_epigenome_directory'],
        nlines = count_number_of_lines,
    run:
        if params.nlines < 100:
            shell("touch {output}")
        else:
            shell("python3 workflow/src/aggregate_epigenomes.py --loci_file {input} --reference_epigenome_dir {params.enformer_epigenome_directory} --output_file {output} --use_multiprocessing")


rule prepare_training_data:
    input:
        p1 = rules.aggregate_epigenomes.output,
        p2 = rules.create_training_set.output.f2
    output:
        p1=os.path.join(AGGREGATION_DIR, f'train_{run}_{config["enformer"]["aggtype"]}.{{tf}}_{{tissue}}.prepared.csv.gz'),
        p2=os.path.join(AGGREGATION_DIR, f'test_{run}_{config["enformer"]["aggtype"]}.{{tf}}_{{tissue}}.prepared.csv.gz')
    message: 
        "preparing {wildcards} training and test data"
    params:
        run = run,
        jobname = '{tf}_{tissue}',
        rscript = config['rscript'],
        aggtype = config['enformer']['aggtype'],
        nlines = count_number_of_lines
    resources:
        mem_mb=24000,
        partition="caslake",
    run:
        if params.nlines < 100:
            shell("touch {output.p1} && touch {output.p2}")
        else:
            shell("{params.rscript} workflow/src/train_test_split.R --data_file {input.p1} --ground_truth_file {input.p2} --aggregation_method {params.aggtype} --train_prepared_file {output.p1} --test_prepared_file {output.p2}")

rule train_TFPred_weights:
    input: rules.prepare_training_data.output.p1
    output: 
        mlogistic=os.path.join(MODELS_DIR, "{tf}_{tissue}", f'{{tf}}_{{tissue}}_{config["date"]}.logistic.rds'),
        mlinear=os.path.join(MODELS_DIR, "{tf}_{tissue}", f'{{tf}}_{{tissue}}_{config["date"]}.linear.rds')
    message:
        "training on {wildcards} training data"
    params:
        run = run,
        jobname = '{tf}_{tissue}',
        rscript = config['rscript'],
        nfolds=5,
        basename=os.path.join(MODELS_DIR, "{tf}_{tissue}", f'{{tf}}_{{tissue}}_{config["date"]}'),
        nlines = count_number_of_lines
    resources:
        #mem_mb=100000,
        mem_mb = helpers.get_mem_mb_allocations,
        partition = helpers.get_cluster_allocation,
        cpu_task=12,
        mem_cpu=12,
        load=50
    run:
        if params.nlines < 100:
            shell("touch {output.mlogistic} && touch {output.mlinear}")
        else:
            shell("{params.rscript} workflow/src/train_enet.R --train_data_file {input} --rds_file {params.basename} --nfolds {params.nfolds}; sleep 12")

rule evaluate_TFPred:
    input: 
        linear_model = rules.train_TFPred_weights.output.mlinear,
        logistic_model = rules.train_TFPred_weights.output.mlogistic,
        train_data = rules.prepare_training_data.output.p1,
        test_data = rules.prepare_training_data.output.p2
    output: 
        os.path.join(MODELS_EVAL_DIR, f'{{tf}}_{{tissue}}_{config["date"]}.logistic.train_eval.txt.gz'),
        os.path.join(MODELS_EVAL_DIR, f'{{tf}}_{{tissue}}_{config["date"]}.logistic.test_eval.txt.gz'),
        os.path.join(MODELS_EVAL_DIR, f'{{tf}}_{{tissue}}_{config["date"]}.linear.train_eval.txt.gz'),
        os.path.join(MODELS_EVAL_DIR, f'{{tf}}_{{tissue}}_{config["date"]}.linear.test_eval.txt.gz')
    message:
        "evaluating on {wildcards} training and test data"
    params:
        run = run,
        jobname = '{tf}_{tissue}',
        rscript = config['rscript'],
        basename=os.path.join(MODELS_EVAL_DIR, f'{{tf}}_{{tissue}}_{config["date"]}'),
        nlines = count_number_of_lines
    resources:
        mem_mb= 100000,
        partition="caslake"
    run:
        if params.nlines < 100:
            shell("touch {output}")
        else:
            shell("{params.rscript} workflow/src/evaluate_enet.R --linear_model {input.linear_model} --logistic_model {input.logistic_model} --train_data_file {input.train_data} --test_data_file {input.test_data} --eval_output {params.basename}")