

def collect_valid_predictor_files(wildcards):
    # tf, tissue, = os.path.join(PREDICTORS_DIR, "{wildcards.tf}_{wildcards.tissue}.predictors.txt")
    chckd = checkpoints.create_training_set.get(**wildcards).output[0]
    AA, BB, = glob_wildcards(os.path.join(chckd, '{tf}_{tissue}.predictors.txt'))
    mmd = expand(os.path.join(chckd, f"{tf}_{tissue}.predictors.txt"), tf = AA, tissue = BB) # for t, tt in zip(tf, tissue)]
    #return(os.path.join(PREDICTORS_DIR, "{tf}_{tissue}.predictors.txt"))
    print(mmd)
    return(mmd)

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

rule find_homer_motifs:
    # input:
    #     os.path.join(config["homer"]['motifs_database'], '{motif_file}')
    output: 
        os.path.join(HOMERFILES_DIR, 'scannedMotifs', 'scanMotifsGenomeWide.{motif_file}.txt')
    params:
        run = run,
        homer_cmd = config['homer']['scanMotifsGenome'],
        genome = config['homer']['genome'],
        mfile = os.path.join(config['homer']['motifs_database'], '{motif_file}'),
        #ofile = lambda output: os.path.join(output, 'scanMotifsGenomeWide_{motif_file}.txt'),
        jobname = '{motif_file}'
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
        lambda wildcards: expand(os.path.join(HOMERFILES_DIR, 'scannedMotifs', 'scanMotifsGenomeWide.{motif_file}.txt'), tf = set(homer_motifs_dict.keys()), motif_file = homer_motifs_dict[wildcards.tf])
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
    input: 
        rules.merge_homer_motifs.output
        #os.path.join(HOMERFILES_DIR, '{tf}', 'merged_motif_file.txt')
    output:
        f1=os.path.join(PREDICTORS_DIR, '{tf}_{tissue}.predictors.txt'),
        f2=os.path.join(PREDICTORS_DIR, '{tf}_{tissue}.ground_truth.txt')
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
        f3=os.path.join(PREDICTORS_DIR, '{tf}_{tissue}.info.txt.gz'),
        f4=os.path.join(PREDICTORS_DIR, '{tf}_{tissue}.summary.txt'),
        mfile = rules.merge_homer_motifs.output
    message: "working on {wildcards}"
    benchmark: os.path.join(f"data/{run}/benchmark/{{tf}}_{{tissue}}.create_training_set.tsv")
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
        {params.rscript} workflow/src/create_training_sets_bedtools.R --transcription_factor {wildcards.tf} --tissue {wildcards.tissue} --predicted_motif_file {params.mfile} --sorted_bedfiles_directory {params.sortedbeds_dir} --bedlinks_directory {params.bedfiles_dir} --predictors_file {output.f1} --ground_truth_file {output.f2} --info_file {params.f3} --peaks_files {params.peaks_files} --test_chromosomes {params.test_chromosomes} --summary_file {params.f4}
        """

checkpoint create_enformer_configuration:
    input: lambda wildcards: rules.create_training_set.output.f1 #lambda wildcards: checkpoints.create_training_set.get(tf=wildcards.tf, tissue=wildcards.tissue).output.f1 
    output: os.path.join(PREDICTION_PARAMS_DIR, f'enformer_config_{run}.{{tf}}_{{tissue}}.json')
    message: "working on {wildcards}"
    benchmark: os.path.join(f"data/{run}/benchmark/{{tf}}_{{tissue}}.create_enformer_configuration.tsv")
    resources:
        partition="caslake"
    params:
        run = run,
        rscript = config['rscript'],
        bdirectives = config['enformer']['base_directives'],
        dset = runname,
        model = config['enformer']['model'],
        fasta_file = config['genome']['fasta'],
        pdir = config['scratch_dir'], #DATA_DIR,
        ddate = rundate,
        jobname = '{tf}_{tissue}',
        personalized_directives = None if 'personalized' not in config.keys() else config['personalized']['directives'] #None if not config['personalized']['directives'] else config['personalized']['directives']
    run:
        if params.personalized_directives is None:
            shell("{params.rscript} workflow/src/create_enformer_config.R --dataset {params.dset} --transcription_factor {wildcards.tf} --tissue {wildcards.tissue} --base_directives {params.bdirectives} --project_directory {params.pdir} --predictors_file {input} --model {params.model} --fasta_file {params.fasta_file} --parameters_file {output} --date {params.ddate}")
        elif params.personalized_directives is not None: # don't delete the outputs
            shell("{params.rscript} workflow/src/create_enformer_config.R --dataset {params.dset} --transcription_factor {wildcards.tf} --tissue {wildcards.tissue} --base_directives {params.bdirectives} --project_directory {params.pdir} --predictors_file {input} --model {params.model} --fasta_file {params.fasta_file} --parameters_file {output} --date {params.ddate} --personalized_directives {params.personalized_directives}")
        #printf 'INFO: This is the input file %s\n' "{input}";

checkpoint predict_with_enformer:
    input:
        #checkpoints.create_enformer_configuration.output
        lambda wildcards: checkpoints.create_enformer_configuration.get(tf=wildcards.tf, tissue=wildcards.tissue).output
    output:
        os.path.join(PREDICTION_PARAMS_DIR, f'aggregation_config_{runname}_{{tf}}_{{tissue}}.json')
    resources:
        partition = "caslake",
        time="12:00:00", 
        mem_cpu=2,
        cpu_task=2
    params:
        run = run,
        jobname = '{tf}_{tissue}',
        enformer_predict_script = config['enformer']['predict'],
        nlines = count_number_of_lines
    message: 
        "working on {params.jobname}"
    benchmark: os.path.join(f"data/{run}/benchmark/{{tf}}_{{tissue}}.predict_with_enformer.tsv")
    # shell:
    #     """
    #         sbatch workflow/src/run_enformer.sbatch {params.enformer_predict_script} {input}
    #     """
    run:
        if params.nlines < 100:
            shell("touch {output}")
        else:
            shell("sbatch workflow/src/run_enformer.sbatch {params.enformer_predict_script} {input}")

rule aggregate_predictions:
    input:
        lambda wildcards: checkpoints.predict_with_enformer.get(tf=wildcards.tf, tissue=wildcards.tissue).output
    output:
        os.path.join(AGGREGATION_DIR, f'{runname}_{config["enformer"]["aggtype"]}_{{tf}}_{{tissue}}.csv.gz')
    message: 
        "working on {wildcards}"
    benchmark: os.path.join(f"data/{run}/benchmark/{{tf}}_{{tissue}}.aggregate_predictions.tsv")
    resources:
        partition="caslake",
        mem_cpu=8,
        cpu_task=8,
        mem_mb=24000
    params:
        run = run,
        jobname = '{tf}_{tissue}',
        aggregation_script = config['enformer']['aggregate'],
        aggtype = config['enformer']['aggtype'],
        output_folder = AGGREGATION_DIR,
        hpc = "caslake",
        parsl_executor = "local",
        delete_enformer_outputs = config["delete_enformer_outputs"],
        nlines = count_number_of_lines
    run:
        if params.nlines < 100:
            shell("touch {output}")
        else:
            if params.delete_enformer_outputs == True:
                shell("python3 {params.aggregation_script} --metadata_file {input} --agg_types {params.aggtype} --output_directory {params.output_folder} --hpc {params.hpc} --parsl_executor {params.parsl_executor} --delete_enformer_outputs")
            elif params.delete_enformer_outputs == False: # don't delete the outputs
                shell("python3 {params.aggregation_script} --metadata_file {input} --agg_types {params.aggtype} --output_directory {params.output_folder} --hpc {params.hpc} --parsl_executor {params.parsl_executor}")

rule prepare_training_data:
    input:
        p1 = rules.aggregate_predictions.output,# os.path.join(AGGREGATION_DIR, f'{runname}_{config["enformer"]["aggtype"]}_{{tf}}_{{tissue}}.csv.gz'), #rules.aggregate_predictions.output,
        p2 = rules.create_training_set.output.f2
    output:
        p1=os.path.join(AGGREGATION_DIR, f'train_{run}_{config["enformer"]["aggtype"]}.{{tf}}_{{tissue}}.prepared.csv.gz'),
        p2=os.path.join(AGGREGATION_DIR, f'test_{run}_{config["enformer"]["aggtype"]}.{{tf}}_{{tissue}}.prepared.csv.gz')
    message: 
        "preparing {wildcards} training and test data"
    benchmark: os.path.join(f"data/{run}/benchmark/{{tf}}_{{tissue}}.prepare_training_data.tsv")
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
    benchmark: os.path.join(f"data/{run}/benchmark/{{tf}}_{{tissue}}.train_TFPred_weights.tsv")
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
        mem_cpu=12
    run:
        if params.nlines < 100:
            shell("touch {output.mlogistic} && touch {output.mlinear}")
        else:
            shell("{params.rscript} workflow/src/train_enet.R --train_data_file {input} --rds_file {params.basename} --nfolds {params.nfolds}; sleep 3")

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
    benchmark: os.path.join(f"data/{run}/benchmark/{{tf}}_{{tissue}}.evaluate_TFPred.tsv")
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

rule compile_statistics:
    input: 
        f1 = expand(os.path.join(MODELS_EVAL_DIR, f'{{tf}}_{{tissue}}_{config["date"]}.logistic.train_eval.txt.gz'), zip, tf = TF_list, tissue = tissue_list),
        f2 = expand(os.path.join(MODELS_EVAL_DIR, f'{{tf}}_{{tissue}}_{config["date"]}.logistic.test_eval.txt.gz'), zip, tf = TF_list, tissue = tissue_list)
        # f1 = expand('{tf}', tf = TF_list),
        # f2 = expand('{tissue}', tissue = tissue_list)
    output:
        os.path.join(STATISTICS_DIR, f'{run}.compiled_stats.txt')
    message:
        "compiling statistics for {wildcards}"
    params:
        run = run,
        jobname = run,
        rscript = config['rscript'],
        input_f1 = ','.join(TF_list),
        input_f2 = ','.join(tissue_list),
        path_pattern = lambda wildcards: os.path.join(MODELS_EVAL_DIR, f'{{1}}_{{2}}_{config["date"]}.logistic.{{3}}_eval.txt.gz'),
        model_path = lambda wildcards: os.path.join(MODELS_DIR, f'{{1}}_{{2}}', f'{{1}}_{{2}}_{config["date"]}.logistic.rds')
    resources:
        mem_mb= 100000,
        partition="caslake",
        time="06:00:00"
    shell:
        """
            {params.rscript} workflow/src/compile_statistics.R --transcription_factors {params.input_f1} --tissues {params.input_f2} --path_pattern {params.path_pattern} --statistics_file {output} --model_path {params.model_path}
        """