process_worker_pool.py  --max_workers=4 -a midway3-login3.rcc.local,192.168.122.1,172.25.0.65,128.135.167.77,127.0.0.1,10.50.247.65,10.50.248.65 -p 0 -c 1.0 -m None --poll 10 --task_port=54259 --result_port=54892 --logdir=/project2/haky/temi/projects/TFPred-snakemake/predict/data/predictions_folder/cistrome_AR_Breast/predictions_2023-04-05/runinfo/001/htex_slurm --block_id=0 --hb_period=30  --hb_threshold=120 --cpu-affinity none --available-accelerators 0 1 2 3