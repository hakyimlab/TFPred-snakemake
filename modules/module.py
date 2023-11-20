

import os, glob, re, sys, json

print_progress = False

# linking motifs
def find_motif_files(tf, homer_path):
    alltfs = [f for f in os.listdir(homer_path) if os.path.isfile(os.path.join(homer_path, f))]
    pattern = f'\w*-{tf}\.motif$|^{tf}.motif$|{tf}-\w*\.motif$|^{tf}\.\w+\.motif$'
    output = [os.path.join(homer_path, i) for i in alltfs if re.match(pattern, i)]
    return(output)

# homer motifs
def link_homer_motifs(TF, tissue, from_db, to_db):
    tissue = str(tissue)
    if ' ' in tissue:
        tissue = tissue.replace(' ', '-')
    if not os.path.isdir(to_db):
        os.makedirs(to_db)
    files_to_link = find_motif_files(TF.lower(), from_db)
    if len(files_to_link) >= 1:
        motif_files_dict = {}
        if print_progress:
            print(f'INFO - Found {len(files_to_link)} motif files to be linked for {TF}')
        for f in files_to_link:
            dst_file = os.path.join(to_db, os.path.basename(f))
            if not os.path.isfile(dst_file):
                os.symlink(f, dst_file)
            if not f'{TF}_{tissue}' in motif_files_dict.keys():
                motif_files_dict[f'{TF}_{tissue}'] = [dst_file]
            else:
                motif_files_dict[f'{TF}_{tissue}'].append(dst_file)
    elif len(files_to_link) == 0:
        if print_progress:
            print(f'INFO - Found no motif files to be linked for {TF}')
        return(None)
    return(motif_files_dict)    



def return_cistrome_ids(TF, tissue, db_df):
    ids = db_df.loc[(db_df['Factor'] == TF) & (db_df['Tissue_type'] == tissue)].loc[:,'DCid'].tolist()
    return(ids)

def return_cistrome_bedfiles(TF_ids, from_db, suffix=None, prefix=None):
    patterns = [f'{os.path.join(from_db, f"{id}_*.bed")}' for id in TF_ids]
    patterns = [glob.glob(p) for p in patterns]
    patterns = [p for p in patterns for p in p]
    return(patterns)


def link_cistrome_bedfiles(TF, tissue, from_db, to_db, db_df):
    ids = return_cistrome_ids(TF, tissue, db_df=db_df)
    if ' ' in tissue:
        tissue = tissue.replace(' ', '-')
    bedfiles_paths = return_cistrome_bedfiles(TF_ids=ids, from_db=from_db)
    if len(bedfiles_paths) > 0:
        if print_progress:
            print(f'INFO - Found {len(bedfiles_paths)} bedfiles to be linked for {TF}-{tissue} pair')
        bedfiles_ids = {f'{TF}_{tissue}': bedfiles_paths}
        if not os.path.isdir(to_db):
            os.makedirs(to_db)
        for key, values in bedfiles_ids.items():
            for value in values:
                dst_file = os.path.join(to_db, os.path.basename(value))
                if not os.path.isfile(dst_file):
                    os.symlink(value, dst_file)
        return(bedfiles_ids)
    else:
        if print_progress:
            print(f'INFO - No bedfiles found for {TF}')
        return(None)
    
def group_tf_motif_files(dts):
    grouping_dict = {}
    for dt in dts:
        if dt[0] not in grouping_dict.keys():
            grouping_dict[dt[0]] = [dt[1]]
        else:
            grouping_dict[dt[0]].append(dt[1])
    return(grouping_dict)


def group_tf_tissues(dts):
    grouping_dict = {}
    for dt in dts:
        if dt[0] not in grouping_dict.keys():
            grouping_dict[dt[0]] = [dt[1]]
        else:
            grouping_dict[dt[0]].append(dt[1])
    return(grouping_dict)



# need to read in the predictions folder
def return_prediction_folder(json_file_path):
    with open(json_file_path) as f:
        data = json.load(f)
        p = data['predictions_folder']
    return(p)

def return_prediction_date(json_file_path):
    with open(json_file_path) as f:
        data = json.load(f)
        p = data['run_date']
    return(p)