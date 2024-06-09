
import os


def get_mem_mb_allocations(wildcards, attempt):
    if attempt > 5:
        return(attempt * 10000)
    else:
        return(attempt * 10000)

def get_cluster_allocation(wildcards, attempt):
    if attempt > 5:
        return('bigmem')
    else:
        return('beagle3')
    


def collectMotifFiles(TFs_list, model_config):
    motif_dict = {tf: model_config[tf]['motifFiles'] for tf in TFs_list}
    return(motif_dict)

def collectPeakFiles(model_details, model_config):
    peak_dict = dict()
    for det in model_details:
        tf = det[0]
        tissue = det[1]
        peak_file = model_config[tf]['peakFiles'][tissue]
        if ' ' in tissue:
            tissue = tissue.replace(' ', '-')
        mname = f'{tf}_{tissue}'
        peak_dict[mname] = peak_file
        return(peak_dict)
    
def createMotifInputs(motifsDict):
    res = []
    for key, values in motifsDict.items():
        if isinstance(values, list):
            res.extend(values)
        elif isinstance(values, str):
            res.append(values)
    return(res)

def createMotifOutputs(motifsDict, hdir):
    res = []
    for key, values in motifsDict.items():
        if isinstance(values, list):
            res.extend([os.path.join(hdir, f'{key}', f'scanMotifsGenomeWide.{v}.txt') for v in values])
        elif isinstance(values, str):
            res.append(os.path.join(hdir, f'{key}', f'scanMotifsGenomeWide.{values}.txt'))
    return(res)